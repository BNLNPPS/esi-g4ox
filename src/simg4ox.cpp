#include <filesystem>
#include <iostream>
#include <string>

#include <argparse/argparse.hpp>

#include "FTFP_BERT.hh"
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalPhysics.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"

#include "G4CX/G4CXOpticks.hh"
#include "SysRap/NP.hh"
#include "SysRap/SEvt.hh"
#include "SysRap/STrackInfo.h"
#include "SysRap/spho.h"
#include "SysRap/sphoton.h"
#include "U4/U4Random.hh"
#include "U4/U4StepPoint.hh"
#include "U4/U4Touchable.h"
#include "U4/U4Track.h"

using namespace std;


SEvt* sev = nullptr;

struct DetectorConstruction : G4VUserDetectorConstruction {

  DetectorConstruction(filesystem::path gdml_file) : gdml_file_(gdml_file) {}

  G4VPhysicalVolume* Construct() override {
    G4GDMLParser parser;
    parser.Read(gdml_file_.string(), false);

    G4VPhysicalVolume* world = parser.GetWorldVolume();

    G4CXOpticks::SetGeometry(world);

    return world;
  }

 private:

  filesystem::path gdml_file_;
};


struct PrimaryGenerator : G4VUserPrimaryGeneratorAction {

  void GeneratePrimaries(G4Event* event) override {
    NP* photons = NP::Make<float>(0, 4, 4);

    cout << "read photons" << endl;
    photons->load("out/photons.npy");
    //photons->dump();

    size_t n_photons = photons->num_items();
    sphoton* sphotons = reinterpret_cast<sphoton*>(photons->bytes());

    for (int i=0; i < n_photons; i++) {
      sphoton &p = sphotons[i];
      //cout << "val: " << i << ": " << p;

      G4ThreeVector position_mm(p.pos.x, p.pos.y, p.pos.z);
      G4double time_ns = p.time;
      G4ThreeVector direction(p.mom.x, p.mom.y, p.mom.z);
      //direction = direction.unit();
      G4double wavelength_nm = p.wavelength ;
      G4ThreeVector polarization(p.pol.x, p.pol.y, p.pol.z);

      G4PrimaryVertex* vertex = new G4PrimaryVertex(position_mm, time_ns);
      G4double kineticEnergy = h_Planck*c_light/(wavelength_nm*nm) ;

      G4PrimaryParticle* particle = new G4PrimaryParticle(G4OpticalPhoton::Definition());
      particle->SetKineticEnergy(kineticEnergy);
      particle->SetMomentumDirection(direction);
      particle->SetPolarization(polarization);

      vertex->SetPrimary(particle);
      event->AddPrimaryVertex(vertex);
    }

    SEvt* sev = SEvt::Get_ECPU();
    sev->SetInputPhoton(photons);
  }
};


struct EventAction : G4UserEventAction {

  void BeginOfEventAction(const G4Event *event) override {
    sev->beginOfEvent(event->GetEventID());
  }

  void EndOfEventAction(const G4Event *event) override {
    int eventID = event->GetEventID();
    sev->addEventConfigArray();
    sev->gather();
    sev->endOfEvent(eventID);

    // GPU-based simulation
    G4CXOpticks* gx = G4CXOpticks::Get() ;
    gx->simulate(eventID, true ) ;
  }
};


void GetLabel(spho& ulabel, const G4Track* track)
{
    spho* label = STrackInfo<spho>::GetRef(track);
    assert( label && label->isDefined() && "all photons are expected to be labelled" );

    std::array<int,spho::N> a_label ;
    label->serialize(a_label) ;

    ulabel.load(a_label);
}


struct SteppingAction : G4UserSteppingAction {

  void UserSteppingAction(const G4Step* step) {
    static int photon_id = 0;

    if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
      return;

    cout << "XXX optical photon: " << photon_id << " - "
         << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << " - ";

    const G4VProcess* process = step->GetPreStepPoint()->GetProcessDefinedStep();

    if (process == nullptr) {
      cout << "no process defined" << endl;
      return;
    } else {
      cout << process->GetProcessName() << endl;
    }

    const G4Track* track = step->GetTrack();
    G4VPhysicalVolume* pv = track->GetVolume();
    const G4VTouchable* touch = track->GetTouchable();

    spho ulabel = {} ;
    GetLabel( ulabel, track );

    const G4StepPoint* pre = step->GetPreStepPoint() ;
    const G4StepPoint* post = step->GetPostStepPoint() ;

    G4ThreeVector delta = step->GetDeltaPosition();
    double step_mm = delta.mag()/mm  ;

    SEvt* sev = SEvt::Get_ECPU();
    sev->checkPhotonLineage(ulabel);

    sphoton& current_photon = sev->current_ctx.p ;

    // first_flag identified by the flagmask having a single bit (all genflag are single bits, set in beginPhoton)
    bool first_flag = current_photon.flagmask_count() == 1 ;
    if(first_flag)
    {
        //LOG(LEVEL) << " first_flag, track " << track ;
        U4StepPoint::Update(current_photon, pre);   // populate current_photon with pos,mom,pol,time,wavelength
        sev->pointPhoton(ulabel);        // sctx::point copying current into buffers
    }

    bool warn = true ;
    bool is_tir = false ;
    unsigned flag = U4StepPoint::Flag<G4OpBoundaryProcess>(post, warn, is_tir ) ;

    bool is_fastsim_flag = flag == DEFER_FSTRACKINFO ;
    bool is_boundary_flag = OpticksPhoton::IsBoundaryFlag(flag) ;  // SD SA DR SR BR BT
    bool is_surface_flag = OpticksPhoton::IsSurfaceDetectOrAbsorbFlag(flag) ;  // SD SA
    bool is_detect_flag = OpticksPhoton::IsSurfaceDetectFlag(flag) ;  // SD

    current_photon.iindex = is_detect_flag ? U4Touchable::ImmediateReplicaNumber(touch) : U4Touchable::AncestorReplicaNumber(touch) ;

    U4StepPoint::Update(current_photon, post);

    current_photon.set_flag( flag );

    sev->pointPhoton(ulabel);     // save SEvt::current_photon/rec/seq/prd into sevent
  }
};


struct TrackingAction : G4UserTrackingAction
{
  const G4Track* transient_fSuspend_track = nullptr;

  void PreUserTrackingAction_Optical_FabricateLabel( const G4Track* track )
  {
      U4Track::SetFabricatedLabel<spho>(track);
      spho* label = STrackInfo<spho>::GetRef(track);
      assert(label) ;
      //saveOrLoadStates(label->id);  // moved here as labelling happens once per torch/input photon
  }


  void PreUserTrackingAction(const G4Track* track) override {
    spho ulabel = {} ;

    spho* label = STrackInfo<spho>::GetRef(track);

    if( label == nullptr )
    {
        PreUserTrackingAction_Optical_FabricateLabel(track) ;
        label = STrackInfo<spho>::GetRef(track);
    }
    assert( label && label->isDefined() );

    std::array<int,spho::N> a_label ;
    label->serialize(a_label) ;

    ulabel.load(a_label);

    U4Random::SetSequenceIndex(ulabel.id);

    bool resume_fSuspend = track == transient_fSuspend_track ;

    SEvt* sev = SEvt::Get_ECPU();
    //LOG_IF(fatal, sev == nullptr) << " SEvt::Get(1) returned nullptr " ;
    assert(sev);


    if(ulabel.gen() == 0)
    {
        if(resume_fSuspend == false)
        {
            sev->beginPhoton(ulabel);  // THIS ZEROS THE SLOT
        }
        else  // resume_fSuspend:true happens following FastSim ModelTrigger:YES, DoIt
        {
            sev->resumePhoton(ulabel);
        }
    }
    else if( ulabel.gen() > 0 )   // HMM: thats going to stick for reemission photons
    {
        if(resume_fSuspend == false)
        {
            sev->rjoinPhoton(ulabel);
        }
        else   // resume_fSuspend:true happens following FastSim ModelTrigger:YES, DoIt
        {
            sev->rjoin_resumePhoton(ulabel);
        }
    }
  }


  void PostUserTrackingAction(const G4Track* track) override {
    G4TrackStatus tstat = track->GetTrackStatus();

    bool is_fStopAndKill = tstat == fStopAndKill;
    bool is_fSuspend     = tstat == fSuspend;
    bool is_fStopAndKill_or_fSuspend = is_fStopAndKill || is_fSuspend;

    assert( is_fStopAndKill_or_fSuspend );

    spho ulabel = {} ;
    GetLabel( ulabel, track );
    //if(!Enabled(ulabel)) return ; // EIDX, GIDX skipping

    if(is_fStopAndKill)
    {
        SEvt* sev = SEvt::Get_ECPU();
        U4Random::SetSequenceIndex(-1);
        sev->finalPhoton(ulabel);
        transient_fSuspend_track = nullptr ;
    }
    else if(is_fSuspend)
    {
        transient_fSuspend_track = track ;
    }
  }
};


struct G4App
{
  G4App(filesystem::path gdml_file) :
    det_cons_(new DetectorConstruction(gdml_file)),
    prim_gen_(new PrimaryGenerator),
    stepping_(new SteppingAction),
    tracking_(new TrackingAction),
    event_act_(new EventAction)
  {
  }

  G4VUserDetectorConstruction*   det_cons_;
  G4VUserPrimaryGeneratorAction* prim_gen_;
  SteppingAction* stepping_;
  TrackingAction* tracking_;
  EventAction*    event_act_;
};


int main(int argc, char **argv)
{
  argparse::ArgumentParser program("simg4ox", "0.0.0");

  string gdml_file, macro_name;
  bool interactive;

  program.add_argument("-g", "--gdml")
    .help("path to GDML file")
    .default_value(string("geom.gdml"))
    .nargs(1)
    .store_into(gdml_file);

  program.add_argument("-m", "--macro")
    .help("path to G4 macro")
    .default_value(string("run.mac"))
    .nargs(1)
    .store_into(macro_name);

  program.add_argument("-i", "--interactive")
    .help("whether to open an interactive window with a viewer")
    .flag()
    .store_into(interactive);

  try {
    program.parse_args(argc, argv);
  }
  catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program;
    exit(EXIT_FAILURE);
  }

  // Create global event
  sev = SEvt::HighLevelCreate(SEvt::ECPU);

  // Configure Geant4
  // The physics list must be instantiated before other user actions
  G4VModularPhysicsList *physics = new FTFP_BERT;
  physics->RegisterPhysics(new G4OpticalPhysics);

  G4RunManager run_mgr;
  run_mgr.SetUserInitialization(physics);

  G4App* g4app = new G4App(gdml_file);
  run_mgr.SetUserInitialization(g4app->det_cons_);
  run_mgr.SetUserAction(g4app->prim_gen_);
  run_mgr.SetUserAction(g4app->event_act_);
  run_mgr.SetUserAction(g4app->tracking_);
  run_mgr.SetUserAction(g4app->stepping_);
  run_mgr.Initialize();

  G4UIExecutive *uix = nullptr;
  G4VisManager  *vis = nullptr;

  if (interactive) {
    uix = new G4UIExecutive(argc, argv);
    vis = new G4VisExecutive;
    vis->Initialize();
  }

  G4UImanager *ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/execute " + macro_name);

  if (interactive) {
    uix->SessionStart();
  }

  delete uix;

  return EXIT_SUCCESS;
}
