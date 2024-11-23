#include <filesystem>

#include "G4BooleanSolid.hh"
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4CX/G4CXOpticks.hh"
#include "G4SDManager.hh"
#include "SysRap/NP.hh"
#include "SysRap/SEvt.hh"
#include "SysRap/STrackInfo.h"
#include "SysRap/spho.h"
#include "SysRap/sphoton.h"
#include "U4/U4Random.hh"
#include "U4/U4StepPoint.hh"
#include "U4/U4Touchable.h"
#include "U4/U4Track.h"
#include "U4.hh"

#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"
#include "G4Profiler.hh"
#include "G4ProfilerMessenger.hh"

bool IsSubtractionSolid(G4VSolid *solid)
{
    if (!solid)
        return false;

    // Check if the solid is directly a G4SubtractionSolid
    if (dynamic_cast<G4SubtractionSolid *>(solid))
        return true;

    // If the solid is a Boolean solid, check its constituent solids
    G4BooleanSolid *booleanSolid = dynamic_cast<G4BooleanSolid *>(solid);
    if (booleanSolid)
    {
        G4VSolid *solidA = booleanSolid->GetConstituentSolid(0);
        G4VSolid *solidB = booleanSolid->GetConstituentSolid(1);

        // Recursively check the constituent solids
        if (IsSubtractionSolid(solidA) || IsSubtractionSolid(solidB))
            return true;
    }

    // For other solid types, return false
    return false;
}

struct PhotonHit : public G4VHit
{
    PhotonHit() = default;

    PhotonHit(unsigned id, G4double energy, G4double time, G4ThreeVector position, G4ThreeVector direction,
              G4ThreeVector polarization)
        : fid(id), fenergy(energy), ftime(time), fposition(position), fdirection(direction), fpolarization(polarization)
    {
    }

    // Copy constructor
    PhotonHit(const PhotonHit &right)
        : G4VHit(right), fid(right.fid), fenergy(right.fenergy), ftime(right.ftime), fposition(right.fposition),
          fdirection(right.fdirection), fpolarization(right.fpolarization)
    {
    }

    // Assignment operator
    const PhotonHit &operator=(const PhotonHit &right)
    {
        if (this != &right)
        {
            G4VHit::operator=(right);
            fid = right.fid;
            fenergy = right.fenergy;
            ftime = right.ftime;
            fposition = right.fposition;
            fdirection = right.fdirection;
            fpolarization = right.fpolarization;
        }
        return *this;
    }

    // Equality operator
    G4bool operator==(const PhotonHit &right) const
    {
        return (this == &right);
    }

    // Print method
    void Print() override
    {
        G4cout << "Detector id: " << fid << " energy: " << fenergy << " nm"
               << " time: " << ftime << " ns"
               << " position: " << fposition << " direction: " << fdirection << " polarization: " << fpolarization
               << G4endl;
    }

    // Member variables
    G4int fid{0};
    G4double fenergy{0};
    G4double ftime{0};
    G4ThreeVector fposition{0, 0, 0};
    G4ThreeVector fdirection{0, 0, 0};
    G4ThreeVector fpolarization{0, 0, 0};
};

using PhotonHitsCollection = G4THitsCollection<PhotonHit>;

struct PhotonSD : public G4VSensitiveDetector
{
    PhotonSD(G4String name) : G4VSensitiveDetector(name), fHCID(-1)
    {
        G4String HCname = name + "_HC";
        collectionName.insert(HCname);
        G4cout << collectionName.size() << "   PhotonSD name:  " << name << " collection Name: " << HCname << G4endl;
    }

    void Initialize(G4HCofThisEvent *hce) override
    {
        fPhotonHitsCollection = new PhotonHitsCollection(SensitiveDetectorName, collectionName[0]);
        if (fHCID < 0)
        {
            G4cout << "PhotonSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0] << G4endl;
            fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        }
        hce->AddHitsCollection(fHCID, fPhotonHitsCollection);
    }

    G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *) override
    {
        G4Track *theTrack = aStep->GetTrack();
        // Only process optical photons
	if (theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
            return false;

        //G4cout << aStep -> GetTrack() -> GetCreatorProcess() -> GetProcessName() << G4endl; # will be relevant later if we have scintillation


	G4double theEnergy = theTrack->GetTotalEnergy() / CLHEP::eV;

        // Create a new hit (CopyNr is set to 0 as DetectorID is omitted)
        PhotonHit *newHit = new PhotonHit(
            0, // CopyNr set to 0
            theEnergy, theTrack->GetGlobalTime(), aStep->GetPostStepPoint()->GetPosition(),
            aStep->GetPostStepPoint()->GetMomentumDirection(), aStep->GetPostStepPoint()->GetPolarization());

        fPhotonHitsCollection->insert(newHit);
        theTrack->SetTrackStatus(fStopAndKill);
        return true;
    }

    void EndOfEvent(G4HCofThisEvent *) override
    {
        G4int NbHits = fPhotonHitsCollection->entries();
        G4cout << "PhotonSD::EndOfEvent Number of PhotonHits: " << NbHits << G4endl;
    }

    void AddOpticksHits()
    {
        SEvt *sev = SEvt::Get_EGPU();
        unsigned int num_hits = sev->GetNumHit(0);

        for (int idx = 0; idx < int(num_hits); idx++)
        {
            sphoton hit;
            sev->getHit(hit, idx);
            G4ThreeVector position = G4ThreeVector(hit.pos.x, hit.pos.y, hit.pos.z);
            G4ThreeVector direction = G4ThreeVector(hit.mom.x, hit.mom.y, hit.mom.z);
            G4ThreeVector polarization = G4ThreeVector(hit.pol.x, hit.pol.y, hit.pol.z);
            int theCreationProcessid;
            if (OpticksPhoton::HasCerenkovFlag(hit.flagmask))
            {
                theCreationProcessid = 0;
            }
            else if (OpticksPhoton::HasScintillationFlag(hit.flagmask))
            {
                theCreationProcessid = 1;
            }
            else
            {
                theCreationProcessid = -1;
            }
            // to be done, first value of hit is detectorid
            // tbd change wavelength to energy
            std::cout << "Adding hit from Opticks:" << hit.wavelength << " " << position << " " << direction << " "
                      << polarization << std::endl;

            PhotonHit *newHit = new PhotonHit(0, hit.wavelength, hit.time, position, direction, polarization);
            fPhotonHitsCollection->insert(newHit);
        }
    }

  private:
    PhotonHitsCollection *fPhotonHitsCollection{nullptr};
    G4int fHCID;
};

struct DetectorConstruction : G4VUserDetectorConstruction
{
    DetectorConstruction(std::filesystem::path gdml_file) : gdml_file_(gdml_file)
    {
    }

    G4VPhysicalVolume *Construct() override
    {
        parser_.Read(gdml_file_.string(), false);
        G4VPhysicalVolume *world = parser_.GetWorldVolume();

        //G4CXOpticks::SetGeometry(world);
        cudaDeviceReset();
        G4CXOpticks::SetGeometry();
	G4LogicalVolumeStore *lvStore = G4LogicalVolumeStore::GetInstance();

        static G4VisAttributes invisibleVisAttr(false);

        // Check if the store is not empty
        if (lvStore && !lvStore->empty())
        {
            // Iterate over all logical volumes in the store
            for (auto &logicalVolume : *lvStore)
            {
                G4VSolid *solid = logicalVolume->GetSolid();

                // Check if the solid uses subtraction
                if (IsSubtractionSolid(solid))
                {
                    // Assign the invisible visual attributes to the logical volume
                    logicalVolume->SetVisAttributes(&invisibleVisAttr);

                    // Optionally, print out the name of the logical volume
                    G4cout << "Hiding logical volume: " << logicalVolume->GetName() << G4endl;
                }
            }
        }

        return world;
    }

    void ConstructSDandField() override
    {
        G4cout << "ConstructSDandField is called." << G4endl;
        G4SDManager *SDman = G4SDManager::GetSDMpointer();

        const G4GDMLAuxMapType *auxmap = parser_.GetAuxMap();
        for (auto const &[logVol, listType] : *auxmap)
        {
            for (auto const &auxtype : listType)
            {
                if (auxtype.type == "SensDet")
                {
                    G4cout << "Attaching sensitive detector to logical volume: " << logVol->GetName() << G4endl;
                    G4String name = logVol->GetName() + "_PhotonDetector";
                    PhotonSD *aPhotonSD = new PhotonSD(name);
                    SDman->AddNewDetector(aPhotonSD);
                    logVol->SetSensitiveDetector(aPhotonSD);
                }
            }
        }
    }

  private:
    std::filesystem::path gdml_file_;
    G4GDMLParser parser_;
};

struct PrimaryGenerator : G4VUserPrimaryGeneratorAction
{
    SEvt *sev;

    PrimaryGenerator(SEvt *sev) : sev(sev)
    {
    }

    void GeneratePrimaries(G4Event *event) override
    {
        
	/*
	NP *photons = NP::Make<float>(0, 4, 4);

        photons->load("out/photons.npy");

        size_t n_photons = photons->num_items();
        sphoton *sphotons = reinterpret_cast<sphoton *>(photons->bytes());

        for (int i = 0; i < n_photons; i++)
        {
            sphoton &p = sphotons[i];

            G4ThreeVector position_mm(p.pos.x, p.pos.y, p.pos.z);
            G4double time_ns = p.time;
            G4ThreeVector direction(p.mom.x, p.mom.y, p.mom.z);
            // direction = direction.unit();
            G4double wavelength_nm = p.wavelength;
            G4ThreeVector polarization(p.pol.x, p.pol.y, p.pol.z);

            G4PrimaryVertex *vertex = new G4PrimaryVertex(position_mm, time_ns);
            G4double kineticEnergy = h_Planck * c_light / (wavelength_nm * nm);

            G4PrimaryParticle *particle = new G4PrimaryParticle(G4OpticalPhoton::Definition());
            particle->SetKineticEnergy(kineticEnergy);
            particle->SetMomentumDirection(direction);
            particle->SetPolarization(polarization);

            vertex->SetPrimary(particle);
            event->AddPrimaryVertex(vertex);
        }

        sev->SetInputPhoton(photons);
    }*/

    	    G4ThreeVector position_mm(0, 0, -1*m);
            G4double time_ns = 0;
            G4ThreeVector direction(0, 0.2, 0.8);
            // direction = direction.unit();
            G4double wavelength_nm = 0.1;

            G4PrimaryVertex *vertex = new G4PrimaryVertex(position_mm, time_ns);
            G4double kineticEnergy = h_Planck * c_light / (wavelength_nm * nm);
  	    G4PrimaryParticle *particle = new G4PrimaryParticle(G4Electron::Definition());
            particle->SetKineticEnergy(0.05*GeV);
            particle->SetMomentumDirection(direction);
            vertex->SetPrimary(particle);
            event->AddPrimaryVertex(vertex);
}
};

struct EventAction : G4UserEventAction
{
    SEvt *sev;

    EventAction(SEvt *sev) : sev(sev)
    {
    }

    void BeginOfEventAction(const G4Event *event) override
    {
        sev->beginOfEvent(event->GetEventID());
    }

    void EndOfEventAction(const G4Event *event) override
    {
        int eventID = event->GetEventID();
        //sev->addEventConfigArray();
        //sev->gather();
        //sev->endOfEvent(eventID);

        // GPU-based simulation
        G4CXOpticks *gx = G4CXOpticks::Get();

        gx->simulate(eventID, false);
        cudaDeviceSynchronize();
        unsigned int num_hits = SEvt::GetNumHit(0);
        std::cout << "Opticks: NumCollected:  " << SEvt::GetNumGenstepFromGenstep(0) << std::endl;

        std::cout << "Opticks: NumCollected:  " << SEvt::GetNumPhotonCollected(0) << std::endl;

        std::cout << "Opticks: NumHits:  " << num_hits << std::endl;
        if (num_hits > 0)
        {
            G4HCtable *hctable = G4SDManager::GetSDMpointer()->GetHCtable();
            for (G4int i = 0; i < hctable->entries(); ++i)
            {
                std::string sdn = hctable->GetSDname(i);
                std::size_t found = sdn.find("Photondetector");
                if (found != std::string::npos)
                {
                    PhotonSD *aSD = (PhotonSD *)G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdn);
                    aSD->AddOpticksHits();
                }
            }
        }
        G4CXOpticks::Get()->reset(eventID);
    }
};

void get_label(spho &ulabel, const G4Track *track)
{
    spho *label = STrackInfo<spho>::GetRef(track);
    assert(label && label->isDefined() && "all photons are expected to be labelled");

    std::array<int, spho::N> a_label;
    label->serialize(a_label);

    ulabel.load(a_label);
}

struct SteppingAction : G4UserSteppingAction
{
    SEvt *sev;

    SteppingAction(SEvt *sev) : sev(sev)
    {
    }

    void UserSteppingAction(const G4Step *aStep)
    {


    if(aStep->GetTrack()->GetTrackStatus() == fAlive){ 
    if( aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName()  !="MirrorPyramid" && aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName()  =="MirrorPyramid")
    {if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
			G4cout << "I am a photon enetering pyramid from outside. Am I coming from the correct way? Not sure" << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName() << G4endl;
aStep->GetTrack()->SetTrackStatus(fStopAndKill);			
    }}



G4int fNumPhotons = 0;  // number of scintillation photons this step
    const G4Track* aTrack;
    G4SteppingManager* fpSteppingManager =
    G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
    G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();
    if(stepStatus != fAtRestDoItProc)
    {
      G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
      size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
      for(size_t i3 = 0; i3 < MAXofPostStepLoops; i3++)
      {
        if((*procPost)[i3]->GetProcessName() == "Cerenkov")
        {
          std::cout << "Nagec"<< std::endl;
          aTrack                             = aStep->GetTrack();
          const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
          G4double charge                    = aParticle->GetDefinition()->GetPDGCharge();
          const G4Material* aMaterial        = aTrack->GetMaterial();
          G4MaterialPropertiesTable* MPT     = aMaterial->GetMaterialPropertiesTable();
          G4MaterialPropertyVector* Rindex   = MPT->GetProperty(kRINDEX);
          G4Cerenkov* proc                   = (G4Cerenkov*) (*procPost)[i3];
          fNumPhotons                        = proc->GetNumPhotons();
    
    //Photoncounter += fNumPhotons;
          if(fNumPhotons > 0)
          {
            G4double Pmin        = Rindex->Energy(0);
            G4double Pmax        = Rindex->GetMaxEnergy();
            G4double nMax        = Rindex->GetMaxValue();
            G4double beta1       = aStep->GetPreStepPoint()->GetBeta();
            G4double beta2       = aStep->GetPostStepPoint()->GetBeta();
            G4double beta        = (beta1 + beta2) * 0.5;
            G4double BetaInverse = 1. / beta;
            G4double maxCos      = BetaInverse / nMax;
            G4double maxSin2     = (1.0 - maxCos) * (1.0 + maxCos);
            G4double MeanNumberOfPhotons1 =
              proc->GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
            G4double MeanNumberOfPhotons2 =
              proc->GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);
            U4::CollectGenstep_G4Cerenkov_modified(aTrack, aStep, fNumPhotons, BetaInverse, Pmin,
                                                   Pmax, maxCos, maxSin2, MeanNumberOfPhotons1,
                                                   MeanNumberOfPhotons2);
           std::cout << "MeanNumberOfPhotons1" << MeanNumberOfPhotons1 << std::endl;

           SEvt::Check(0);
        cudaDeviceSynchronize();
           G4RunManager* rm     = G4RunManager::GetRunManager();
        const G4Event* event = rm->GetCurrentEvent();
           G4int eventid        = event->GetEventID();
           G4CXOpticks::Get()->simulate(eventid, false);
           cudaDeviceSynchronize();
           unsigned int num_hits = SEvt::GetNumHit(0);
           std::cout <<"gpl" << num_hits << std::endl;

		
	    if(num_hits > 0)
      {
        G4HCtable* hctable = G4SDManager::GetSDMpointer()->GetHCtable();
        for(G4int i = 0; i < hctable->entries(); ++i)
        {
          std::string sdn   = hctable->GetSDname(i);
          std::size_t found = sdn.find("PhotonDetector");
          if(found != std::string::npos)
          {
                      std::cout << "PhotonDetector: " << sdn << std::endl;
            PhotonSD* aSD = (PhotonSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdn);
            aSD->AddOpticksHits();
          }
        }
      }



           G4CXOpticks::Get()->reset(eventid);
            //std << GenStepcounter++;
          }
        }


	        if((*procPost)[i3]->GetProcessName() == "Scintillation")
        {
                std::cout << "Nagec2"<< std::endl;
          G4Scintillation* proc1 = (G4Scintillation*) (*procPost)[i3];
          fNumPhotons            = proc1->GetNumPhotons();
          ////Photoncounter += fNumPhotons;
          std::cout << "fNumPhotons" << fNumPhotons << std::endl;

          G4double timeconst = 0.0;
          if(fNumPhotons > 0)
          {
            aTrack                         = aStep->GetTrack();
            const G4Material* aMaterial    = aTrack->GetMaterial();
            G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
            timeconst                      = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
            const G4Track* aTrack          = aStep->GetTrack();
            //U4::CollectGenstep_DsG4Scintillation_r4695(aTrack, aStep, fNumPhotons, 1, timeconst);
            //GenStepcounter++;
            //
            //

          }
        }
      }
    }

	    
	 return;

	 // std::cout<<  step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << std::endl;
        if (aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
            return;

        const G4VProcess *process = aStep->GetPreStepPoint()->GetProcessDefinedStep();

        if (process == nullptr)
            return;

        const G4Track *track = aStep->GetTrack();
        G4VPhysicalVolume *pv = track->GetVolume();
        const G4VTouchable *touch = track->GetTouchable();

        spho ulabel = {};
        get_label(ulabel, track);

        const G4StepPoint *pre = aStep->GetPreStepPoint();
        const G4StepPoint *post = aStep->GetPostStepPoint();

        G4ThreeVector delta = aStep->GetDeltaPosition();
        double step_mm = delta.mag() / mm;

        sev->checkPhotonLineage(ulabel);

        sphoton &current_photon = sev->current_ctx.p;

        if (current_photon.flagmask_count() == 1)
        {
            U4StepPoint::Update(current_photon, pre); // populate current_photon with pos, mom, pol, time, wavelength
            sev->pointPhoton(ulabel);                 // copying current into buffers
        }

        bool tir;
        unsigned flag = U4StepPoint::Flag<G4OpBoundaryProcess>(post, true, tir);
        bool is_detect_flag = OpticksPhoton::IsSurfaceDetectFlag(flag);

        current_photon.iindex =
            is_detect_flag ? U4Touchable::ImmediateReplicaNumber(touch) : U4Touchable::AncestorReplicaNumber(touch);

        U4StepPoint::Update(current_photon, post);

        current_photon.set_flag(flag);

        sev->pointPhoton(ulabel);
    }
};

struct TrackingAction : G4UserTrackingAction
{
    const G4Track *transient_fSuspend_track = nullptr;
    SEvt *sev;

    TrackingAction(SEvt *sev) : sev(sev)
    {
    }

    void PreUserTrackingAction_Optical_FabricateLabel(const G4Track *track)
    {
	return;
        U4Track::SetFabricatedLabel<spho>(track);
        spho *label = STrackInfo<spho>::GetRef(track);
        assert(label);
    }

    void PreUserTrackingAction(const G4Track *track) override
    {
        
	    return;
	    spho *label = STrackInfo<spho>::GetRef(track);

        if (label == nullptr)
        {
            PreUserTrackingAction_Optical_FabricateLabel(track);
            label = STrackInfo<spho>::GetRef(track);
        }

        assert(label && label->isDefined());

        std::array<int, spho::N> a_label;
        label->serialize(a_label);

        spho ulabel = {};
        ulabel.load(a_label);

        U4Random::SetSequenceIndex(ulabel.id);

        bool resume_fSuspend = track == transient_fSuspend_track;

        if (ulabel.gen() == 0)
        {
            if (resume_fSuspend == false)
                sev->beginPhoton(ulabel);
            else
                sev->resumePhoton(ulabel);
        }
        else if (ulabel.gen() > 0)
        {
            if (resume_fSuspend == false)
                sev->rjoinPhoton(ulabel);
            else
                sev->rjoin_resumePhoton(ulabel);
        }
    }

    void PostUserTrackingAction(const G4Track *track) override
    {
	   return;
        G4TrackStatus tstat = track->GetTrackStatus();

        bool is_fStopAndKill = tstat == fStopAndKill;
        bool is_fSuspend = tstat == fSuspend;
        bool is_fStopAndKill_or_fSuspend = is_fStopAndKill || is_fSuspend;

        assert(is_fStopAndKill_or_fSuspend);

        spho ulabel = {};
        get_label(ulabel, track);

        if (is_fStopAndKill)
        {
            U4Random::SetSequenceIndex(-1);
            sev->finalPhoton(ulabel);
            transient_fSuspend_track = nullptr;
        }
        else if (is_fSuspend)
        {
            transient_fSuspend_track = track;
        }
    }
};

struct G4App
{
    G4App(std::filesystem::path gdml_file)
        : sev(SEvt::CreateOrReuse_EGPU()), det_cons_(new DetectorConstruction(gdml_file)),
          prim_gen_(new PrimaryGenerator(sev)), event_act_(new EventAction(sev)), stepping_(new SteppingAction(sev)),
          tracking_(new TrackingAction(sev))
    {
    }

    //~G4App(){ G4CXOpticks::Finalize();}

    // Create "global" event
    SEvt *sev;

    G4VUserDetectorConstruction *det_cons_;
    G4VUserPrimaryGeneratorAction *prim_gen_;
    EventAction *event_act_;
    SteppingAction *stepping_;
    TrackingAction *tracking_;
};
