#include <filesystem>

#include "G4BooleanSolid.hh"
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
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
#include "SysRap/NP.hh"
#include "SysRap/SEvt.hh"
#include "SysRap/STrackInfo.h"
#include "SysRap/spho.h"
#include "SysRap/sphoton.h"
#include "U4/U4Random.hh"
#include "U4/U4StepPoint.hh"
#include "U4/U4Touchable.h"
#include "U4/U4Track.h"

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

struct DetectorConstruction : G4VUserDetectorConstruction
{

    DetectorConstruction(std::filesystem::path gdml_file) : gdml_file_(gdml_file)
    {
    }

    G4VPhysicalVolume *Construct() override
    {
        G4GDMLParser parser;
        parser.Read(gdml_file_.string(), false);

        G4VPhysicalVolume *world = parser.GetWorldVolume();

        G4CXOpticks::SetGeometry(world);

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

  private:
    std::filesystem::path gdml_file_;
};

struct PrimaryGenerator : G4VUserPrimaryGeneratorAction
{
    SEvt *sev;

    PrimaryGenerator(SEvt *sev) : sev(sev)
    {
    }

    void GeneratePrimaries(G4Event *event) override
    {
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
        sev->addEventConfigArray();
        sev->gather();
        sev->endOfEvent(eventID);

        // GPU-based simulation
        G4CXOpticks *gx = G4CXOpticks::Get();
        gx->simulate(eventID, true);
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

    void UserSteppingAction(const G4Step *step)
    {

        if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
            return;

        const G4VProcess *process = step->GetPreStepPoint()->GetProcessDefinedStep();

        if (process == nullptr)
            return;

        const G4Track *track = step->GetTrack();
        G4VPhysicalVolume *pv = track->GetVolume();
        const G4VTouchable *touch = track->GetTouchable();

        spho ulabel = {};
        get_label(ulabel, track);

        const G4StepPoint *pre = step->GetPreStepPoint();
        const G4StepPoint *post = step->GetPostStepPoint();

        G4ThreeVector delta = step->GetDeltaPosition();
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
        U4Track::SetFabricatedLabel<spho>(track);
        spho *label = STrackInfo<spho>::GetRef(track);
        assert(label);
    }

    void PreUserTrackingAction(const G4Track *track) override
    {
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
        : sev(SEvt::HighLevelCreate(SEvt::ECPU)), det_cons_(new DetectorConstruction(gdml_file)),
          prim_gen_(new PrimaryGenerator(sev)), event_act_(new EventAction(sev)), stepping_(new SteppingAction(sev)),
          tracking_(new TrackingAction(sev))
    {
    }

    // Create "global" event
    SEvt *sev;

    G4VUserDetectorConstruction *det_cons_;
    G4VUserPrimaryGeneratorAction *prim_gen_;
    EventAction *event_act_;
    SteppingAction *stepping_;
    TrackingAction *tracking_;
};
