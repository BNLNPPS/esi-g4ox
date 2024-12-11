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
#include "G4SDManager.hh"
#include "SysRap/NP.hh"
#include "SysRap/SEvt.hh"
#include "SysRap/SGenerate.h"
#include "SysRap/STrackInfo.h"
#include "SysRap/spho.h"
#include "SysRap/sphoton.h"
#include "U4/U4Random.hh"
#include "U4/U4StepPoint.hh"
#include "U4/U4Touchable.h"
#include "U4/U4Track.h"
#include "U4/U4VPrimaryGenerator.h"

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
        SEvt *sev = SEvt::Get_ECPU();
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
        int idx_arg = event->GetEventID();
        NP *gs = SEvent::MakeTorchGenstep(idx_arg);
        NP *ph = SGenerate::GeneratePhotons(gs);
        U4VPrimaryGenerator::GeneratePrimaries_From_Photons(event, ph);
        delete ph;
        SEvent::SetGENSTEP(gs);
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

struct SteppingAction : G4UserSteppingAction
{
    SEvt *sev;

    SteppingAction(SEvt *sev) : sev(sev)
    {
    }

    void UserSteppingAction(const G4Step *step)
    {
        // std::cout<<  step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << std::endl;
        if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
            return;
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
    }

    void PreUserTrackingAction(const G4Track *track) override
    {
    }

    void PostUserTrackingAction(const G4Track *track) override
    {
    }
};

struct G4App
{
    G4App(std::filesystem::path gdml_file)
        : sev(SEvt::CreateOrReuse_ECPU()), det_cons_(new DetectorConstruction(gdml_file)),
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
