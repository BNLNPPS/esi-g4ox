#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <atomic>
#include <vector>

#include "G4TaskGroup.hh"

#include "G4BooleanSolid.hh"
#include "G4CX/G4CXOpticks.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "SysRap/NP.hh"
#include "SysRap/SEvt.hh"
#include "SysRap/STrackInfo.h"
#include "SysRap/spho.h"
#include "SysRap/sphoton.h"
#include "SysRap/scerenkov.h"
#include "SysRap/OpticksGenstep.h"
#include "SysRap/sgs.h"
#include "SysRap/NP.hh"
#include "U4.hh"
#include "U4/U4Random.hh"
#include "U4/U4StepPoint.hh"
#include "U4/U4Touchable.h"
#include "U4/U4Track.h"
#include "G4RunManagerFactory.hh"
#include "G4AutoLock.hh"


namespace { G4Mutex genstep_mutex = G4MUTEX_INITIALIZER; }

// =============================================================================
// ASYNC GPU SIMULATION WITH DOUBLE-BUFFERING
//
// Architecture:
// - GenstepBuffer: Thread-safe buffer that collects gensteps from CPU
// - GPUTaskManager: Manages async GPU processing with worker thread
//
// Flow:
// 1. CPU collects gensteps into "active" buffer
// 2. When threshold reached, swap buffers and queue GPU task
// 3. GPU thread processes queued buffer while CPU fills the other
// 4. Double-buffering allows CPU/GPU overlap
//
// Backpressure:
// - If GPU is slow, queue fills up
// - When queue is full, CPU blocks (prevents OOM)
// - Queue size configurable via GPU_MAX_QUEUE_SIZE env var
// =============================================================================

// Thread-safe genstep buffer for double-buffering
struct GenstepBuffer {
    std::vector<quad6> gensteps;
    std::vector<sgs> labels;
    int64_t photon_count = 0;
    int64_t genstep_count = 0;
    int event_id = 0;

    void clear() {
        gensteps.clear();
        labels.clear();
        photon_count = 0;
        genstep_count = 0;
    }

    void addGenstep(const quad6& gs, int64_t numphotons) {
        sgs label;
        label.index = gensteps.size();
        label.photons = numphotons;
        label.offset = photon_count;
        label.gentype = gs.gentype();

        gensteps.push_back(gs);
        labels.push_back(label);
        photon_count += numphotons;
        genstep_count++;
    }

    bool empty() const { return gensteps.empty(); }
    size_t size() const { return gensteps.size(); }
};

// =============================================================================
// GPUTaskManager using G4TaskGroup for async GPU processing
//
// Uses Geant4's built-in tasking infrastructure instead of custom threading.
// G4TaskGroup handles: thread pool, task queue, synchronization, backpressure.
// =============================================================================

class GPUTaskManager {
public:
    static constexpr int64_t DEFAULT_PHOTON_THRESHOLD = 10000000;  // 10M photons

private:
    int64_t photon_threshold_;

    // Double-buffer: CPU writes to active, GPU processes submitted
    std::shared_ptr<GenstepBuffer> active_buffer_;
    G4Mutex buffer_mutex_;

    // G4TaskGroup manages async tasks - no custom threading needed
    std::unique_ptr<G4TaskGroup<void, void>> task_group_;

    // GPU mutex - only ONE task can use GPU at a time
    G4Mutex gpu_mutex_;

    // Statistics
    std::atomic<int> batch_counter_{0};
    std::atomic<int> completed_batches_{0};
    std::atomic<uint64_t> total_hits_processed_{0};
    std::atomic<uint64_t> total_photons_processed_{0};
    std::atomic<uint64_t> total_gpu_time_us_{0};  // microseconds for atomic

public:
    GPUTaskManager(G4Mutex* /* unused */ = nullptr,
                   int64_t threshold = DEFAULT_PHOTON_THRESHOLD)
        : photon_threshold_(threshold)
        , active_buffer_(std::make_shared<GenstepBuffer>())
        , task_group_(nullptr)
    {
        const char* env_thresh = std::getenv("GPU_PHOTON_FLUSH_THRESHOLD");
        if (env_thresh) {
            photon_threshold_ = std::atoll(env_thresh);
        }
    }

    ~GPUTaskManager() {
        shutdown();
    }

    void start() {
        // Create task group - uses Geant4's thread pool
        task_group_ = std::make_unique<G4TaskGroup<void, void>>();

        G4cout << "GPUTaskManager [G4TaskGroup]: Started" << G4endl;
        G4cout << "  Photon threshold: " << photon_threshold_ << G4endl;
    }

    void shutdown() {
        // Flush any remaining gensteps
        {
            G4AutoLock lock(&buffer_mutex_);
            if (active_buffer_ && !active_buffer_->empty()) {
                submitBuffer(active_buffer_);
                active_buffer_ = std::make_shared<GenstepBuffer>();
            }
        }

        // Wait for all tasks to complete
        waitForCompletion();

        task_group_.reset();

        G4cout << "GPUTaskManager [G4TaskGroup]: Shutdown complete" << G4endl;
        G4cout << "  Total batches: " << completed_batches_.load() << G4endl;
        G4cout << "  Total photons: " << total_photons_processed_.load() << G4endl;
        G4cout << "  Total hits: " << total_hits_processed_.load() << G4endl;
        G4cout << "  Total GPU time: " << (total_gpu_time_us_.load() / 1e6) << "s" << G4endl;
    }

    // Called from SteppingAction to add a genstep
    // This is the HOT PATH - must be fast
    void addGenstep(const quad6& gs, int64_t numphotons, int eventID) {
        std::shared_ptr<GenstepBuffer> buffer_to_submit;

        {
            G4AutoLock lock(&buffer_mutex_);

            active_buffer_->event_id = eventID;
            active_buffer_->addGenstep(gs, numphotons);

            // Check threshold
            if (active_buffer_->photon_count >= photon_threshold_) {
                // Swap buffer and prepare to submit
                buffer_to_submit = active_buffer_;
                active_buffer_ = std::make_shared<GenstepBuffer>();
            }
        }

        // Submit OUTSIDE the lock to avoid blocking CPU threads
        if (buffer_to_submit) {
            submitBuffer(buffer_to_submit);
        }
    }

    // Force flush (called at end of run)
    void flushRemaining(int eventID) {
        std::shared_ptr<GenstepBuffer> buffer_to_submit;

        {
            G4AutoLock lock(&buffer_mutex_);
            if (active_buffer_ && !active_buffer_->empty()) {
                active_buffer_->event_id = eventID;
                buffer_to_submit = active_buffer_;
                active_buffer_ = std::make_shared<GenstepBuffer>();
            }
        }

        if (buffer_to_submit) {
            G4cout << "GPUTaskManager: Final flush of " << buffer_to_submit->photon_count << " photons" << G4endl;
            submitBuffer(buffer_to_submit);
        }

        // Wait for all tasks to complete
        waitForCompletion();
    }

    void waitForCompletion() {
        if (task_group_) {
            task_group_->join();  // G4TaskGroup handles waiting
        }
    }

    // Statistics
    int getCompletedBatches() const { return completed_batches_.load(); }
    uint64_t getTotalHitsProcessed() const { return total_hits_processed_.load(); }
    int64_t getThreshold() const { return photon_threshold_; }
    double getTotalGPUTime() const { return total_gpu_time_us_.load() / 1e6; }

private:
    void submitBuffer(std::shared_ptr<GenstepBuffer> buffer) {
        if (!buffer || buffer->empty()) return;
        if (!task_group_) return;

        int batch_id = batch_counter_.fetch_add(1);
        int event_id = buffer->event_id;

        G4cout << "GPUTaskManager: Queued batch " << batch_id
               << " (" << buffer->photon_count << " photons, "
               << buffer->genstep_count << " gensteps)" << G4endl;

        // Submit task to G4TaskGroup - it handles all the threading
        task_group_->exec([this, batch_id, event_id, buffer]() {
            processGPUTask(batch_id, event_id, buffer);
        });
    }

    void processGPUTask(int batch_id, int event_id, std::shared_ptr<GenstepBuffer> buffer) {
        // CRITICAL: Only one GPU task at a time - G4CXOpticks/SEvt are not thread-safe
        G4AutoLock gpu_lock(&gpu_mutex_);

        G4cout << "=== GPU Batch " << batch_id << " Processing ===" << G4endl;
        G4cout << "  Photons: " << buffer->photon_count << G4endl;
        G4cout << "  Gensteps: " << buffer->genstep_count << G4endl;

        // Get Opticks instances
        G4CXOpticks* gx = G4CXOpticks::Get();
        SEvt* sev = SEvt::Get_EGPU();

        if (!gx || !sev) {
            G4cerr << "GPUTaskManager: G4CXOpticks or SEvt not available!" << G4endl;
            return;
        }

        // Clear any existing gensteps in SEvt
        sev->clear_genstep();

        // Load our buffered gensteps into SEvt using bulk method
        // Create NP array directly from vector data (single memcpy)
        NP* gs_array = NP::Make<float>(buffer->gensteps.size(), 6, 4);
        memcpy(gs_array->values<float>(), buffer->gensteps.data(), buffer->gensteps.size() * sizeof(quad6));
        sev->addGenstep(gs_array);

        // Run GPU simulation
        auto start = std::chrono::high_resolution_clock::now();
        gx->simulate(event_id, false);
        cudaDeviceSynchronize();
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        // Get results
        unsigned int num_hits = sev->GetNumHit(0);

        // Update stats (all atomic)
        total_gpu_time_us_ += elapsed_us;
        total_hits_processed_ += num_hits;
        total_photons_processed_ += buffer->photon_count;

        G4cout << "  GPU time: " << (elapsed_us / 1e6) << "s" << G4endl;
        G4cout << "  Hits: " << num_hits << G4endl;

        // Write hits to file
        if (num_hits > 0) {
            writeHitsToFile(batch_id, sev, num_hits);
        }

        // Reset for next batch
        gx->reset(event_id);

        completed_batches_++;
        G4cout << "=== GPU Batch " << batch_id << " Complete ===" << G4endl;
    }

    void writeHitsToFile(int batchID, SEvt* sev, unsigned int num_hits) {
        std::ostringstream fname;
        fname << "gpu_hits_batch_" << batchID << ".txt";
        std::ofstream outFile(fname.str());

        if (!outFile.is_open()) {
            G4cerr << "GPUTaskManager: Failed to open " << fname.str() << G4endl;
            return;
        }

        outFile << "# Batch " << batchID << ", Hits: " << num_hits << "\n";
        outFile << "# time wavelength pos_x pos_y pos_z mom_x mom_y mom_z pol_x pol_y pol_z process_flag\n";

        for (unsigned int idx = 0; idx < num_hits; idx++) {
            sphoton hit;
            sev->getHit(hit, idx);

            int process_flag = -1;
            if (OpticksPhoton::HasCerenkovFlag(hit.flagmask))
                process_flag = 0;
            else if (OpticksPhoton::HasScintillationFlag(hit.flagmask))
                process_flag = 1;

            outFile << hit.time << " " << hit.wavelength << " "
                    << hit.pos.x << " " << hit.pos.y << " " << hit.pos.z << " "
                    << hit.mom.x << " " << hit.mom.y << " " << hit.mom.z << " "
                    << hit.pol.x << " " << hit.pol.y << " " << hit.pol.z << " "
                    << process_flag << "\n";
        }

        outFile.close();
        G4cout << "GPUTaskManager: Saved " << num_hits << " hits to " << fname.str() << G4endl;
    }
};


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

std::string str_tolower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
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
            //G4cout << "PhotonSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0] << G4endl;
            fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        }
        hce->AddHitsCollection(fHCID, fPhotonHitsCollection);
    }

    G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *) override
    {
        G4Track *theTrack = aStep->GetTrack();
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
        //G4cout << "PhotonSD::EndOfEvent Number of PhotonHits: " << NbHits << G4endl;

        // Open an output file (text mode)

	/*    
        int tid = G4Threading::G4GetThreadId();
	std::ostringstream fname;
	fname << "g4_photon_hits_thread" << tid << ".txt";
	std::ofstream outFile(fname.str().c_str(), std::ios::app);

	
	if (!outFile.is_open())
        {
            G4cerr << "Error opening output file g4_photon_hits.txt!" << G4endl;
            return;
        }

        // Loop over all recorded hits TOBEDONE: move this to endofrunaction
        for (G4int i = 0; i < NbHits; i++)
        {
            PhotonHit *hit = (*fPhotonHitsCollection)[i];

            G4int id = hit->fid;
            G4double energy = hit->fenergy;
            G4double time = hit->ftime;
            G4ThreeVector position = hit->fposition;
            G4ThreeVector direction = hit->fdirection;
            G4ThreeVector pol = hit->fpolarization;

            // Write out info in a style similar to Opticks hits
            
	    outFile << "Adding hit from Geant4: " << energy << " eV  "
                    << "(" << position.x() << ", " << position.y() << ", " << position.z() << ")  "
                    << "(" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")  "
                    << "(" << pol.x() << ", " << pol.y() << ", " << pol.z() << ")  "
                    << "Time=" << time << " "
                    << "ID=" << id << G4endl;
	    
        }

        // Close the file
        outFile.close(); */
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
            std::cout << hit.wavelength << " " << position << " " << direction << " " << polarization << std::endl;

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
        G4ThreeVector position_mm(-0.4 * m, -0.3 * m, -0.3 * m);
	G4double time_ns = 0;
        G4ThreeVector direction(0, 0.2, 0.8);
        G4double wavelength_nm = 0.1;

        G4PrimaryVertex *vertex = new G4PrimaryVertex(position_mm, time_ns);
        G4PrimaryParticle *particle = new G4PrimaryParticle(G4Electron::Definition());
        particle->SetKineticEnergy(0.005 * GeV);
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
    }

    void EndOfEventAction(const G4Event *event) override
    {
    }
};

struct RunAction : G4UserRunAction
{
    GPUTaskManager *gpu_task_mgr_;  // Pointer to GPU task manager (owned by G4App)

    RunAction(GPUTaskManager *gpu_mgr = nullptr) : gpu_task_mgr_(gpu_mgr)
    {
    }

    void setGPUTaskManager(GPUTaskManager *mgr) {
        gpu_task_mgr_ = mgr;
    }

    void BeginOfRunAction(const G4Run *run) override
    {
        // Start the GPU task manager thread
        if (G4Threading::IsMasterThread() && gpu_task_mgr_) {
            gpu_task_mgr_->start();
        }
    }

    void EndOfRunAction(const G4Run *run) override
    {
        if (G4Threading::IsMasterThread())
        {
            if (gpu_task_mgr_) {
                // Flush any remaining photons and wait for GPU to finish
                gpu_task_mgr_->flushRemaining(0);
                gpu_task_mgr_->shutdown();

                G4cout << "=== GPU Task Manager Summary ===" << G4endl;
                G4cout << "Total batches processed: " << gpu_task_mgr_->getCompletedBatches() << G4endl;
                G4cout << "Total hits from GPU: " << gpu_task_mgr_->getTotalHitsProcessed() << G4endl;
                G4cout << "Photon threshold per batch: " << gpu_task_mgr_->getThreshold() << G4endl;
            } else {
                // Fallback to original behavior if no task manager
                G4CXOpticks *gx = G4CXOpticks::Get();

                auto start = std::chrono::high_resolution_clock::now();
                gx->simulate(0, false);
                cudaDeviceSynchronize();
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "Simulation time: " << elapsed.count() << " seconds" << std::endl;

                SEvt *sev = SEvt::Get_EGPU();
                unsigned int num_hits = sev->GetNumHit(0);
                std::cout << "Opticks: NumGensteps:  " << sev->GetNumGenstepFromGenstep(0) << std::endl;
                std::cout << "Opticks: NumPhotonsCollected:  " << sev->GetNumPhotonCollected(0) << std::endl;
                std::cout << "Opticks: NumHits:  " << num_hits << std::endl;
            }
        }
    }
};

struct SteppingAction : G4UserSteppingAction
{
    SEvt *sev;
    GPUTaskManager *gpu_task_mgr_;  // Pointer to GPU task manager (owned by G4App)

    SteppingAction(SEvt *sev, GPUTaskManager *gpu_mgr = nullptr)
        : sev(sev), gpu_task_mgr_(gpu_mgr)
    {
    }

    void setGPUTaskManager(GPUTaskManager *mgr) {
        gpu_task_mgr_ = mgr;
    }

    // Create a Cerenkov genstep quad6 (same as U4::MakeGenstep_G4Cerenkov_modified)
    static quad6 MakeGenstep_Cerenkov(
        const G4Track* aTrack,
        const G4Step* aStep,
        G4int numPhotons,
        G4double betaInverse,
        G4double pmin,
        G4double pmax,
        G4double maxCos,
        G4double maxSin2,
        G4double meanNumberOfPhotons1,
        G4double meanNumberOfPhotons2)
    {
        G4StepPoint* pPreStepPoint  = aStep->GetPreStepPoint();
        G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();

        G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4double      t0 = pPreStepPoint->GetGlobalTime();
        G4ThreeVector deltaPosition = aStep->GetDeltaPosition();

        const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
        const G4Material* aMaterial = aTrack->GetMaterial();

        G4double Wmin_nm = h_Planck * c_light / pmax / nm;
        G4double Wmax_nm = h_Planck * c_light / pmin / nm;

        quad6 gs;
        gs.zero();

        // Cast to scerenkov structure for easier field access
        scerenkov* ck = (scerenkov*)(&gs);

        ck->gentype = OpticksGenstep_G4Cerenkov_modified;
        ck->trackid = aTrack->GetTrackID();
        ck->matline = aMaterial->GetIndex() + SEvt::G4_INDEX_OFFSET;
        ck->numphoton = numPhotons;

        ck->pos.x = x0.x();
        ck->pos.y = x0.y();
        ck->pos.z = x0.z();
        ck->time = t0;

        ck->DeltaPosition.x = deltaPosition.x();
        ck->DeltaPosition.y = deltaPosition.y();
        ck->DeltaPosition.z = deltaPosition.z();
        ck->step_length = aStep->GetStepLength();

        ck->code = aParticle->GetDefinition()->GetPDGEncoding();
        ck->charge = aParticle->GetDefinition()->GetPDGCharge();
        ck->weight = aTrack->GetWeight();
        ck->preVelocity = pPreStepPoint->GetVelocity();

        ck->BetaInverse = betaInverse;
        ck->Wmin = Wmin_nm;
        ck->Wmax = Wmax_nm;
        ck->maxCos = maxCos;

        ck->maxSin2 = maxSin2;
        ck->MeanNumberOfPhotons1 = meanNumberOfPhotons1;
        ck->MeanNumberOfPhotons2 = meanNumberOfPhotons2;
        ck->postVelocity = pPostStepPoint->GetVelocity();

        return gs;
    }

    void UserSteppingAction(const G4Step *aStep)
    {
        G4Track *aTrack;
        G4int fNumPhotons = 0;

        G4StepPoint *preStep = aStep->GetPostStepPoint();
        G4VPhysicalVolume *volume = preStep->GetPhysicalVolume();

        if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
            // Kill if step count exceeds 10000 to avoid reflection forever
            if (aStep->GetTrack()->GetCurrentStepNumber() > 10000) {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            }
        }

        if (volume && std::strstr(volume->GetName().c_str(), "MirrorPyramid") != nullptr)
        {
            aTrack = aStep->GetTrack();
            if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
            {
                aTrack->SetTrackStatus(fStopAndKill);
            }
        }

        G4SteppingManager *fpSteppingManager =
            G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
        G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();

        if (stepStatus != fAtRestDoItProc)
        {
            G4ProcessVector *procPost = fpSteppingManager->GetfPostStepDoItVector();
            size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();

            for (size_t i3 = 0; i3 < MAXofPostStepLoops; i3++)
            {
                if ((*procPost)[i3]->GetProcessName() == "Cerenkov")
                {
                    aTrack = aStep->GetTrack();
                    const G4DynamicParticle *aParticle = aTrack->GetDynamicParticle();
                    G4double charge = aParticle->GetDefinition()->GetPDGCharge();
                    const G4Material *aMaterial = aTrack->GetMaterial();
                    G4MaterialPropertiesTable *MPT = aMaterial->GetMaterialPropertiesTable();

                    G4MaterialPropertyVector *Rindex = MPT->GetProperty(kRINDEX);
                    if (!Rindex || Rindex->GetVectorLength() == 0)
                    {
                        G4cout << "WARNING: Material has no valid RINDEX data. Skipping Cerenkov calculation." << G4endl;
                        return;
                    }

                    G4Cerenkov *proc = (G4Cerenkov *)(*procPost)[i3];
                    fNumPhotons = proc->GetNumPhotons();

                    if (fNumPhotons > 0)
                    {
                        G4double Pmin = Rindex->Energy(0);
                        G4double Pmax = Rindex->GetMaxEnergy();
                        G4double nMax = Rindex->GetMaxValue();
                        G4double beta1 = aStep->GetPreStepPoint()->GetBeta();
                        G4double beta2 = aStep->GetPostStepPoint()->GetBeta();
                        G4double beta = (beta1 + beta2) * 0.5;
                        G4double BetaInverse = 1. / beta;
                        G4double maxCos = BetaInverse / nMax;
                        G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
                        G4double MeanNumberOfPhotons1 =
                            proc->GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
                        G4double MeanNumberOfPhotons2 =
                            proc->GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);

                        const G4Event* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
                        if (!event) return;
                        G4int eventid = event->GetEventID();

                        // ASYNC MODE: Add genstep directly to GPUTaskManager's buffer
                        // This bypasses SEvt allowing CPU to continue while GPU processes
                        if (gpu_task_mgr_) {
                            quad6 gs = MakeGenstep_Cerenkov(aTrack, aStep, fNumPhotons,
                                BetaInverse, Pmin, Pmax, maxCos, maxSin2,
                                MeanNumberOfPhotons1, MeanNumberOfPhotons2);
                            gpu_task_mgr_->addGenstep(gs, fNumPhotons, eventid);
                        } else {
                            // SYNC MODE: Use original U4 path (adds to SEvt)
                            G4AutoLock lock(&genstep_mutex);
                            U4::CollectGenstep_G4Cerenkov_modified(aTrack, aStep, fNumPhotons,
                                BetaInverse, Pmin, Pmax, maxCos, maxSin2,
                                MeanNumberOfPhotons1, MeanNumberOfPhotons2);
                        }
                    }
                }
            }
        }
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
    G4App(std::filesystem::path gdml_file, bool enable_async_gpu = true)
        : sev(SEvt::CreateOrReuse_EGPU()),
          gpu_task_mgr_(enable_async_gpu ? new GPUTaskManager(&genstep_mutex) : nullptr),
          det_cons_(new DetectorConstruction(gdml_file)),
          prim_gen_(new PrimaryGenerator(sev)),
          event_act_(new EventAction(sev)),
          run_act_(new RunAction(gpu_task_mgr_)),
          stepping_(new SteppingAction(sev, gpu_task_mgr_)),
          tracking_(new TrackingAction(sev))
    {
        if (gpu_task_mgr_) {
            G4cout << "G4App: Async GPU simulation enabled (threshold="
                   << gpu_task_mgr_->getThreshold() << " photons)" << G4endl;
        } else {
            G4cout << "G4App: Using synchronous GPU simulation (end of run)" << G4endl;
        }
    }

    ~G4App() {
        // GPUTaskManager destructor handles shutdown
        delete gpu_task_mgr_;
        // G4CXOpticks::Finalize();
    }

    // Create "global" event
    SEvt *sev;

    // GPU task manager for async processing (owned by G4App)
    GPUTaskManager *gpu_task_mgr_;

    G4VUserDetectorConstruction *det_cons_;
    G4VUserPrimaryGeneratorAction *prim_gen_;
    EventAction *event_act_;
    RunAction *run_act_;
    SteppingAction *stepping_;
    TrackingAction *tracking_;
};
