/**
CSGOptiXSimTest
======================

Canonically used by cxsim.sh

**/

#include <cuda_runtime.h>

#include "CSGFoundry.h"
#include "CSGOptiX.h"
#include "G4GDMLParser.hh"
#include "OPTICKS_LOG.hh"
#include "QSim.hh"
#include "SEventConfig.hh"
#include "SEvt.hh"
#include "SSim.hh"
#include "U4GDML.h"
#include "U4Tree.h"

int main(int argc, char **argv)
{
    OPTICKS_LOG(argc, argv);
    // SEventConfig::SetRGModeSimulate();
    // SEventConfig::SetDebugLite();

    SEvt *evt = SEvt::Create(SEvt::EGPU);

    // Build the full path to the gensteps file

    const char *gdmlpath = "/esi/esi-g4ox/geom/pfrich_min_added_parameters.gdml";
    G4GDMLParser parser_;

    parser_.Read(gdmlpath, false);
    G4VPhysicalVolume *world = parser_.GetWorldVolume();
    assert(world);

    SSim *sim = SSim::CreateOrReuse();
    assert(sim && "sim instance should have been grabbed/created in ctor");

    stree *st = sim->get_tree();
    st->get_num_triangulated();
    const U4Tree *tr;
    static U4SensorIdentifier *SensorIdentifier;
    tr = U4Tree::Create(st, world, SensorIdentifier);
    sim->initSceneFromTree(); // not so easy to do at lower level as do not want to change to SSim arg to U4Tree::Create
                              // for headeronly testing
    CSGFoundry *fd_ = CSGFoundry::CreateFromSim(); // adopts SSim::INSTANCE
    float4 ce = make_float4(0.f, 0.f, 0.f, 100.f); // TODO: this should come from the geometry

    SEventConfig::SetMaxExtent(
        ce.w); // must do this config before QEvent::init which happens with CSGOptiX instanciation
    SEventConfig::SetMaxTime(10.f);

    CSGOptiX *cx = CSGOptiX::Create(fd_); // encumbent SSim used for QSim setup in here
    QSim *qs = cx->sim;
    // if(!SEvt::HasInputPhoton(SEvt::EGPU)) SEvt::AddTorchGenstep();

    int eventID = 0;
    bool end = false;
    auto start = std::chrono::high_resolution_clock::now();
    qs->simulate(eventID, end);
    cudaDeviceSynchronize();
    auto endtime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endtime - start;
    std::cout << "Simulation time: " << elapsed.count() << " seconds" << std::endl;

    SEvt *sev = SEvt::Get_EGPU();
    unsigned int num_hits2 = sev->GetNumHit(0);
    std::cout << "Opticks: NumHits:  " << num_hits2 << std::endl;

    return 0;
}
