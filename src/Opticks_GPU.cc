/**
CSGOptiXSimTest
======================

Canonically used by cxsim.sh 

**/

#include <cuda_runtime.h>

#include "OPTICKS_LOG.hh"
#include "SEventConfig.hh"
#include "SEvt.hh"
#include "CSGFoundry.h"
#include "CSGOptiX.h"
#include "QSim.hh"
#include "G4GDMLParser.hh"
#include "U4GDML.h"
#include "U4Tree.h"
#include "SSim.hh"


int main(int argc, char** argv)
{
    OPTICKS_LOG(argc, argv); 
    //SEventConfig::SetRGModeSimulate();
    //SEventConfig::SetDebugLite();

    SEvt* evt = SEvt::Create(SEvt::EGPU) ;

    // Build the full path to the gensteps file

    const char* gdmlpath = "/esi/esi-g4ox/geom/basic_detector_diff_physics.gdml";  
    G4GDMLParser parser_;

    parser_.Read(gdmlpath, false);
    G4VPhysicalVolume *world = parser_.GetWorldVolume();
    assert(world);

    SSim* sim = SSim::CreateOrReuse();
    assert(sim && "sim instance should have been grabbed/created in ctor" ); 

    stree* st = sim->get_tree();
    st->get_num_triangulated();
    const U4Tree*  tr;
    static U4SensorIdentifier* SensorIdentifier ;
    tr = U4Tree::Create(st, world, SensorIdentifier ) ;
    sim->initSceneFromTree(); // not so easy to do at lower level as do not want to change to SSim arg to U4Tree::Create for headeronly testing
    CSGFoundry* fd_ = CSGFoundry::CreateFromSim() ; // adopts SSim::INSTANCE
    float4 ce = make_float4( 0.f, 0.f, 0.f, 100.f );  // TODO: this should come from the geometry 

    SEventConfig::SetMaxExtent( ce.w );  // must do this config before QEvent::init which happens with CSGOptiX instanciation
    SEventConfig::SetMaxTime( 10.f ); 

    CSGOptiX* cx = CSGOptiX::Create(fd_); // encumbent SSim used for QSim setup in here 
    QSim* qs = cx->sim ; 
    //if(!SEvt::HasInputPhoton(SEvt::EGPU)) SEvt::AddTorchGenstep();      

    int eventID = 0 ; 
    bool end = false ; 
    qs->simulate(eventID, end);  
    cudaDeviceSynchronize(); 
    evt->save(); 

    return 0 ; 
}
