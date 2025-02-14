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

    const char* comp = "genstep,photon,hit,domain,record,rec,seq" ;  // NB no "simtrace" here
    SEventConfig::SetGatherComp(comp); 
    SEventConfig::SetSaveComp(comp); 
    // TODO: this now automated ? check that 


    SEvt* evt = SEvt::Create(SEvt::EGPU)  ;  // holds gensteps and output NPFold of component arrays
    SEvt::AddCarrierGenstep();   // normally gensteps added after geometry setup, but can be before in this simple test 

    // TODO: this is missing setFrame

    const char* gdmlpath = "/esi/esi-g4ox/geom/basic_detector_diff_physics.gdml";  
    G4GDMLParser parser_;

    
    parser_.Read(gdmlpath, false);
    G4VPhysicalVolume *world = parser_.GetWorldVolume();

    SSim* sim = SSim::CreateOrReuse();
    stree* st = sim->get_tree();
    st->get_num_triangulated();
    const U4Tree*  tr;
    static U4SensorIdentifier* SensorIdentifier ;
    tr = U4Tree::Create(st, world, SensorIdentifier ) ;
    sim->initSceneFromTree(); // not so easy to do at lower level as do not want to change to SSim arg to U4Tree::Create for headeronly testing
    CSGFoundry* fd_ = CSGFoundry::CreateFromSim() ; // adopts SSim::INSTANCE


    CSGOptiX* cx = CSGOptiX::Create(fd_);   // uploads geometry, instanciates QSim 

    QSim* qs = cx->sim ; 







    int eventID = 0 ; 
    bool end = true ; 

    qs->simulate(eventID, end);  // internally calls CSGOptiX::simulate_launch following genstep uploading by QSim

    cudaDeviceSynchronize(); 

    evt->save();  // uses SGeo::LastUploadCFBase_OutDir to place outputs into CFBase/ExecutableName folder sibling to CSGFoundry  
 
    return 0 ; 
}
