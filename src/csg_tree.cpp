#include <string>

#include <plog/Log.h>

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;


void from_gdml(string gdmlpath)
{
  G4GDMLParser parser;
  parser.Read(gdmlpath, false);

  const G4VPhysicalVolume* world = parser.GetWorldVolume();

  if (!world) {
    LOG_ERROR << "Failed creatng G4 volume from GDML " << gdmlpath << endl;
  }
  else {
    LOG_INFO << "Created G4 volume " << world->GetName() << " from " << gdmlpath << endl;

    G4CXOpticks* g4cx = G4CXOpticks::SetGeometry(world);
    g4cx->saveGeometry("./csg_tree");

    delete g4cx;
  }
}
