#include <filesystem>

#include <plog/Log.h>

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;


void from_gdml(filesystem::path gdmlpath)
{
  G4GDMLParser parser;
  parser.Read(gdmlpath.string(), false);

  const G4VPhysicalVolume* world = parser.GetWorldVolume();

  if (!world) {
    LOG_ERROR << "Failed creatng G4 volume from GDML " << gdmlpath << endl;
  }
  else {
    LOG_INFO << "Created G4 volume " << world->GetName() << " from " << gdmlpath << endl;

    G4CXOpticks* g4cx = G4CXOpticks::SetGeometry(world);
    g4cx->saveGeometry("./out/csg");

    delete g4cx;
  }
}
