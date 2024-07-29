#include <string>

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;


void from_gdml(string gdmlpath)
{
  G4GDMLParser parser;
  parser.Read(gdmlpath, false);
  const G4VPhysicalVolume& world = *parser.GetWorldVolume(); 

  cout << "from_gdml: " << gdmlpath << ": " << world.GetName() << endl;

  G4CXOpticks* g4cx = G4CXOpticks::SetGeometry(&world);
  g4cx->saveGeometry("./csgeom");

  delete g4cx;
}
