#include <algorithm>
#include <string>
#include <iostream>

#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4GenericPhysicsList.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;


struct G4App : public G4VUserDetectorConstruction
{
  G4VPhysicalVolume* Construct();
  virtual ~G4App() {}
};


int main(int argc, char **argv)
{
  G4App* g4app = new G4App;

  G4RunManager runMgr;
  runMgr.SetUserInitialization(new G4GenericPhysicsList);
  runMgr.SetUserInitialization(static_cast<G4VUserDetectorConstruction*>(g4app));
  runMgr.Initialize();
  //runMgr.BeamOn(2);

  return EXIT_SUCCESS;
}


G4VPhysicalVolume* G4App::Construct()
{
  G4GDMLParser parser;
  string gdmlpath("geom.gdml");
  parser.Read(gdmlpath, false);

  return parser.GetWorldVolume();
}
