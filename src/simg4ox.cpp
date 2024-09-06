#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

#include <argparse/argparse.hpp>

#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4GenericPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;


struct DetectorConstruction : G4VUserDetectorConstruction
{
  DetectorConstruction(filesystem::path gdml_file) : gdml_file_(gdml_file) {}

  G4VPhysicalVolume* Construct()
  {
    G4GDMLParser parser;
    parser.Read(gdml_file_.string(), false);

    return parser.GetWorldVolume();
  }

 private:

  filesystem::path gdml_file_;
};


struct PrimaryGenerator : G4VUserPrimaryGeneratorAction
{
  void GeneratePrimaries(G4Event* evt) {}
};


struct G4App
{
  G4App(filesystem::path gdml_file) :
    det_cons_(new DetectorConstruction(gdml_file)),
    prim_gen_(new PrimaryGenerator)
  {
  }

  G4VUserDetectorConstruction*   det_cons_;
  G4VUserPrimaryGeneratorAction* prim_gen_;
};


int main(int argc, char **argv)
{
  argparse::ArgumentParser program("simg4ox", "0.0.0");

  string gdml_file;

  program.add_argument("-g", "--gdml")
    .help("path to GDML file")
    .default_value(string("geom.gdml"))
    .nargs(1)
    .store_into(gdml_file);

  try {
    program.parse_args(argc, argv);
  }
  catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program;
    exit(EXIT_FAILURE);
  }

  // Configure Geant4
  G4RunManager run_mgr;
  run_mgr.SetUserInitialization(new G4GenericPhysicsList);

  G4App* g4app = new G4App(gdml_file);

  run_mgr.SetUserInitialization(g4app->det_cons_);
  run_mgr.SetUserAction(g4app->prim_gen_);
  run_mgr.Initialize();
  run_mgr.BeamOn(2);

  return EXIT_SUCCESS;
}
