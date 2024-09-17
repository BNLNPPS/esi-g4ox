#include <filesystem>
#include <iostream>
#include <string>

#include <argparse/argparse.hpp>

#include "FTFP_BERT.hh"
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalPhysics.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"

#include "G4CX/G4CXOpticks.hh"
#include "SysRap/NP.hh"
#include "SysRap/sphoton.h"

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
  void GeneratePrimaries(G4Event* event) {

    NP* photons = NP::Make<float>(0, 4, 4);

    cout << "read photons" << endl;
    photons->load("out/photons.npy");
    //photons->dump();

    size_t n_photons = photons->num_items();
    sphoton* sphotons = reinterpret_cast<sphoton*>(photons->bytes());

    for (int i=0; i < n_photons; i++)
    {
      sphoton &p = sphotons[i];
      //cout << "val: " << i << ": " << p;

      G4ThreeVector position_mm(p.pos.x, p.pos.y, p.pos.z);
      G4double time_ns = p.time;
      G4ThreeVector direction(p.mom.x, p.mom.y, p.mom.z);
      //direction = direction.unit();
      G4double wavelength_nm = p.wavelength ;
      G4ThreeVector polarization(p.pol.x, p.pol.y, p.pol.z);

      G4PrimaryVertex* vertex = new G4PrimaryVertex(position_mm, time_ns);
      G4double kineticEnergy = h_Planck*c_light/(wavelength_nm*nm) ;

      G4PrimaryParticle* particle = new G4PrimaryParticle(G4OpticalPhoton::Definition());
      particle->SetKineticEnergy(kineticEnergy);
      particle->SetMomentumDirection(direction);
      particle->SetPolarization(polarization);

      vertex->SetPrimary(particle);
      event->AddPrimaryVertex(vertex);
    }
  }
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

  string gdml_file, macro_name;
  bool interactive;

  program.add_argument("-g", "--gdml")
    .help("path to GDML file")
    .default_value(string("geom.gdml"))
    .nargs(1)
    .store_into(gdml_file);

  program.add_argument("-m", "--macro")
    .help("path to G4 macro")
    .default_value(string("run.mac"))
    .nargs(1)
    .store_into(macro_name);

  program.add_argument("-i", "--interactive")
    .help("whether to open an interactive window with a viewer")
    .flag()
    .store_into(interactive);

  try {
    program.parse_args(argc, argv);
  }
  catch (const exception& err) {
    cerr << err.what() << endl;
    cerr << program;
    exit(EXIT_FAILURE);
  }

  // Configure Geant4
  // The physics list must be instantiated before other user actions
  G4VModularPhysicsList *physics = new FTFP_BERT;
  physics->RegisterPhysics(new G4OpticalPhysics);

  G4RunManager run_mgr;
  run_mgr.SetUserInitialization(physics);

  G4App* g4app = new G4App(gdml_file);
  run_mgr.SetUserInitialization(g4app->det_cons_);
  run_mgr.SetUserAction(g4app->prim_gen_);
  run_mgr.Initialize();

  G4UIExecutive *uix = nullptr;
  G4VisManager  *vis = nullptr;

  if (interactive) {
    uix = new G4UIExecutive(argc, argv);
    vis = new G4VisExecutive;
    vis->Initialize();
  }

  G4UImanager *ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/execute " + macro_name);

  if (interactive) {
    uix->SessionStart();
  }

  delete uix;

  return EXIT_SUCCESS;
}
