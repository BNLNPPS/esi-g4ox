#include <string>

#include <argparse/argparse.hpp>

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4VModularPhysicsList.hh"

#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "SysRap/OPTICKS_LOG.hh"

#include "g4app.h"

#include "G4VUserActionInitialization.hh"
#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"


using namespace std;


struct ActionInitialization : public G4VUserActionInitialization
{
private:
    G4App* fG4App; // Store the pointer to G4App

public:
    // Note the signature: now we take a pointer to the G4App itself
    ActionInitialization(G4App* app)
        : G4VUserActionInitialization(), fG4App(app)
    {
    }

    virtual void BuildForMaster() const override
    {
        // Possibly define master actions if needed
    }

    virtual void Build() const override
    {
        // Now you can safely refer to fG4App here
        // because it is a member of ActionInitialization.
        SetUserAction(fG4App->prim_gen_);
        SetUserAction(fG4App->run_act_);
        SetUserAction(fG4App->event_act_);
        SetUserAction(fG4App->tracking_);
        SetUserAction(fG4App->stepping_);
    }
};

int main(int argc, char **argv)
{
    OPTICKS_LOG(argc, argv);

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

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const exception &err)
    {
        cerr << err.what() << endl;
        cerr << program;
        exit(EXIT_FAILURE);
    }

    // Configure Geant4
    // The physics list must be instantiated before other user actions
    G4VModularPhysicsList *physics = new FTFP_BERT;
    physics->RegisterPhysics(new G4OpticalPhysics);

    auto *run_mgr = G4RunManagerFactory::CreateRunManager();
    run_mgr->SetUserInitialization(physics);


    G4App *g4app = new G4App(gdml_file);

    ActionInitialization* ActionInit = new ActionInitialization(g4app);
    run_mgr->SetUserInitialization(ActionInit);
    run_mgr->SetUserInitialization(g4app->det_cons_);
    run_mgr->Initialize();

    G4UIExecutive *uix = nullptr;
    G4VisManager *vis = nullptr;

    if (interactive)
    {
        uix = new G4UIExecutive(argc, argv);
        vis = new G4VisExecutive;
        vis->Initialize();
    }

    G4UImanager *ui = G4UImanager::GetUIpointer();
    ui->ApplyCommand("/control/execute " + macro_name);

    if (interactive)
    {
        uix->SessionStart();
    }

    delete uix;

    return EXIT_SUCCESS;
}
