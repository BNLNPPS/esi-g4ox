#include <filesystem>

#include <plog/Log.h>

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include "G4CX/G4CXOpticks.hh"

using namespace std;

void from_gdml(filesystem::path gdml_file, filesystem::path out_prefix)
{
    G4GDMLParser parser;
    parser.Read(gdml_file.string(), false);

    const G4VPhysicalVolume *world = parser.GetWorldVolume();

    if (!world)
    {
        LOG_ERROR << "Failed creatng G4 volume from GDML " << gdml_file << endl;
    }
    else
    {
        G4CXOpticks *g4cx = G4CXOpticks::SetGeometry(world);

        auto outpath = out_prefix / gdml_file.stem();

        g4cx->saveGeometry(outpath.string().c_str());

        LOG_INFO << "Created G4 volume " << world->GetName() << " from " << gdml_file << endl;
        LOG_INFO << "Saved CSG tree to " << outpath << endl;

        delete g4cx;
    }
}
