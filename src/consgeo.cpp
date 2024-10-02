#include <algorithm>
#include <iostream>
#include <string>

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>

#include <argparse/argparse.hpp>

#include "g4ox.h"

using namespace std;

int main(int argc, char **argv)
{
    using PLogFormat = plog::TxtFormatter;
    static plog::ColorConsoleAppender<PLogFormat> consoleAppender;
    plog::init(plog::debug, &consoleAppender);

    argparse::ArgumentParser program("consgeo", "0.0.0");

    string gdml_file;

    program.add_argument("-g", "--gdml")
        .help("path to GDML file")
        .default_value(string("geom.gdml"))
        .nargs(1)
        .store_into(gdml_file);

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

    cout << "gdml_file: " << gdml_file << endl;

    from_gdml(gdml_file);

    return EXIT_SUCCESS;
}
