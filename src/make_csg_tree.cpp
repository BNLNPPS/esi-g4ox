#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include <plog/Init.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ColorConsoleAppender.h>

#include "g4ox.h"

using namespace std;


class ArgParser
{
 public:
  ArgParser(int &argc, char **argv);
  string get_value(const string &option) const;
 private:
   vector<string> args;
};


int main(int argc, char **argv)
{
  ArgParser arg_parser(argc, argv);

  using PLogFormat = plog::TxtFormatter;
  static plog::ColorConsoleAppender<PLogFormat> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

  string gdmlpath = arg_parser.get_value("-f");

  from_gdml(gdmlpath);

  return EXIT_SUCCESS;
} 


ArgParser::ArgParser(int &argc, char **argv)
{
  for (int i=1; i < argc; ++i)
    args.push_back( string(argv[i]) );
}

string ArgParser::get_value(const string &option) const
{
  auto itr = find(args.begin(), args.end(), option);

  string value{""};

  if (itr != args.end() && ++itr != args.end()) {
    value = *itr;
  }

  // Default values
  if (value.empty() && option == "-f") return "geom.gdml";

  return value;
}
