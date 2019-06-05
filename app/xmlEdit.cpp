#include <pugixml.hpp>
#include <sstream>
#include <iostream>
#include <string>

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "No input file given\nExiting...\n";
    return -1;
  }

  std::stringstream line(argv[2]);
  std::stringstream objfolder, datafolder;
  int newResolution = 64;
  line >> newResolution;
  objfolder << "results/" << newResolution << "_p2_a2";
  datafolder << "results/data_" << newResolution << "_p2_a2";

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(argv[1]);

  pugi::xml_node node;
  node = xmlDoc.child("simulation");
  node.child("resolution").attribute("xyz").set_value(newResolution);
  node.remove_child("objFolder");
  node.append_child("objFolder")
      .append_child(pugi::node_pcdata)
      .set_value(objfolder.str().c_str());
  node.remove_child("dataFolder");
  node.append_child("dataFolder")
      .append_child(pugi::node_pcdata)
      .set_value(datafolder.str().c_str());

  xmlDoc.save_file(argv[1]);

  return 0;
}
