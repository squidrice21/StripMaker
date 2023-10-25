#include <string>

#include "Util.h"
#include "feature/FeatureIO.h"

#include <CLI/CLI.hpp>

int main(int argc, char *argv[]) {
  CLI::App app{"Input preprocessing."};

  std::string scap_filename;
  std::string output_scap_filename;

  app.add_option("-i,--input", scap_filename, "The ground truth scap input.");
  app.add_option("-o,--output", output_scap_filename, "The output csv file.");

  CLI11_PARSE(app, argc, argv);

  // 1. Preprocess
  int width, height;
  bool to_preprocess = true;
  Capture capture;
  read_input(scap_filename, capture, width, height, to_preprocess);

  // 2. Write preprocessed results
  std::ofstream scap_ofs(output_scap_filename);
  std::string out_buffer = capture.to_string();
  scap_ofs << "#" << width << "\t" << height << std::endl;
  scap_ofs.write(out_buffer.c_str(), out_buffer.size());
  scap_ofs.close();

  return 0;
}
