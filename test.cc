#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "fit-sphere.hh"

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <filename> <radius>" << std::endl;
    return 1;
  }

  char *end;
  double radius = std::strtod(argv[2], &end);
  if (end == argv[2] || radius <= 0) {
    std::cerr << "Not a valid radius: " << argv[2] << std::endl;
    return 2;
  }

  std::vector<double> data;
  try {
    std::ifstream f(argv[1]);
    f.exceptions(std::ios::failbit | std::ios::badbit);
    while (!f.eof()) {
      for (size_t i = 0; i < 3; ++i) {
        double x;
        f >> x >> std::ws;
        data.push_back(x);
      }
    }
  } catch(std::ifstream::failure &) {
    std::cerr << "Error reading file \"" << argv[1] << "\"" << std::endl;
    return 3;
  }

  std::cout << data.size() / 3 << " points read." << std::endl;

  std::vector<double> result(3);
  auto found = fit_sphere(&data[0], data.size() / 3, radius, &result[0]);

  if (!found)
    std::cout << "Warning: iteration exceeded max. count before a good solution was found.\n";

  std::cout << result[0] << ", " << result[1] << ", " << result[2] << std::endl;

  return 0;
}
