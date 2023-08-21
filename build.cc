#include <iostream>
#include <string>

#include "include/lcde/testunit.h"
#include "include/lcde/cxxopts.hpp"
#include "test/rdata.h"

using namespace lcde;

template <typename T> 
std::vector<T> load_data(const std::string& filename) {
  std::vector<T> data;
  std::ifstream in("./data/" + filename, std::ios::binary);
  if (!in.is_open()) {
    std::cerr << "unable to open " << filename << std::endl;
    exit(EXIT_FAILURE);
  }
  // Read size.
  uint64_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
  data.resize(size);
  // Read values.
  in.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
  in.close();
  return data;
}

// This file tests the operations of a single LCDE object. 
int main(int argc, char** argv) {

  cxxopts::Options options("plot", "Plot cumulative density using LCDE");
  options.add_options()
    ("data", "Data file with keys", cxxopts::value<std::string>())
    ("r,sampling-rate", "Sampling rate for building LCDE", 
                                 cxxopts::value<double>())
    ("f,fanout", "Fanout of the root linear model", 
                                 cxxopts::value<int>())
    ;
  options.parse_positional({"data"});

  auto result = options.parse(argc, argv);

  const std::string dataset = result["data"].as<std::string>();
  const double sampling_rate = result["sampling-rate"].as<double>();
  const int fanout = result["fanout"].as<int>();
  
  TestUnit<mpf> tu;
  std::vector<double> params;
  params.push_back(sampling_rate);
  params.push_back(static_cast<double>(fanout));

  std::vector<uint64_t> uint_data = load_data<uint64_t>(dataset);
  std::vector<mpf> data(uint_data.begin(), uint_data.end());

  tu.build(data, params);
  tu.printJSON(dataset);
  
  return 0;
}