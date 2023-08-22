#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "include/lcde/testunit.h"
#include "include/lcde/cxxopts.hpp"
#include "test/rdata.h"

using namespace lcde;

template <typename T> 
std::vector<T> load_data(const std::string& filename) {
  std::vector<T> data;
  std::ifstream in("../data/" + filename, std::ios::binary);
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
  TestUnit<mpf> tu;

  double sampling_rate = 0.01;
  double fanout = 20000;

  std::vector<double> params;
  params.push_back(sampling_rate);
  params.push_back(static_cast<double>(fanout));

  std::vector<std::string> datasets = {
    "books_200M_uint64",
    "fb_200M_uint64",
    "lognormal_200M_uint64",
    "normal_200M_uint64",
    "osm_cellids_200M_uint64",
    "uniform_dense_200M_uint64",
    "uniform_sparse_200M_uint64",
    "wiki_ts_200M_uint64"
  };

  long data_idx = 4;
  std::vector<uint64_t> uint_data = load_data<uint64_t>(datasets[data_idx]);
  std::vector<mpf> data(uint_data.begin(), uint_data.end());

  tu.build(data, params);

  std::ofstream file;
  std::string filename = datasets[data_idx] + "_log.json";
  file.open(filename);

  int max_err = 0;
  int max_i = 0;
  int max_res = 0;

  file << "{\"errors\":[";
  mpf err;

  for (int i = 0; i < data.size(); ++i) {
    if (i > 0) file << ",";
    int res = tu.find(data[i]);
    file << "{";
    file << "\"x\":" << data[i] << ",";
    
    err = res - i;
    if (err < 0) err = -err;

    file << "\"err\":" << err;
    file << "}";
  }

  file << "]}";
  file.close();
  
  return 0;
}