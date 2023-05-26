// code template for tesing the timing of functions

#include <chrono>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>

#include "../src/cnmlcd.h"

#include "header.h"
#include "rdata.h"

using namespace std;

template <class T>
void printTime(T time) {
  cout << endl;
  cout << "\033[1;36;47m                \033[m";
  cout << "\033[1;36;47m 😂 Elapsed time(ns): " << static_cast<double>(time) / 1e+9 << "s 😂 \033[m";
  cout << "\033[1;36;47m                \033[m\n";
  cout << endl;
  return;
}

int main() {

  vector<string> sizes = {
    "32",
    "64"
  };

  vector<string> datasets = {
    "books",
    "fb",
    "lognormal",
    "osm_cellids",
    "normal",
    "uniform_dense",
    "uniform_sparse",
    "wiki_ts"
  };

  vector<int> samplingRates = {
    10, 
    100, 
    1000,
    10000,
    100000
  };

  int idx;
  int size;
  int doSample;
  int samplingRate;

  cout << "[1] books   [2] fb   [3] lognormal   [4] osm   [5] normal" << endl;
  cout << "[6] uniform(dense)   [7] uniform(sparse)   [8] wiki" << endl;

  while (1) {
    cout << "Enter an index for dataset selection : ";
    cin >> idx;
    if (1 <= idx && idx <= 8) break;
    cout << "Please enter a valid index..." << endl;
  }

  cout << endl << "[1] uint32   [2] uint64" << endl;

  while (1) {
    cout << "Select the data size : ";
    cin >> size;
    if (1 <= size && size <= 2) break;
    cout << "Please enter a valid data size..." << endl;
  }

  string dataset = "data/" + datasets[idx - 1] + "_200M_uint" + sizes[size - 1];
  vector<uint64_t> original_data = load_data<uint64_t>(dataset);

  cout << "Allocating memory for data...";
  vector<cnmlcd::mpf> data(original_data.size());
  cout << "Done!\n";


  cout << "Inserting data...";
  for (size_t i = 0; i < original_data.size(); ++i) {
    data[i] = original_data[i];
  }
  cout << "Done!\n";

  while (1) {
    cout << endl << "[1] Sample data   [2] Original data" << endl;
    cout << "Select whether to sample the data or use the original data : ";

    cin >> doSample;
    if (doSample == 1) {
      cout << "Data will be sampled for faster computation." << endl;
      break;
    } else if (doSample == 2) {
      cout << "Data will be not be sampled." << endl;
      break;
    }
    cout << "Please enter a valid prompt." << endl;
  }

  if (doSample == 1) {
    while (1) {
      cout << endl;
      cout << "[1] 10   [2] 100   [3] 1000   [4] 10000   [5] 100000" << endl;
      cout << "Select the sampling rate : ";
      cin >> samplingRate;
      if (1 <= samplingRate && samplingRate <= 5) break;
      cout << "Please enter a valid sampling rate..." << endl;
    }
    cout << "Sampling data...";
    int ratio = samplingRates[samplingRate - 1];
    for (int i = 0, size = data.size(); i < size / ratio - 1; ++i) {
      data[i] = data[i * ratio];
    }
    data.resize(data.size() / ratio);
    cout << "Done!\n";

    dataset += ("_sr_" + to_string(ratio));
  }
  
  cnmlcd::Builder c;

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  auto start = chrono::high_resolution_clock::now();
  ///////////////////////////////////////////////////////////////////////
  ///////////////////// code for testing goes here. /////////////////////
  ///////////////////////////////////////////////////////////////////////


  c.solve(data);


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  auto end = chrono::high_resolution_clock::now();
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  chrono::nanoseconds elapsedNS = end - start;
  printTime(elapsedNS.count());

  c.printLCD();

  c.printMemory();

  c.setCanvasPath(dataset, "pdf");
  c.plot("all", true);

  return 0;
}