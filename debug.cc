#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>

#include "include/lcde/builder.h"
#include "test/header.h"

using namespace std;

template <class T>
void printTime(T time) {
  cout << endl;
  cout << "\033[1;36;47m                \033[m";
  cout << "\033[1;36;47m 😂 Elapsed time(ns): " << static_cast<double>(time) << "ns 😂 \033[m";
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

  int idx;
  int size;

  cout << "[1] books   [2] fb   [3] lognormal   [4] osm   [5] normal" << endl;
  cout << "[6] uniform(dense)   [7] uniform(sparse)   [8] wiki" << endl;

  while (1) {
    cout << "Enter an index for dataset selection : ";
    cin >> idx;
    if (1 <= idx && idx <= 8) break;
    cout << "Please enter a valid index..." << endl;
  }

  cout << endl << "[1] uint32   [2] uint64" << endl;

  // while (1) {
  //   cout << "Select the data size : ";
  //   cin >> size;
  //   if (1 <= size && size <= 2) break;
  //   cout << "Please enter a valid data size..." << endl;
  // }
  size = 2;

  string dataset = "../test/data/" + datasets[idx - 1] + "_200M_uint" + sizes[size - 1];
  vector<uint64_t> data = load_data<uint64_t>(dataset);

  cout << data[0] << " <-----> " << data[199999999] << endl;
  lcde::Builder<uint64_t> b;

  b.build(data);
  

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  auto start = chrono::high_resolution_clock::now();
  ///////////////////////////////////////////////////////////////////////
  ///////////////////// code for testing goes here. /////////////////////
  ///////////////////////////////////////////////////////////////////////

  b.find(data[0], 0);

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  auto end = chrono::high_resolution_clock::now();
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  chrono::nanoseconds elapsedNS = end - start;
  printTime(elapsedNS.count());

  // b.printLCD();

  // for (auto& e : b.)

  uint64_t total_error = 0;
  uint64_t max_error = 0;
  uint64_t min_error = 100000;

  uint64_t error_count = 0;

  for (size_t i = 1; i < data.size(); i += 1000) {
    if (b.find(data[i], i)) continue;
    else error_count += 1;
    // uint64_t foundidx = b.find(data[i]);
    // uint64_t uerror = foundidx > i ? foundidx - i : i - foundidx;
    // if (uerror > max_error) {
    //   max_error = uerror;
    //   // cout << "Max error updated : " << tmpidx << " //// " << foundidx << " to " << i << endl;
    // }
    // if (uerror < min_error) min_error = uerror;
    // total_error += uerror;
  }

  cout << endl << "Total error : " << error_count << "   Max error : " << max_error << "   Min error : " << min_error << endl;



  // b.build(rdata);

}