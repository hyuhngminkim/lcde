#include <iostream>
#include <vector>
#include <string>
#include <sciplot/sciplot.hpp>

#include "header.h"

using namespace std;
using namespace sciplot;

int main() {
  vector<string> files = {
    "data/books_200M_uint64",
    "data/fb_200M_uint64",
    "data/lognormal_200M_uint64",
    "data/normal_200M_uint64",
    "data/osm_cellids_200M_uint64",
    "data/uniform_dense_200M_uint64",
    "data/uniform_sparse_200M_uint64",
    "data/wiki_ts_200M_uint64"
  };

  vector<uint64_t> data = load_data<uint64_t>(files[0]);
  vector<int> x,y;

  // cout << data[1000000] << endl;
  // double a = static_cast<double>(data[1000000]);
  // cout << a << endl;

  Plot2D plot;
  plot.xlabel("Data");
  plot.ylabel("Cumulative Distribution");

  

  for (int i = 0, size = data.size(); i < size; i += 2) {
    x.push_back(data[i]);
    y.push_back(i);
    // double a = static_cast<double>(data[i]);
    // plot.drawPoints(a, i);
  }

  plot.drawPoints(data, y);

  Figure fig = {{plot}};
  Canvas canvas = {{fig}};
  canvas.size(600, 400);
  canvas.show();

  // for (auto f : files) {
  //   Plot2D plot;
  //   plot.xlabel("Data");
  //   plot.ylabel("Cumulative Distribution");
  //   plot.legend().hide();
  //   vector<uint64_t> data = load_data<uint64_t>(f);
  //   // plot.xrange(data[0], data[199999999]);
  //   // for (int i = 0; i < data.size(); ++i) {
  //   //   double x = data[i] / 10;
  //   //   plot.drawPoints(x, i).lineColor("red").pointType(7);
  //   // }
  //   plot.drawPoints(data, data);
  //   Figure fig = {{plot}};
  //   Canvas canvas = {{fig}};
  //   canvas.size(600, 400);
  //   canvas.show();
  // }

  return 0;
}