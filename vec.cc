#include <iostream>
#include "include/lcde/utils.h"

using namespace lcde;
using namespace std;

struct point {
  mpf x;
  mpf y;
};

int main() {

  vector<point> v = {{1,2},{1,2},{2,2},{3,3},{3,3},{3,3},{4,5},{5,6},{6,6},{6,6},{7,8},{9,9},{9,9}};
  VectorT<mpf> x(v.size());
  VectorT<int> w(v.size());

  weight(v, x, w);

  cout << x << endl;
  cout << w << endl;

  return 0;
}