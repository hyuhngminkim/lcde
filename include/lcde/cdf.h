#ifndef LCDE_CDF_H_
#define LCDE_CDF_H_

#include "vector.h"
#include "utils.h"
#include "lcd.h"

namespace lcde {

class LCD;

class CDF {
  public:
  VectorT<mpf> knots;
  VectorT<mpf> cpk;
  mpf cap_;
  mpf base_;
  const uint64_t DIST = 100;

  // Creates a null CDF object
  CDF(mpf cap) {
    cap_ = cap;
    base_ = cap;
    cpk.resize(1);
    cpk(0) = cap;
  }

  CDF(const LCD& lcd) {
    knots = lcd.knots;
    cpk = lcd.cpk.row(0);
  }

  CDF(const LCD& lcd, mpf base, mpf ratio) {
    knots = lcd.knots;
    cpk = (lcd.cpk.row(0) * ratio).array() + base;
    cap_ = base + ratio;
    base_ = base;
  }

  void normalize(mpf base, mpf ratio) {
    cpk = (cpk * ratio).array() + base;
    cap_ = base + ratio;
    base_ = base;
  }

  size_t size() {
    size_t mpf_size = sizeof(mpf);
    return mpf_size * knots.cols() + mpf_size * cpk.cols() + mpf_size * 2;
  }
};

}

#endif      // LCDE_CDF_H_