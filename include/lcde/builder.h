#ifndef LCDE_BUILDER_H_
#define LCDE_BUILDER_H_

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <map>
#include <utility>
#include <iterator>
#include <cmath>
#include <Eigen/Eigenvalues>

#include "utils.h"
#include "lcd.h"
#include "loglik.h"
#include "gradient.h"
#include "mm.h"
#include "pnnls.h"
#include "line.h"
#include "cdf.h"

namespace lcde {

template <typename KeyType>
class Model;

template <typename KeyType>
class Builder {
  private:
  // The sampling ratio
  double sampling_rate;
  // Minimum number of elements for an estimation
  unsigned long min_size = 3;
  // The number of second layer lcds
  long fanout;
  // Slope and intercept for the linear regression of the first layer
  mpf slope;
  mpf intercept;

  mpf current_base = 0;
  mpf current_ratio = 0;
  size_t data_size_;

  // Temporary search boundaries
  int LRANGE, RRANGE;

  // Parameters ll(log-likelihood) and C(cumulative density) are commented as
  // "Remove" because these parameters can be removed when deemed unnecessary. 
  struct cdfPoint {
    KeyType knot;
    size_t cdf;
    mpf ll;                                                   // Remove
    mpf C;                                                    // Remove

    cdfPoint() {
      knot = 0;
      cdf = 0;
      ll = 1;                                                 // Remove
      C = 1;                                                  // Remove
    }

    cdfPoint(const Builder* b, mpf k, mpf c, mpf l, mpf normc)
    : knot(static_cast<KeyType>(k)),
      cdf(static_cast<size_t>(c * b->data_size_)),
      ll(l), C(normc) {}                                      // Remove
    
    cdfPoint(const cdfPoint& c) {
      knot = c.knot;
      cdf = c.cdf;
      ll = c.ll;                                              // Remove
      C = c.C;                                                // Remove
    }
  };

  // Auxiliary data structure for storing sampled data and its sample density
  struct point {
    mpf x;
    mpf y;
  };

  // For storing sampled data and their sample density
  std::vector<point> sample_data;

  // Map first layer to the vector on next layer
  std::vector<int> rootMap; 
  std::vector<int> sizePerInterval;

  // Second layer vector that stores all necessary data
  std::vector<cdfPoint> globalLayout;

  void LCDtoVector(const LCD& lcd) {
    auto cpk = (lcd.cpk.row(0) * current_ratio).array() + current_base;

    // If the resulting LCD has a single knot
    if (lcd.knots.cols() == 1) {
      globalLayout.push_back(cdfPoint(
        this,
        lcd.knots(0),
        current_ratio + current_base,
        lcd.ll,                                                  // Remove
        lcd.C                                                    // Remove
      ));
      return;
    }

    // If the resulting LCD has > 1 knots
    for (long i = 0; i < lcd.knots.cols(); ++i) {
      globalLayout.push_back(cdfPoint(
        this, 
        lcd.knots(i),
        cpk(i),
        lcd.ll,                                                  // Remove
        lcd.C                                                    // Remove
      ));
    }
  }

  void solve(const std::vector<point>& data) {
    int max_iter = 100;
    mpf tol = 1e-6;      // experimental changes
    auto data_size = data.size();
    VectorT<mpf> x(data_size);
    VectorT<int> w(data_size);

    weight(data, x, w);
    double lambda = 1e-15;
    int n = data_size;
    int nx = x.cols();

    if (nx <= 2) {
      std::cerr << "Not enough elements. Aborting...\n";
      return;
    }

    mpf lower = x(0);
    mpf upper = x(nx - 1);
    VectorT<mpf> theta;
    VectorT<mpf> pi;

    LCD lcd_ = LCD(0, lower, upper, theta, pi); 
    if (lower < lcd_.lower) lcd_.lower = lower;
    if (upper > lcd_.upper) lcd_.upper = upper;

    VectorT<mpf> xx(x.cols());
    xsq(xx, x, w);
    loglik(lcd_, x, w);
    mpf ll_old = std::numeric_limits<mpf>::lowest();

    // main loop
    for (int i = 0; i < max_iter; ++i) {
      if (lcd_.ll <= ll_old + tol) {
        break;
      }
      ll_old = lcd_.ll;

      VectorT<mpf> g_theta = maxima_gradient(lcd_, x, w, xx);

      if (g_theta.cols() >= 1) {
        VectorT<mpf> g_pi = VectorT<mpf>::Zero(g_theta.cols());
        int oi = 0;     // index of theta of original lcd
        int ni = 0;     // index of new theta
        while (oi < lcd_.theta.cols()) {
          if (lcd_.theta(oi) == g_theta(ni)) g_pi(ni) = lcd_.pi(oi++);
          ++ni;
        }
        lcd_ = LCD(lcd_.alpha, lcd_.lower, lcd_.upper, g_theta, g_pi);
      }

      VectorT<mpf>& knots = lcd_.knots;
      MatrixT<mpf>& cpk = lcd_.cpk;

      int nk = knots.cols();

      MatrixT<mpf> cpkr = getElementFromIndex(cpk, {nk - 1}).replicate(1, nk)
                          - cbind(0.0, dropElementByIndex(cpk, nk - 1));

      // E{(X - theta)_+}
      VectorT<mpf> mu = cpkr.row(1).array() - knots.array() * cpkr.row(0).array();
      // index of knots against x
      std::vector<int> knots_idx = getIndexOnly(knots, x);
      // gradient vector
      VectorT<mpf> grad = (mu * n) + getElementFromIndex(xx, knots_idx);
      // E{(X - theta_j)_+ (X - theta_k)_+}
      DynamicMatrixT<mpf> mm = computeDiff(knots, cpkr, nk);
      // negative Hessian matrix
      mm.triangularView<Eigen::Upper>() = mm.transpose();
      DynamicMatrixT<mpf> H = mm - tcrossprod(mu);

      // Catalogue: 
      // https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
      // https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
      Eigen::SelfAdjointEigenSolver<DynamicMatrixT<mpf>> es(H);

      std::vector<mpf> v2;
      std::vector<int> R_idx;
      v2.reserve(nk);
      R_idx.reserve(nk);
      // Note that the eigen() function of R sorts the eigenvalues in decreasing
      // order, but SelfAdjointEigenSolver of Eigen sorts the eigenvalues in
      // increasing order.
      mpf lower_bound = es.eigenvalues()(nk - 1) * lambda;
      for (auto i = 0; i < nk; ++i) {
        if (es.eigenvalues()(i) >= lower_bound) {
          R_idx.push_back(i);
          v2.push_back(sqrt(es.eigenvalues()(i)));
        }
      }

      int kr = v2.size();
      DynamicMatrixT<mpf> first_kr_transposed = 
                es.eigenvectors()(Eigen::placeholders::all, Eigen::seqN(0, kr))
                  .transpose();
      DynamicMatrixT<mpf> R = columnWiseMult(first_kr_transposed, v2);
      VectorT<mpf> p = grad / n + concatenate(-lcd_.alpha, lcd_.pi) * H;
      VectorT<mpf> b = 
      auxiliaryDivision((first_kr_transposed * p.transpose()).transpose(), v2);

      VectorT<double> nnls = pnnls(R.cast<double>(), b.cast<double>(), 1);
      nnls(0) = -nnls(0);

      // perform line search
      line_lcd(lcd_, x, w, xx, nnls, ll_old);
      // remove zero slope changes if exists
      if ((lcd_.pi.array() == 0).any()) lcd_.simplify();
    }

    LCDtoVector(lcd_);
    
    return;
  }

  friend class Model<KeyType>;

  public:
  Builder(const Builder&) = delete;
  Builder() = default;

  void build(const std::vector<KeyType>& data, int lrange, int rrange,
             std::vector<double> params) {
    // Set parameters
    LRANGE = lrange;
    RRANGE = rrange;
    data_size_ = data.size();

    sampling_rate = params[0];
    // fanout = 2
    fanout = static_cast<long>(params[1]);

    const long input_size = data.size();
    const long sample_size = std::min<long>(
      input_size, std::max<long>(sampling_rate * input_size, min_size)
    );

    sample_data.reserve(sample_size);

    long offset = static_cast<long>(1. * input_size / sample_size);
    for (long i = 0; i < input_size; i += offset) {
      sample_data.push_back({static_cast<mpf>(data[i]), 
                              static_cast<mpf>(1. * i / input_size)});
    }

    // Train first layer 
    point min = sample_data.front();
    point max = sample_data.back();

    slope = 1. / (max.x - min.x);
    intercept = -slope * min.x;

    slope *= fanout - 1;
    intercept *= fanout - 1;

    // Allocate memory for second layer
    std::vector<std::vector<point>> training_data(fanout);
    for (const auto& d : sample_data) {
      long rank = static_cast<long>( slope * d.x + intercept);
      rank = std::max(0L, std::min(fanout - 1, rank));
      training_data[rank].push_back(d);
    }

    // Train each subdata 
    rootMap.reserve(fanout);
    sizePerInterval.reserve(fanout);                                                  // Remove
    globalLayout.push_back(cdfPoint());

    for (long model_idx = 0; model_idx < fanout; ++model_idx) {
      std::vector<point>& current_training_data = training_data[model_idx];
      if (model_idx == 0) {                       // First model
        rootMap.push_back(0);
        if (current_training_data.size() < min_size) {
          point tp;
          tp.x = 0;
          tp.y = 0;
          current_training_data.push_back(tp);
        } else {
          min = current_training_data.front();
          max = current_training_data.back();

          current_ratio = max.y;
          solve(current_training_data);
          current_base += current_ratio;
          // Maybe remove
          // Copy ll and C of next point
          globalLayout[0].ll = globalLayout[1].ll;
          globalLayout[0].C = globalLayout[1].C;
        }
        sizePerInterval.push_back(current_training_data.size());                                                  // Remove
      } else if (model_idx == fanout - 1) {      // Last model
        if (current_training_data.size() < min_size) {
          rootMap.push_back(rootMap[rootMap.size() - 1]);
          sizePerInterval.push_back(sizePerInterval[sizePerInterval.size() - 1]);            // Remove
        } else {
          rootMap.push_back(globalLayout.size());
          min = training_data[model_idx - 1].back();

          current_ratio = 1 - min.y;
          solve(current_training_data);
          sizePerInterval.push_back(current_training_data.size());                    // Remove
        }
        globalLayout.push_back(cdfPoint(
          this, 
          std::numeric_limits<KeyType>::max(),
          1.,    // dummy max cdf
          1.,                                                  // Remove
          1.
        ));
      } else {                                    // Intermediate models
        if (current_training_data.size() < min_size) {
          rootMap.push_back(rootMap[rootMap.size() - 1]);

          point tp;
          tp.x = training_data[model_idx - 1].back().x;
          tp.y = training_data[model_idx - 1].back().y;
          current_training_data.push_back(tp);
          sizePerInterval.push_back(sizePerInterval[sizePerInterval.size() - 1]);            // Remove
        } else {
          rootMap.push_back(globalLayout.size());
          min = training_data[model_idx - 1].back();
          max = current_training_data.back();

          current_ratio = max.y - min.y;
          solve(current_training_data);
          current_base += current_ratio;
          sizePerInterval.push_back(current_training_data.size());                    // Remove
        }
        
      }
    }
  }

  std::pair<size_t, size_t> find(const KeyType& key) const {
    long rank = slope * key + intercept;
    rank = std::max(0L, std::min(static_cast<long>(fanout - 1), rank));
    long nrank = std::min(rank + 1, static_cast<long>(fanout - 1));
  
    auto it = std::upper_bound(globalLayout.begin() + rootMap[rank], 
                               globalLayout.begin() + rootMap[nrank], key,
    [](const KeyType& k, const cdfPoint& c) {
      return k < c.knot;
    });
    
    auto LITER = std::max(std::prev(it, LRANGE), globalLayout.begin());;
    auto RITER = std::min(std::next(it, RRANGE), globalLayout.end() - 1);

    return std::make_pair(LITER->cdf, RITER->cdf);
  }

};

}         // namespace lcde

#endif    // LCDE_BUILDER_H_