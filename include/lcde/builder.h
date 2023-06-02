#ifndef LCDE_BUILDER_H_
#define LCDE_BUILDER_H_

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include <sciplot/sciplot.hpp>

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
class Builder {
  private:
  mpf slope_;
  mpf intercept_;

  static constexpr mpf sampling_rate_ = .01;
  static constexpr long fanout_ = 10;
  static constexpr long min_size_ = 3;
  long input_size;
  std::vector<CDF> estimates_;

  struct point {
    mpf x;
    mpf y;
  };

  std::vector<point> sample_data_;

  public:
  Builder(const Builder&) = delete;
  Builder() = default;

  // CDF solve(std::vector<mpf>::iterator begin, std::vector<mpf>::iterator end) {
  CDF solve(const std::vector<point>& data) {
    int max_iter_ = 100;
    mpf tol_ = 1e-6;
    auto data_size = data.size();
    VectorT<mpf> x(data_size);
    VectorT<int> w(data_size);

    weight(data, x, w);
    double lambda = 1e-15;
    int n = data_size;
    int nx = x.cols();

    if (nx <= 2) {
      std::cerr << "Not enough elements. Aborting...\n";
      LCD l;
      return l;
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
    for (int i = 0; i < max_iter_; ++i) {
      if (lcd_.ll <= ll_old + tol_) {
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

    return CDF(lcd_);
  }

  void build(const std::vector<KeyType>& data) {
    estimates_.reserve(fanout_);

    // Sample Data
    auto begin = data.begin();
    auto end = data.end();

    input_size = data.size();
    const long sample_size = std::min<long>(
      input_size, std::max<long>(sampling_rate_ * input_size, min_size_)
    );

    std::cout << "sample size : " << sample_size << std::endl;

    sample_data_.reserve(sample_size);

    long offset = static_cast<long>(1. * input_size / sample_size);
    for (long i = 0; i < input_size; i += offset) {
      sample_data_.push_back({static_cast<mpf>(data[i]), 
                              static_cast<mpf>(1. * i / input_size)});
    }

    // Train first layer
    point min = sample_data_.front();
    point max = sample_data_.back();

    slope_ = 1. / (max.x - min.x);
    intercept_ = -slope_ * min.x;

    slope_ *= fanout_ - 1;
    intercept_ *= fanout_ - 1;

    // Allocate memory for second layer
    std::vector<std::vector<point>> training_data(fanout_);
    for (const auto& d : sample_data_) {
      long rank = static_cast<long>( slope_ * d.x + intercept_);

      rank = std::max(0L, std::min(fanout_ - 1, rank));

      training_data[rank].push_back(d);
    }

    for (long model_idx = 0; model_idx < fanout_; ++model_idx) {
      std::vector<point>& current_training_data = training_data[model_idx];

      if (model_idx == 0) {
        if (current_training_data.size() < 2) {
          estimates_.push_back(CDF(0.));

          point tp;
          tp.x = 0;
          tp.y = 0;
          current_training_data.push_back(tp);
        } else {
          min = current_training_data.front();
          max = current_training_data.back();

          CDF c = solve(current_training_data);
          c.normalize(0, max.y);

          estimates_.push_back(c);
        }
      } else if (model_idx == fanout_ - 1) {
        if (current_training_data.size() <= 2) {
          estimates_.push_back(CDF(1.));
        } else {
          min = training_data[model_idx - 1].back();
          max = current_training_data.back();

          CDF c = solve(current_training_data);
          c.normalize(estimates_[model_idx - 1].cap_, 1. - min.y);

          estimates_.push_back(c);
        }
      } else {
        if (current_training_data.size() <= 2) {
          estimates_.push_back(CDF(estimates_[model_idx - 1].cap_));

          point tp;
          tp.x = training_data[model_idx - 1].back().x;
          tp.y = training_data[model_idx - 1].back().y;
          current_training_data.push_back(tp);
        } else {
          min = training_data[model_idx - 1].back();
          max = current_training_data.back();

          CDF c = solve(current_training_data);
          c.normalize(estimates_[model_idx - 1].cap_, max.y - min.y);

          estimates_.push_back(c);
        }
      }
    }
  }

  bool find(KeyType key, int ground_truth) {
    long rank = slope_ * key + intercept_;
    rank = std::max(0L, std::min(static_cast<long>(fanout_ - 1), rank));

    // std::cout << "!!!" << key << " : " << rank << std::endl;

    CDF& c = estimates_[rank];
    long lhs, rhs;

    // std::cout << c.knots << std::endl;
    // std::cout << c.cpk << std::endl;
    // return true;

    if (c.cpk.cols() == 1) {
      lhs = c.base_ * 200000000;
      rhs = c.cap_ * 200000000;
      if (lhs <= ground_truth && ground_truth <= rhs) return true;
      else return false;
    } else {
      auto it = std::upper_bound(c.knots.data(), c.knots.data() + c.knots.cols(), key);
      int idx = std::distance(c.knots.data(), it);

      lhs = idx > 3 ? c.cpk[idx - 4] * 200000000 : c.base_ * 200000000;
      rhs = idx <= c.cpk.cols() - 5 ? c.cpk[idx + 4] * 200000000 : c.cap_ * 200000000;
      if (lhs <= ground_truth && ground_truth <= rhs) return true;
      else {
        std::cout << key << std::endl;
        std::cout << c.knots << std::endl;
        std::cout << c.cpk << std::endl;
        std::cout << estimates_[rank + 1].cpk << std::endl;
        std::cout << estimates_[rank + 1].knots << std::endl;
        return false;
      }
    }
    
    // if (c.cpk.cols() >= 2) {
    //   std::cout << c.knots << std::endl;
    //   std::cout << c.cpk << std::endl;
    //   std::cout << c.cap_ * 200000000 - c.DIST << std::endl;
    //   auto it = std::upper_bound(c.knots.data(), c.knots.data() + c.knots.cols(), key);
    //   int idx = std::distance(c.knots.data(), it);


      


    //   std::cout << "[" << ground_truth << "] " << c.cpk(idx - 1) * 200000000 << " , " << c.cpk(idx) * 200000000 << std::endl;
      
    //   int lhs = c.cpk(idx - 1) * 200000000;
    //   int rhs = c.cpk(idx) * 200000000;
    //   if (lhs <= ground_truth && ground_truth <= rhs) return true;
    //   else return false;
    // } else {

    //   std::cout << "[" << ground_truth << "] " << c.cap_ * 200000000 - c.DIST << " , " << c.cap_ * 200000000 + c.DIST << std::endl;
    //   // int lhs = c.cap_ * 200000000 - c.DIST;
    //   // int rhs = c.cap_ * 200000000 + c.DIST;
    //   int rhs = c.cap_ * 200000000;
    //   int lhs = rank > 0 ? estimates_[rank - 1].cap_ * 200000000 : 0;
    //   if (lhs <= ground_truth && ground_truth <= rhs) return true;
    //   else return false;
    // }
    // mpf cdf = get_cdf(estimates_[rank], key);
    // uint64_t res = static_cast<uint64_t>(input_size * (cdf + rank) / fanout_);
    // return 22;
  }

  // void printLCD() {
  //   int idx = 0; 
  //   for (auto i : estimates_) {
  //     // i.print_lcd();
  //     ++idx;
  //     std::cout << "[" << idx << "] " << i.lower << ", " << i.upper << std::endl;
  //   }
  // }

};

}         // namespace lcde

#endif    // LCDE_BUILDER_H_