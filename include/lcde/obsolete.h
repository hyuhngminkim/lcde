#ifndef LCDE_OBSOLETE_H_
#define LCDE_OBSOLETE_H_

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

namespace lcde {

class OLCDE {
  private:
  LCD lcd_;
  std::string canvasPath_;
  std::vector<mpf> data_;
  size_t size_;

  sciplot::Vec get_range(mpf lhs, mpf rhs) {
    double dlhs = static_cast<double>(lhs), drhs = static_cast<double>(rhs);
    int interval = 2000;

    return sciplot::linspace(dlhs, drhs, interval);
  }

  sciplot::Vec get_interval(mpf lhs, mpf rhs) {
    double dlhs = cast<double>(lhs), drhs = cast<double>(rhs);

    return sciplot::linspace(dlhs, drhs, 1);
  }

  void plot_subroutine(std::vector<sciplot::Vec>& ranges, 
                       std::vector<sciplot::Vec>& intervals) {
    if (lcd_.theta.cols() == 0) {
      ranges.push_back(get_range(lcd_.lower, lcd_.upper));
      intervals.push_back(get_range(lcd_.lower, lcd_.upper));
      return;
    }
    ranges.push_back(get_range(lcd_.lower, lcd_.theta[0]));
    intervals.push_back(get_interval(lcd_.lower, lcd_.theta[0]));
    for (int i = 1; i < lcd_.theta.size(); ++i) {
      ranges.push_back(get_range(lcd_.theta[i - 1], lcd_.theta[i]));
      intervals.push_back(get_interval(lcd_.theta[i - 1], lcd_.theta[i]));
    }
    ranges.push_back(get_range(lcd_.theta[lcd_.theta.size() - 1], lcd_.upper));
    intervals.push_back(get_interval(lcd_.theta[lcd_.theta.size() - 1], lcd_.upper));

    return;
  }

  void save_canvas(sciplot::Canvas canvas) {
    if (canvasPath_.empty()) {
      std::cout << "Filename is not specified. Aborting..." << std::endl;
      return;
    } else {
      std::cout << "Saving file to " << canvasPath_ << std::endl;
      canvas.save(canvasPath_);
      return;
    }
  }

  public:
  OLCDE(const OLCDE&) = delete;
  OLCDE() = default;
  OLCDE(std::string s) : canvasPath_(s) {}
  OLCDE(std::vector<mpf> v) {
    ingest(v);
  }

  void ingest(std::vector<mpf>& x_in) {
    data_ = x_in;
    size_ = data_.size();
  }

  void build(int max_iter_=100, mpf tol_=1e-6) {
    VectorT<mpf> x(data_.size());
    VectorT<int> w(data_.size());

    weight(data_, x, w);
    double lambda = 1e-15;
    int n = data_.size();
    int nx = x.cols();

    if (nx <= 2) {
      std::cerr << "Not enough elements. Aborting...\n";
      return;
    }

    mpf lower = x(0);
    mpf upper = x(nx - 1);
    VectorT<mpf> theta;
    VectorT<mpf> pi;

    lcd_ = LCD(0, lower, upper, theta, pi); 
    
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

    return;
  }


  mpf find(mpf key) {
    auto it = std::upper_bound(lcd_.knots.data(), lcd_.knots.data() + lcd_.knots.cols(), key);
    int idx = std::distance(lcd_.knots.data(), it) - 1;
    return (idx > 0 ? lcd_.cpk(0,idx - 1) : 0) 
          + ((std::exp(lcd_.slope[idx] * key + lcd_.intercept[idx]) / lcd_.C) 
          - lcd_.fk[idx]) / lcd_.slope[idx];
  }

  void printLCD() const {
    lcd_.print_lcd();
    return;
  }

  void printMemory() const {
    lcd_.calculate_memory();
    return;
  }

  // If no format is specified, save the resulting image as png
  void setCanvasPath(std::string path) {
    canvasPath_ = path + ".png";
  }

  void setCanvasPath(std::string path, std::string format) {
    canvasPath_ = path + "." + format;
  }

  void plot(std::string prompt, bool save=false, int width=500, int height=300) {
    using namespace sciplot;

    std::vector<Vec> ranges;
    std::vector<Vec> intervals;

    plot_subroutine(ranges, intervals);

    VectorT<double> slope = lcd_.slope.cast<double>();
    VectorT<double> intercept = lcd_.intercept.cast<double>();
    VectorT<double> fk = lcd_.fk.cast<double>();
    double C = static_cast<double>(lcd_.C);

    // density
    Plot2D plot0;

    plot0.xlabel("Data");
    plot0.ylabel("Density");

    plot0.xrange(cast<double>(lcd_.lower), cast<double>(lcd_.upper));

    for (size_t i = 0; i < ranges.size(); ++i) {
      plot0.drawCurve(ranges[i], exp(slope[i] * ranges[i] + intercept[i]) /C)
           .lineColor("red");
      plot0.drawPoints(intervals[i], 
                      exp(slope[i] * intervals[i] + intercept[i]) / C)
           .lineColor("red")
           .pointType(7);
    }
    plot0.legend().hide();

    // logdensity
    Plot2D plot1;

    plot1.xlabel("Data");
    plot1.ylabel("Log-density");

    plot1.xrange(cast<double>(lcd_.lower), cast<double>(lcd_.upper));

    for (size_t i = 0; i < ranges.size(); ++i) {
      plot1.drawCurve(ranges[i], slope[i] * ranges[i] + intercept[i])
           .lineColor("red");
      plot1.drawPoints(intervals[i], slope[i] * intervals[i] + intercept[i])
           .lineColor("red")
           .pointType(7);
    }
    plot1.legend().hide();

    // cdf
    Plot2D plot2;

    plot2.xlabel("Data");
    plot2.ylabel("Cumulative Distribution");

    plot2.xrange(cast<double>(lcd_.lower), cast<double>(lcd_.upper));

    plot2.drawCurve(ranges[0], (exp(slope[0] * ranges[0] + intercept[0])) / (slope[0] * C))
         .lineColor("red");
    for (size_t i = 1; i < ranges.size(); ++i) {
      plot2.drawCurve(ranges[i], 
        (i > 0 ? lcd_.cpk(0, i - 1) : 0)
        + ((std::exp(slope[i] * ranges[i] + intercept[i]) / C) - fk[i]) / slope[i]
      ).lineColor("red");
      plot2.drawPoints(intervals[i], 
        (i > 0 ? lcd_.cpk(0, i - 1) : 0)
        + ((std::exp(slope[i] * intervals[i] + intercept[i]) / C) - fk[i]) / slope[i]
      ).lineColor("red")
       .pointType(7);
    }
    plot2.legend().hide();

    if (prompt == "pdf" || prompt == "density") {
      Figure fig = {{plot0}};
      Canvas canvas = {{fig}};
      canvas.size(width, height);
      canvas.show();
      if (save) save_canvas(canvas);
    } 
    else if (prompt == "logdensity") {
      Figure fig = {{plot1}};
      Canvas canvas = {{fig}};
      canvas.size(width, height);
      canvas.show();
      if (save) save_canvas(canvas);
    } 
    else if (prompt == "cdf") {
      Figure fig = {{plot2}};
      Canvas canvas = {{fig}};
      canvas.size(width, height);
      canvas.show();
    } 
    else if (prompt == "all") {
      Figure fig = {{plot1, plot0, plot2}};
      Canvas canvas = {{fig}};
      canvas.size(width * 3, height);
      canvas.show();
      if (save) save_canvas(canvas);
    } else {
      std::cout << "Not a valid prompt. Aborting..." << std::endl;
      return;
    }
  }

};        // class OLCDE

}         // namespace lcde

#endif    // LCDE_OBSOLETE_H_