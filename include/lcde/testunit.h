#ifndef LCDE_TEST_UNIT_H_
#define LCDE_TEST_UNIT_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <map>
#include <utility>
#include <iterator>
#include <cmath>
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

using namespace sciplot;

// The purpose of class TestUnit is twofold:
// 1. Serve as a wrapper class for the builder such that we can plot the cdf of
//    the resulting LCDE object. 
// 2. Serve as a testing unit that we can test individual operations of LCDE.
//    Current LCDE implementation is built for the SOSD benchmark, and therefore
//    lacks the function of testing individual LCDE objects. 
// As such, class TestUnit has many similarities with its original class Builder.
template <typename KeyType>
class TestUnit {
  public:
  TestUnit() {
    testObject.setOuter(*this);
  }

  private:
  // The sampling ratio
  double sampling_rate;
  // Minimum number of elements for an estimation => Under test
  unsigned long min_size = 3;
  // The number of second layer lcds => Under test
  long fanout;
  // Slope and intercept for the linear regression of the first layer
  mpf slope;
  mpf intercept;
  // Auxiliary data structure for storing sampled data and its sample density
  struct point {
    mpf x;
    mpf y;

    point(mpf lhs, mpf rhs) : x(lhs), y(rhs) {}
  };
  // For storing sampled data and their sample density
  std::vector<point> sample_data;

  mpf current_base = 0;
  mpf current_ratio = 0;
  size_t data_size;

  // Struct cdfPoint is used to reduce the memory consumption.
  // However, memory consumption is not a matter in this testing environment.
  // We therefore define a new struct that stores all necessary elements for
  // plotting. 
  class {
    private:
    // Reference to the outer class TestUnit
    // Required for storing current_base and current_ratio
    TestUnit* tu = nullptr;

    public:
    // Constructor for testObject
    void setOuter(TestUnit& testUnitInstance) {
      tu = &testUnitInstance;
    }

    // LCD attributes required for plotting
    std::vector<double> C;
    std::vector<std::vector<double>> slopes;
    std::vector<std::vector<double>> intercepts;
    std::vector<std::vector<double>> cpk;
    std::vector<std::vector<double>> fk;
    std::vector<std::vector<double>> theta;

    // attributes required for calculating the proportional base and ratio of
    // the CDF of each data slice
    std::vector<double> bases;
    std::vector<double> ratios;

    // sciplot attributes required for plotting
    std::vector<std::vector<Vec>> lines;
    std::vector<std::vector<Vec>> knots;

    // Parameters for plotting with Python
    // The original equation for calculating the CDF is
    /*
      (base + ratio * (cpk - (fk / slope))) + (ratio / (C * slope)) * e^{slope * x + intercept}
    = (base + ratio * (cpk - (fk / slope))) + e^{slope * x + intercept + ln(ratio / (C * slope))}
    */
    // which requires 6 parameters. 
    // We reduce the number of parameters into 3, by precomputing each term.
    // The resulting equation is
    /*
      addend + e^{slope * x + newIntercept}
    */
    std::vector<std::vector<mpf>> mpf_slopes;
    std::vector<std::vector<mpf>> mpf_thetas;
    std::vector<std::vector<mpf>> addends;
    std::vector<std::vector<mpf>> newIntercepts;

    double lower;
    double upper = 0;

    void append(const LCD& lcd) {
      C.push_back(static_cast<double>(lcd.C));
      slopes.push_back(convertToSTL(lcd.slope));
      intercepts.push_back(convertToSTL(lcd.intercept));
      VectorT<mpf> cpk_first_row = lcd.cpk.row(0);
      std::vector<double> tmp_cpk_first_row = convertToSTL(cpk_first_row);
      tmp_cpk_first_row.insert(tmp_cpk_first_row.begin(), 0.);
      cpk.push_back(tmp_cpk_first_row);
      fk.push_back(convertToSTL(lcd.fk));

      std::vector<double> temp_theta = d_concatenate3(lcd.lower, lcd.theta, lcd.upper);
      theta.push_back(temp_theta);

      upper = lcd.upper > upper ? lcd.upper : upper;

      bases.push_back(static_cast<double>(tu->current_base));
      ratios.push_back(static_cast<double>(tu->current_ratio));
    }

    void mpfAppend(const LCD& lcd) {
      // Append C
      const mpf t_C = lcd.C;

      // Append slopes
      mpf_slopes.push_back(convert(lcd.slope));
      const auto& t_slope = lcd.slope;

      // Append intercepts
      const auto& t_intercept = lcd.intercept;

      // Append cpks
      VectorT<mpf> cpk_first_row = lcd.cpk.row(0);
      std::vector<mpf> t_cpk = convert(cpk_first_row);
      t_cpk.insert(t_cpk.begin(), 0.);

      // Append fk
      const auto& t_fk = lcd.fk;

      // Append theta
      std::vector<mpf> t_theta = concatenate3(lcd.lower, lcd.theta, lcd.upper);
      mpf_thetas.push_back(t_theta);

      std::vector<mpf> t_addend;
      std::vector<mpf> t_newIntercept;
      const mpf& t_base = tu->current_base;
      const mpf& t_ratio = tu->current_ratio;

      const size_t interval_size = t_slope.size();
      for (size_t i = 0; i < interval_size; ++i) {
        mpf addend;
        mpf newIntercept;
        if (t_slope[i] != 0) {
          addend = (t_base + t_ratio * (t_cpk[i] - (t_fk[i] / t_slope[i])));
          mpf t_ln;
          if (t_slope[i] < 0)  {
            t_ln = std::log(t_ratio / (t_C * -t_slope[i]));
          } else {
            t_ln = std::log(t_ratio / (t_C * t_slope[i]));
          }
          newIntercept = t_intercept[i] + t_ln;
        } else {
          // When slope = 0, we cannot divide the integral by slope. 
          // We calclate the cdf by addend + intercept * x for two reasons:
          // 1. we must keep the slope variable to check if it is zero
          // 2. we do not use the intercept variable
          // We therefore store the value multiplied to x in intercept not slope
          addend = t_base + t_ratio * (t_cpk[i] - (std::exp(t_intercept[i]) * t_theta[i]) / t_C);
          newIntercept = (t_ratio * std::exp(t_intercept[i])) / t_C;
        }
        t_addend.push_back(addend);
        t_newIntercept.push_back(newIntercept);
      }
      addends.push_back(t_addend);
      newIntercepts.push_back(t_newIntercept);
    }

    void printJSON(const std::string& dataset) {
      std::ofstream file;
      std::string filename = "parameters/" + dataset + "_plot_parameters.json";
      file.open(filename);

      file << std::fixed << std::setprecision(30);
      
      file << "{\"parameters\":[";

      for (size_t i = 0; i < mpf_thetas.size(); ++i) {
        if (i > 0) file << ",";
        file << "{";
        // print slopes
        file << "\"slope\":[";
        for (size_t j = 0; j < mpf_slopes[i].size(); ++j) {
          if (j > 0) file << ",";
          file << mpf_slopes[i][j];
        }
        file << "],";
        // print thetas
        file << "\"theta\":[";
        for (size_t j = 0; j < mpf_thetas[i].size(); ++j) {
          if (j > 0) file << ",";
          file << mpf_thetas[i][j];
        }
        file << "],";
        // print addends
        file << "\"addend\":[";
        for (size_t j = 0; j < addends[i].size(); ++j) {
          if (j > 0) file << ",";
          file << addends[i][j];
        }
        file << "],";
        // print intercepts
        file << "\"intercept\":[";
        for (size_t j = 0; j < newIntercepts[i].size(); ++j) {
          if (j > 0) file << ",";
          file << newIntercepts[i][j];
        }
        file << "]";
        file << "}";
      }

      file << "]}";

      file.close();

      return;
    }

    void mpfPlot() {
      Plot2D plot0;

      plot0.xlabel("Data");
      plot0.ylabel("Cumulative Density");

      plot0.xrange(theta[0][0], upper);

      for (size_t i = 0; i < C.size(); ++i) {
        for (size_t j = 0; j < lines[i].size(); ++j) {
          plot0.drawCurve(lines[i][j], static_cast<double>(addends[i][j]) + (mpf_slopes[i][j] >= 0 ? 1 : -1) * exp(static_cast<double>(mpf_slopes[i][j]) * lines[i][j] + static_cast<double>(newIntercepts[i][j])))
               .lineColor("red");
        }
      }

      plot0.legend().hide();

      Figure fig = {{plot0}};
      Canvas canvas = {{fig}};
      canvas.size(600, 450);
      canvas.show();
    }

    void finalize() {
      for (size_t i = 0; i < C.size(); ++i) {
        preparePlot(lines, knots, theta[i]);
      }
    }

    void plot(const bool doSave = false, const std::string saveDirectory = "") {
      // cdf
      Plot2D plot0;

      plot0.xlabel("Data");
      plot0.ylabel("Cumulative Density");

      plot0.xrange(theta[0][0], upper);

      for (size_t i = 0; i < C.size(); ++i) {
        for (size_t j = 0; j < lines[i].size(); ++j) {
          plot0.drawCurve(lines[i][j], bases[i] + ratios[i] * (cpk[i][j] + ((exp(slopes[i][j] * lines[i][j] + intercepts[i][j]) / C[i]) - fk[i][j]) / slopes[i][j]))
               .lineColor("red");
          // plot0.drawPoints(knots[i][j], bases[i] + ratios[i] * (cpk[i][j] + ((exp(slopes[i][j] * knots[i][j] + intercepts[i][j]) / C[i]) - fk[i][j]) / slopes[i][j]))
          //      .lineColor("red")
          //      .pointType(7);
        }
      }

      plot0.legend().hide();

      Figure fig = {{plot0}};
      Canvas canvas = {{fig}};
      canvas.size(600, 450);
      canvas.show();

      if (doSave) {
        canvas.save(saveDirectory);
      }
    }

    private:
    // Auxiliary function to append a mpf vector to a double vector
    inline void appendVectorT(VectorT<double>& v1, const VectorT<mpf>& v2) {
      VectorT<double> v2_d = v2.cast<double>();
      VectorT<double> temp(v1.size() + v2_d.size());
      temp.head(v1.size()) = v1;
      temp.tail(v2_d.size()) = v2_d;
      v1 = temp;
    }

    inline void appendVectorT(VectorT<double>& v1, const VectorT<double>& v2) {
      VectorT<double> temp(v1.size() + v2.size());
      temp.head(v1.size()) = v1;
      temp.tail(v2.size()) = v2;
      v1 = temp;
    }

    // Auxiliary function to append a single mpf element to a double vector
    inline void appendVectorT(VectorT<double>& v, const mpf& elem) {
      VectorT<double> temp(v.size() + 1);
      temp.head(v.size()) = v;
      temp[v.size()] = static_cast<double>(elem);
      v = temp;
    }

    // Auxiliary function to convert Eigen vector to std::vector
    // Assumes the input Eigen vector has type mpf
    inline std::vector<double> convertToSTL(const VectorT<mpf>& v) {
      VectorT<double> d_v = v.cast<double>();
      std::vector<double> res(d_v.data(), d_v.data() + d_v.size());
      return res;
    }

    // Auxiliary function to concatenate lower, theta, and upper
    inline std::vector<double> d_concatenate3(mpf lower, const VectorT<mpf> t, mpf upper) {
      VectorT<mpf> temp = concatenate(concatenate(lower, t), upper);
      return convertToSTL(temp);
    }

    inline std::vector<mpf> convert(const VectorT<mpf>& v) {
      std::vector<mpf> res(v.data(), v.data() + v.size());
      return res;
    }

    // Auxiliary function to concatenate lower, theta, and upper
    inline std::vector<mpf> concatenate3(mpf lower, const VectorT<mpf> t, mpf upper) {
      VectorT<mpf> temp = concatenate(concatenate(lower, t), upper);
      return convert(temp);
    }

    // Prepare the points for line plotting
    inline Vec prepareLine(double lhs, double rhs) {
      int interval = 100;    // plot 2000 points
      return linspace(lhs, rhs, interval);
    }

    // Prepare the points for knot plotting
    inline Vec prepareKnot(double lhs, double rhs) {
      return linspace(lhs, rhs, 1);
    }

    void preparePlot(std::vector<std::vector<Vec>>& l,
                     std::vector<std::vector<Vec>>& k,
                     const std::vector<double>& t) {
      assert(t.size() >= 2);
      std::vector<Vec> res_line;
      std::vector<Vec> res_knot;
      if (t.size() == 2) {  // theta has no elements other than lower / upper
        res_line.push_back(prepareLine(t[0], t[1]));
        res_knot.push_back(prepareKnot(t[0], t[1]));
      } else {
        for (size_t i = 1; i < t.size(); ++i) {
          res_line.push_back(prepareLine(t[i - 1], t[i]));
          res_knot.push_back(prepareKnot(t[i - 1], t[i]));
        }
      }
      l.push_back(res_line);
      k.push_back(res_knot);
    }

    void saveCanvas(Canvas canvas, std::string dir) {
      if (dir.empty()) {
        std::cerr << "Dataset name is not specified. File is not saved. " << std::endl;
        return;
      } else {
        canvas.save(dir);
        return;
      }
    }

  } testObject;

  public:
  // Build LCDE into a single object. 
  void buildSingle(const std::vector<point>& data) {
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
      // std::cout << "Loop\n";
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

    testObject.append(lcd_);
    testObject.mpfAppend(lcd_);
    
    return;
  }

  // Build LCDE using the standard linear model - fanout structure. 
  // Parameters are brought in as a vector for scalability and testing purposes.
  void build(const std::vector<KeyType>& data, std::vector<double> params) {
    // Set parameters
    sampling_rate = params[0];
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
      long rank = static_cast<long>(slope * d.x + intercept);
      rank = std::max(0L, std::min(fanout - 1, rank));
      training_data[rank].push_back(d);
    }

    // Train each subdata 
    for (long model_idx = 0; model_idx < fanout; ++model_idx) {
      std::vector<point>& current_training_data = training_data[model_idx];

      // The case for when the current_training_data.size() < 0 is a 
      // safety measure for when a slot of training_data has no element.
      // Bring the largest element of the previous training_data that has at least
      // a single element. If the first training_data slot has no element, insert
      // 0 as its element.
      if (model_idx == 0) {                      // First model
        if (current_training_data.size() < min_size) {
          current_training_data.push_back(point(0, 0));
          continue;
        } else {
          max = current_training_data.back();

          current_ratio = max.y;
          buildSingle(current_training_data);
          current_base += current_ratio;
        }
      } else if (model_idx == fanout - 1) {      // Last model
        if (current_training_data.size() < min_size) {
          continue;
        } else {
          min = training_data[model_idx - 1].back();

          current_ratio = 1 - min.y;
          buildSingle(current_training_data);
        }
      } else {                                    // Intermediate models
        if (current_training_data.size() < min_size) {
          current_training_data.push_back(training_data[model_idx - 1].back());
          continue;
        } else {
          min = training_data[model_idx - 1].back();
          max = current_training_data.back();

          current_ratio = max.y - min.y;
          buildSingle(current_training_data);
          current_base += current_ratio;
        }
      }
    }
  }

  void countThetas() const {
    std::cout << std::fixed;
    std::cout.precision(15);
    for (size_t i = 0; i < testObject.slopes.size(); ++i) {
      std::cout << std::endl;
      auto& vec = testObject.slopes[i];
      std::cout << "Slopes : [" << vec[0] << " , " << vec[vec.size() - 1] << "]" \
                << " : " << vec.size() << std::endl;
      auto& vec2 = testObject.intercepts[i];
      std::cout << "Intercepts [" << vec2[0] << " , " << vec2[vec2.size() - 1] << "]" \
                << " : " << vec2.size() << std::endl;
      auto& vec3 = testObject.theta[i];
      std::cout << "Thetas [" << vec3[0] << " , " << vec3[vec3.size() - 1] << "]" \
                << " : " << vec3.size() << std::endl;
      auto& vec4 = testObject.C[i];
      std::cout << "C [" << vec4 << "] " << std::endl;
      auto& vec5 = testObject.fk[i];
      std::cout << "fks [" << vec5[0] << " , " << vec5[vec5.size() - 1] << "]" \
                << " : " << vec5.size() << std::endl;
    }
  }

  // Find the appropriate index from the theta vector
  void findFromVector(const KeyType& key) const {
    long rank = static_cast<long>(slope * key + intercept);
    rank = std::max(0L, std::min(static_cast<long>(fanout - 1), rank));

    const auto& targetTheta = testObject.theta[rank];

    for (auto element : targetTheta) {
      std::cout << element << " ";
    }
    std::cout << std::endl;

    auto it = std::upper_bound(targetTheta.begin(), targetTheta.end(), key);

    std::cout << *it << std::endl;
  }

  // Find the index of a given key
  void find(const KeyType& key) const {
    long rank = static_cast<long>(slope * key + intercept);
    rank = std::max(0L, std::min(static_cast<long>(fanout - 1), rank));
    
    // auto iter = std::upper_bound(testObject.theta[rank].begin(),
    //                              testObject.theta[rank].end(), key);


    std::cout << rank << std::endl;
    std::cout << testObject.theta.size() << std::endl;
    // for (auto i : testObject.theta[rank]) {
    //   std::cout << i << " ";
    // }
    // std::cout << std::endl;
    // long idx = std::distance(testObject.theta[rank].begin(), iter - 1);

    // std::cout << testObject.theta[rank].size() << " , " << idx << std::endl;


    // auto s = testObject.slopes[rank][idx];
    // auto i = testObject.intercepts[rank][idx];
    // auto cpk = testObject.cpk[rank][idx];
    // auto fk = testObject.fk[rank][idx];
    // auto C = testObject.C[rank];
    // auto base = testObject.bases[rank];
    // auto ratio = testObject.ratios[rank];

    // std::cout << base + ratio * (cpk + ((std::exp(s * key + i) / C) - fk) / s);
  }

  // Compare the found index of a key with the ground truth index
  void compare(const::std::vector<KeyType>& data) {
    // for (size_t i = 0; i < data.size(); ++i) {

    // }
    size_t i = 100000000;
    std::cout << data[i] << std::endl;
    KeyType key = data[i];
    find(key);
  }
  /*
  https://quick-bench.com/q/738WZAS9aPv4Spoo8KVXNZ6t3S0
  https://codereview.stackexchange.com/questions/252427/c-performance-linear-regression-in-other-way
  http://spfrnd.de/posts/2018-03-10-fast-exponential.html
  */

  void plot(bool doSave = false, std::string dataset = "plot_result", 
            std::string extension = "pdf") {
    testObject.finalize();
    std::string saveFileName = "./plot/" + dataset + "." + extension;
    testObject.plot(doSave, saveFileName);
  }

  void printJSON(const std::string& dataset) {
    testObject.printJSON(dataset);
  }
};

}

#endif // LCDE_TEST_UNIT_H_