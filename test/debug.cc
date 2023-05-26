#include <numeric>
#include <algorithm>

// #include "../lcd.h"
#include "header.h"
#include "../utils.h"
// #include "../lcd.h"
// #include "../gradient.h"
// #include "../mm.h"
// #include "../pnnls.h"

#include "../../gn/gnrandom.h"

using namespace std;
// using mpf = boost::multiprecision::mpf_float_50;


// extern "C" {
//   void pnnls_(double** fa_, int* fmda_, int* fm_, int* fn_, double* fb_, 
//               double* fx_, double* frnorm_, double* fw_, double* fzz_,
//               int* findex_, int* fmode_, int* fk_);
// }

// /*!
//  * @function    pnnls
//  * @abstract    Computes the non-negative least squares
//  * @discussion  A C wrapper function for the pnnls fortran subroutine by 
//  *              Yong Wang, Charles L. Lawson, and Richard J. Hanson. 
//  * @param   a   The M by N matrix
//  * @param   b   The M vector
//  * @param   k   The first k variables are not NN-restricted
//  */
// template <class T>
// void pnnls(DynamicMatrixT<T> a, VectorT<T> b, int k) {
//   int m = a.rows();
//   int n = a.cols();
//   assert(n == b.cols());

//   // https://stackoverflow.com/questions/54282732/c-convert-array-of-elements-into-2-d-matrix
//   double (*fa)[m] = (double(*)[m])(a.data());
//   double* fb = b.data();

//   cout << fixed;
//   cout.precision(12);
//   for (int i = 0; i < 2; ++i) {
//     for (int j = 0; j < 2; ++j) {
//       cout << fa[i][j] << " ";
//     }
//     cout << endl;
//   }

//   double x[n];            // only for output
//   double rnorm[1];        // only for output

//   double w[n];            // n-vector of working space
//   double zz[m];           // m-vector of working space
//   int index[n];           // n-vector index, only for output
//   int mode[1];            // success/failure flag, 1 = success

//   pnnls_((double**)fa, &m, &m, &n, fb, x, rnorm, w, zz, index, mode, &k);

//   for (int i = 0; i < 2; ++i) {
//     cout << x[i] << " ";
//   }
//   cout << endl;
  
// }


int main() {

  // VectorT<double> target(10);
  // target << 199, 299, 399, 499, 599, 699, 799, 899, 999, 1099;
  
  // VectorT<double> x(5);
  // x << 1, 7, 14, 21, 29;

  // VectorT<double> v(10);
  // v << -1, 2, 5, 6, 10, 15, 20, 22, 24, 30;

  // 0316: What happens when you append a scalar to a null vector?
  // VectorT<double> test;
  // double val = 9.9;
  // VectorT<double> res = cumsum(test, true);
  // cout << res << endl;

  // 0316: Debugging dropByIndex()
  // VectorT<double> test(1);
  // test << 1;
  // VectorT<double> res = dropByIndex(test, 0);
  // cout << res << endl;

  // 0316: Debugging getElementFromIndex
  // MatrixT<int> m(3, 3);
  // m << 1, 2, 3,
  //      4, 5, 6,
  //      7, 8, 9;

  // vector<int> question = {0, 1, 1, 0, 1};
  // Eigen::Block<MatrixT<int>> b = m(Eigen::seq(0, 1), Eigen::seq(1, 2));
  // // MatrixT<int> temp = m(Eigen::seq(0, 1), Eigen::seq(1, 2));

  // DynamicMatrixT<int> res = getElementFromIndex(b, question);

  // cout << res << endl;

  // 0317: Debugging cbind for Eigen::Block
  // VectorT<double> x(10);
  // x << 100, 105, 110, 130, 141, 152, 158, 189, 194, 200;
  // VectorT<double> theta;
  // VectorT<double> pi;
  // LCD<double> lcd(0, 100, 200, theta, pi);

  // lcd.print_lcd();

  // auto cpk = lcd.cpk;
  // auto top = cpk.topRows(1);
  // cout << cbind(0.0, lcd.cpk.topRows(3)) << endl;
  // cout << lcd.cpk(Eigen::seq(0, 2), Eigen::all) << std::endl << std::endl;
  // cout << cbind(0.0, lcd.cpk(Eigen::seq(0, 2), Eigen::all)) << endl;

  // auto cpx_val = cpx(lcd, x, 2);


  // 0320: Debugging getIndexOnlyFromSorted
  // vector<int> idx = getIndexOnlyFromSorted(x, v);
  // for (auto i : idx) cout << i << endl;


  // 0320: Debugging getMaxIndexFromRange
  // VectorT<double> rnds = VectorT<double>::Random(10);
  // cout<< getMaxIndexFromRange(rnds, 4, 9) << endl;
  // int t = getMaxIndexFromRange(v, 2, 4);
  // cout<< t << endl;

  // 0320: Debugging getIndexOnlyFromIntervals
  // std::vector<int> test_range = {0, 4, 6, 9};
  // print_vector(rnds);
  // print_vector(getIndexOnlyFromIntervals(rnds, test_range));

  // 0321: Debugging getIndexOnlyFromIntervals - interleaved version
  // std::vector<int> test_range = {0, 4, 6, 9};
  // print_vector(rnds);
  // print_vector(getIndexOnlyFromIntervals(rnds, test_range, true));

  // 0321: Debugging class copying
  // SomeClass<double> s(4);
  // s = SomeClass<double>(5);
  // print_vector(s.knots);

  // 0321: Debugging dropElementByIndex for matrix
  // MatrixT<int> m(3, 3);
  // m << 1, 2, 3,
  //      4, 5, 6,
  //      7, 8, 9;

  // cout << dropElementByIndex(m, 2);


  // 0322: Debugging tcrossprod
  // VectorT<double> x(5);
  // x << 1, 2, 3, 4, 5;
  // DynamicMatrixT<double> d = tcrossprod(x);

  // cout << d << endl;

  // 0322: Debugging repeatelementN
  // VectorT<double> v(10);
  // v << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;

  // cout << repeatN(v, 5) << endl;
  // cout << repeatElementwiseN(v, 5) << endl;

  // 0322: Debugging auxiliaryMatrixSubtraction
  // int nk = 3;
  // DynamicMatrixT<double> d = DynamicMatrixT<double>::Random(3,3);
  // VectorT<double> v = VectorT<double>::Random(9);

  // cout << d << endl << endl << v << endl << endl;
  // cout << auxiliaryMatrixSubtraction(d, v) << endl;
  // cout << auxiliaryMatrixAddition(d, v) << endl;


  // 0323: Debugging copyDiagonal
  // DynamicMatrixT<double> d = DynamicMatrixT<double>::Random(5, 5);
  // cout << d << endl << endl;
  // copyDiagonal(d);
  // cout << d << endl;


  // 0324: Debugging SelfAdjoingEigenSolver
  // DynamicMatrixT<double> d = DynamicMatrixT<double>::Random(5, 5);
  // copyDiagonal(d);

  // cout << d << endl << endl;

  // Eigen::SelfAdjointEigenSolver<DynamicMatrixT<double>> saes(d);

  // cout << saes.eigenvalues().transpose() << endl << endl;
  // cout << saes.eigenvectors() << endl << endl;

  // Eigen::EigenSolver<DynamicMatrixT<double>> es(d);

  // cout << es.eigenvalues().transpose() << endl << endl;
  // cout << es.eigenvectors() << endl << endl;


  // 0324: Testing SelfAdjoingEigenSolver iteration
  // DynamicMatrixT<double> d = DynamicMatrixT<double>::Random(5, 5);
  // copyDiagonal(d);
  // cout << d << endl << endl;
  // Eigen::SelfAdjointEigenSolver<DynamicMatrixT<double>> saes(d);
  // cout << saes.eigenvalues()(4) << " ";
  

  // 0325: Testing push_back on reserved vector
  // vector<int> r;
  // r.reserve(15);
  // for (int i = 0; i < 15; ++i) r.push_back(i);
  // print_vector(r);
  // cout << r.capacity() << endl;


  // 0406: Testing auxiliarydivision
  // VectorT<double> lhs(5);
  // lhs << 5, 12, 5, 6, 98;
  // vector<double> rhs = {4, 8, 12, 42, 21};

  // print_vector(auxiliaryDivision(lhs, rhs));

  // 0509: Testing vector element multiplication
  // VectorT<double> a(2);
  // a << 99.9, 10.8;
  // cout << a << endl;
  // a(0) = -a(0);
  // cout << a << endl;

  // 0510: Testing pnnls
  // DynamicMatrixT<double> a(2, 2);
  // VectorT<double> b(2);

  // a <<  -285.5224651, -6.875762,
  //        0.3336661, -13.855795;
  // b << 0.1881488, -0.1976159;
  // stored as [ 1, 4, 7, 2, 5, 8 ]
  // m = rows = 3
  // n = cols = 2
  // for the array to be converted to a matrix, it must be sliced by the row
  // number (m)
  // which probably changes the orientation to the following
  // 1, 4, 7
  // 2, 5, 8
  // VectorT<double> res = pnnls(a, b, 1);
  // cout << res << endl;

  // for (int i = 0; i < 50; ++i) {
  //   VectorT<int> r(10);
  // }


  // 0521: Test LCD
  // VectorT<double> r(3);
  // r << 10.99, 11.99, 13.99;

  // DynamicMatrixT<double> a(3, 3);
  // a << 1, 2, 3,
  //      4, 5, 6, 
  //      7, 8, 9;

  // VectorT<double> q = pnnls(a, r, 1);

  // cout << q << endl;
  // LCD lcd();

  // 0522: Test boost exp
  // mpf x = 3;
  // mpf res = boost::multiprecision::exp(x);
  // cout << res << endl;


  // 0522: Test type casting
  // DynamicMatrixT<mpf> a(3,3);
  // a << 1, 2, 3,
  //      4, 5, 6,
  //      7, 8, 9;
  // cout << a << endl << endl;
  // VectorT<mpf> v(3);
  // v << 10, 12, 14;

  vector<double> v1 = gn::generate_random(1000, 1000);
  sort(v1.begin(), v1.end());
  
  // vector<double> v2 = gn::generate_random(1000, 1000);

  VectorT<mpf> w1(1000);
  // VectorT<double> w2(1000);

  for (int i = 0; i < 1000; ++i) {
    w1(i) = v1[i];
    // w2(i) = v2[i];
  }

  // uint64_t ns = timing([&] {
  //   w2 - w1;
  // });

  // cout << ns << "ns, " << static_cast<double>(ns) / 1e+3 << "us\n";

  // ns = timing([&] {
  //   w2.array() - 10.9;
  // });

  // cout << ns << "ns, " << static_cast<double>(ns) / 1e+3 << "us\n";


  // 0522: Iteration time
  // size_t size = 1000;
  // mpf a = 0;
  // uint64_t ns = timing([&]{
  //   for (size_t i = 0; i < size; ++i) {
  //     a += *(w1.data() + i);
  //   }
  // });
  
  // cout << a << endl;
  // cout << ns << "ns, " << static_cast<double>(ns) / 1e+3 << "us\n";
  

  // 0522: partial sum
  // dropelementbyindex

  // vector<long double> v2;
  // for (auto i : v1) v2.push_back(i);


  // VectorT<long double> x(v1.size());
  // VectorT<int> w(v1.size());

  // weight(v2, x, w);

  // VectorT<long double> xx(x.cols());
  // cout << sizeof(long double) << endl;


  // 0523
  // VectorT<int> a(3);
  // a << 1, 2, 3;

  // MatrixT<int> b(3, 3);
  // b << 1, 2, 3,
  //      4, 5, 6,
  //      7, 8, 9;
  
  // b.matrix().colwise() += a.matrix();
  
  // cout << b << endl;

  // Eigen::MatrixXf mat1(2,4);
  // Eigen::VectorXf v1(2);
  // DynamicMatrixT<int> mat(3,3);
  // VectorT<int> v(3);
  
  // mat << 1, 2, 3, 
  //        5, 6, 7, 
  //        9, 10, 11;
         
  // v << 0,
  //      1,
  //      2;
  
  // // cout << columnWiseMult(mat, v); around 4us
       
  // cout << mat.array() * v.transpose().replicate(1, mat.cols()).array() << endl;


  // 0524: multiplication time
  VectorT<double> r5(1000);
  VectorT<int> w5(1000);

  w5.cast<double>();

  uint64_t ns = timing([&]{
    r5 * w5;
  });
  cout << ns << "ns, " << static_cast<double>(ns) / 1e+3 << "us\n";
  
  return 0;
}