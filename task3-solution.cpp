#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include "matplotlibcpp.h"

#define print(Matrix) ( std::cout << "OUTPUT " << #Matrix << ":\n" << Matrix << "\n\n" )

#define printTask(name, err)                        \
(                                                   \
  std::cout << "--- TASK " << name << ") " <<       \
  (err<1E-8 ? "solved correctly" : "IS INCORRECT")<<\
  " (err=" << err << ")\n"                          \
)


class ClosedQuadraticSplineCurve {
 public:
  // Constructor
  explicit ClosedQuadraticSplineCurve(const Eigen::Matrix2Xd &p);
  // ~ClosedQuadraticSplineCurve() = default;
  // ClosedQuadraticSplineCurve() = delete;
  // ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &); // = delete;
  // ClosedQuadraticSplineCurve(const ClosedQuadraticSplineCurve &&); // = delete;
  // ClosedQuadraticSplineCurve &operator=(const ClosedQuadraticSplineCurve &) = delete;
  // Point evaluation operator for sorted parameter arguments
  [[nodiscard]] Eigen::Matrix2Xd curve_points(const Eigen::VectorXd &v) const;
  [[nodiscard]] Eigen::Matrix2Xd curve_points_MASTERS(const Eigen::VectorXd &v) const;
  // Curvature evaluation operator for sorted parameter arguments
  // [[nodiscard]] Eigen::VectorXd local_curvatures(
      // const Eigen::VectorXd &v) const;
  // Approximate computation of the length of the curve
  // [[nodiscard]] double length(double rtol = 1E-6) const;

 private:
  // [[nodiscard]] bool checkC1(double tol = 1.0E-4) const;
  // Number of points to be interpolated
  unsigned int n_;
  // Coordinates of interpolated points, duplicate: $\cob{\Vp^0=\Vp^n}$
  Eigen::Matrix2Xd p_;
  // Knot sequence
  Eigen::VectorXd t_;
  // Coefficients in local representation
  Eigen::Matrix2Xd x_;
};

void plot(Eigen::Matrix2Xd p, Eigen::Matrix2Xd curve) {
  namespace plt = matplotlibcpp;
  
  // if we plot the rows directly, it somehow doesn't work
  Eigen::VectorXd x = p.row(0);
  Eigen::VectorXd y = p.row(1);
  plt::plot(x,y, "r*");
  
  x = curve.row(0);
  y = curve.row(1);
  plt::plot(x,y);

  
  plt::savefig("./cx_out/PLOT.png");
}


void task3() {
  using namespace Eigen;
  
  // SMALL EXAMPLE
  int n=5;
  Matrix2Xd p(2, n);
  p << 0.5, 0, 0, 1, 1,
       1.4, 1, 0, 0, 1;
  
  // GENERATE CQSC (CONSTRUCTOR)
  ClosedQuadraticSplineCurve cqsc = ClosedQuadraticSplineCurve(p);
  
  // EVALUATE AT POINTS
  int points = 200;
  VectorXd v = VectorXd::LinSpaced(points, 0, 1);
  Matrix2Xd mySol = cqsc.curve_points(v);
  
  // CHECK (MATHEMATICALLY)
  Matrix2Xd correctSol = cqsc.curve_points_MASTERS(v);
  double err = (correctSol-mySol).norm();
  printTask("3c", err);
  
  // CHECK (GRAPHICALLY)
  plot(p, mySol);
  
}

ClosedQuadraticSplineCurve::ClosedQuadraticSplineCurve(const Eigen::Matrix2Xd &p) : n_(p.cols()), p_(2, p.cols() + 1), t_(p.cols() + 1), x_(2, p.cols()) {
  assert((n_ > 2) && "At least three points have to be supplied");
  assert((n_ % 2 == 1) && "Number of points must be odd!");
  // Save point coordinates, duplicate endpoints
  p_.col(0) = p.col(n_ - 1);
  p_.rightCols(n_) = p;
  
  // START (CONSTRUCTOR)
  
  // INITIALIZE t_
  double sum = 0; // lin. interpolation distance
  for (int k=1; k<=n_; k++) {
    sum += (p_.col(k)-p_.col(k-1)).norm();
  }
  
  t_(0) = 0;
  double partialSum = 0;
  for (int k=1; k<=n_; k++) {
    partialSum += (p_.col(k)-p_.col(k-1)).norm();
    t_(k) = partialSum / sum;
  }
  
  // INITIALIZE x_
  Eigen::MatrixXd twoDiag = Eigen::MatrixXd::Zero(n_, n_);
  twoDiag.diagonal() = Eigen::VectorXd::Ones(n_);
  twoDiag.diagonal(+1) = Eigen::VectorXd::Ones(n_-1);
  twoDiag(n_-1,0) = 1;
  
  auto hj = [&](int j) { return t_(j) - t_(j-1); };
  
  Eigen::Matrix2Xd rhs(2,n_);
  for (int i=1; i<n_; i++) {
    rhs.col(i-1) = (p_.col(i)-p_.col(i-1))/hj(i) - (p_.col(i+1)-p_.col(i))/hj(i+1);
  }
  rhs.col(n_-1) = (p_.col(n_)-p_.col(n_-1))/hj(n_) - (p_.col(1)-p_.col(0))/hj(1);
  
  x_ = twoDiag.lu().solve(rhs.transpose()).transpose();
  
  // END
}



Eigen::Matrix2Xd ClosedQuadraticSplineCurve::curve_points(const Eigen::VectorXd &v) const {
  Eigen::Matrix2Xd sol(2, v.size());
  
  // --- TODO: START ---
  
  // NOTE:
  // t_ goes from 0..n (inclusive)
  // t_(0) = 0; t_(n) = 1
  
  int j=1; // meaning between t_0 and t_1
  for (int vVec=0; vVec<v.size(); vVec++) {
    
    // search right segment
    while (!(t_(j-1) <= v(vVec) && v(vVec) <= t_(j))) {
      j++;
      if (j>n_) std::cout << "Error, matching segment wasn't found\n\n";
    }
    
    // calculate s(t)
    double myH   = t_(j) - t_(j-1);
    double myTau = (v(vVec) - t_(j-1)) / myH;
    
    sol.col(vVec) = p_.col(j-1)*(1.-myTau) + x_.col(j-1)*myH*myTau*(1.-myTau) + p_.col(j)*myTau;
  }
  
  // ---  TODO: END  ---
  
  return sol;
}



// MASTERS SOLUTION
Eigen::Matrix2Xd ClosedQuadraticSplineCurve::curve_points_MASTERS(
    const Eigen::VectorXd &v) const {
  unsigned int N = v.size();
  assert(N > 0);
  // Matrix containing points to be computed
  Eigen::Matrix2Xd s(2, N);
  // Lengths of knot intervals
  const Eigen::VectorXd h = t_.tail(n_) - t_.head(n_);
  // TODO: (5-15.d) Compute the points on the curve for the sorted parameter
  // values v.
  // START

  // Run through all the provided parameter values
  unsigned int knot_idx = 0;
  for (unsigned int k = 0; k < N; ++k) {
    assert(((v[k] >= 0.0) && (v[k] <= 1.0)) && "Out of range!");
    if (k > 0) {
      assert(v[k] >= v[k - 1] && "Parameter values not sorted!");
    }
    // Find right knot interval: knot\_idx stores index of knot to the right of
    // current parameter value
    while ((knot_idx < n_) && (v[k] >= t_[knot_idx])) {
      knot_idx++;
    }
    assert(knot_idx > 0);
    const double tau = (v[k] - t_[knot_idx - 1]) / h[knot_idx - 1];
    s.col(k) = ((1.0 - tau) * p_.col(knot_idx - 1)) +
               (h[knot_idx - 1] * tau * (1 - tau)) * x_.col(knot_idx - 1) +
               (tau * p_.col(knot_idx));
  }
  // END
  return s;
}
