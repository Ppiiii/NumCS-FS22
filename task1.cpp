#include <iostream>
#include <Eigen/Core>

#define print(Matrix) ( std::cout << "OUTPUT " << #Matrix << ":\n" << Matrix << "\n\n" )

#define printTask(name, err)                        \
(                                                   \
  std::cout << "--- TASK " << name << ") " <<       \
  (err<1E-8 ? "solved correctly" : "IS INCORRECT")<<\
  " (err=" << err << ")\n"                          \
)


class CompactStorageQR {
  public:
    // Constructor taking raw data array, which has to be persistent
    CompactStorageQR (const double *data, unsigned int n) : n_(n), M_(data, n, n) {
      for (unsigned int k = 0; k < n_-1; ++k ) {
        double sn { 0 };
        for (unsigned int j = k + 1; j < n_; ++j ) {
          const double t {M_(j, k)};
          sn += t*t ;
        }
        if ( sn > 1.0 ) {
          throw std::runtime_error(
            "CompactStorageQR: Illegal subdiagonal column norm! ");
        }
      }
    }
    // Determinant of the matrix stored in this object
    double det() const;
    // Right multiplication with another matrix
    Eigen::MatrixXd matmult(const Eigen::MatrixXd &X) const;
    // Solution of linear systems of equations
    Eigen::MatrixXd solve(const Eigen::MatrixXd &B) const;
    
  private:
    unsigned int n_; // Matrix dimensions
    // Raw data wrapped into a matrix, see [Lecture → ??]
    Eigen::Map<const Eigen::MatrixXd> M_;
    
};



// own (wikipedia) implementation of QR-decomp, can be ignored
Eigen::MatrixXd decomp(const Eigen::MatrixXd& A) {
  int n = A.cols();
  
  Eigen::MatrixXd vVectors(n,n);
  Eigen::MatrixXd R = A;
  for (int i=0; i<n-1; i++) {
    Eigen::VectorXd z = R.col(i).tail(n-i);
    
    Eigen::VectorXd uk = z;
    uk(0) += (z(0)>0 ? 1 : -1)*z.norm();
    uk = uk/uk.norm();
    
    Eigen::VectorXd vk = Eigen::VectorXd::Zero(n);
    vk.tail(n-i) = uk;
    
    R = R - 2*vk*(vk.transpose()*R);
    vVectors.col(i) = vk;
  }
  vVectors.triangularView<Eigen::Upper>() = R;
  
  return vVectors;
}


void task1() {
  using namespace Eigen;
  
  // SMALL TEST
  int n = 3;
  int k = 3;
  
  MatrixXd A(n,n);
  MatrixXd X(n,k);
  // A <<  12,-51,  4,
  //       6,167,-68,
  //       -4, 24,-41;
  A << 1,0,0,
      0,3,0,
      0,4,1;
  X << 1,2,3,
      4,5,6,
      7,8,9;
      
  // ------------
  
  // PREPARATION OF FORTRAN77-STYLE DATA
  Eigen::Matrix<double,Dynamic,Dynamic,ColMajor> M = decomp(A);
  double* raw = M.data();
  
  // CALL OF COMPACTSTORAGEQR FUNCTION
  CompactStorageQR obj = CompactStorageQR(raw, n);
  MatrixXd mySol1a = obj.matmult(X);
  MatrixXd mySol1b = obj.solve(X);
  double mySol1c   = obj.det();
  
  // EIGEN CORRECT SOLUTION
  MatrixXd correctSol1a = A*X;
  MatrixXd correctSol1b = A.lu().solve(X);
  double correctSol1c   = A.determinant();
  
  // CHECK IF SOLUTION IS CORRECT
  double err1a = (mySol1a-correctSol1a).norm();
  double err1b = (mySol1b-correctSol1b).norm();
  double err1c = (mySol1c-correctSol1c);
  
  printTask("1a", err1a);
  printTask("1b", err1b);
  printTask("1c", err1c);
}

Eigen::MatrixXd CompactStorageQR::matmult(const Eigen::MatrixXd &X) const {
  assert((X.rows() == n_) && "Wrong size of matrix factor");
  Eigen::MatrixXd sol(n_,X.rows());
  
  // --- TODO: START ---
  
  // extract R matrix from the upper right half of the matrix M
  Eigen::MatrixXd R = M_.template triangularView<Eigen::Upper>();
  
  // calculate sol as Q * R * X, where Q is multiplied in the next step
  sol = R*X;
  
  // multiply Q_0 * Q_1 * ... * Q_(n-2) * sol, starting with the rightmost multiplication
  for (int i=n_-2; i>=0; i--) {
    int vecLength = n_ - i;
    Eigen::VectorXd v(n_);
    
    // each v_k is calculated using the lower part of the k-th column of M
    v.tail(vecLength-1) = M_.col(i).tail(vecLength-1);
    v(i) = std::sqrt(1 - v.tail(vecLength-1).squaredNorm());
    
    // each Q_k * sol is calculated as (I - 2*v*v^T)*sol
    // using associativity and distributivity, we get sol - 2*v*(v^T*sol)
    // (v*v^T)*sol is O(n^3) and thus way slower than v*(v^T*sol), which is O(n^2)
    sol = sol - 2 * (v*(v.transpose()*sol));
  }
  
  // ---  TODO: END  ---
  
  return sol;
}

Eigen::MatrixXd CompactStorageQR::solve(const Eigen::MatrixXd &B) const {
  assert((B.rows() == n_) && "Wrong size of matrix factor");
  Eigen::MatrixXd sol(n_,B.rows());
  
  // --- TODO: START ---
  
  // calculate sol as R.solve(Q_(n-2)*...*Q_0*B)
  sol = B;
  
  // similar to (a), but starting with the multiplication Q_0*B instead of Q_(n-2)*B
  // this is because Q^(-1) = Q^T = (Q_0*Q_1*...*Q_(n-2))^T = Q_(n-2)*...*Q_0
  for (int i=0; i<n_-1; i++) {
    Eigen::VectorXd v(n_);
    
    v.tail(n_-i-1) = M_.col(i).tail(n_-i-1);
    v(i) = std::sqrt(1 - v.tail(n_-i-1).squaredNorm());
    
    sol = sol - 2*v*(v.transpose()*sol);
  }
  
  // we are able to use solve because Eigen knows that this is a triangular matrix
  sol = M_.template triangularView<Eigen::Upper>().solve(sol);
  
  // ---  TODO: END  ---
  
  return sol;
}

double CompactStorageQR::det() const {
  
  // --- TODO: START ---
  
  // det(A) = det(QR) = det(Q)*det(R) = +/- det(R) = +/- prod(diagonalEntries(R))
  return M_.diagonal().prod() * (n_%2==0 ? -1 : 1);
  
  // ---  TODO: END  ---
}




    
    
    
