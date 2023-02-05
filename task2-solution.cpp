#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#define print(Matrix) ( std::cout << "OUTPUT " << #Matrix << ":\n" << Matrix << "\n\n" )

#define printTask(name, err)                        \
(                                                   \
  std::cout << "--- TASK " << name << ") " <<       \
  (err<1E-8 ? "solved correctly" : "IS INCORRECT")<<\
  " (err=" << err << ")\n"                          \
)

Eigen::VectorXd eval_uN_Stupid(const Eigen::VectorXd &x, unsigned int M);
Eigen::VectorXd eval_uN(const Eigen::VectorXd &x, unsigned int M);
Eigen::VectorXd eval_F_Master(const Eigen::VectorXd &x);
Eigen::VectorXd eval_F(const Eigen::VectorXd &x);
Eigen::MatrixXd eval_DF_Master(const Eigen::VectorXd &x);
Eigen::MatrixXd eval_DF(const Eigen::VectorXd &x);

void task2() {
  using namespace Eigen;
  
  // SMALL TEST
  int N=3;
  int M=5;
  VectorXd x(N+1);
  x << 1, 12, 17, 42;
  
  // CALL (SLOW) CORRECT SOLUTION
  VectorXd correctSol2a = eval_uN_Stupid(x, M);
  VectorXd correctSol2c = eval_F_Master(x);
  MatrixXd correctSol2d = eval_DF_Master(x);
  
  // CALL FUNCTION OF CORRESPONDING EXERCISE
  VectorXd mySol2a = eval_uN(x, M);
  VectorXd mySol2c = eval_F(x);
  MatrixXd mySol2d = eval_DF(x);
  
  // CHECK IF SOLUTION IS CORRECT
  double err2a = (mySol2a-correctSol2a).norm();
  double err2c = (mySol2c-correctSol2c).norm();
  double err2d = (mySol2d-correctSol2d).norm();
  
  printTask("2a", err2a);
  printTask("2c", err2c);
  printTask("2d", err2d);
  
}

// to check for part (a), this is an O(n^2) solution
Eigen::VectorXd eval_uN_Stupid(const Eigen::VectorXd &x, unsigned int M) {
  Eigen::VectorXd sol(M);
  
  for (int i=0; i<M; i++) {
    double t = ((double) i) / M;
    double sum = 0;
    for (int j=0; j<x.size(); j++) {
      sum += x(j) * std::cos(2*M_PI*j*t);
    }
    sol(i) = sum;
  }
  
  return sol;
}

Eigen::VectorXd eval_uN(const Eigen::VectorXd &x, unsigned int M) {
  Eigen::VectorXd sol(M);
  
  // --- TODO: START ---
  
  // zero-extend our vector (faster solution when padding s.t. M' = 2^m)
  Eigen::VectorXcd expanded(M);
  expanded << x, Eigen::VectorXcd::Zero(M-x.size());
  
  // fast fourier transform to obtain the solution
  Eigen::FFT<double> fft;
  sol = fft.fwd(expanded).real();
  
  // ---  TODO: END  ---
  
  return sol;
}

// to check for part (c), this is the masters solution, but written as ugly as possible
Eigen::VectorXd eval_F_Master(const Eigen::VectorXd &x) {
  using A = Eigen::VectorXd; using D = int; using B = Eigen::ArrayXd;
  auto g = [&](A a) { return a.array().pow(3).matrix(); };
  auto G = [&](A a, D d) { return eval_uN_Stupid(a, d); };
  auto N = x.size()+1; auto C = M_PI*M_PI;
  B s = B::LinSpaced(N-1,0,N-2); A t = ((std::sqrt(C)/(N-1))*s).sin();
  return G(x.array()*4*C*s*s,N-1)+g(G(x,N-1))-t;
}

Eigen::VectorXd eval_F(const Eigen::VectorXd &x) {
  int n = x.size(); int N = n-1;
  Eigen::VectorXd sol(n);
  
  // --- TODO: START ---
  
  Eigen::VectorXd newX(n);
  for (int i=0; i<n; i++) newX(i) = 4*M_PI*M_PI*i*i*x(i);
  
  Eigen::VectorXd sinPart = Eigen::VectorXd::LinSpaced(n,0,M_PI*N/n);
  sinPart = sinPart.array().sin();
  
  // apply 0-2.b)
  sol = eval_uN(newX, n) + eval_uN(x, n).array().cube().matrix() - sinPart;
  
  // ---  TODO: END  ---
  
  return sol;
}

// to check for part (d), this is the masters solution, but written as ugly as possible
Eigen::MatrixXd eval_DF_Master(const Eigen::VectorXd &x) {
  using A = Eigen::VectorXd; using B = Eigen::MatrixXd;
  auto N=x.size()-1; auto E=[](int j) { return 4*M_PI*M_PI*j*j; };
  A D{eval_uN_Stupid(x, N+1)}; B C(N+1,N+1);
  for (auto k=0; k<=N; k++) for (auto j=0; j<=N; j++) 
  C(k, j) = (E(j)+3.0*D[k]*D[k]) * std::cos((2*M_PI*k)/(N+1)*j);
  return C;
}

Eigen::MatrixXd eval_DF(const Eigen::VectorXd &x) {
  const int n = x.size();
  Eigen::MatrixXd sol(n,n);
  
  // --- TODO: START ---
  
  auto tk = [&](int k) {return (double)k / n; };
  
  // part 1: F_k = \sum_{j=0}^N 4*\pi^2*j^2*x_j * \cos(2*\pi*j*t_k)
  // where the partial derivative \partial_{x_i}F_k = 4*\pi^2*i^2 * \cos(2*\pi*i*t_k)
  for (int k=0; k<n; k++) {
    for (int i=0; i<n; i++) {
      // k-th function of F with partial derivative x_i
      sol(k,i) = 4*M_PI*M_PI*i*i*std::cos(2*M_PI*i*tk(k));
    }
  }
  
  // part 2: F_k = (\sum_{j=0}^N x_j * \cos(2*\pi*j*t_k))^3
  // where the partial derivative \partial_{x_i}F_k = 3 * (\sum_{j=0}^N x_j*\cos(2*\pi*j*t_k))^2 * \cos(2*\pi*i*t_k)
  Eigen::VectorXd part2 = 3*eval_uN(x,n).array().square();
  
  for (int k=0; k<n; k++) {
    for (int i=0; i<n; i++) {
      // k-th function of F with partial derivative x_i
      sol(k,i) += part2(k) * std::cos(2*M_PI*i*tk(k));
    }
  }
  
  
  // part 3: F_k = -\sin(\pi * t_k)
  // where the partial derivative is 0.
  ;
  
  // ---  TODO: END  ---
  
  return sol;
}
