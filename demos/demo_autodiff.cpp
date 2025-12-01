#include <iostream>
#include <vector>
#include <matrix.hpp>
#include <fstream>
#include <nonlinfunc.hpp>
#include <autodiff.hpp>
#include <timestepper.hpp>


using namespace ASC_ode;


template <typename T>
T func1 (T x, T y)
{
  return x * sin(y);
  // return 1e6 + y;
}

//Exercise 18.4 Add Legendre Polynomials
template <typename T>
void LegendrePolynomials(int n, T x, std::vector<T>& P) {
  if (n < 0) {
      P.clear();
      return;
  }
  P.resize(n + 1);
  P[0] = T(1);
  if (n == 0) return;
  P[1] = x;
  for (int k = 2; k <= n; ++k) {
      P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
  }
}


//Exercise 18.5 Pendulum ODE with AutoDiff
class PendulumAD : public NonlinearFunction
{
private:
  double m_length;
  double m_gravity;

public:
  PendulumAD(double length, double gravity=9.81) : m_length(length), m_gravity(gravity) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    T_evaluate<double>(x, f);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    Vector<AutoDiff<2>> x_ad(2);
    Vector<AutoDiff<2>> f_ad(2);

    x_ad(0) = Variable<0>(x(0));
    x_ad(1) = Variable<1>(x(1));
    T_evaluate<AutoDiff<2>>(x_ad, f_ad);

    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
         df(i,j) = f_ad(i).deriv()[j];
  }

  template <typename T>
  void T_evaluate (VectorView<T> x, VectorView<T> f) const

  {
    f(0) = x(1);
    f(1) = -m_gravity/m_length*sin(x(0));
  }
};


int main()
{
  double x = 1, y = 2;
  AutoDiff<2> adx = Variable<0>(x);
  AutoDiff<2> ady = Variable<1>(y);

  std::cout << "adx = " << adx << std::endl;
  std::cout << "ady = " << ady << std::endl;

  AutoDiff<2> prod = adx * ady;
  std::cout << "prod = " << prod << std::endl;

  std::cout << "func1(adx, ady) = " << func1(adx, ady) << std::endl;

  double eps = 1e-8;
  std::cout << "numdiff df/dx = " << (func1(x + eps, y) - func1(x-eps, y)) / (2*eps) << std::endl;
  std::cout << "numdiff df/dy = " << (func1(x, y + eps) - func1(x, y-eps)) / (2*eps) << std::endl;


  {
    // we can do second derivatives:
    AutoDiff<1, AutoDiff<1>> addx{Variable<0>(2)};
    std::cout << "addx = " << addx << std::endl;
    // func = x*x
    // func' = 2*x
    // func'' = 2
    std::cout << "addx*addx = " << addx * addx << std::endl;

    // std::cout << "sin(addx) = " << sin(addx) << std::endl;
  }

  {
    //TODO: Exercise 18.4 Evaluate Legendre Polynomials
    //Evaluate Legendre Polynimials up to Order 5 in Interval [-1,1]

    //Create Output File for the Polynomials

    std::ofstream outfile("../outputs/Legendre/legendre_polynomials.txt");
    for (double x = -1.0; x <= 1.0; x += 0.01) {
        AutoDiff<1> adx = Variable<0>(x);
        std::vector<AutoDiff<1>> P;
        LegendrePolynomials(5, adx, P);
        outfile << x;
        for (size_t i=0; i<P.size(); ++i) {
            outfile << "," << P[i].value() << "," << P[i].deriv()[0];
        }
        outfile << "\n";
    }
    outfile.close();
  }

  {
    //TODO: 18.5 Test AutoDiff class for the pendulum
    // ODE: a'' = -g/l * sin(a); with a(0) = a0; a'(0) = a0p

    PendulumAD pendulum(1.0); // length = 1.0 m

    Vector<> x(2);
    x(0) = M_PI/4; // initial angle 45 degrees
    x(1) = 0.0;    // initial angular velocity 0

    Vector<> f(2);
    pendulum.evaluate(x, f);
    std::cout << "Pendulum ODE evaluation at a=pi/4, a'=0: " << std::endl;
    std::cout << "a' = " << f(0) << ", a'' = " << f(1) << std::endl;

    Matrix<> df(2,2);
    pendulum.evaluateDeriv(x, df);
    std::cout << "Pendulum ODE Jacobian at a=pi/4, a'=0: " << std::endl;
    std::cout << df << std::endl;
    
  }
  
  return 0;
}