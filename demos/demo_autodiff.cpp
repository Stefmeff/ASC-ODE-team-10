#include <iostream>
#include <autodiff.hpp>
#include <vector>
#include <fstream>


using namespace ASC_ode;


template <typename T>
T func1 (T x, T y)
{
  return x * sin(y);
  // return 1e6 + y;
}

// TODO: Exercise 18.4 Add Legengre Polynomials
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
  
  return 0;
}