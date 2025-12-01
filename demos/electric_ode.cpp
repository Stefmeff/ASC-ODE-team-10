#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


class RCCircuit : public NonlinearFunction
{
private:
  double Resistance;
  double Capacitance;

public:
  RCCircuit(double R, double C) : Resistance(R), Capacitance(C) {}
  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};