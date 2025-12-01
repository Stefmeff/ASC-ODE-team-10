

#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

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


int main(int argc, char** argv)
{
  //Read command line argument: ExplicitEuler, ImplicitEuler, ImprovedEuler, CrankNicolson

  std::string method;
  if (argc > 1) {
      method = argv[1];
  } else {
      method = "ExplicitEuler"; // default method
  }

  // create the ODE right-hand side
  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

  std::unique_ptr<TimeStepper> stepper;
  std::string path;

  Vector<> Radau(3), RadauWeight(3);
  GaussRadau(Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  std::cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  std::endl;

  // Initialize Butcher tableaux for Gauss methods
  Vector<> Gauss2c(2), Gauss3c(3);

  // arbitrary order Gauss-LegendreR
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussLegendre(c, b1);
  auto [a, b] = ComputeABfromC(c);


  if (method=="ImprovedEuler") 
  {
    std::cout << "Using Improved Euler Method" << std::endl;
    stepper = std::make_unique<ImprovedEuler>(rhs);
    path = "../outputs/ImprovedEuler/";
  } else if (method=="CrankNicolson") {
    std::cout << "Using Crank-Nicolson Method" << std::endl;
    stepper = std::make_unique<CrankNicolson>(rhs);
    path = "../outputs/CrankNicolson/";
  } else if (method=="ImplicitEuler") {     
    std::cout << "Using Implicit Euler Method" << std::endl;
    stepper = std::make_unique<ImplicitEuler>(rhs);
    path = "../outputs/ImplicitEuler/";
  } else if(method == "RungeKutta") {
    std::cout << "Using Runge Kutta Method: " << path << std::endl;
    stepper = std::make_unique<RungeKutta>(rhs, Gauss2a, Gauss2b, Gauss2c);
    path = "../outputs/RungeKutta/";
  } else if(method=="ImplicitRungeKutta") {
    std::cout << "Using Implicit Runge Kutta: " << path << std::endl;
    auto [Gauss3a,Gauss3b] = ComputeABfromC (Gauss3c);
    stepper = std::make_unique<ImplicitRungeKutta>(rhs, Gauss3a, Gauss3b, Gauss3c);
    path = "../outputs/ImplicitRungeKutta/";
  } else {
    std::cout << "Using Explicit Euler Method" << std::endl;
    stepper = std::make_unique<ExplicitEuler>(rhs);
    path = "../outputs/ExplicitEuler/";
  } 

  double tend = 4*M_PI;
  int steps[] = {50, 100, 150, 200};

  // Gauss3c .. points tabulated, compute a,b:
  //ImplicitRungeKutta stepper(rhs, Gauss3a, Gauss3b, Gauss3c);


  /*
  // arbitrary order Gauss-Legendre
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussLegendre(c, b1);

  auto [a, b] = ComputeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */

  /* 
  // arbitrary order Radau
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussRadau(c, b1);

  auto [a, b] = ComputeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */


  for (int s : steps){

    double tau = tend/s;

    Vector<> y = { 1, 0 };  // initializer list
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
    std::ofstream outfile (path + "steps" + std::to_string(s) + ".txt");


    



    //std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

    for (int i = 0; i < s; i++)
    {
      stepper->DoStep(tau, y);

      //std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }


}
