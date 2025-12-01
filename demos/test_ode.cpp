

#include <iostream>
#include <fstream> 

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

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
    } else {
      std::cout << "Using Explicit Euler Method" << std::endl;
      stepper = std::make_unique<ExplicitEuler>(rhs);
      path = "../outputs/ExplicitEuler/";
    }

  double tend = 4*M_PI;
  int steps[] = {50, 100, 150, 200};


  for (int s : steps){

      double tau = tend/s;

      // initializer list
      Vector<> y = { 1, 0 };  
      auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

      std::ofstream outfile (path + "steps" + std::to_string(s) + ".txt");

      std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
      outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

      for (int i = 0; i < s; i++)
      {
        stepper->DoStep(tau, y);

        std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
        outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      }
  }


}
