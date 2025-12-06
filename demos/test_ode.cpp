

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

//TODO: Exercise 17.4.1
class ElectricNetwork : public NonlinearFunction
{
private:
    //  variables
    double resistance;
    double capacitance;
    double angularFrequency; 

public:
    // Constructor
    ElectricNetwork(double R, double C) 
        : resistance(R), capacitance(C) 
    {
        // 100 * PI is the angular frequency (omega)
        angularFrequency = 100.0 * M_PI;
    }

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }
    
    // f(x)
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
        double Uc = x(0); // Voltage
        double t  = x(1); // Time
        
        // dUc/dt = (1/RC) * (U0(t) - Uc)
        
        f(0) = (1.0 / (resistance * capacitance)) * (std::cos(angularFrequency * t) - Uc);
        
        // dt/dt = 1
        f(1) = 1.0; 
    }
    
    // Jacobian J = df/dx
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
        double t = x(1);
        
        // Use new names for calculation
        double rc_inv = 1.0 / (resistance * capacitance);

        df = 0.0;

        // d(f0)/d(x0) -> d(U'c)/dUc = -1/RC
        df(0,0) = -rc_inv;
        
        // d(f0)/d(x1) -> d(U'c)/dt  = (1/RC) * (-omega * sin(omega*t))
        df(0,1) = rc_inv * (-angularFrequency * std::sin(angularFrequency * t));

        // d(f1)/... -> derivatives of 1 are 0
        df(1,0) = 0.0;
        df(1,1) = 0.0;
    }
};

int main(int argc, char** argv)
{
  // Usage: ./test_ode [System] [Method] [num_stages]
  // System: MassSpring, ElectricNetwork
  // Method: ExplicitEuler, ImplicitEuler, ImprovedEuler, CrankNicolson, RungeKutta, ImplicitRungeKutta
  // For ImplicitRungeKutta, num_stages specifies the number of stages

  //Read command line argument: ExplicitEuler, ImplicitEuler, ImprovedEuler, CrankNicolson
  std::string system;
  std::string method;
  int num_stages = 2; // default number of stages for RK methods

  if (argc == 4 ) {
      system = argv[1];
      method = argv[2];
      num_stages = std::stoi(argv[3]);
  } else if (argc = 3 ) {
      system = argv[1];
      method = argv[2];
  } else if (argc == 2 ) {
      system = argv[1];
      method = "ExplicitEuler"; // default method
  } else {
      system = "MassSpring"; // default system
      method = "ExplicitEuler"; // default method
  }

  //Select System: Electrical vs Mechanical
  std::shared_ptr<NonlinearFunction> rhs;
  if(system == "ElectricNetwork") {
    std::cout << "Using Electric Network System" << std::endl;
    rhs = std::make_shared<ElectricNetwork>(1.0,1.0); // R=1kOhm, C=1uF
  }else{
    std::cout << "Using Mass-Spring System" << std::endl;
    rhs = std::make_shared<MassSpring>(1.0, 1.0);
  }

  std::unique_ptr<TimeStepper> stepper;
  std::string path;

  Vector<> Radau(3), RadauWeight(3);
  GaussRadau(Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  std::cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  std::endl;

  // Initialize Butcher tableaux for Gauss methods
  //Vector<> Gauss2c(2), Gauss3c(3);

  Matrix<> a(num_stages, num_stages);
  Vector<> b(num_stages), c(num_stages);
  if(num_stages == 2) {
    auto [Gauss_2a, Gauss_2b] = ComputeABfromC(Gauss2c);
    std::cout << "Gauss2c = " << Gauss2c << ", b = " << Gauss_2b << ", a = " << Gauss_2a << std::endl;
    a = Gauss_2a;
    b = Gauss_2b;
    c = Gauss2c;

  } else if(num_stages == 3) {
    auto [Gauss_3a, Gauss_3b] = ComputeABfromC(Gauss3c);
    std::cout << "Gauss3c = " << Gauss3c << ", b = " << Gauss_3b << ", a = " << Gauss_3a << std::endl;
    a = Gauss_3a;
    b = Gauss_3b;
    c = Gauss3c;
  } else {
    //arbitrary order Gauss-Legendre
    Vector<> b1(num_stages);
    GaussLegendre(c, b1);
    auto [a_temp, b_temp] = ComputeABfromC(c);
    std::cout << "GaussLegendre c = " << c << ", b = " << b_temp << ", a = " << a_temp << std::endl;
    a = a_temp;
    b = b_temp;
  }
 



  // Select Time Stepping Method
  if (method=="ImprovedEuler") 
  {
    std::cout << "Using Improved Euler Method" << std::endl;
    stepper = std::make_unique<ImprovedEuler>(rhs);
    path = "../outputs/" + system + "/ImprovedEuler/";
  } else if (method=="CrankNicolson") {
    std::cout << "Using Crank-Nicolson Method" << std::endl;
    stepper = std::make_unique<CrankNicolson>(rhs);
    path = "../outputs/" + system + "/CrankNicolson/";
  } else if (method=="ImplicitEuler") {     
    std::cout << "Using Implicit Euler Method" << std::endl;
    stepper = std::make_unique<ImplicitEuler>(rhs);
    path = "../outputs/" + system + "/ImplicitEuler/";
  } else if(method == "ExpRungeKutta") {
    std::cout << "Using Runge Kutta Method: " << path << std::endl;
    stepper = std::make_unique<RungeKutta>(rhs, a, b, c);
    path = "../outputs/" + system + "/ExplicitRungeKutta/stages:" + std::to_string(num_stages);
  } else if(method=="ImpRungeKutta") {
    std::cout << "Using Implicit Runge Kutta: " << path << std::endl;
    auto [Gauss3a,Gauss3b] = ComputeABfromC (Gauss3c);
    stepper = std::make_unique<ImplicitRungeKutta>(rhs, a, b, c);
    path = "../outputs/" + system + "/ImplicitRungeKutta/stages:" + std::to_string(num_stages);
  } else {
    std::cout << "Using Explicit Euler Method" << std::endl;
    stepper = std::make_unique<ExplicitEuler>(rhs);
    path = "../outputs/" + system + "/ExplicitEuler/";
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

  auto [a, b] = computeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */

  /* 
  // arbitrary order Radau
  int stages = 5;
  Vector<> c(stages), b1(stages);
  GaussRadau(c, b1);

  auto [a, b] = computeABfromC(c);
  ImplicitRungeKutta stepper(rhs, a, b, c);
  */


  for (int s : steps){

    double tau = tend/s;

    Vector<> y = { 1, 0 };  // initializer list
    std::ofstream outfile (path + "steps:" + std::to_string(s) + ".txt");


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
