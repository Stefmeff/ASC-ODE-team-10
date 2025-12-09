

#include <iostream>
#include <fstream> 
#include <getopt.h>
#include <stdlib.h>
#include <cstdlib>
#include <filesystem>


#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>
#include <filesystem>
#include <iostream>

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
        angularFrequency = 100 * M_PI;
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
  // Usage: ./test_ode [--system <type>] [--method <name>] [--stages N]
  // -s: MassSpring, ElectricNetwork
  // -m: ExplicitEuler, ImplicitEuler, ImprovedEuler, CrankNicolson, ExpRungeKutta, ImpRungeKutta
  // -k: Number of stages for Runge-Kutta methods

  //arguments:
  std::string system;
  std::string method;
  int num_stages = 2; // default number of stages for RK methods

  //PARSE COMMAND LINE OPTIONS:
  static struct option long_options[] = {
      {"system", required_argument, 0, 's'},
      {"method", required_argument, 0, 'm'},
      {"stages", required_argument, 0, 'k'}, 
      {0, 0, 0, 0}
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "s:m:k:", long_options, NULL)) != -1) {
      switch (opt) {
      case 's':
          system = optarg;
          break;
      case 'm':
          method = optarg;
          break;
      
      case 'k':
          num_stages = atoi(optarg);
          break;
      default:
          fprintf(stderr,
              "Usage: %s --system <type> --method <name> [--stages N]\n", argv[0]);
          exit(EXIT_FAILURE);
      }   
  }

  //Initial Conditions
  Vector<> y = {1.0, 0.0};  // initializer list
  double tend = 4*M_PI;
  int steps[] = {5, 50, 100, 150, 200, 20000};

  //Select System: Electrical vs Mechanical
  std::shared_ptr<NonlinearFunction> rhs;
  if(system == "ElectricNetwork") {
    //ELECTRIC NETWORK
    std::cout << "Using Electric Network System" << std::endl;
    rhs = std::make_shared<ElectricNetwork>(1.0,1.0); // R=1kOhm, C=1uF
    tend = 1;
  }else{
    //MASS SPRING SYSTEM
    std::cout << "Using Mass-Spring System" << std::endl;
    rhs = std::make_shared<MassSpring>(1.0, 1.0);
    system = "MassSpring";
  }

  std::unique_ptr<TimeStepper> stepper;
  

  Vector<> Radau(3), RadauWeight(3);
  GaussRadau(Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  std::cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  std::endl;

  // Initialize Butcher tableaux for Gauss methods
  //Vector<> Gauss2c(2), Gauss3c(3);

  Matrix<> a(num_stages, num_stages);
  Vector<> b(num_stages), c(num_stages);

  

  if(num_stages == 2) {
    auto [Gauss_2a, Gauss_2b] = computeABfromC(Gauss2c);
    std::cout << "Gauss2c = " << Gauss2c << ", b = " << Gauss_2b << ", a = " << Gauss_2a << std::endl;
    a = Gauss_2a;
    b = Gauss_2b;
    c = Gauss2c;

  } else if(num_stages == 3) {
    auto [Gauss_3a, Gauss_3b] = computeABfromC(Gauss3c);
    std::cout << "Gauss3c = " << Gauss3c << ", b = " << Gauss_3b << ", a = " << Gauss_3a << std::endl;
    a = Gauss_3a;
    b = Gauss_3b;
    c = Gauss3c;
  } else {
    //arbitrary order Gauss-Legendre
    Vector<> b1(num_stages);
    GaussLegendre(c, b1);
    auto [a_temp, b_temp] = computeABfromC(c);
    std::cout << "GaussLegendre c = " << c << ", b = " << b_temp << ", a = " << a_temp << std::endl;
    a = a_temp;
    b = b_temp;
  }
  std::string path;
 
  //SELECT TIME STEPPING METHOD:
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
    auto [Gauss3a,Gauss3b] = computeABfromC (Gauss3c);
    stepper = std::make_unique<ImplicitRungeKutta>(rhs, a, b, c);
    path = "../outputs/" + system + "/ImplicitRungeKutta/stages:" + std::to_string(num_stages);
  } else {
    std::cout << "Using Explicit Euler Method" << std::endl;
    stepper = std::make_unique<ExplicitEuler>(rhs);
    path = "../outputs/" + system + "/ExplicitEuler/";
  } 


  //clear output directory at path

  
  //SIMULATION LOOP:
  for (int s : steps){

    double tau = tend/s;

    y = (system == "ElectricNetwork") ? Vector<>({0.0, 0.0}) : Vector<>({1.0, 0.0}); // reset initial conditions

    printf("Simulating with %d steps, tau=%f\n", s, tau);
    std::ofstream outfile (path + "steps:" + std::to_string(s) + ".txt");


    //std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

    for (int i = 0; i < s; i++)
    {
      stepper->doStep(tau, y);

      //std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }
}
