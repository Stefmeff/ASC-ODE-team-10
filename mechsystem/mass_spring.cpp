#include "mass_spring.hpp"
#include "Newmark.hpp"

int main()
{
  //erstelle ein MassSpringSystem in 2D
  MassSpringSystem<2> mss;
  mss.setGravity( {0,-9.81} );
  auto fA = mss.addFix( { { 0.0, 0.0 } } );
  auto mA = mss.addMass( { 1, { 1.0, 0.0 } } );
  auto mB = mss.addMass( { 1, { 2.0, 0.0 } } );
  auto mC = mss.addMass( { 3, { 3.0, 0.0 } } );
  auto mD = mss.addMass( { 4, { 4.0, 0.0 } } );

  mss.addSpring ( { 1, 200, { fA, mA }}  );
  mss.addSpring ( { 1, 100, { mA, mB }} );

  //mss.addDistanceConstraint({1,{fA,mA}});
  mss.addSpring ( { 1, 200, { mB, mC }}  );
  mss.addSpring ( { 1, 100, { mC, mD }} );

  std::cout << "mss: " << std::endl << mss << std::endl;


  double tend = 10;
  double steps = 1000;

  Vector<> x(2*mss.masses().size());
  Vector<> dx(2*mss.masses().size());
  Vector<> ddx(2*mss.masses().size());

  auto mss_func = std::make_shared<MSS_Function<2>> (mss);
  auto mass = std::make_shared<IdentityFunction> (x.size());

  mss.getState (x, dx, ddx);

  std::cout << "Initial state: x = " << x
            << ", dx = " << dx << std::endl;
  
  SolveODE_Newmark(tend, steps, x, dx,  mss_func, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t
                                                             << ", x = " << Vec<4>(x) << std::endl; });
}
