#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <autodiff.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


/**
* Link für website zu github:
* Animation einfügen in jupyer notebooks:
* Wie bearbeite ich jupyter Reports
* Interface erklären:
* Alle abbiludungen in reports hochladen

*/

template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

class DistanceConstraint
{
public:
  double length;
  std::array<Connector,2> connectors;
};


template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<DistanceConstraint> m_constraints;  // distance constraints
  Vec<D> m_gravity=0.0;

//TODO: add distance constraints
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  size_t addDistanceConstraint (DistanceConstraint c)
  {
    m_constraints.push_back (c);
    return m_constraints.size()-1;
  }


  

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & constraints() { return m_constraints; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  //add additional dimensions for constraint Lagrange multipliers
  virtual size_t dimX() const override { return D*mss.masses().size() + mss.constraints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size() + mss.constraints().size(); }

  virtual void evaluate(VectorView<double> x, VectorView<double> f) const
  {

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto fmat = f.asMatrix(mss.masses().size(), D);

    //calculate forces on each mass

    //1. For each mass, calculate gravitational force
    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    //2. For each spring, calculate spring forces
    for (auto spring : mss.springs())
    {
      //get positions of connected masses/fixes
      auto [c1,c2] = spring.connectors;
      Vec<D> p1, p2;
      
      if (c1.type == Connector::FIX)
        // Position of fixed connector
        p1 = mss.fixes()[c1.nr].pos;
      else
        // Position of mass connector
        p1 = xmat.row(c1.nr);
      if (c2.type == Connector::FIX)
        p2 = mss.fixes()[c2.nr].pos;
      else
        p2 = xmat.row(c2.nr);

      // compute spring force magnitude
      double force = spring.stiffness * (norm(p1-p2)-spring.length);
      // direction from p1 to p2
      Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);

      // apply forces to masses in respective directions
      if (c1.type == Connector::MASS)
        fmat.row(c1.nr) += force*dir12;
      if (c2.type == Connector::MASS)
        fmat.row(c2.nr) -= force*dir12;
    }

    
    //3. Add Distance Constraints 
    for (size_t c = 0; c < mss.constraints().size(); c++)
    {
      DistanceConstraint distConst = mss.constraints()[c];
      auto [conn1, conn2] = distConst.connectors;
      Vec<D> p1, p2;
      
      // Get positions
      if (conn1.type == Connector::MASS)
        p1 = xmat.row(conn1.nr);
      else
        p1 = mss.fixes()[conn1.nr].pos;
        
      if (conn2.type == Connector::MASS)
        p2 = xmat.row(conn2.nr);
      else
        p2 = mss.fixes()[conn2.nr].pos;
      
      //Add Constraint forces: λ ∇g(x)
      double dist = norm(p1 - p2);
      Vec<D> dir12 = (dist > 1e-12) ? 1/dist * (p2-p1) : Vec<D>(0.0);
      double lambda = x(D * mss.masses().size() + c);
      
      if (conn1.type == Connector::MASS)
        fmat.row(conn1.nr) += lambda * dir12;
      if (conn2.type == Connector::MASS)
        fmat.row(conn2.nr) -= lambda * dir12;

      //Add Constraint equation g(x) = |p2 - p1| - length = 0
      f(D * mss.masses().size() + c) = dist - mss.constraints()[c].length;
    }

    //Divide forces by mass to get accelerations
    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;
}

virtual void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
{
    size_t nMass = mss.masses().size();
    size_t nConstr = mss.constraints().size();

    auto xmat = x.asMatrix(nMass, D);

    // Zero Jacobian
    df = 0.0;

    // 1. Spring contributions
    for (const auto& spring : mss.springs())
    {
        auto [c1, c2] = spring.connectors;

        Vec<D> p1 = (c1.type == Connector::MASS) ? xmat.row(c1.nr) : mss.fixes()[c1.nr].pos;
        Vec<D> p2 = (c2.type == Connector::MASS) ? xmat.row(c2.nr) : mss.fixes()[c2.nr].pos;

        Vec<D> diff = p2 - p1;
        double dist = norm(diff);
        if (dist < 1e-12) continue;

        Vec<D> dir = (1.0 / dist) * diff;

        double k = spring.stiffness;

        // Compute stiffness matrix K (D x D)
        double K[D][D];
        for (int a = 0; a < D; ++a)
            for (int b = 0; b < D; ++b)
                K[a][b] = k * ((dist - spring.length) / dist * (a == b ? 1.0 : 0.0) + dir(a) * dir(b));

        if (c1.type == Connector::MASS)
        {
            size_t i = c1.nr;
            if (c2.type == Connector::MASS)
            {
                size_t j = c2.nr;
                for (int a = 0; a < D; ++a)
                    for (int b = 0; b < D; ++b)
                    {
                        df(i*D + a, j*D + b) -= K[a][b] / mss.masses()[i].mass;
                        df(i*D + a, i*D + b) += K[a][b] / mss.masses()[i].mass;
                        df(j*D + a, i*D + b) -= K[a][b] / mss.masses()[j].mass;
                        df(j*D + a, j*D + b) += K[a][b] / mss.masses()[j].mass;
                    }
            }
            else
            {
                for (int a = 0; a < D; ++a)
                    for (int b = 0; b < D; ++b)
                        df(i*D + a, i*D + b) += K[a][b] / mss.masses()[i].mass;
            }
        }
        else if (c2.type == Connector::MASS)
        {
            size_t j = c2.nr;
            for (int a = 0; a < D; ++a)
                for (int b = 0; b < D; ++b)
                    df(j*D + a, j*D + b) += K[a][b] / mss.masses()[j].mass;
        }
    }

    // 2. Distance constraints contributions
    for (size_t c = 0; c < nConstr; ++c)
    {
        auto distConst = mss.constraints()[c];
        auto [conn1, conn2] = distConst.connectors;

        Vec<D> p1 = (conn1.type == Connector::MASS) ? xmat.row(conn1.nr) : mss.fixes()[conn1.nr].pos;
        Vec<D> p2 = (conn2.type == Connector::MASS) ? xmat.row(conn2.nr) : mss.fixes()[conn2.nr].pos;

        Vec<D> diff = p2 - p1;
        double dist = norm(diff);
        if (dist < 1e-12) continue;

        Vec<D> dir = (1.0 / dist) * diff;

        // df/dx for constraints
        if (conn1.type == Connector::MASS)
            for (int a = 0; a < D; ++a){
                df(nMass*D + c, conn1.nr*D + a) = -dir(a);               
                df(conn1.nr*D + a, nMass*D + c) = dir(a) / mss.masses()[conn1.nr].mass;
            }

        if (conn2.type == Connector::MASS)
            for (int a = 0; a < D; ++a){
                df(nMass*D + c, conn2.nr*D + a) = dir(a);                
                df(conn2.nr*D + a, nMass*D + c) = -dir(a) / mss.masses()[conn2.nr].mass;
            }
    }
}

};

#endif