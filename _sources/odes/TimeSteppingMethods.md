# 2. Time Stepping Methods

The application implements various time-stepping methods used for finding approximate solutions to ODE's found in `src/timestepper.hpp`. Time stepping advances an ODE solution from time steps $t_n$ to $t_{n+1}=t_n+\tau$ by approximating the RHS $f(y)$. Our steppers share a common interface `TimeStepper`, which stores the RHS as a `NonlinearFunction` and provides `doStep` function to calculate the next time step.

```cpp
class TimeStepper
{
protected:
  std::shared_ptr<NonlinearFunction> m_rhs;
public:
  TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
  virtual ~TimeStepper() = default;
  virtual void doStep(double tau, VectorView<double> y) = 0;
};
```
The following subchapters explain the idea and concrete implementation of various time-stepping methods like Explicit/Implicit Euler, Crank Nicolson, Runge-Kutta etc.

## 2.1 Explicit Euler Method

The Explicit Euler method is a numerical scheme for solving ordinary differential equations (ODEs) of the form:

$$\frac{dy}{dt} = f(y(t)), \quad y(0) = y_0, \quad \forall t \in [0,T]$$

The method advances the solution from time step $n$ to $n+1$ via:

$$y_{n+1} = y_n + \tau \, f(y_n)$$

where:
- $y_n$ is the solution at time $t_n$
- $\tau$ = T/n is the time step size 
- $f(y_n)$ is the right-hand side (RHS) function evaluated at the current state

The `ExplicitEuler` class is implemented in `src/timestepper.hpp` and inherits from the abstract `TimeStepper` base class. The do-step method is implemented as 

```cpp
void doStep(double tau, VectorView<double> y) override
{
    //evaluate f at y(n)
    this->m_rhs->evaluate(y, m_vecf);
    //update y(n+1) = y(n) + tau*f(y(n))
    y += tau * m_vecf;
}
```
## 2.2 Improved Euler Method

Improved Euler is based on Explicit Euler but tries to achieve better accuracy by approximating the midpoint:

$$ \tilde{y} = y_n+\tfrac{\tau}{2}f(y_n)$$

And using this additional information to compute a better approximation for $y_{n+1}$ using:

$$y_{n+1} = y_n + \tau \, f(\tilde{y})$$

We implement this method in the The `ImprovedEuler` class:

```cpp

void doStep(double tau, VectorView<double> y) override
{
    //evaluate midpoint
    this->m_rhs->evaluate(y, m_vecf);
    Vector<> y_tmp = y + (tau/2) * m_vecf;
    //update y_n+1
    this->m_rhs->evaluate(y_tmp, m_vecf);
    y += (tau) * m_vecf;
}
```

## 2.3 Implicit Euler Method

Implicit Euler uses the future state in the slope to achieve more accurate results:

$$y_{n+1} = y_n + \tau \, f(y_{n+1})$$

Since $f(y_{n+1})$ is not given explicitly as before, but implicitly, we need to solve the resulting system of equations using Newtons method. The method is implemented in the class `ImplicitEuler`. 
```cpp
ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
: TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
{
  m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
  auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
  m_equ = ynew - m_yold - m_tau * m_rhs;  // r(y_{n+1})
}

void doStep(double tau, VectorView<double> y) override
{
  //set known values
  m_yold->set(y);   
  m_tau->set(tau); 
  //solve for y_n+1 
  NewtonSolver(m_equ, y); 
}
```

## 2.4 Crank Nicolson

Crank–Nicolson is another implicit time-stepping method, averaging the derivatives at the current next time step:

$$y_{n+1} = y_n + \frac{\tau}{2}\bigl(f(y_n) + f(y_{n+1})\bigr)$$

It requires evaluating $f(y_n)$ first, then solving for $y_{n+1}$ with Newton just like Implicit Euler.

The `CrankNicolson` class builds the residual with both slopes and stores $f(y_n)$ as a constant:

```cpp
CrankNicolson(std::shared_ptr<NonlinearFunction> rhs) 
: TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)),
  f_yold(rhs->dimF()), m_vecf(std::make_shared<ConstantFunction>(rhs->dimF()))
{
  m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
  auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
  m_equ = ynew - m_yold - m_tau * (m_vecf + m_rhs);  // τ/2 set at runtime
}

void doStep(double tau, VectorView<double> y) override
{
  //set known values
  m_yold->set(y);
  m_tau->set(tau*0.5);
  //evaluate f(y_n)            
  this->m_rhs->evaluate(y, f_yold);
  m_vecf->set(f_yold);
  //solve for y_{n+1}
  NewtonSolver(m_equ, y);
}
```

## 2.5 Explicit Runge-Kutta for Arbitrary Butcher Tableaus

Explicit Runge–Kutta generalizes single-step methods using a Butcher tableau with stages $s$, coefficients $a_{ij}$, $b_j$, and nodes $c_j$. Each stage computes a slope $k_j$ using weighted combinations of prior stages:

$$k_j = f\left(y_n + \tau\sum_{l=1}^{j-1}a_{jl}k_l\right), \quad y_{n+1} = y_n + \tau\sum_{j=1}^{s}b_j k_j$$

This allows constructing higher-order methods. The `RungeKutta` class in `src/implicitRK.hpp` stores the tableau and iteratively builds all stage slopes before advancing:

```cpp
RungeKutta(std::shared_ptr<NonlinearFunction> rhs,
  const Matrix<> &a, const Vector<> &b, const Vector<> &c) 
: TimeStepper(rhs), m_a(a), m_b(b), m_c(c),
m_stages(c.size()), m_n(rhs->dimX()), m_k(m_stages*m_n), m_y(m_stages*m_n)
{ }   

void doStep(double tau, VectorView<double> y) override
{
  // Initialize stage values with current y
  for (int j = 0; j < m_stages; j++)
    m_y.range(j*m_n, (j+1)*m_n) = y;

  // Compute stage derivatives k_j
  m_k = 0.0;
  for (int j = 0; j < m_stages; j++)
  {
    VectorView<> ystage = m_y.range(j*m_n, (j+1)*m_n);
    for (int l = 0; l < j; l++)
      ystage += tau * m_a(j,l) * m_k.range(l*m_n, (l+1)*m_n);
    this->m_rhs->evaluate(ystage, m_k.range(j*m_n, (j+1)*m_n));
  }

  // Combine stages: y_{n+1} = y_n + tau * sum(b_j * k_j)
  for (int j = 0; j < m_stages; j++)
    y += tau * m_b(j) * m_k.range(j*m_n, (j+1)*m_n);
}
```
