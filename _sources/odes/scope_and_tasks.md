# Project Scope and Tasks

The goal of this project is to implement and analyze numerical methods for ODEs
as described in the course’s JupyterBook:

- [Overview and theory](https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html)

## Corresponding JupyterBook sections

Our work corresponds to the following concrete JupyterBook sections:

- **Explicit and Improved Euler – Exercise 1**  
  Mass–spring system with explicit and improved Euler:
  - [Theory & implementation](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ee.html)  
  - [Exercises](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ee.html#exercise)

- **Implicit Euler and Crank–Nicolson – Exercise 1**  
  Implicit methods applied to the mass–spring system:
  - [Theory & implementation](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ie.html)  
  - [Exercises (mass–spring)](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ie.html#excercises)

- **Electric Network – Exercise 2**  
  Modeling an electric network (RC/RLC-type system) as an ODE:
  - [Electric network exercise (last task in block)](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ie.html#excercises)

- **Automatic Differentiation – Exercise 2**  
  AD and its application to Legendre polynomials and the pendulum:
  - [AD theory & implementation](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ad.html)  
  - [AD exercises (Legendre polynomials, pendulum)](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ad.html#exercises)  
  - [Pendulum AD exercise](https://jschoeberl.github.io/IntroSC/ODEs/implementation_ad.html#exercise-test-the-autodiff-class-for-the-pendulum)

- **Runge–Kutta methods – Exercise 2**  
  Higher-order (implicit and explicit) RK schemes:
  - [RK theory & implementation](https://jschoeberl.github.io/IntroSC/ODEs/RungeKutta.html)  
  - [RK exercises](https://jschoeberl.github.io/IntroSC/ODEs/RungeKutta.html#exercises)

## Focus of Team 10 implementation

Our Team 10 implementation focuses on:

1. **Time-Stepping for Linear ODEs (mass–spring system)**  
   Based on the explicit and implicit Euler sections:
   - Using different step counts (step sizes) to study stability and accuracy.
   - Generating phase and position–time plots.

2. **Modeling an Electric Network (RC system)**  
   Based on the electric network exercise in the implicit Euler / Crank–Nicolson chapter:
   - Implementing and solving a first-order ODE describing an electric network.
   - Comparing different steppers and producing time-domain voltage plots.

3. **Automatic Differentiation (AD)**  
   Based on the AD chapter and exercises:
   - Extending and using an `AutoDiff` class.
   - Computing Legendre polynomials and their derivatives, including visualization.

4. **Implicit Runge–Kutta Methods**  
   Based on the Runge–Kutta chapter and exercises:
   - Implementing implicit RK schemes using quadrature information.
   - Using them as time-stepping methods via a generic timestepper framework.