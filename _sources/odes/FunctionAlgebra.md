# 1. Function Algebra

In order to simulate the behaviour of various systems and solve their underlying ODE's numerically in C++, we first need to define a representation of nonlinear functions, that our programs can efficiently work with, in order apply the various time-stepping techniques for solving ODE's discussed in the subsequent chapter.

## 1.1 Representation of Functions

We introduce the `NonlinearFunction` abstract base class in `src/nonlinfunc.hpp`, which defines an interface for evaluating functions and their Jacobians. Every function stores its domain and codomain dimensions and must implement evaluation and differentiation:

$$f: \mathbb{R}^n \to \mathbb{R}^m$$

```cpp
class NonlinearFunction
{
public:
  virtual ~NonlinearFunction() = default;
  virtual size_t dimX() const = 0;
  virtual size_t dimF() const = 0;
  virtual void evaluate (VectorView<double> x, VectorView<double> f) const = 0;
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
};
```

Concrete implementations include `IdentityFunction`, `ConstantFunction` (stores a fixed vector), and function composition/arithmetic built via operators. This design enables Newton's method and ODE solvers to work uniformly with any function expression.

## 1.2 Differentiation via AutoDiff

The `AutoDiff` template in `src/autodiff.hpp` tracks both the value and all directional derivatives (up to $N$ directions) of an expression:

```cpp
template <size_t N, typename T = double>
class AutoDiff
{
private:
  T m_val;
  std::array<T, N> m_deriv;
public: 
  AutoDiff () : m_val(0), m_deriv{} {}
  AutoDiff (T v) : m_val(v), m_deriv{} { }
  
  T value() const { return m_val; }
  std::array<T, N>& deriv() { return m_deriv; }
  const std::array<T, N>& deriv() const { return m_deriv; }
};
```

The basic rules for differentiation are implemented via operators (`+`, `*`, `-`, `/`) and elementary functions (`sin`, `cos`, `exp`, `log`, `pow`, `sqrt`), allowing an easy composition of functions. For instance, multiplication applies the product rule for the derivative of the result:

```cpp
template <size_t N, typename T = double>
AutoDiff<N, T> operator* (const AutoDiff<N, T>& a, const AutoDiff<N, T>& b)
{
  AutoDiff<N, T> result(a.value() * b.value());
  for (size_t i = 0; i < N; i++)
    result.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
  return result;
}
```
The same goes concept applies for using elementary functions, e.g. sine:
```cpp
template <size_t N, typename T = double>
AutoDiff<N, T> cos(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> result(cos(a.value()));
    for (size_t i = 0; i < N; i++)
        result.deriv()[i] = -sin(a.value()) * a.deriv()[i];
    return result;
}
```
This avoids hand-coding Jacobians and enables rapid prototyping of complex models.