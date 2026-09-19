# RootForge

**RootForge** is a C++ library of classical numerical and mathematical
algorithms, built to be dropped into other software as a lightweight,
dependency-free computation engine. It also ships with an optional
**ImGui-based desktop GUI** that lets you explore and run every function in
the library interactively — no need to write test code to see how a method
behaves.

The project is organized as a set of small, focused modules under the
`RootForge` namespace: general-purpose math helpers, numerical calculus
(differentiation, integration, ODE solving), root-finding, interpolation,
and a full symbolic-ish `Polynomial` class with FFT-based multiplication.

---

## Table of Contents

- [Project Idea](#project-idea)
- [Features at a Glance](#features-at-a-glance)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Building](#building)
- [The GUI: RootForge Library Explorer](#the-gui-rootforge-library-explorer)
- [Library Reference](#library-reference)
  - [Core Types (`RootForgeCore.h`)](#core-types-rootforgecoreh)
  - [Mathematical Functions](#mathematical-functions)
  - [Differentiation](#differentiation)
  - [Integration](#integration)
  - [Root Finding](#root-finding)
  - [Interpolation](#interpolation)
  - [ODE Solvers](#ode-solvers)
  - [Polynomial](#polynomial)
- [Using RootForge as a Library in Your Own Project](#using-rootforge-as-a-library-in-your-own-project)
- [Known Limitations](#known-limitations)
- [Roadmap Ideas](#roadmap-ideas)
- [License](#license)

---

## Project Idea

Most small-to-medium C++ projects that need numerical methods end up either
pulling in a heavyweight dependency (Eigen, GSL, Boost.Math) or reinventing
the same handful of algorithms badly. RootForge aims to sit in between:

- **Self-contained** — no external math dependency, just the STL.
- **Educational and inspectable** — every algorithm is implemented from
  first principles (Newton-Raphson, Simpson's rule, RK4, Lagrange
  interpolation, etc.) so it's easy to read, learn from, and modify.
- **Reusable** — organized into clean headers/namespaces so any subset can
  be copied into another project.
- **Interactive** — the included GUI turns the whole library into a
  point-and-click calculator, useful for testing, demos, teaching, or just
  sanity-checking a method before wiring it into real code.

## Features at a Glance

| Area | What's inside |
|---|---|
| **Mathematical Functions** | sign, fractional part, custom `power`, primality/sieve, factorial, roots of unity, Gamma function (exact + Spouge approximation), log-Gamma, binomial coefficients, `erf`/`erfc`, Taylor-series `exp`/`sin`/`cos`/`ln`, simplified square roots |
| **Differentiation** | forward, backward, centered, second-order centered, and Richardson-extrapolated derivatives |
| **Integration** | trapezoidal rule, midpoint rule, Simpson's 1/3 rule, Richardson-extrapolated Simpson, 3-point Gaussian quadrature |
| **Root Finding** | Newton's method for nth roots, Bisection, Newton-Raphson, Secant, Regula Falsi |
| **Interpolation** | linear interpolation, numerical (Lagrange-form) polynomial interpolation, explicit Lagrange polynomial construction |
| **ODE Solvers** | Euler's method, RK4, explicit Midpoint method, RK2 (Heun's), and a vector/system RK4 solver for coupled first-order ODEs |
| **Polynomial** | evaluation (real & complex), addition, subtraction, scalar multiplication, long division, FFT-based multiplication, differentiation, integration, numerical derivative/integral, deflation-based root finding (Newton-Horner), and the quadratic formula |
| **GUI** | A full ImGui desktop app — "RootForge Library Explorer" — that exposes nearly every function above through menus and input fields |

## Project Structure

```
RootForge/
├── CMakeLists.txt              # Build configuration (fetches GLFW, links OpenGL + ImGui)
├── RootForge.h                 # Umbrella header (includes every module)
├── RootForgeCore.h             # Shared types: Function, DFunction, NFunction, Status, NumericalResult
├── MathematicalFunctions.h/.cpp# General-purpose math utilities
├── NumericalAnalysis.h/.cpp    # Differentiation, Integration, Root, Interpolation, ODE namespaces
├── Polynomial.h/.cpp           # Polynomial class
├── NumberTheory.h              # (referenced by RootForge.h — add your own implementation)
├── GUI_Test.cpp                # Entry point for the ImGui "Library Explorer" application
└── imgui/                      # Dear ImGui sources + GLFW/OpenGL3 backends (not vendored here — see Requirements)
```

> **Note:** `RootForge.h` includes `NumberTheory.h`, which isn't part of
> this snapshot of the project. If you're only using the computation
> modules (not the umbrella header), include the specific headers you need
> instead (`RootForgeCore.h`, `MathematicalFunctions.h`,
> `NumericalAnalysis.h`, `Polynomial.h`) to avoid that dependency.

## Requirements

To build the **GUI application**:

- CMake 3.14+
- A C++14-capable compiler
- OpenGL (provided by your OS)
- [GLFW](https://www.glfw.org/) — fetched automatically by `CMakeLists.txt` via `FetchContent`
- [Dear ImGui](https://github.com/ocornut/imgui) — **not fetched automatically**; place the `imgui` folder (including `backends/imgui_impl_glfw.cpp` and `backends/imgui_impl_opengl3.cpp`) at the project root before building

To use the **math library alone** (no GUI): just a C++14 compiler and the
STL — no third-party dependencies at all.

## Building

```bash
# 1. Clone Dear ImGui into the project root (one-time setup)
git clone https://github.com/ocornut/imgui.git

# 2. Configure and build
mkdir build && cd build
cmake ..
cmake --build .

# 3. Run
./RootForge
```

CMake will fetch and build GLFW automatically the first time you configure.

## The GUI: RootForge Library Explorer

`GUI_Test.cpp` builds a single-window desktop app with a simple, repeatable
workflow:

1. **Pick a Library** from the first dropdown — Mathematical Functions,
   Differentiation, Integration, Root Finding, Interpolation, ODE Solvers,
   or Polynomial.
2. **Pick a Function** from the second dropdown — options update
   automatically based on the chosen library.
3. **Fill in the inputs** that appear for that specific function (numbers,
   or comma-separated lists for vector inputs like `1, 2, 3`).
4. Press **Calculate** — the result appears immediately below, and is also
   appended to a scrollable **History** panel so you can compare runs.
5. Repeat as many times as you like, switching libraries/functions freely.

**On functions that need a callable `f(x)` or `dy/dx = f(x, y)`:** the
library itself has no expression parser, so instead of free-text math
input, the GUI offers a small dropdown of predefined test functions for
Differentiation, Integration, and Root Finding (e.g. `x² - 4`, `sin(x)`,
`eˣ - 2`, `x³ - x - 2`), and predefined `dy/dx` forms for the ODE solvers
(e.g. `dy/dx = y`, `dy/dx = x + y`). `SystemRK4` runs a fixed illustrative
2-equation system (a simple harmonic oscillator: `y₀' = y₁`, `y₁' = -y₀`).

Some functions (`taylorExp`, `taylorSin`, `taylorCos`, `taylorLn`,
`NewtonMethod`) also print a detailed iteration/error report to the
console via `NumericalResult::ShowResult`, in addition to the value shown
in the GUI.

## Library Reference

### Core Types (`RootForgeCore.h`)

```cpp
using Function  = std::function<double(double)>;                 // f(x)
using DFunction = std::function<double(double, double)>;         // f(x, y)
using NFunction = std::function<double(double, vector<double>)>; // f(x, y-vector)

enum class Status { Success, InvalidInput, DivisionByZero, MaxIterations, NotConverged };

struct NumericalResult {
    double value, absError, relativeError, percentError;
    std::size_t iterations;
    Status status;
    static void ShowResult(const NumericalResult result); // prints a formatted report
};
```

### Mathematical Functions

Namespace: `RootForge::MathematicalFunctions`

| Function | Signature | Notes |
|---|---|---|
| `sgn` | `int sgn(double x)` | -1, 0, or 1 |
| `FractionalPart` | `double FractionalPart(double x)` | `x - floor(x)` |
| `power` | `double power(double a, double n)` | Fast exponentiation by squaring for integer `n`; falls back to `exp(n·ln(a))` for fractional `n` |
| `is_integer` | `bool is_integer(double num)` | |
| `sieve_of_eratosthenes` | `vector<int> sieve_of_eratosthenes(int n)` | Returns all primes ≤ n |
| `factorial` | `long double factorial(long n)` | Computed via prime-factorization for large-n stability |
| `nthRoot` | `vector<complex<double>> nthRoot(int n, bool invert = false)` | The n-th roots of unity |
| `Gamma` | `double Gamma(double n)` | Exact for integers/half-integers, `NAN` otherwise |
| `SpougeGamma` | `double SpougeGamma(double n)` | Spouge's approximation, works for general real `n` |
| `LogGamma` | `double LogGamma(double n)` | `log(Gamma(n))` |
| `BinomialCoeff` | `long BinomialCoeff(long n, long k)` | |
| `erf` / `erfc` | `double erf(double x)` / `double erfc(double x)` | Numerically integrated via Simpson's rule |
| `taylorExp` / `taylorSin` / `taylorCos` / `taylorLn` | `double taylorX(double a, double x)` | Taylor-series evaluation with a built-in error report |
| `simplify_sqrt` | `string simplify_sqrt(int n)` | Simplifies `√n` into `a√b` form |
| `isPrime` | `bool isPrime(int n)` | |

### Differentiation

Namespace: `RootForge::Differentiation`

- `forwardDifference(f, x, h)`
- `backwardDifference(f, x, h)`
- `centeredDifference(f, x, h)`
- `D2centeredDifference(f, x, h)` — second derivative
- `RichardsonExtrapolationDiff(f, x, h)` — 4th-order accurate via Richardson extrapolation of the centered difference

### Integration

Namespace: `RootForge::Integration`

- `trapezoidal(f, start, end, n)`
- `midpoint(f, start, end, n)`
- `simpsons_13(f, start, end, n)` — requires an even `n`
- `RichardsonExtrapolationInt(f, start, end, n)` — Richardson-extrapolated Simpson's rule
- `GaussianQuadrature(f, a, b)` — 3-point Gauss-Legendre quadrature

### Root Finding

Namespace: `RootForge::Root`

- `NewtonMethod(num, n)` — computes the n-th root of `num`
- `Bisection(f, start, end)` — requires a sign change over `[start, end]`
- `NewtonRaphson(f)` — starts from `x = 1`, uses the centered difference as the derivative
- `Secant(f)` — automatically searches `[-100, 100]` for a bracketing sign change
- `RegulaFalsi(f)` — same automatic bracketing as `Secant`

### Interpolation

Namespace: `RootForge::Interpolation`

- `LinearInterpolation(x1, y1, x2, y2, x_i)`
- `NumericalPolyInterpolation(x_axis, y_axis, x)` — evaluates the Lagrange interpolant at a point without building it explicitly
- `LagrangePolyInterpolation(x_axis, y_axis)` — returns the explicit interpolating polynomial's coefficients

### ODE Solvers

Namespace: `RootForge::ODE`

- `EulerMethod(f, x0, y0, h, target)`
- `RK4(f, x0, y0, h, target)` — classic 4th-order Runge-Kutta
- `Midpoint(f, x0, y0, h, target)` — explicit midpoint method
- `RK2(f, x0, y0, h, target)` — Heun's method
- `SystemRK4(functions, x0, y, h, target)` — RK4 for a system of coupled first-order ODEs, where `y` is the initial state vector and `functions[i]` computes `dy_i/dx`

### Polynomial

Class: `Polynomial` (top-level, not inside `RootForge`)

Coefficients are stored **low-to-high degree**: `coeffs[0]` is the constant
term, `coeffs[coeffs.size()-1]` is the leading coefficient.

| Method | What it does |
|---|---|
| `Polynomial(coefficients)` | Constructs from a coefficient vector; also builds a callable `RootForge::Function` |
| `degree()` | `coeffs.size() - 1` |
| `evaluate(x)` / `evaluate(complex<double> x)` | Horner's method evaluation, real or complex |
| `add(p1, p2)` | Static; polynomial addition |
| `negate(p1, p2)` | Static; computes `p1 - p2` |
| `ScalarProd(p, x)` | Static; scales all coefficients by `x` |
| `divide(p1, p2)` | Static; returns `{quotient, remainder}` |
| `multiply(p1, p2)` | Static; **FFT-based** multiplication using `nthRoot` |
| `NewtonHorner(p)` | Static; finds all real roots by repeated Newton-Raphson + deflation |
| `Differentiate(p)` | Static; term-by-term derivative |
| `Integrate(p)` | Static; term-by-term antiderivative (constant of integration = 0) |
| `NumericalDiff(x)` | Centered-difference derivative at `x` |
| `NumericalInt(x_start, x_end)` | Simpson's-rule definite integral |
| `quadratic_formula(a, b, c)` | Static; prints the real or complex roots of `ax² + bx + c` to the console |

## Using RootForge as a Library in Your Own Project

You don't need the GUI to use RootForge — the computational modules have
no dependency on ImGui/GLFW/OpenGL at all. To embed it in another project:

```cpp
#include "RootForgeCore.h"
#include "MathematicalFunctions.h"
#include "NumericalAnalysis.h"
#include "Polynomial.h"

int main() {
    // Root-find sin(x) = 0 near x = 3
    double root = RootForge::Root::Bisection([](double x){ return sin(x); }, 3.0, 4.0);

    // Integrate x^2 from 0 to 1
    double area = RootForge::Integration::simpsons_13(
        [](double x){ return x * x; }, 0.0, 1.0, 100);

    // Work with a polynomial: 3x^2 + 2x - 5
    Polynomial p({-5, 2, 3});
    double y = p.evaluate(2.0);
}
```

Add the corresponding `.cpp` files to your build (`MathematicalFunctions.cpp`,
`NumericalAnalysis.cpp`, `Polynomial.cpp`) and you're set — no GUI-related
sources or dependencies required.

## Known Limitations

- `NumberTheory.h`, referenced by the umbrella `RootForge.h` and by
  `CMakeLists.txt`, is not implemented in this snapshot — include the
  specific headers you need instead of `RootForge.h` if you hit a missing
  header error.
- `MathematicalFunctions::FFT` is declared in the header but has no
  implementation yet.
- Several root-finding and ODE routines (`Secant`, `RegulaFalsi`,
  `NewtonRaphson`) don't take convergence/iteration caps, so a poorly
  chosen test function can loop for a long time or diverge.
- `Gamma()` only returns exact results for integer and half-integer
  arguments; use `SpougeGamma()` for a general real-valued approximation.
- The GUI's differentiation/integration/root-finding tools use a fixed set
  of predefined test functions rather than a free-form expression parser.

## Roadmap Ideas

- Implement `NumberTheory.h` and `FFT` to complete the umbrella header.
- Add a lightweight expression parser so the GUI can accept arbitrary
  `f(x)` input instead of a preset list.
- Add convergence/iteration limits and `Status`/`NumericalResult` returns
  consistently across all root-finding and ODE methods.
- Add unit tests (e.g. with Catch2 or GoogleTest) for each module.

## License

Add your preferred license here (e.g. MIT, Apache 2.0) before publishing
the repository publicly.
