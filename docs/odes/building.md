# Building the Project

We briefly summarize how to build and run the Team 10 implementation.

Implementation repository:

- <https://github.com/Stefmeff/ASC-ODE-team-10>

## Requirements

- C++ compiler with C++20 support (e.g. `g++` or `clang++`),
- CMake (≥ 3.x),
- Python 3 with `numpy` and `matplotlib` (for plotting).

## Clone and Configure

```bash
git clone https://github.com/Stefmeff/ASC-ODE-team-10.git
cd ASC-ODE-team-10

# Initialize submodules (e.g. nanoblas, if present)
git submodule update --init --recursive 
```
## Building
```bash
mkdir build
cd build
cmake ..
make 
```
