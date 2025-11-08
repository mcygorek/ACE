
# ACE: 
Code for numerically-exact simulations of open quantum systems
using the Automated Compression of Environments (ACE) method

Original author: Moritz Cygorek

# Installation:

## Linux: 

Go to ACE directory and type "make" 

The Eigen library has to available. The Makefile tries "/usr/include/eigen3/". If it's located elsewhere, please set the environment variable "EIGEN_HOME" during compilation.

The Intel MKL is automatically used if set up correctly on the system.

Successful compilation produces a binary "ACE" in the "bin" subdirectory

## Windows: 

A precompiled ACE.exe is included for easy usage, which however comes with limited functionality and suboptimal performance, and it may be outdated.
It has been tested to run under Windows with MSYS2:

- install and run MSYS2
- pacman -S mingw-w64-ucrt-x86_64-toolchain
- pacman -S mingw-w64-ucrt-x86_64-eigen3
- pacman -S mingw-w64-ucrt-x86_64-pybind11
- pacman -S git
- git clone https://github.com/mcygorek/ACE
- export EIGEN_HOME=/ucrt64/include/eigen3
- cd ACE; make -f Makefile_static

Successful compilation produces a binary "ACE.exe" in the "bin" subdirectory
