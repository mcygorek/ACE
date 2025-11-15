
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

If "pybind11" is installed, python bindings can be generated with "make pybind". See examples in pybind subdirectory.

## Windows: 

### WSL2:

First, install WSL2 with Ubuntu and update the repository. To this end, open the PowerShell and type:

> wsl --install Ubuntu

> sudo apt-get update

> sudo apt-get upgrade

Next, install C++ and Python compilers as well as required libraries

> sudo apt install g++ autotools libeigen3-dev

> sudo apt install python3-dev python3-matplotlib python3-pybind11

Then, clone the github repository and compile (including python bindings)

> git clone https://github.com/mcygorek/ACE

> cd ACE

> make

> make pybind

Now, there should be a binary file "ACE" in the "bin" subdirectory. Ideally, set the path to this directory by appending to your local .bashrc
> export PATH=$PATH:/PATH_TO/ACE/bin
where PATH_TO is to be replaced by the path to the ACE directory, typically /home/USERNAME/ with the name of the user USERNAME

Try running examples, e.g.
> ACE ACE/examples/JCP2024/Fig5/phonon_assisted.param

or the python examples in the pybind subdirectories, e.g., 
> python3 /PATH_TO/ACE/pybind/06_cQED.py
but be sure to edit the second line and replace the /.../ by /PATH_TO/


### MSYS2:

A precompiled ACE.exe is included for easy usage, which however comes with limited functionality and suboptimal performance, and it may be outdated.
It has been tested to run under Windows with MSYS2:

First, install and run MSYS2. Then, in the UCRT version of MSYS2, run
> pacman -S mingw-w64-ucrt-x86_64-toolchain

> pacman -S mingw-w64-ucrt-x86_64-eigen3

> pacman -S git

> pacman -S make

> git clone https://github.com/mcygorek/ACE
 
> export EIGEN_HOME=/ucrt64/include/eigen3

> cd ACE; make -f Makefile_static

Successful compilation produces a binary "ACE.exe" in the "bin" subdirectory

Python bindings on MSYS2 have not been tested yet.
