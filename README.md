
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

> sudo apt install g++ make libeigen3-dev

> sudo apt install python3-dev python3-matplotlib

Then, clone the github repository and compile (including python bindings)

> cd

> git clone https://github.com/mcygorek/ACE

> cd ACE

> make

Now, there should be a binary file "ACE" in the "bin" subdirectory. Ideally, set the path to this directory by appending to your local .bashrc
> export PATH=$PATH:/PATH_TO/ACE/bin
where PATH_TO is to be replaced by the path to the ACE directory, typically /home/USERNAME/ with the name of the user USERNAME

Try running examples, e.g.
> ACE ACE/examples/JCP2024/Fig5/phonon_assisted.param

#### Python bindings

The ACE library provides Python bindings via pybind11. To install this, run 

> sudo apt install python3-pybind11

and in the ACE directory

> make pybind

You can check out the python examples in the ACE/pybind subdirectory, e.g., 
> python3 /PATH_TO/ACE/pybind/06_cQED.py
but be sure to edit the second line and replace the /.../ by /PATH_TO/

#### Jupyter notebooks

If you want to use .ipynb notebooks, you can install jupyter on WSL, then access the notebooks via your regular browser (outside of WSL):

> sudo apt install python3-notebook

If you run 

> jupyter-notebook

you will see a link (starting with "localhost:8888/"). Copy this long URL into your browser and go there. To try out the examples in the ACE/pybind subdirectory, make sure to copy the content of the .py files into newly created .ipynb notebook files. Note: the ".../ACE/pybind" in the second line in the examples should be replaced by "/home/WSLUSER/ACE/pybind", where WSLUSER is the username provided when installing WSL.

#### VS Code with WSL extension

There is an extension to VS Code that lets you access WSL. Just make sure you have 

> sudo apt install python3-ipython

installed on WSL and you should be able to access .ipynb notebooks on WSL using VS Code.

### MSYS2:

A precompiled ACE.exe is included for easy usage, which however comes with limited functionality and suboptimal performance, and it may be outdated.
A new Windows executable can be created with MSYS2:

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
