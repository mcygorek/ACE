
# ACE: 
Code for numerically-exact simulations of open quantum systems. 

Original author: Moritz Cygorek

Contributors: Thomas Bracht

## Scope and Methods

Its name derives from the Automated Compression of Environments (ACE) algorithm [1] based on real-time path integrals and tensor networks. Over the years a number of other methods have been implemented so that it grew to a quite versatile and general open quantum systems simulator.
These methods include:

| | Authors | Title | Journal |
|-|---------|-------|---------|
|[1]| M. Cygorek, M. Cosacchi, A. Vagov, V.M. Axt, B.W. Lovett, J.Keeling & E.M. Gauger| *Simulation of open quantum systems by automated compression of arbitrary environments*| [Nat. Phys. 18, 662–668 (2022).](https://doi.org/10.1038/s41567-022-01544-9)|
|[2]| M.R. Jørgensen, F.A. Pollock| *Exploiting the Causal Tensor Network Structure of Quantum Processes to Efficiently Simulate Non-Markovian Path Integrals* | [Phys. Rev. Lett. 123, 240602 (2019).](https://doi.org/10.1103/PhysRevLett.123.240602)|
|[3]| M. Cygorek, J. Keeling, B.W. Lovett, E.M. Gauger| *Sublinear Scaling in Non-Markovian Open Quantum Systems Simulations*| [Phys. Rev. X 14, 011010 (2024).](https://doi.org/10.1103/PhysRevX.14.011010)| 
|[4]| M. Cygorek, B.W. Lovett, J. Keeling, E.M. Gauger| *Treelike process tensor contraction for automated compression of environments*|[Phys. Rev. Research 6, 043203 (2024).](https://doi.org/10.1103/PhysRevResearch.6.043203)|
|[5]| V. Link, H.-H. Tu, W.T. Strunz | *Open Quantum System Dynamics from Infinite Tensor Network Contraction*| [Phys. Rev. Lett. 132, 200403 (2024)](https://doi.org/10.1103/PhysRevLett.132.200403)|

The code is written in C/C++ but offers extensive Python bindings. Usage of the C/C++-Code based on Parameter files is described in:

| | Authors | Title | Journal |
|-|---------|-------|---------|
|[6]| M. Cygorek, E.M. Gauger| *ACE: A general-purpose non-Markovian open quantum systems simulation toolkit based on process tensors* | [J. Chem. Phys. 161, 074111 (2024).](https://doi.org/10.1063/5.0221182)|

More about the concept of process tensors can be found in 

| | Authors | Title | Journal |
|-|---------|-------|---------|
|[7]| M. Cygorek, E.M. Gauger| *Understanding and utilizing the inner bonds of process tensors* | [SciPost Phys. 18, 024 (2025).](https://doi.org/10.21468/SciPostPhys.18.1.024)|


# Installation:

The ACE code can be used in two ways: 

- The C/C++ binaries simulate open quantum systems based on parameter files.
- Python bindings access the same internal C/C++ functions but allow easy integration into common Python-based data processing frameworks

The binaries can be directly downloaded from the releases. Simulation are then started from the command line by providing a parameter file as, e.g., 

`ACE.exe example.param`

Note that the binary must be chosen according to the operating system, and it must be placed in a PATH that can be resolved by the terminal.


The standalone Python package can be straightforwardly installed via

`pip install ACE-OQS`

Once installed, methods can be called after `import ACE` in Python scripts (see examples in the pybind/examples/ subdirectory).

It should be noted that these binaries are not linked against BLAS/LAPACK functions, e.g., provided by the Intel MKL. Linking against BLAS/LAPACK may lead to a sizable speed-up. To this end we recommend to download the source using `git clone https://github.com/mcygorek/ACE` and compile it. This process is described for different operating systems below.

## Compilation from source in Linux: 

Go to ACE directory and type "make" 

The Eigen library has to available. The Makefile tries "/usr/include/eigen3/". If it's located elsewhere, please set the environment variable "EIGEN_HOME" during compilation.

The Intel MKL is automatically used if set up correctly on the system.

Successful compilation produces a binary "ACE" in the "bin" subdirectory

If "pybind11" is installed (if not, run `python3 -m pip install build pybind11` first), Python bindings can be generated with "make pybind". This creates a .whl file in the pybind/ subdirectory that can be installed using `pip install pybind/*.whl`.

## Compilation from source in Windows: 

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
