###############################
# CTEMPO and CCTEMPO examples #
###############################

This directory contains scripts to reproduce the examples in the manuscript on the CTEMPO and CCTEMPO algorithms.


-------------
Prerequisites
-------------

These algorithm are implemented within the C++-Code ACE. Python bindings for CTEMPO can be found in the example pybind/examples/08_QUAPI.py. Python bindings for CCTEMPO are still under construction.

To obtain the C++ binaries, download the ACE code using 

> git clone https://github.com/mcygorek/ACE

and compile the binaries using 

> make ACE TEMPO CTEMPO CCTEMPO

Binaries with the corresponding name should now be found in the ACE/bin/ subdirectory. For the remainder, it is assumed that this subdirectory is added to the PATH, or that the names of the binaries are replaced by their name including their absolute path.


----------------
Spin-Boson model
----------------

The subdirectory "spinboson" contains exemplary parameter files that can be used to simulate the spin-boson dynamics as indicated in the manuscript for coupling strength alpha=0.1 (or 1), total propagation time te=20.48, time step dt=0.004, and compression threshold epsilon=10^{-8}. Please change these values according to your need.
The same parameter file can be executed with different binaries to employ the respective method. We suggest to specify the name of the output file to indicate hthe method used, for example:

> CTEMPO ohmic_alpha0.1_te20.48_dt0.004_thr1e-8.param -outfile ohmic_alpha0.1_te20.48_dt0.004_thr1e-8_CTEMPO.out
> TEMPO ohmic_alpha0.1_te20.48_dt0.004_thr1e-8.param -outfile ohmic_alpha0.1_te20.48_dt0.004_thr1e-8_TEMPO.out
> ACE ohmic_alpha0.1_te20.48_dt0.004_thr1e-8.param -outfile ohmic_alpha0.1_te20.48_dt0.004_thr1e-8_PTMPO.out

"CTEMPO" and "TEMPO" implement the algorithm with the same name. The binary "ACE" is the general binary for PT-MPO simulations. The Jorgensen-Pollock algorithm is selected by the line "use_Gaussian true" in the parameter file.

The corresponding output file ("*.out") contains the time in the first column and the value <sigma_z> in the second column.



-----------
FMO complex
-----------

The subdirectory "FMO" contains exemplary parameter files for the Fenna-Matthews-Olson complex. To execute them, call

> CCTEMPO T77_te1_dt0.001_thr1e-5.param 

which calculates the dynamics for the situation when the first site (index 0) is initially excited while all other sites are in the ground state. The Hamiltoninan is described by on-site and two-site terms. 
The time evolution is then given in the output file, where the first column contains the time and column nr. 2,4,6,8,10,12,14 contain the occupations of site 0,1,2,3,4,5,6, respectively.

Moreover, the example
> CCTEMPO E0_T77_te1_dt0.001_thr1e-5.param 
calculates the dynamics, where the initial state is the 0-th eigenstate of the system Hamiltonian in the single-excitation manifold. Note that the time evolution of the eigenstates will be written in "E0_T77_te1_dt0.001_thr1e-5.eigs" 

To obtain convergence, the compression threshold should be decreased to about 10^{-6}


--------------------
Superradiance of QDs
--------------------

The subdirectory "superradiance" contains exemplary parameter files for the superradiant emission from several semiconductor quantum dots (QDs).

For example, executing
> CCTEMPO N8_T4_Tdecay500_te1000_dt0.1_thr1e-7.param
produces an output file "N8_T4_Tdecay500_te1000_dt0.1_thr1e-7.out", whose second column contains the occupation of the first of 8 superradiantly emitting QDs when each QD is coupled to a super-Ohmic phonon environment (all emitters are identical and so the total occupation is N times the second column of the output file).

A Bash script for creation of parameters files is also provided, and can be called like, e.g.,

> N=8 T=4 Tdecay=500 te=1000 dt=0.1 thr=1e-7 ./superrad_generate_param.sh 

N is the numer of QDs, T is the temperature in Kelvin (a value of -1 means phonon coupling is switched off), Tdecay is the life time that a single emitter would have in isolation, te is the total propagation time (in ps), dt is the time step (also in ps), and thr is the compression threshold.

Note that simulations with large N and small thresholds require sizable memory on the computer. 



