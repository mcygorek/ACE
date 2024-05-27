
This example shows how to calculate superradiant emission from two quantum dots, each coupled to a local phonon bath:

We use the reverse rotating wave approximation (RWA) as outline in Phys. Rev. X 14, 011010 (2024). We also work in the rotated basis where the Rabi model becomes the spin-boson model (coupling via sigma_z Operators) to keep the memory and disk space footprint smaller. This is why, e.g., occupation operators |1><1|_2 are mapped to 0.5*(Id_2+sigma_x), which is the 2x2 matrix with all 0.5 entries.

To calculate the phonon and photon PT-MPOs, run

ACE generate_phonon_PT.param  

and

ACE generate_photon_PT.param  

After this is done, run 

ACE propagate.param 

which takes the photon PT-MPO and two copies of the phonon PT-MPO for the calculation. Note that the propagation can take lots of memory, both in RAM (please reserve 8BG) and on hard disk (please reserve 20GB). 

Once this is done, one can take the derivative of the occupations using the tool:

differentiate_outfile -infile propagate.out -col 2 > propagate.diff

and plot the results, e.g., using gnuplot:

plot "propagate.diff" u 1:(-$4*64) w l, 2*(1+2*x/64)*exp(-2*x/64), 2*exp(-x/64) 


