
To reproduce Figure S.3.1:
-------------------------

To acquire data, run:

./script_ACE.sh && ./script_Gaussian.sh &&  ./script_TEMPO.sh && M=3 ./script_ACE.sh && ./script_TEMPO_tmem.sh 

To process data, run:

./analyze_ACE.sh && ./analyze_Gaussian.sh && ./analyze_TEMPO_tmem.sh && ./analyze_ACE_M3.sh && ./analyze_TEMPO_tmem.sh 

_manually_ check with the "max_err.sh" script that the TEMPO results with tmem=2.5 lead to relative errors >1% when \Delta E < 1.5 meV.
Then, split up the results "analyze_TEMPO_tmem2.5.dat" into two files name
"analyze_TEMPO_tmem2.5_inaccurate.dat" and "analyze_TEMPO_tmem2.5_accurate.dat", where the former contains the entries \Delta E smaller or equal 1.5 meV and the latter contains the rest.

To generate the data for the spectral density, run

ACE print_J_dE0.5.param 
and
ACE print_J_dE10.param


To produce plot, run:

gnuplot figS3_1.plt  && pdflatex figS3_1.tex

Output:

figS3_1.pdf



