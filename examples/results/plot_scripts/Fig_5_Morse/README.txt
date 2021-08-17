
To reproduce Figure 5:
-------------------------

To acquire data:
Run the example "04_morse" and copy the *.out file to this directory. 

The Morse potential eigenstates are calculated by running (~ 5 min):
ACE printEV.param

To produce plot, run:

gnuplot fig5_resub.plt && pdflatex fig5_resub.tex

Output:

fig5_resub.pdf

