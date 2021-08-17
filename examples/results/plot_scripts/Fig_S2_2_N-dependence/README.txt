
To reproduce Figure S.2.2:
-------------------------

To acquire data, run:

./script_te2.5_dt0.05_thr1e-6.sh && ./script_te5_dt0.05_thr1e-6.sh


To process data, run:

./analyze_te2.5.sh &&  ./analyze_te5.sh    

To produce plot, run:

gnuplot figS2_2.plt && pdflatex figS2_2.tex

Output:
figS2_2.pdf


