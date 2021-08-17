
To reproduce Figure S.2.1:
-------------------------

To acquire data, run (Takes 1-2 days!):


./script_te2.5_dt0.1.sh && ./script_te2.5_dt0.05.sh && ./script_te2.5_dt0.01.sh  && ./script_te2.5_dt0.005.sh  &&  ./script_te2.5_dt0.001.sh  && ./script_Trotter.sh 


To process data, run:

./analyze_dt0.1.sh && ./analyze_dt0.05.sh && ./analyze_dt0.01.sh && ./analyze_dt0.005.sh && ./analyze_dt0.001.sh  && ./Trotter_analyze.sh


To produce plot, run:

gnuplot figS2_1.plt && pdflatex figS2_1.tex

Output:
figS2_1.pdf


