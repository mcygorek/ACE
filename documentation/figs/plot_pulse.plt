set term epslatex color standalone size 12cm,6cm 
set out "plot_pulse.tex"


set yrange [0:1.18]
set ylabel 'occupations'
set xlabel 'time'

plot "ACE2.out" u 1:2 w l t '\texttt{-driver driver2.param}'

