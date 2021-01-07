set term epslatex color standalone size 12cm,6cm 
set out "plot_rabi.tex"


set yrange [0:1.18]
set ylabel 'occupations'
set xlabel 'time'

plot "ACE1.out" u 1:2 w l t '\texttt{-dt 0.001 -te 20 -TLS\_const\_driving 1}'

