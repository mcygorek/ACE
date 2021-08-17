set term epslatex standalone color size 12cm,5cm
set out "figS5_1_data.tex"


set tmargin screen 0.98
set bmargin screen 0.22
set lmargin screen 0.11
set rmargin screen 0.98

LW=2.5
PS=1.1
set xlabel '$\kappa t$' offset screen 0,0.02
set ylabel 'Total occupation' offset screen 0.015
set xrange [0:2.5]
set yrange [0:2.2]

set key samp 1.3

plot "02_delta0.out" u 1:($2+$4) every 50 w p ps PS lw LW t '$\delta=0$' \
, "03_delta1.out" u 1:($2+$4) every 50 w p ps PS lw LW t '$\delta=\kappa$' \
, "04_delta10.out" u 1:($2+$4) every 50 w p ps PS pt 6 lw LW t '$\delta=10\,\kappa$' \
, 2*exp(-x) lc rgb "gray" lw LW  t 'independent emission' \
, 2*(1+1*x)*exp(-2*x) lc rgb "red" lw LW t 'coherent emission'



