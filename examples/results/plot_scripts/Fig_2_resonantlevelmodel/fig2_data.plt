set term epslatex standalone color size 24cm, 7.5cm
set out "fig2_data.tex"


set multiplot layout 1,3


LM=0.0525
X1=0.26
DX=0.20 #0.19
RM=0.88
TM=0.977
BM=0.125



DXLAB=-0.05
DYLAB=0
DYKEY=0.01  #DYLAB-0.05
DXKEY1=0.065
DXKEY2=DXKEY1-0.005

LW=3

set tmargin screen TM
set bmargin screen BM
set lmargin screen LM
set rmargin screen LM+X1

set label 1 '\textbf{a}' at screen LM+DXLAB,TM-DYLAB
set yrange [0:1]
set xrange [0:5]
set ytics (0, 0.2, 0.4, 0.6, 0.8, 1)
set xtics nomirror
set ylabel 'Occupation' offset screen 0.01
set xlabel '$g t$' offset screen 0,0.03
set key samplen 1.2 center at screen LM+X1+DXKEY1,(TM+BM)/2  spacing 1.2

plot  "01_resonantlevelmodel_N2_degenerate.out" every 10 w p pt 6 lw LW t '$N_E=2$' \
, (sin(sqrt(2)*x)**2) w l lw LW t '$\sin^2(\sqrt{2}gt)$'


set lmargin screen LM+X1+DX
set rmargin screen RM
set label 1 '\textbf{b}' at screen LM+X1+DX+DXLAB,TM-DYLAB
set key center at screen RM+DXKEY2,(TM+BM)/2  #TM-DYKEY 


set xrange [0:2.5]
set ylabel 'Occupation' offset screen 0.01

set xlabel '$\gamma t$' offset screen 0,0.03
plot "02_resonantlevelmodel_N4.out" every 5 w p lw LW t '$N_E=4$' \
, "03_resonantlevelmodel_N10.out" every 5 w p lw LW t '$N_E=10$' \
, "04_resonantlevelmodel_N100.out" every 5 w p pt 6 lw LW t '$N_E=100$' \
, (1-exp(-x)) lw LW t 'Markov' \
, x<0.4 ? 10/(2.*3.14159)*x*x : 1/0 lw LW lc rgb "red" t 'quadratic'


set lmargin screen RM-0.12
set rmargin screen RM-0.015
set tmargin screen BM+0.45
set bmargin screen BM+0.13

unset label 1
#set key left
unset key
set ytics 100 format '%g' offset screen 0.006,0
set xtics 50 offset screen 0, 0.01
set ylabel '$d_\textrm{max}$' offset screen 0.015,-0.01
set xlabel '$N_E$' offset screen 0,0.035
set xrange [0:108]
set yrange [0:229]

plot "maxdim_DOS.dat" w lp pt 7 lw LW lc rgb "black" t '$d_\textrm{max}$'


unset multiplot

