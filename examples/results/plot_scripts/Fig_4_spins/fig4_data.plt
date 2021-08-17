set term epslatex standalone color size 12cm,11.5cm
set out "fig4_data.tex"


set multiplot layout 3,1

LW=2.5
LW2=3


analyt(x)=x<0 ? 1/0 : 0.5*cos(x/2./0.6582119569)

LM=0.13
RM=0.98
TM=0.87 #0.85
YMAX1=0.6 #95
YMIN=0.6
FAC=(YMAX1+YMIN)/(2*YMIN)
H2=0.223
H=H2*FAC
DH=0.05
DXLAB=0
DYLAB=0 #0.01 
DYKEY=0.992

set lmargin screen LM
set rmargin screen RM
set tmargin screen TM
set bmargin screen TM-H

#set ytics (-0.4, -0.2, 0, 0.2, 0.4)
set ytics (-0.5, -0.25, 0, 0.25, 0.5)
set xtics format ''
set yrange [-YMIN:YMAX1]
set xrange [-4.9:20]
YLABELX=0.035 #0.025
YLABELY=0 # 0.01
set ylabel '$S_x/ \hbar$' offset screen YLABELX,YLABELY #-0.035

unset key

set label 1 '\textbf{a}' at screen DXLAB, TM-DYLAB

plot "01_polarised_N10_thr1e-10.out" u 1:2 every 30 w p lw LW pt 2 t '$N=10$'\
, "01_polarised_N100_thr1e-10.out" u 1:2 every 30 w p lw LW pt 2 lc 4 t '$N=100$'\
, "01_polarised_N1000_thr1e-10.out" u 1:2 every 30 w p lw LW pt 2 lc 2 t '$N=1000$'\
, analyt(x)  w l lw LW lc rgb "red" t 'analytic'


set yrange [-YMIN:YMIN]
#set ylabel '$S_x/ \hbar$' offset screen YLABELX,YLABELY


set tmargin screen TM-H-DH
set bmargin screen TM-H-DH-H2

set label 1 '\textbf{b}' at screen DXLAB, TM-H-DH-DYLAB
set label 2 '' at screen LM+DXLAB, TM-H-DH-DYLAB
unset key

plot "02_b20_N10_thr1e-10.out" u 1:2 every 50 w p lw LW pt 2 lc 1 t '$\epsilon=10^{-10}$' \
, "02_b20_N100_thr1e-10.out" u 1:2 every 50 w p lw LW pt 2 lc 4 t '$\epsilon=10^{-10}$' \
, "02_b20_N10_thr1e-13.out" u 1:2 every 50 w p lw LW pt 6 lc 1 t '$\epsilon=10^{-13}$' \
, "02_b20_N100_thr1e-13.out" u 1:2 every 50 w p lw LW pt 6 lc 4 t '$\epsilon=10^{-13}$' \
, "02_b20_N10_thr1e-16.out" u 1:2 every 50 w p lw LW pt 8 lc 1 t '$\epsilon=10^{-16}$' \
, "02_b20_N100_thr1e-16.out" u 1:2 every 50 w p lw LW pt 8 lc 4 t '$\epsilon=10^{-16}$'




set xtics format '%g'
set tmargin screen TM-H-DH*2-H2
set bmargin screen TM-H-DH*2-H2*2

set xlabel '$t J/\hbar $'
set label 1 '\textbf{c}' at screen DXLAB, TM-H-DH*2-H2-DYLAB
set key samplen 1.25 maxrow 3 at screen RM+0.01,DYKEY  spac 1.05



plot  "03_unpolarised_N10_thr1e-10.out" u 1:2 every 50 w p lw LW pt 2 lc 1 t '$\phantom{mmmmm}$' \
, "03_unpolarised_N100_thr1e-10.out" u 1:2 every 50 w p lw LW pt 2 lc 4 t '$\phantom{.}$' \
,-1 w p lw LW pt 2 lc 2 t '\phantom{.}' \
, "03_unpolarised_N10_thr1e-13.out" u 1:2 every 50 w p lw LW pt 6 lc 1 t '$\phantom{.}$' \
, "03_unpolarised_N100_thr1e-13.out" u 1:2 every 50 w p lw LW pt 6 lc 4 t '$\phantom{.}$' \
, -1 lc rgb "white" t '\phantom{.}' \
, "03_unpolarised_N10_thr1e-16.out" u 1:2 every 50 w p lw LW pt 8 lc 1 t '$\phantom{.}$' \
, "03_unpolarised_N100_thr1e-16.out" u 1:2 every 50 w p lw LW pt 8 lc 4 t '$\phantom{.}$' \
, -1  w l lw LW lc rgb "red" t '\phantom{.}'



unset multiplot


