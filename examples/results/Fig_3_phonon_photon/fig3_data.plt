set term epslatex color standalone size 18cm,8.7cm
set out 'fig3_data.tex'


set multiplot layout 2,1

LW=4
LW2=1.5*LW

DXLAB=-0.064
DYLAB=0 #0.0375

LM=0.08 #0.065
TM=0.985
BM=0.125
X1=0.6
X2=0.475
DY=0.17
Y=(TM-BM-DY)/2

DXKEY1=0.28 #0.23
#DXKEY2=0.15

DBUFF=0.05+0.06
DBUFF2=DBUFF-0.05

set lmargin screen LM
set rmargin screen LM+X1
set tmargin screen TM
set bmargin screen TM-Y

set yrange [0:0.22]
set xtics nomirror 

set label 1 '\textbf{a}' at screen 0,TM-DYLAB
set ylabel 'Occupation' offset screen 0.018
set key samp 1.25 center at screen LM+X1+DXKEY1, TM-Y/2

set xlabel 'Time (ps)' offset screen 0,0.01

plot -1 lc rgb "white" t '\phantom{.}' \
, "01_none_Lindblad.out" w l lc 4 lw LW2 dt 3 t '\phantom{.}'\
, "03_calculate_phonon_PT.out" w l lw LW lc rgb "pink" t '\phantom{.}'\
, "06_combine_phonon_photon.out" w l lw LW lc 1 t '\phantom{.}'\
, "02_iQUAPI.out" every 5 w p pt 6 lc 3 lw LW t '\phantom{.}'\



set lmargin screen LM
set rmargin screen LM+X2
set tmargin screen TM-Y-DY
set bmargin screen TM-2*Y-DY

set key samp 1.5 center at screen LM+X1+DXKEY1, TM-Y-DY-Y/2

set xtics nomirror format '%g'
set yrange [0:1.05]
set ytics (0, 0.2, 0.4, 0.6, 0.8, 1)

set label 1 '\textbf{b}' at screen 0, TM-Y-DY-DYLAB
set ylabel 'Occupation' offset screen 0.005
set key samp 1.5 

plot -1 lc rgb "white" t '\phantom{.}' \
, "11_initial_iQUAPI.out" every 5:1:5 w p pt 6 lc 3 lw LW t '\phantom{.}'\
, "09_initial_nophonon_photon.out" w l dt 3 lc 1 lw 2*LW  t '\phantom{.}'\
, "10_initial_phonon_photon.out" w l lw LW lc 2 t '\phantom{.}'\
, "07_calculate_photon_PT_BW0.4.out" w l dt 3 lc rgb "red" lw LW2 t '\phantom{.}'\
, "08_combine_phonon_photon_BW0.4.out" w l lw LW lc 4 t '\phantom{.}'\


unset multiplot

