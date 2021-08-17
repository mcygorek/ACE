set term epslatex standalone color size 22cm, 12cm
set out 'figS2_2.tex'

set multiplot layout 6,1
TMA=0.9
BMA=0.6
LMA=0.09 #0.08
RMA=0.49
set tmargin screen TMA
set bmargin screen BMA
set lmargin screen LMA
set rmargin screen RMA


set ylabel 'Occupation' offset screen -0.01
set xlabel '$\gamma t$'

set label 1 '\textbf{a}' at screen LMA-0.073, TMA+0.077  #0.08
#set label 1 '\textbf{a}' at screen LMA-0.062, TMA+0.05

set key center top outside maxrow 2 samplen 1.5
LW=3

set yrange [0:1.1]
set xtics 0.5


plot \
 "N2_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW t '$\,N_E=2$' \
, "N4_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW t '$\,N_E=4$' \
, "N6_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW t '$\,N_E=6$' \
, "N8_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW t '$\,N_E=8$' \
, "N100_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW lc rgb "black" t '$\,N_E=100$' \
, "N8_te5_dt0.05_thr1e-6.out" u 1:2 w l lw LW lt 4  notit \



#-----------------------
set xrange [0:5]
unset xtics
set yrange [0:300]
set ytics 100
unset ylabel
#set ylabel '$d_i$'

TM=0.977
LM=0.6
RM=0.88

Y=0.1 #0.095
DY=0.05

DXKEY=1-RM
DYKEY=0
SPL=1.

unset xlabel
unset key
#set key outside samplen 1.5 maxcolumns 1 
set lmargin screen LM
set rmargin screen RM

set label 1 '\textbf{b}' at screen LM-0.066, TM

LW=2.5

scale=0.05
BW=1*scale

set style fill solid 0.8 border -1
set boxwidth 1 relative

array A[1]=[500]

n=0
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

plot "N2_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW t '$\,N_E=2$'  \
#, A w p lc rgb "white" t '$t_\textrm{calc.}<1 s$'

n=1
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

plot "N4_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW lt 2 t '$\,N_E=4$' \
#, A w p lc rgb "white" t '$t_\textrm{calc.}=$'

n=2
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

plot "N6_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW lt 3 t '$\,N_E=6$' 

n=3
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

plot "N8_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW lt 4 t '$\,N_E=8$' 

n=4
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

plot "N10_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW lc rgb "red" t '$\,N_E=10$' 

n=5
THISTM=TM-n*(Y+DY)
THISBM=TM-n*(Y+DY)-Y
set tmargin screen THISTM
set bmargin screen THISBM
set key at screen RM+DXKEY,THISTM+DYKEY samplen SPL right

set label 1 'Inner dimension $d_i$' center rotate at screen LM-0.06, TM-2.5*(Y+DY)-Y/2.
set xlabel '$\gamma t$'
set xtics 0.5

plot "N100_te5_dt0.05_thr1e-6.dims2" u (($0+0.5)*scale):1 w boxes lw LW lc rgb "black" t '$\,N_E=100$' 



#==================================
TMC=0.33 #0.36
BMC=0.1
RMC=RMA #-0.07
LMC=LMA
set tmargin screen TMC
set bmargin screen BMC
set lmargin screen LMC
set rmargin screen RMC

set label 1 '\textbf{c}' norotate at screen LMC-0.071, TMC+0.13 #+0.08

set ylabel 'Computation time' offset screen 0.02 
set xlabel '$N_E$
set xtics auto
set log x
set xrange [1:1000]

set ytics auto # nomirror
set ytics format '$10^{%T}$'
set yrange [0.1:1.5e4] #[0.1:5e4]  

set log y 

LW=3
DT2=2
DT3=3

timeunit=1

alph=1e-8
g(x,nmax)=alph*nmax
fit log(g(x,100)) "analyze_te5.out" u ($1):( log($2/($1*$3*$3*$3))) via alph

#alph=10
#f(x,nmax)=alph*nmax*x
#set xrange restore
#fit f(x, 100)  "analyze_te5.out" u ($1*$3*$3*$3):2 via alph

DYLABC=0.05
DYLABC_OFFSET=0.04
set key left at screen LMC+0.12,TMC+0.1 spac 1.3 samplen 1.5 

set label 2 '$t_\textnormal{CPU}$' at screen LMC,TMC+DYLABC_OFFSET+DYLABC 
set label 3 '$\alpha\, n_\textrm{max} \,N_E \,d_\textrm{max}^3$' at screen LMC,TMC+DYLABC_OFFSET
set label 4 '$\gamma t_e=5$' at screen LMC+0.149,TMC+DYLABC_OFFSET+DYLABC*2 -0.02
set label 5 '$\gamma t_e=2.5$' at screen LMC+0.231,TMC+DYLABC_OFFSET+DYLABC*2 -0.02


set ytics ("1 s" 1, "10 s" 10, "1 min" 60, "10 min" 600, "1 h" 3600, "6 h" 6*3600) offset screen 0.005

plot "analyze_te5.out" u 1:(($2)/timeunit) w lp lw LW t '\phantom{tes}' \
, "" u 1:(alph*100*$1*$3*$3*$3) w lp lw LW t '\phantom{.}' \
, "analyze_te2.5.out" u 1:(($2)/timeunit) w lp lw LW lc 1 dt DT2 t '\phantom{.}' \
, "" u 1:(alph*50*$1*$3*$3*$3) w lp lw LW lc 2 dt DT2 t '\phantom{.}' \

print "alph: ", alph

unset multiplot
