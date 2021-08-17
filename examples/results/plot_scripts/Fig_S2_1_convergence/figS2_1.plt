set term epslatex standalone color size 24cm, 14cm
set out 'figS2_1.tex'


LM=0.08
TM=0.92
X=0.395
DX=0.11
Y=0.3
DY=0.20
DYKEY=0.04

set multiplot

nx=0; ny=0
set lmargin screen LM+nx*(X+DX)
set rmargin screen LM+nx*(X+DX)+X
set tmargin screen TM-ny*(Y+DY)
set bmargin screen TM-ny*(Y+DY)-Y
set ylabel 'Compression error'
set xlabel 'SVD truncation threshold $\epsilon$'
set label 1 '\textbf{a}' at screen LM+nx*(X+DX)-0.067, TM-ny*(Y+DY)+0.055

set log xy 
set xtics format '$10^{%T}$'
set ytics format '$10^{%T}$'
set key center at screen LM+nx*(X+DX)+X/2,TM-ny*(Y+DY)+DYKEY maxrow 2 samplen 1.5
LW=2

set xrange [1e-3:3e-12]
set yrange [1e-6:1e0]

plot \
 "analyze_dt0.001.out" u ($1):5 w lp lw LW t '$\gamma\Delta t=0.001$' \
, "analyze_dt0.005.out" u ($1):5 w lp lw LW t '$\gamma\Delta t=0.005$' \
, "analyze_dt0.01.out" u ($1):5 w lp lw LW t '$\gamma\Delta t=0.01$' \
,"analyze_dt0.05.out" u ($1):5 w lp lw LW t '$\gamma\Delta t=0.05$' \
,"analyze_dt0.1.out" u ($1):5 w lp lw LW t '$\gamma\Delta t=0.1$' \


nx=1; ny=0
set lmargin screen LM+nx*(X+DX)
set rmargin screen LM+nx*(X+DX)+X
set tmargin screen TM-ny*(Y+DY)
set bmargin screen TM-ny*(Y+DY)-Y
set label 1 '\textbf{b}' at screen LM+nx*(X+DX)-0.067, TM-ny*(Y+DY)+0.055
set xlabel 'Maximal inner dimension $d_\textrm{max}$'
#set xtics auto

set log xy 
#set xtics format '$10^{%T}$'
set xtics format '%g'
set xtics ("4" 4, "6" 6, "10" 10, "16" 16, "30" 30, "50" 50, "100" 100, "200" 200,  "400" 400)
#set xtics add (4) add (30) add (500)
set ytics format '$10^{%T}$'
set key center at screen LM+nx*(X+DX)+X/2,TM-ny*(Y+DY)+DYKEY maxrow 2 samplen 1.5

set xrange [4:400]
set yrange [1e-6:1e0]

plot \
 "analyze_dt0.001.out" u ($4):5 w lp lw LW t '$\gamma\Delta t=0.001$' \
, "analyze_dt0.005.out" u ($4):5 w lp lw LW t '$\gamma\Delta t=0.005$' \
, "analyze_dt0.01.out" u ($4):5 w lp lw LW t '$\gamma\Delta t=0.01$' \
,"analyze_dt0.05.out" u ($4):5 w lp lw LW t '$\gamma\Delta t=0.05$' \
,"analyze_dt0.1.out" u ($4):5 w lp lw LW t '$\gamma\Delta t=0.1$' \


nx=0; ny=1
set lmargin screen LM+nx*(X+DX)
set rmargin screen LM+nx*(X+DX)+X
set tmargin screen TM-ny*(Y+DY)
set bmargin screen TM-ny*(Y+DY)-Y
set label 1 '\textbf{c}' at screen LM+nx*(X+DX)-0.067, TM-ny*(Y+DY)+0.055
set xlabel 'CPU time'

set log xy 
set xtics ("1 s" 1, "10 s" 10, "1 min" 60, "10 min" 600, "1 h" 3600)
#set xtics format '$10^{%T}$'
set ytics format '$10^{%T}$'
set key center at screen LM+nx*(X+DX)+X/2,TM-ny*(Y+DY)+DYKEY maxrow 2 samplen 1.5

set xrange [1:1*3600]
set yrange [1e-6:1e0]

scale=1.

plot \
 "analyze_dt0.001.out" u (scale*($2)):5 w lp lw LW t '$\gamma\Delta t=0.001$' \
, "analyze_dt0.005.out" u (scale*($2)):5 w lp lw LW t '$\gamma\Delta t=0.005$' \
, "analyze_dt0.01.out" u (scale*($2)):5 w lp lw LW t '$\gamma\Delta t=0.01$' \
,"analyze_dt0.05.out" u (scale*($2)):5 w lp lw LW t '$\gamma\Delta t=0.05$' \
, "analyze_dt0.1.out" u (scale*($2)):5 w lp lw LW t '$\gamma\Delta t=0.1$' \


#---------------------------------------------
nx=1; ny=1
set lmargin screen LM+nx*(X+DX)
set rmargin screen LM+nx*(X+DX)+X
set tmargin screen TM-ny*(Y+DY)
set bmargin screen TM-ny*(Y+DY)-Y
set label 1 '\textbf{d}' at screen LM+nx*(X+DX)-0.067, TM-ny*(Y+DY)+0.055

set log xy 
unset xtics 
unset ytics
#set xtics format '$10^{%T}$'
set xtics auto
set ytics format '$10^{%T}$'
#unset key
set key center at screen LM+nx*(X+DX)+X/2+0.1,TM-ny*(Y+DY)+DYKEY maxrow 2 samplen 1.5

set xlabel 'Time step width $\gamma\Delta t$'
set ylabel 'Trotter error'
set xrange [1e-2:1e-1]
set yrange [1e-6:1e-3]

g(x)= a*(x**n)
n=2
a=0.01
fit [log(0.00001):log(10000)] log(g(x)) "Trotter_analyze.out" u 1:(log($2)) via a

plot \
  "Trotter_analyze.out" u 1:2 w p pt 7 t 'ACE' \
, [0.1:0.01] g(x) t 'Fit to $\alpha \Delta t^2$'

#print a

unset multiplot
