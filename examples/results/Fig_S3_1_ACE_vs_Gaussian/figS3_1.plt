set term epslatex color standalone size 18cm, 7cm
set out 'figS3_1.tex'

TM=0.99
H=0.82
LM=0.1
X1=0.3
DX=0.15
X2=0.42
DXLAB=-0.1
DXLAB2=-0.12
DYLAB=0.02

set multiplot

set tmargin screen TM
set bmargin screen TM-H
set lmargin screen LM
set rmargin screen LM+X1

set yrange [0:0.13]
set xrange [0:4]
set xlabel 'Energy $\hbar\omega$ (meV)'
set ylabel '$J(\omega)$ (ps$^{-1}$)'

set label 1 '\textbf{a}' at screen LM+DXLAB, TM-DYLAB
set arrow 2 from 1.25,0.109 to 1.75,0.109 heads size 0.05,90,45
set label 3 '$\Delta E$' center at 1.5,0.119

hbar=0.6582119569
LW=2
unset key

set style fill solid 1.00 border rgb "#406090"

plot "print_J_dE0.5.dat" u ($1*hbar):2 w boxes  lt rgb "#80C0FF" \
, "print_J_dE10.dat" u ($1*hbar):2 w l lw LW lc rgb "black"\


#-------------------

set lmargin screen LM+X1+DX
set rmargin screen LM+X1+DX+X2

set xrange [0:3]
set log y
set yrange [1:3600]
set key bottom samplen 1.25

LW=2
LW2=3

set label 1 '\textbf{b}' at screen LM+X1+DX+DXLAB2, TM-DYLAB
unset label 3
unset arrow 2

set xlabel 'Spectral density width $\Delta E$ (meV)'
set ylabel 'Computation time'
#set ytics ("1 s" 1, "10 s" 10, "30 s" 30, "1 min" 60, "10 min" 600, "30 min" 1800)
set ytics ("1 s" 1, "3 s" 3, "10 s" 10, "30 s" 30, "1 min" 60, "3 min" 3*60, "10 min" 600, "30 min" 1800) #, "1 h" 3600)

plot "analyze_ACE.dat" u 2:3 w lp lw LW t 'ACE $M=2$' \
, "analyze_ACE_M3.dat" u 2:3 w lp lw LW t 'ACE $M=3$' \
, "analyze_Gaussian.dat" u 2:3 w lp lw LW pt 10 ps 1.25 t 'Gaussian PT' \
, "analyze_TEMPO.dat" u 2:3 w lp lw LW t 'TEMPO (full)' \
, "analyze_TEMPO_tmem2.5_inaccurate.dat" u 2:3 w p ps 0.6 lc rgb "0xFFAAAA" notit \
, "analyze_TEMPO_tmem2.5_accurate.dat" u 2:3 w lp pt 5 lc rgb "0xFFAAAA" lw LW t 'TEMPO ($\tau_\textrm{mem}=2.5$ ps)'


unset multiplot

