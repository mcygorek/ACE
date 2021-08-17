set term epslatex color standalone size 12cm,9.5cm
set out "fig5_resub.tex"


set multiplot layout 2,2 

LW=2

set palette rgb 33,13,10; #rainbow

set cbrange [0:1.0]
set style line 1 lw LW lc rgb "black"
set style line 2 lw LW lc palette cb 0
set style line 3 lw LW lc palette cb 0.2
set style line 4 lw LW lc palette cb 0.4
set style line 5 lw LW lc palette cb 0.6
set style line 6 lw LW lc palette cb 0.8
set style line 7 lw LW lc palette cb 1

B=0.125
L=0.125
WA=0.34  
WB=WA
WCD=0.34

DX1=0.18
DX2=0.18

H1=0.25
H2=0.32
DY=0.28

DXLAB=-0.12
DYLAB=0.0
DXLAB2=DXLAB
DYLAB2=DYLAB
LABX3=0
LABY3=0.12

THIST=B+H1+H2+DY
THISB=B+H2+DY
THISL=L
THISR=L+WA

set tmargin screen THIST
set bmargin screen THISB
set lmargin screen THISL
set rmargin screen THISR

set label 1 '\textbf{a}' at screen THISL+DXLAB, THIST+DYLAB

set ylabel '$v(x)$' #'$E/(2\epsilon)$'
set xlabel '$x$'    #'$r/a_0$'

set xrange [-1:7.5]
set yrange [-25:8]

scale=30
LWL=2
LW=2

# from printEV.E and printEV.X files:
array E[5]=[-20.2452,-12.2343,-6.22501,-2.22649,-0.240453]
array X[5]=[0.515171,1.76693,3.52585,6.36541,13.3156]

xcoupl=sqrt(2*5)

#set key samplen 3 spac 1.5 at screen THISL+KEYX, THIST
unset key

plot "printEV.EV" u 1:($2) w l lw 3 lc rgb "black" t 'v(x)' \
, "" u 1:(scale*$3+E[1]):(E[1]) w filledcurves lw LW lc 1 t '$\psi_0(x)$' \
, "" u 1:(scale*-$4+E[2]):(E[2]) w filledcurves lw LW lc 2 t '$\psi_1(x)$' \
, "" u 1:(scale*-$5+E[3]):(E[3]) w filledcurves lw LW lc 3 t '$\psi_2(x)$' \
, "" u 1:(scale*$6+E[4]):(E[4]) w filledcurves lw LW lc 4 t '$\psi_3(x)$' \
, "" u 1:(scale*$7+E[5]):(E[5]) w filledcurves lw LW lc 5 t '$\psi_4(x)$' \
, "" u 1:($2) w l lw 2 lc rgb "black" notit \
, [i=1:5:1] '+' u (X[i]/xcoupl):(E[i]) w p ps 1.7 lw 4 pt 2 lc rgb "0x000000" t '$\langle i|\hat{x}| i\rangle$'



THISL=L+WA+DX1
THISR=L+WA+DX1+WB

set lmargin screen THISL
set rmargin screen THISR

set label 1 '\textbf{b}' at screen THISL+DXLAB, THIST+DYLAB

set xrange [0:7.6]
set yrange [0:0.16]
set ytics 0.05

scale=2.
LWL=2
LW=2

dE=2*0.0244386
hbar=0.6582119569

set xlabel '$\omega_k/\Omega$' 
set ylabel '$g_k/\Omega$' offset screen 0.015
set boxwidth dE/hbar
set style fill solid 1.00 border rgb "#406090"

unset key

plot "printEV.E_g" u ($1/hbar):2 w boxes lt rgb "#80C0FF" t '$\gamma_k$'


THIST=B+H2
THISB=B
THISL=L
THISR=L+WCD

set tmargin screen THIST
set bmargin screen THISB
set lmargin screen THISL
set rmargin screen THISR

unset colorbox
set key samplen 2 spac 1.3 right at screen L+2*WCD+DX2, THIST+LABY3 maxrows 2 width 1.8
set yrange [0:1]
set ytics 0.2

set xlabel '$\Omega t$'
set ylabel 'Occupation $n_e$' offset screen 0

set label 1 '\textbf{c}' at screen THISL+DXLAB2, THIST+DYLAB2

plot  "01_independent_boson_T0.5_M5_Jscale0.1.out" every 3 w p lc rgb "#000000" pt 6 ps 1.3 t 'SBM' \
, "02_harmonic_oscillator_T0.5_M5_Jscale0.1_HO.out" every 3 w p lc rgb "#000000" pt 2  t 'HO' \
, "08_morse_T0.5_M5_Jscale0.1_L100.out" w l ls 2 t '$\Lambda=100$' \
, "07_morse_T0.5_M5_Jscale0.1_L10.out" w l ls 3 t '$\Lambda=10$' \
, "06_morse_T0.5_M5_Jscale0.1_L5.out" w l ls 4 t '$\Lambda=5$' \
, "05_morse_T0.5_M4_Jscale0.1_L4.out" w l ls 5 t '$\Lambda=4$' \
, "04_morse_T0.5_M3_Jscale0.1_L3.out" w l ls 6 t '$\Lambda=3$' \
, "03_morse_T0.5_M2_Jscale0.1_L2.out" w l ls 7 t '$\Lambda=2$' \

unset key

THISL=L+WCD+DX2
THISR=L+2*WCD+DX2

set lmargin screen THISL
set rmargin screen THISR

set label 1 '\textbf{d}' at screen THISL+DXLAB2, THIST+DYLAB2

plot  "01_independent_boson_T0.5_M5_Jscale0.1.out" every 3 w p lc rgb "#000000" pt 6 ps 1.3 t 'Independent Boson' \
, "02_harmonic_oscillator_T0.5_M5_Jscale0.1_HO.out" every 3 w p lc rgb "#000000" pt 2  t 'Harmonic Oscillator' \
, "14_morse_T0.5_M5_Jscale0.1_L100_renorm.out" w l ls 2 t '$\Lambda=100$' \
, "13_morse_T0.5_M5_Jscale0.1_L10_renorm.out" w l ls 3 t '$\Lambda=10$' \
, "12_morse_T0.5_M5_Jscale0.1_L5_renorm.out" w l ls 4 t '$\Lambda=5$' \
, "11_morse_T0.5_M4_Jscale0.1_L4_renorm.out" w l ls 5 t '$\Lambda=4$' \
, "10_morse_T0.5_M3_Jscale0.1_L3_renorm.out" w l ls 6 t '$\Lambda=3$' \
, "09_morse_T0.5_M2_Jscale0.1_L2_renorm.out" w l ls 7 t '$\Lambda=2$' \


unset multiplot


