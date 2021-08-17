set term epslatex color standalone size 20cm,16cm
set out "figS6_1_data.tex"

set multiplot

TM=0.97
Y=0.163
DY=0.08   #0.075
RM=0.995
LM=0.07
DXKEY=0
DYKEY=0.03
DXLAB=0.
DYLAB=0.0178

LW=3
LW2=2
set ylabel '$n_e$'
set yrange [0:1]  #[-0.1:1.1]

zw=0.02
gw=0.2/2
stdcolor="#000000"

set xtics format ''
set lmargin screen LM
set rmargin screen RM
n=0
set tmargin screen TM-n*(Y+DY)
set bmargin screen TM-n*(Y+DY)-Y
set object 2 rectangle from 10-zw,0 to 10+zw,1 fc rgb "#000080" fs solid noborder
set object 3 rectangle from 20-zw,0 to 20+zw,1 fc rgb "#0000ff" fs solid noborder
set object 4 rectangle from 30-zw,0 to 30+zw,1 fc rgb "#0066ff" fs solid noborder
set object 5 rectangle from 40-zw,0 to 40+zw,1 fc rgb "#80b3ff" fs solid noborder
set key samplen 1.25 right at screen RM+DXKEY,TM-n*(Y+DY)+DYKEY
plot "02_instant_Fock.out" w l lw LW lc rgb stdcolor t 'Instantaneous Fock state preparation' 

n=1
set tmargin screen TM-n*(Y+DY)
set bmargin screen TM-n*(Y+DY)-Y

set object 2 rectangle from 10-gw,0 to 10+gw,1 fc rgb "#000080" fs solid noborder
set object 3 rectangle from 20-gw,0 to 20+gw,1 fc rgb "#0000ff" fs solid noborder
set object 4 rectangle from 30-gw,0 to 30+gw,1 fc rgb "#0066ff" fs solid noborder
set object 5 rectangle from 40-gw,0 to 40+gw,1 fc rgb "#80b3ff" fs solid noborder
set key samplen 1.25 right at screen RM+DXKEY,TM-n*(Y+DY)+DYKEY
plot "03_dispersive_kappa0.out" w l lw LW lc rgb stdcolor t 'Pulsed excitation' 

n=2
set tmargin screen TM-n*(Y+DY)
set bmargin screen TM-n*(Y+DY)-Y
set object 2 rectangle from 10-gw,0 to 10+gw,1 fc rgb "#000080" fs solid noborder
unset object 3 
unset object 4 
unset object 5 
set key samplen 1.25 right at screen RM+DXKEY,TM-n*(Y+DY)+DYKEY*2 maxrow 2
set label 1 'Single mode:' right at screen RM+DXKEY-0.2,TM-n*(Y+DY)+DYKEY+DYLAB

plot "01_singlemode.out" u 1:2 w l lw LW lc rgb stdcolor t 'all $n$' \
, "" u 1:4 w l lw LW2 lc 1 t '$n=0$' \
, "" u 1:6 w l lw LW2 lc 2  t '$n=1$' \
, "" u 1:8 w l lw LW2 lc 4  t '$n=2$' \
, "" u 1:2 w l lw LW2 lc rgb stdcolor notit \

unset label 1


n=3
set tmargin screen TM-n*(Y+DY)
set bmargin screen TM-n*(Y+DY)-Y
set object 2 rectangle from 10-gw,0 to 10+gw,1 fc rgb "#000080" fs solid noborder
set object 3 rectangle from 20-gw,0 to 20+gw,1 fc rgb "#0000ff" fs solid noborder
set object 4 rectangle from 30-gw,0 to 30+gw,1 fc rgb "#0066ff" fs solid noborder
set object 5 rectangle from 40-gw,0 to 40+gw,1 fc rgb "#80b3ff" fs solid noborder
set key samplen 1.25 right at screen RM+DXKEY,TM-n*(Y+DY)+DYKEY

#set xtics format '%g'
set xtics ('$\tau_1$' 10, '$\tau_2$' 20, '$\tau_3$' 30, '$\tau_4$' 40)
set xlabel '$gt$'
plot "04_dispersive_kappa0.1.out" w l lw LW lc rgb stdcolor t 'Pulsed excitation with cavity losses' 

unset multiplot
