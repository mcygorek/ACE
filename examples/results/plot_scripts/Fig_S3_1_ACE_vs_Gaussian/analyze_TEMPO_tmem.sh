#!/bin/bash


#!/bin/bash

if [ -z $tmem ]; then
  tmem=2.5
fi
if [ -z $M ]; then
  M=2
fi
if [ -z $T ]; then
  T=0
fi
if [ -z $dt ]; then
  dt=0.1
fi
if [ -z $thr ]; then
  thr=1e-7
fi

outfile=analyze_TEMPO_tmem${tmem}.dat
rm -f $outfile

for N in $(seq 60 -2 2); do
  dE=$(echo $N | awk '{print $1/60.*3.}' )
  Emin=$(echo $dE | awk '{print 1.5-$1/2}' )
  Emax=$(echo $dE | awk '{print 1.5+$1/2}' )

  prefix=TEMPO_dE${dE}_T${T}_tmem${tmem}_t${dt}_thr${thr}

  debugfile=${prefix}.debug
  string=$(grep "Elapsed" $debugfile  |awk '{print $8}' )

  elapsed=$(echo $string |  awk -F ":"  '{if(NF==3){print (($1*60)+$2)*60+$3} else {print ($1*60)+$2}}')
  echo $N $dE $elapsed   >> $outfile

done

