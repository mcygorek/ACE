#!/bin/bash

if [ -z $gamma ]; then 
  gamma=0.1
fi

if [ -z $c ]; then 
  c=1.
fi

if [ -z $outfile ]; then
  outfile=00_lorentzian_c1_gamma0.1.dat
fi

rm -f $outfile

for x in $(seq -5 0.001 5); do
  echo $x $(echo "pi=4*a(1); 1./pi* $gamma /( ($x -$c)^2 + $gamma^2 )" | bc -l) >> $outfile
done 

