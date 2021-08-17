#!/bin/bash


te=2.5
dt=0.05
thr=1e-6

outfile=analyze_te${te}.out
rm -f $outfile

for N in 2 4 6 8 10 14 20 40 60 80 100 140 200 400; do

prefix=N${N}_te${te}_dt${dt}_thr${thr}

elapsed=$(./extract_elapsed.sh ${prefix}.debug)
maxdim=$(cat ${prefix}.debug | grep "Maxdim: after sweep" |tail -1 | awk '{print $4}')

echo $N $elapsed $maxdim>> $outfile

cat ${prefix}.dims |awk '{print $2}' > ${prefix}.dims2
echo 1 >> ${prefix}.dims2

done

