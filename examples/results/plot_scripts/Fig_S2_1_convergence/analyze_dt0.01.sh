#!/bin/bash


N=10
te=2.5
dt=0.01
thr_ref=1e-11

outfile="analyze_dt${dt}.out"
touch $outfile; rm $outfile

for thr in 1e-3 6e-4 3e-4 1e-4 6e-5 3e-5 1e-5 6e-6 3e-6 1e-6 \
6e-7 3e-7 1e-7 6e-8 3e-8 1e-8 6e-9 3e-9 1e-9 6e-10 3e-10 1e-10 \
6e-11 3e-11 1e-11; do

prefix=N${N}_te${te}_dt${dt}_thr${thr}


file=${prefix}.debug
ofile=${prefix}.out
referencefile=N${N}_te${te}_dt${dt}_thr${thr_ref}.out
echo $file

elapsed=$(./extract_elapsed.sh ${file})
mem=$(cat $file |grep "Maximum resident" | awk '{print  $6}')
maxdim=$(cat $file | grep "Maxdim: after sweep" |tail -1 | awk '{print $4}')
err=$(./max_err.sh $ofile $referencefile)

echo $thr $elapsed $mem $maxdim $err>> $outfile

done 

TMP=$(mktemp)
cp $outfile $TMP
cat $TMP | sort -gr >$outfile


