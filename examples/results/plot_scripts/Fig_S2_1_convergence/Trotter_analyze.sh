#!/bin/bash


ref=N10_te2.5_dt0.005_thr3e-11.out

outfile="Trotter_analyze.out"
rm -f $outfile

cat $ref | awk '{if(NR%2==1){print $0}}' > ref_every2.out
cat $ref | awk '{if(NR%4==1){print $0}}' > ref_every4.out
cat $ref | awk '{if(NR%6==1){print $0}}' > ref_every6.out
cat $ref | awk '{if(NR%8==1){print $0}}' > ref_every8.out
cat $ref | awk '{if(NR%10==1){print $0}}' > ref_every10.out
cat $ref | awk '{if(NR%20==1){print $0}}' > ref_every20.out

echo 0.1 $(./max_err.sh N10_te2.5_dt0.1_thr6e-9.out ref_every20.out) >> $outfile
echo 0.05 $(./max_err.sh N10_te2.5_dt0.05_thr3e-10.out ref_every10.out) >> $outfile
echo 0.03 $(./max_err.sh N10_te2.5_dt0.03_thr1e-10.out ref_every6.out) >> $outfile
echo 0.02 $(./max_err.sh N10_te2.5_dt0.02_thr1e-10.out ref_every4.out) >> $outfile
echo 0.01 $(./max_err.sh N10_te2.5_dt0.01_thr1e-11.out ref_every2.out) >> $outfile




