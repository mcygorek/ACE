#!/bin/bash


for pfile in \
N10_te2.5_dt0.02_thr1e-10.param \
N10_te2.5_dt0.03_thr1e-10.param \
N10_te2.5_dt0.04_thr1e-10.param \
; do
#N10_te2.5_dt0.06_thr3e-10.param \
#N10_te2.5_dt0.07_thr3e-10.param \
#N10_te2.5_dt0.08_thr6e-10.param \
#N10_te2.5_dt0.09_thr1e-9.param \
#; do

prefix=$(basename $pfile ".param")
 
/usr/bin/time -v ACE ${prefix}.param 2>&1 |tee ${prefix}.debug


done 
