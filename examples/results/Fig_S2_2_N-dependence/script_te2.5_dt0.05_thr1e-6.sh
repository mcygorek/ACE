#!/bin/bash

N=10
te=2.5
dt=0.05
thr=1e-6

for N in 2 4 6 8 10 14 20 40 60 80 100 140 200 400; do

prefix=N${N}_te${te}_dt${dt}_thr${thr}

cat > ${prefix}.param  << EOF 
te                  ${te}
dt                  ${dt}
threshold           ${thr}

IF_print_timesteps  true
Fermion_N_modes       $N
Fermion_rate          1
Fermion_omega_min    -5
Fermion_omega_max     5
Fermion_EFermi        1e4
temperature         0

print_IF_dims_to_file  ${prefix}.dims
outfile             ${prefix}.out 
EOF

echo "parameter file: $(ls ${prefix}.param)"

/usr/bin/time -v ACE ${prefix}.param 2>&1 | tee ${prefix}.debug

done

