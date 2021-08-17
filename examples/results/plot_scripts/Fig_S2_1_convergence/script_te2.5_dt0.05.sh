#!/bin/bash

N=10
te=2.5
dt=0.05

for thr in 1e-3 6e-4 3e-4 1e-4 6e-5 3e-5 1e-5 6e-6 3e-6 1e-6 \
6e-7 3e-7 1e-7 6e-8 3e-8 1e-8 6e-9 3e-9 1e-9 6e-10 3e-10; do

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

#write_PT            ${prefix}.pt
outfile             ${prefix}.out 
EOF

echo "parameter file: $(ls ${prefix}.param)"

/usr/bin/time -v ACE ${prefix}.param 2>&1 | tee ${prefix}.debug

done

