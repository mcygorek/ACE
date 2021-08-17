#!/bin/bash

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

for N in $(seq 60 -2 2); do
  dE=$(echo $N | awk '{print $1/60.*3.}' )
  Emin=$(echo $dE | awk '{print 1.5-$1/2}' )
  Emax=$(echo $dE | awk '{print 1.5+$1/2}' )

  prefix=ACE_N${N}_dE${dE}_M${M}_T${T}_dt${dt}_thr${thr}

cat > ${prefix}.param <<EOF
ta                    0
dt                    $dt
te                   20
Nintermediate        20

threshold             $thr
dict_zero             1e-12


add_Hamiltonian            {-1.5*|1><1|_2}
add_Pulse Gauss  7 5 {3*pi} 0   {|1><0|_2}
add_Lindblad     0.1       {|0><1|_2}

Boson_N_modes                      $N
Boson_subtract_polaron_shift       false
Boson_M                            $M
Boson_E_min                        $Emin
Boson_E_max                        $Emax
Boson_J_type                       QDPhonon
Boson_sort                         none
temperature                        $T

outfile                            ${prefix}.out
EOF

/usr/bin/time -v ACE ${prefix}.param 2>&1 |tee ${prefix}.debug

done

