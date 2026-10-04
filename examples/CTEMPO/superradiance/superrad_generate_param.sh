#!/bin/bash
LANG=US
trap 'echo "Ctrl + C detected"; exit 1' INT


CCTEMPO="CCTEMPO"    #use explicit path if CCTEMPO binary not in $PATH

if [ -z "$Tdecay" ]; then Tdecay=500; fi
if [ -z "$te" ];     then te=1000; fi
if [ -z "$tmem" ];   then tmem=25.6; fi
if [ -z "$dt" ];     then dt=0.1; fi
if [ -z "$T" ];      then T=4; fi
if [ -z "$thr" ];    then thr=1e-7; fi
if [ -z "$N" ];      then N=4; fi
if [ -z "$bs" ];     then bs=-1; fi

prefix=N${N}_T${T}_Tdecay${Tdecay}_te${te}_dt${dt}_thr${thr}  #_bs${bs}

cat << EOF >${prefix}.param
te         $te 
t_mem      $tmem
dt         $dt
threshold  $thr
set_precision 12
buffer_blocksize $bs
outfile    ${prefix}.out
N_sites    $N
add_collective_decay {1/${Tdecay}}
S0_add_Output  {|1><1|_2}
EOF
for i in $(seq 0 1 $(bc -l <<< $N-1)); do 
  echo "S${i}_initial {|1><1|_2}" >> ${prefix}.param
done
if [ "$T" != "-1" ]; then
  echo "S0_Boson_J_type       QDPhonon" >> ${prefix}.param
  echo "S0_Boson_omega_max   10" >> ${prefix}.param
  echo "S0_Boson_temperature  $T" >> ${prefix}.param
  for i in $(seq 1 1 $(bc -l <<< $N-1)); do
    echo "S${i}_env_sameas_S 0" >> ${prefix}.param
  done
fi

