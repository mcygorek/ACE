#!/bin/bash

for i in $(seq 1 1 4); do 
file=02_instant_Fock.mparam${i}
fac=$(echo "10+$i*1"|bc -l)
echo "
add_Hamiltonian  { hbar * 1* ((sigma_z) otimes n_3)  + hbar*${fac}*(Id_2 otimes n_3)}  
apply_Operator_left  $(echo "$i*10"|bc -l)  {(Id_2 otimes bdagger_3)}
apply_Operator_right $(echo "$i*10"|bc -l)  {(Id_2 otimes b_3)}
" > $file

done
