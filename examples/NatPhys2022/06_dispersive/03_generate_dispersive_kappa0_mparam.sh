#!/bin/bash

for i in $(seq 1 1 4); do 
file=03_dispersive_kappa0.mparam${i}
fac=$(echo "10+$i*1"|bc -l)
echo "
Nintermediate 10
add_Hamiltonian  { hbar * 1* ((sigma_z) otimes n_5)  + hbar*${fac}*(Id_2 otimes n_5)}  
add_Pulse  Gauss $(echo "$i*10"|bc -l) 0.2 {2.} {hbar*${fac}}  { Id_2 otimes bdagger_5}
" > $file

done
