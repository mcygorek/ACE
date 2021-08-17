#!/bin/bash

debugfile=$1

string=$(grep "Elapsed" $debugfile  |awk '{print $8}' )

echo $string |  awk -F ":"  '{if(NF==3){print (($1*60)+$2)*60+$3} else {print ($1*60)+$2}}'
