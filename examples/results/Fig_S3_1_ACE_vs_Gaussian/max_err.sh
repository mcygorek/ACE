#!/bin/bash

field=2

#prints max_t |f(x)-g(x)|
if [[ $# -lt 2 ]]; then
    echo "Usage: max_error file1 file2 [field]" >&2
    exit 1
fi

if [[ $# -gt 2 ]]; then 
    field=$3
fi

if [ -n "$field" ] && [ "$field" -eq "$field" ] 2>/dev/null; then
  echo "test" >/dev/null
else
  echo "Third parameter '$field' not a number!" >&2
  exit 1
fi


file1=$1
file2=$2


maxfield=$(tail -n 1 $file1 | awk '{print NF}')
if [ -n "$maxfield" ] && [ "$maxfield" -eq "$maxfield" ] 2>/dev/null; then
  echo "test" >/dev/null
else
  echo "maxfield: '$maxfield' not a number!" >&2
  exit 1
fi

field2=$(echo "$maxfield + $field" | bc -l)

#echo "file1='$file1' file2='$file2' field=$field field2=$field2"  >&2

#paste $file1 $file2 | awk '{print $'$field',$'$field2'}' 

max_err=$(paste $file1 $file2 | awk '{print sqrt(($'$field'-$'$field2')*($'$field'-$'$field2'))}' |sort -g | tail -1 )

echo $max_err

