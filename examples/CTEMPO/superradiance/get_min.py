import sys
import numpy as np

col=2
if len(sys.argv)<2:
  raise ValueError(f"usage: python3 get_max.py FILENAME [COL=2]")
if len(sys.argv)>2:
  col=int(sys.argv[2])

#print(sys.argv)
#print(col)

filename=sys.argv[1]

a = np.loadtxt(filename)
i = np.argmin(a[:,col-1])
print(f'{a[i,0]} {a[i,col-1]}')

