import numpy as np

for w in np.arange(0, 30+0.01/2., 0.01):
  print(w, 0.2*w*np.exp(-w/3))

