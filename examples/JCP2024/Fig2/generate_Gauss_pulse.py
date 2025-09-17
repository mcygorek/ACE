import numpy as np
A = 3*np.pi
t0 = 10
tau = 4

def f(t, t0, tau):
  sigma = tau/(2*np.sqrt(2*np.log(2)))
  return np.exp(-0.5*(t-t0)**2/sigma**2)/sigma/np.sqrt(2*np.pi)

for t in np.arange(0, 20.005, 0.01):
  print(t, A*f(t,t0,tau), 0)

