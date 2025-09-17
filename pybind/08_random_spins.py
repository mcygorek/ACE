#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('/.../ACE/pybind/') #<---plug in your directory

from ACEutils import *


# Parameters
ta = 0      # start time
t_op = 0    # time at which first operator is applied
te = 2.5     # end time 
dt = 0.01   # time step
thr = 1e-8
tr = 1
N = 10       # number of spins

np.random.seed(12345)

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
sigma_y = np.array([[0,1j],[-1j,0]], dtype=complex)
sigma_z = np.array([[-1,0],[0,1]], dtype=complex)
Id2 = np.array([[1,0],[0,1]], dtype=complex)

#check commutation relations
print(sigma_x@sigma_y-sigma_y@sigma_x -2j*sigma_z)
print(sigma_y@sigma_z-sigma_z@sigma_y -2j*sigma_x)
print(sigma_z@sigma_x-sigma_x@sigma_z -2j*sigma_y)


rho0 = np.array([[0.5,0.5],[0.5,0.5]], dtype=complex)
initial = InitialState( rho0 )


fprop = FreePropagator()
fprop.add_Hamiltonian(np.zeros((2,2), dtype=complex))

plist =  [f'dt {dt}']
plist += [f'ta {ta}']
plist += [f'te {te}']
plist += [f'threshold {thr}']
plist += [f'threshold_range_factor {tr}']
param =Parameters(plist)
print(param)


#given distribution f(x), apply F^{-1}(uniform(0,1))
def sample_sin(x): 
    return np.arccos(x)


Heisenberg = hbar/4*( np.kron(sigma_x, sigma_x)+np.kron(sigma_y, sigma_y)+np.kron(sigma_z, sigma_z)) 


mode_list = []
S = np.array([0,0,0])
for i in range(N):
    m = ModePropagator()
    J = np.random.normal(0,1)  #coupling strength
    m.add_Hamiltonian( J*Heisenberg )
    phi = np.random.uniform(0,2*np.pi)
    theta = sample_sin(np.random.uniform(0,1))
#    if beta>0:
#        while np.random.uniform(0,1)>np.exp(beta*np.cos(theta)):
#            theta = sample_sin(np.random.uniform(0,1))
    b=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
    S=S+b
    m.initial=1/2*(Id2 + b[0]*sigma_x + b[1]*sigma_y + b[2]*sigma_z)
    mode_list.append(m)
    print(f'J={J} b={b}')
   
print(f'Total Environment Spin S={S}')
    
PT = ProcessTensors(param, mode_list)   


outfile = "06_random_spins.out"
outp  = OutputPrinter( outfile,
                      [sigma_x, sigma_y, sigma_z])

tgrid = TimeGrid(ta, te, dt)

sim   = Simulation()
# Run:
sim.run(fprop, PT, initial, tgrid, outp)

# Read:
(times, data) = read_outfile(outfile)
#print(data)


# Now, plot the results:
import matplotlib.pyplot as plt
plt.xlabel("Time (ps)")
plt.ylabel("Spin (hbar/2)")
plt.plot(times, data[:,0].real, label='S_x')
plt.plot(times, data[:,1].real, label='S_y')
plt.plot(times, data[:,2].real, label='S_z')
plt.legend(loc="upper right")
plt.show()


