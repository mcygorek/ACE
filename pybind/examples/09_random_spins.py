# if not installed via pip, comment out the following and plug in the correct path
# import sys
# sys.path.append('.../ACE/pybind/') #<---plug in correct path
from ACE import *
import numpy as np


# Parameters
ta = 0      # start time
t_op = 0    # time at which first operator is applied
te = 2.5     # end time 
dt = 0.01   # time step
thr = 1e-8
tr = 1
N = 10       # number of spins

np.random.seed(12345)


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


Heisenberg = hbar/4*( np.kron(ACE_sigma_x, ACE_sigma_x)+np.kron(ACE_sigma_y, ACE_sigma_y)+np.kron(ACE_sigma_z, ACE_sigma_z)) 


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
    m.initial=1/2*(np.eye(2,2) + b[0]*ACE_sigma_x + b[1]*ACE_sigma_y + b[2]*ACE_sigma_z)
    mode_list.append(m)
    print(f'J={J} b={b}')
   
print(f'Total Environment Spin S={S}')
    
PT = ProcessTensors(param, mode_list)   


outfile = "09_random_spins.out"
outp  = OutputPrinter( outfile,
                      [ACE_sigma_x, ACE_sigma_y, ACE_sigma_z])

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


