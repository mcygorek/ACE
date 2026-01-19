#!/usr/bin/env python
# coding: utf-8


import sys
sys.path.append('/home/.../ACE/pybind/') #<---plug in your directory
from ACEutils import *

# This example demonstrate the (subtle) effects of phonons on superradiance of two QDs.
# Both QDs are coupled to phonons. We precompute the PT-MPO for a single QD first.
# The result is written to the file QDPhonon.pt

dt = 0.1
plist  = [f'dt {dt}']
plist += [f't_mem 6.4'] # Should be a power of (2^x)*dt for some integer x
plist += [f'te 12.8']   # Should be 2*t_mem
plist += [f'threshold 1e-7'] # <- Might not be fully converged but good enough
plist += [f'threshold_range_factor 10']
plist += [f'use_Gaussian_periodic true'] # Utilize periodic PT-MPOs
plist += [f'Boson_J_type QDPhonon']      # Define spectral density
plist += [f'Boson_omega_max 10']
plist += [f'Boson_temperature 4.2']
plist += [f'Boson_SysOp {{|1><1|_2}}']
plist += [f'write_PT QDPhonon.pt']

PTparam = Parameters(plist)
PT = ProcessTensors(PTparam)  # Calculates the PT-MPO
del PT # PT-MPO is written to file when the object is destroyed


# Now, simulate results

tgrid = TimeGrid(-3000/2, 3000/2, 0.1)
initial = InitialState(np.kron(KetBra(0,0,2),KetBra(0,0,2)))

fprop = FreePropagator()
fprop.add_Lindblad(0.002, np.kron(KetBra(1,0,2),np.eye(2,2)))
fprop.add_Lindblad(0.002, np.kron(np.eye(2,2),KetBra(1,0,2)))
fprop.add_Lindblad(0.002, np.kron(KetBra(0,1,2),np.eye(2,2))                         +np.kron(np.eye(2,2),KetBra(0,1,2)))
fprop.apply_Operator_left(0, np.kron(KetBra(0,1,2),np.eye(2,2))                         +np.kron(np.eye(2,2),KetBra(0,1,2)))
fprop.apply_Operator_right(0, np.kron(KetBra(1,0,2),np.eye(2,2))                         +np.kron(np.eye(2,2),KetBra(1,0,2)))

PT    = ProcessTensors()

outp  = OutputPrinter([(np.kron(KetBra(1,0,2),np.eye(2))+np.kron(np.eye(2),KetBra(1,0,2)))                      @(np.kron(KetBra(0,1,2),np.eye(2))+np.kron(np.eye(2),KetBra(0,1,2)))])

# First: calculate without phonons:
Simulation(fprop, PT, initial, tgrid, outp)
(times, data0) = outp.extract()

# Second: load phonon PT-MPOs. 
# The PT-MPOs is caluculate assuming a single two-level system. We can temporarily 
# extend the system Hilbert space. 
# For the first QD, we add a Hilbert space of dimension 2 behind the TLS Hilbert space 
# upon which the PT-MPO acts: 
PT.add_PT("QDPhonon.pt", 0, 2)
#For the second QD, we add the two-dimensional Hilbert space in front:
PT.add_PT("QDPhonon.pt", 2, 0)

outp.clear_results() #remove old results from outp
Simulation(fprop, PT, initial, tgrid, outp) #run again
(times, data1) = outp.extract()


# Now, plot the results
import matplotlib.pyplot as plt
plt.xlabel("Time")
plt.ylabel("Observable")
plt.title("Rabi rotations")
plt.plot(times, data0[:,0].real, label='G2 without phonons')
plt.plot(times, data1[:,0].real, label='G2 with phonons')
plt.legend(loc="upper right")
plt.show()



