import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *


# This time, we additionally couple the TLS to a spin-boson environment with Ohmic spectral density


te = 102.4  # end time 
dt =   0.1  # time step
A = 3 # pulse area (pi)

alpha, omega_c = 0.05, 10  # parameters of spectral density 
T = 0                      # temperature
thr = 1e-7                 # PT-MPO compression threshold


#Functions:
def myPulse(t) -> complex:
    return 1/(np.sqrt(2.*np.pi)*10)*np.exp(-0.5*((t-30)/10)**2)

def J(w):
    return alpha*w*np.exp(-w/omega_c)
    

pulseTimes = np.arange(0, te, dt/2, dtype=float)
pulseShape = ( pulseTimes, np.multiply(A*np.pi, myPulse(pulseTimes)))
dipole = hbar/2*KetBra(1,0,2)

fprop = FreePropagator()
fprop.add_Pulse( pulseShape, dipole )


# Discretize the spectral density and write it to a file of the proper ACE format
omega_list = np.arange(0, 200, 0.01, dtype=float)
Jfile='04_J.dat'
write_outfile(Jfile, (omega_list, J(omega_list))) 

# PT-MPO generation parameters are still based on ACE parameter files 
PT  = ProcessTensors(Parameters([f'dt {dt}', 
                                 f'te {te}', 
                                 'use_Gaussian_divide_and_conquer true',
                                 f'threshold {thr}', 
                                 f'Boson_J_from_file {Jfile}', 
                                 f'Boson_omega_max {200}', 
                                 f'Boson_temperature {T}'])) 

outfile = "04_pulse_ohmic.out"
# When a string is the first parameter, OutputPrinter will directly print to file 
# (function outp.extract() would return empty)
outp  = OutputPrinter( outfile, [KetBra(1,1,2), KetBra(1,0,2)])


Simulation(fprop, PT, KetBra(0,0,2), TimeGrid(0, te, dt), outp)
(times, data) = read_outfile(outfile)

# Now, plot the results:
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3,  figsize=(15, 4))
ax[0].set_xlabel(r'$\omega$ (ps)$^{-1}$')
ax[0].plot(omega_list, J(omega_list), color='black', label=r'$J(\omega)$')
ax[0].legend(loc="upper right")
ax[1].set_xlabel("Time (ps)")
ax[1].plot( pulseShape[0], pulseShape[1].real, 'r', label='Pulse')
ax[1].legend(loc="upper right")
ax[2].set_xlabel("Time (ps)")
ax[2].set_ylim([0,1])
ax[2].plot(times, data[:,0].real, label='Occupation')
ax[2].plot(times, data[:,1].imag, label='Coherence')
ax[2].legend(loc="upper right")
plt.show()

