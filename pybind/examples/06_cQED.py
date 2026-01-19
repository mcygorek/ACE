import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *

# This example considers phonon-assisted state preparation of a QD strongly coupled to a microcavity


# Parameters
te = 100    # end time 
dt = 0.1   # time step
t_mem = 6.4  
thr = 1e-8
T = 4.2
a_e = 3 #<- Electron radius (used in spectral density)
bdim = 3 #<- Truncated dimension of cavity mode

g     = 0.05/hbar  # QD-cavity coupling (meV -> ps^-1)
gamma = 0.02/hbar  # radiative decay
kappa = 0.05/hbar  # cavity loss

# Gaussian pulse 
def myPulse(t) -> complex:
    A     = 13       # pulse area (pi)
    tc    = 30       # pulse center (ps)
    FWHM  = 7        # pulse duration (ps)
    delta = 1.1/hbar # Detuning (meV -> ps^-1)
    sigma = FWHM/(2*np.sqrt(2*np.log(2)))
    return (A*np.pi)/(np.sqrt(2.*np.pi)*sigma)*np.exp(-0.5*((t-tc)/sigma)**2) * np.exp(-1j*delta*(t-tc))


pulseTimes = np.arange(0, te, dt/2, dtype=float)
pulseShape = ( pulseTimes, myPulse(pulseTimes) )
dipole = hbar/2*np.kron( KetBra(1,0,2), np.eye(bdim) )   # <- Kronecker product with identity matrix

H_JC = hbar*g*(np.kron(KetBra(0,1,2), Boson_create(bdim)) +   # <- Jaynes-Cummings coupling
               np.kron(KetBra(1,0,2), Boson_destroy(bdim)) )
    
# Define free system propagator and OutputPrinter
fprop   = FreePropagator()
fprop.add_Pulse( pulseShape, dipole )
fprop.add_Hamiltonian( H_JC )
fprop.add_Lindblad( gamma, np.kron(KetBra(0,1,2), np.eye(bdim)) )
fprop.add_Lindblad( kappa, np.kron(np.eye(2), Boson_destroy(bdim)) )

# Calculate PT-MPO
plist =  [f'dt {dt}']
plist += [f'te {te}']
plist += [f't_mem {t_mem}']
plist += [f'threshold {thr}']
plist += [f'use_Gaussian_repeat true']
plist += [f'Boson_J_type QDPhonon']
plist += [f'Boson_omega_max 10']
plist += [f'Boson_temperature {T}']
plist += [f'Boson_J_a_e {a_e}']
plist += [f'Boson_SysOp {{|1><1|_2 otimes Id_{bdim} }}'] # <- coupling Hamiltonian has to match system dimension
param = Parameters(plist)
print(param)
PT = ProcessTensors(param)

outp  = OutputPrinter( [np.kron( KetBra(1,1,2), KetBra(0,0,bdim) ),    # Observables: |e,0> occupations
                        np.kron( KetBra(0,0,2), KetBra(1,1,bdim) ),    # |g,1> occupations
                        np.kron( KetBra(1,1,2), KetBra(1,1,bdim) ),    # |e,1> occupations
                        np.kron( KetBra(0,0,2), KetBra(2,2,bdim) )])   # |g,2> occupations

Simulation(fprop, PT, np.kron( KetBra(0,0,2), KetBra(0,0,bdim) ), TimeGrid(0, te, dt), outp)
(times, data) = outp.extract()



# Now, plot the results:
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1,figsize=(12,8))
ax.set(xlabel="Time (ps)")
ax.set(ylabel="Occupations")
ax.fill(pulseShape[0], np.abs(pulseShape[1])/(np.pi), 'pink', label='Pulse')
ax.plot(times, data[:,0].real, label='|e,0>')
ax.plot(times, data[:,1].real, label='|g,1>')
ax.plot(times, data[:,2].real, label='|e,1>')
ax.plot(times, data[:,3].real, label='|g,2>')
ax.legend(loc="upper right")
plt.show()


