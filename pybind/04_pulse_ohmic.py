import sys
sys.path.append('/home/.../ACE/pybind/') #<---plug in your directory

from ACEutils import *

#Prepare 

ta =   0    # start time
te = 102.4  # end time 
dt =   0.1  # time step

def myPulse(t) -> complex:
    return 1*np.pi/(np.sqrt(2.*np.pi)*10)*np.exp(-0.5*((t-30)/10)**2)

def J(w):
    return 0.05*w*np.exp(-w/10)
    

rho0 = np.array([[1,0],[0,0]], dtype=complex)
initial = InitialState( rho0 )

pulseTimes = np.arange(ta, te, dt/2, dtype=float)
pulseShape = ( pulseTimes, myPulse(pulseTimes))
dipole = hbar/2*np.array([[0,0],[1,0]], dtype=complex)

fprop = FreePropagator()
fprop.add_Pulse( pulseShape, dipole )


omega_list = np.arange(0, 200, 0.01, dtype=float)
write_outfile("04_J.dat", (omega_list, J(omega_list))) 
PT  = ProcessTensors(Parameters([f"dt {dt}", f"te {te}", "use_Gaussian_divide_and_conquer true", "threshold 1e-7", "Boson_J_from_file 04_J.dat", "Boson_omega_max 200", "Boson_temperature 0"])) 

outfile = "04_pulse_ohmic.out"
outp  = OutputPrinter( outfile,
                      [np.array([[0,0],[0,1]],dtype=complex),
                       np.array([[0,0],[1,0]],dtype=complex)])

tgrid = TimeGrid(ta, te, dt)

sim   = Simulation()

# Run:
sim.run(fprop, PT, initial, tgrid, outp)

# Read:
(times, data) = read_outfile(outfile)
print(data)

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

