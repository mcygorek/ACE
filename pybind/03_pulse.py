import sys
sys.path.append('/home/.../ACE/pybind/') #<---plug in your directory

from ACEutils import *

#Prepare 

rho0 = np.array([[1,0],[0,0]], dtype=complex)
initial = InitialState( rho0 )

pulseTimes = np.arange(0, 100, 0.005, dtype=float)
pulseShape = ( pulseTimes, np.array([
               3*np.pi/(np.sqrt(2.*np.pi)*5)*np.exp(-0.5*((t-30)/5)**2)*(1+0j) 
                                                   for t in pulseTimes]))
dipole = hbar/2*np.array([[0,0],[1,0]], dtype=complex)

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
fprop = FreePropagator()
fprop.add_Pulse( pulseShape, dipole )

PT  = ProcessTensors()

outfile = "03_pulse.out"
outp  = OutputPrinter( outfile,
                      [np.array([[0,0],[0,1]],dtype=complex),
                       np.array([[0,0],[1,0]],dtype=complex)])

tgrid = TimeGrid(0, 100, 0.01)

sim   = Simulation()

# Run:
sim.run(fprop, PT, initial, tgrid, outp)

# Read:
(times, data) = read_outfile(outfile)
print(data)

# Now, plot the results:
import matplotlib.pyplot as plt
plt.xlabel("Time")
plt.ylabel("Observable")
plt.title("Gaussian 3*pi pulse")
plt.plot( pulseShape[0], pulseShape[1].real, label='Pulse')
plt.plot(times, data[:,0].real, label='Occupation')
plt.plot(times, data[:,1].imag, label='Coherence')
plt.legend(loc="upper right")
plt.show()

