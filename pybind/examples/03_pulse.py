import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *

# We now drive the TLS using a laser with a pulse shape defined using Python

# Parameters:
te = 100    # end time 
dt = 0.01   # time step
A = 3 # pulse area (pi)


#Pulse function. Set output type to complex to generate the correct type of numpy array
def myPulse(t) -> complex:
    return 1/(np.sqrt(2.*np.pi)*5)*np.exp(-0.5*((t-30)/5)**2)

# We discretize the pulse function on a regular grid
pulseTimes = np.arange(0, te, dt/2, dtype=float)
# The pulse shape consists of time points and function values. We rescale the latter by A*np.pi
pulseShape = ( pulseTimes, np.multiply(A*np.pi, myPulse(pulseTimes)))
# Dipole operator d. 
dipole = hbar/2*KetBra(1,0,2)
# Note that we specify only f(t)*d. The hermitian conjugate f^*(t)*d^\dagger is added automatically
    
# Define free system propagator and OutputPrinter
fprop   = FreePropagator()
fprop.add_Pulse( pulseShape, dipole )

outp    = OutputPrinter([KetBra(1,1,2), KetBra(1,0,2)])

# Finally, we define a simulation and run it
Simulation(fprop, ProcessTensors(), KetBra(0,0,2) , TimeGrid(0, te, dt), outp)
# Extract data from OutputPrinter without using an output file
(times, data) = outp.extract()

print(data)

# Now, plot the results:
import matplotlib.pyplot as plt
plt.gcf().set_size_inches((12,5))
plt.xlabel("Time")
plt.ylabel("Observable")
plt.title("Gaussian 3*pi pulse")
plt.plot( pulseShape[0], pulseShape[1].real, label='Pulse')
plt.plot(times, data[:,0].real, label='Occupation')
plt.plot(times, data[:,1].imag, label='Coherence')
plt.legend(loc="upper right")
plt.show()

