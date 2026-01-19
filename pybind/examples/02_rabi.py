import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *


# We calculate Rabi oscillations in a two-level model. 
# Here, we make heavy use of python bindings

dt = 0.01
te = 100
Omega = 0.1

sigma_x = np.array([[0,1],[1,0]], dtype=complex)


# We need the free system propagator, a set of PT-MPOs, an initial state, a time grid 
# and an "OutputPrinter", which handles what observables are extracted and how they are written/returned
fprop   = FreePropagator()
fprop.add_Hamiltonian(hbar/2*Omega*sigma_x)

PT      = ProcessTensors() # Empty PT-MPO for now
initial = KetBra(0,0,2)    # <- KetBra(i,j,k) yields a numpy array representing {|i><j|_k}
outp    = OutputPrinter([KetBra(1,1,2), KetBra(1,0,2)])

# Finally, we define a simulation and run it
Simulation(fprop, PT, initial, TimeGrid(0, te, dt), outp)
# Extract data from OutputPrinter without using an output file
(times, data) = outp.extract()

print(data)

# Now, plot the results
import matplotlib.pyplot as plt
plt.gcf().set_size_inches((12,5))
plt.xlabel("Time")
plt.ylabel("Observable")
plt.title("Rabi rotations")
plt.plot(times, data[:,0].real, label='Occupation')
plt.plot(times, data[:,1].imag, label='Coherence')
plt.legend(loc="upper right")
plt.show()


