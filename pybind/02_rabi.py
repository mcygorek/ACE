import sys
sys.path.append('/home/.../ACE/pybind/') #<---plug in your directory

from ACEutils import *

#Prepare 
plist = []
param = Parameters(plist)

rho0 = np.array([[1,0],[0,0]], dtype=complex)
initial = InitialState( rho0 )

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
fprop = FreePropagator()
fprop.add_Hamiltonian(hbar/2*0.1*sigma_x)

PT  = ProcessTensors(param)

outfile = "02_rabi.out"
outp  = OutputPrinter( outfile,
                      [np.array([[0,0],[0,1]],dtype=complex),
                       np.array([[0,0],[1,0]],dtype=complex)])

tgrid = TimeGrid(0, 100, 0.01)

sim   = Simulation(param)
# Run:
sim.run(fprop, PT, initial, tgrid, outp)


(times, data) = read_outfile(outfile)
print(data)

# Now, plot the results
import matplotlib.pyplot as plt
plt.xlabel("Time")
plt.ylabel("Observable")
plt.title("Rabi rotations")
plt.plot(times, data[:,0].real, label='Occupation')
plt.plot(times, data[:,1].imag, label='Coherence')
plt.legend(loc="upper right")
plt.show()

