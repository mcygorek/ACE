import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *


# We calculate Rabi oscillations in a two-level model. 
# Here, we do this in a non-pythonic way by specifying the parameters 
# (as would be set in parameter files) as a list of strings

dt = 0.01
te = 100
Omega = 0.1
outfile = "01_rabi.out"

plist = [f'dt {dt}', f'te {te}', f'add_Hamiltonian {{hbar/2*{Omega}*sigma_x}}']
plist += ['add_Output {|1><1|_2}', 'add_Output {|1><0|_2}']
plist += [f'outfile {outfile}']
print('Parameters:')
[print(p) for p in plist]
print('')
param = Parameters(plist)

# We need the free system propagator, a set of PT-MPOs, an initial state, a time grid 
# and an "OutputPrinter", which handles what observables are extracted and how they are written/returned
fprop   = FreePropagator(param)
PT      = ProcessTensors(param)
initial = InitialState(param)
outp    = OutputPrinter(param)
tgrid   = TimeGrid(param)

# Finally, we define a simulation and run it
sim     = Simulation(param)
sim.run(fprop, PT, initial, tgrid, outp)

# Read data from outfile again using tools in ACEutils
(times, data) = read_outfile(outfile)

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

