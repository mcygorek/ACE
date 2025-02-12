import sys
#sys.path.append('/home/.../ACE/pybind/') #<---plug in your directory

from ACEutils import *

#This function will process the parameter list and run ACE:
def run_ACE(parameter_list, do_extract_densmat=False):
    param = Parameters(parameter_list)
    fprop = FreePropagator(param)
    PT    = ProcessTensors(param)
    initial = InitialState(param)
    outp  = OutputPrinter(param)
    if do_extract_densmat:
      outp.do_extract=True

    tgrid = TimeGrid(param)
    sim   = Simulation(param)
    sim.run(fprop, PT, initial, tgrid, outp)

    times = tgrid.get_all()
    densmats= np.array(outp.extract())
    return (times, densmats)

# Let's define the parameters

plist = ["dt 0.01", "te 100", "add_Hamiltonian {hbar/2*0.1*sigma_x}"]
plist += ["add_Output |1><1|_2", "add_Output |1><0|_2"]

outfile = "01_rabi.out"
plist += [f"outfile {outfile}"]

# Now run; write output to file
#(times, densmats) = run_ACE(plist, True) #<-This would also extract the complete density matrix
run_ACE(plist)

# Read data from outfile again
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

