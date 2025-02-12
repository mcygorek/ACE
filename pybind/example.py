import sys
#sys.path.append('/.../ACE/pybind/') #set to the corresponding directory
from ACE import *

def run_ACE(parameter_list):
    param = Parameters(parameter_list)
    fprop = FreePropagator(param)
    PT    = ProcessTensors(param)
    initial = InitialState(param)
    outp  = OutputPrinter(param)
    tgrid = TimeGrid(param)
    sim   = Simulation(param)
    sim.run(fprop, PT, initial, tgrid, outp)


plist = ["dt   0.01", "te 100", "add_Hamiltonian {hbar/2*0.1*sigma_x}"]
run_ACE(plist)


