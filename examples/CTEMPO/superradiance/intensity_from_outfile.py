import sys
import numpy as np
from ACE import read_outfile,write_outfile
_, fname, Tdecay = sys.argv


# This Python script calculates the normalized emitted intensity 
# from an outfile of CCTEMPO. To this end, run, e.g.
#
# > python3 intensity_from_outfile.py N8_T-1_Tdecay500_te1000_dt0.1_thr1e-7.out 500
#
# This generates a file "intensity_N8_T-1_Tdecay500_te1000_dt0.1_thr1e-7.out"
# With time and intensity as first and second columns

times, occ = read_outfile(fname)
def compute_Intensity(times, occ, Tdecay):
    return [np.array(times[:-1]), 
            np.array([-Tdecay * (occ[i+1].real-occ[i].real)/(times[i+1]-times[i]) for i in range(len(times)-1)])]


write_outfile('intensity_'+fname, compute_Intensity(times, occ, float(Tdecay)) )
