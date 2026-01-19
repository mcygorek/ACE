import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *

# The example calculates emission spectra from a strongly driven quantum dot coupled to phonons 
# by Fourier transforming the first-order coherence <sigma^+(t_op+tau) sigma^-(t_op)> 

# Parameters
ta = -2000      # start time
t_op = 0        # time at which first operator is applied
te = 2000       # end time 
dt = 0.1        # time step

t_mem = 12.8    # memory time for PT-MPO calculation
thr = 1e-10     # PT-MPO compression threshold
T = 4.2         # temperature (K)

gamma = 0.002        # radiative decay rate (ps^-1)
hbarOmega = 0.5      # cw driving: hbar times Rabi frequency (meV)


# Generate PT-MPO
plist =  [f'dt {dt}']
plist += [f'ta {ta}']
plist += [f'te {te}']
plist += [f't_mem {t_mem}']
plist += [f'threshold {thr}']
plist += [f'use_Gaussian_periodic true']
plist += [f'Boson_J_type QDPhonon']
plist += [f'Boson_omega_max 10']
plist += [f'Boson_temperature {T}']

PT = ProcessTensors(Parameters(plist))
#PT = ProcessTensors()   #<--- Use this instead to test phonon-free spectra

sigma_minus = KetBra(0,1,2)
sigma_plus = KetBra(1,0,2)
sigma_x = KetBra(0,1,2)+KetBra(1,0,2)

fprop = FreePropagator()
fprop.add_Hamiltonian(0.5*hbarOmega*sigma_x)        # cw driving
fprop.add_Lindblad(gamma, sigma_minus )             # rad. decay
fprop.apply_Operator_left(t_op, sigma_minus, True)  # apply sigma^- at time t_op

outp  = OutputPrinter( [sigma_plus] )      # evaluate <sigma^+ ...> at time t

# Run simulation:
Simulation(fprop, PT, KetBra(0,0,2), TimeGrid(ta, te, dt), outp)
(times, data) = outp.extract()


#Consider only coherences starting from t=t_op and subtract stationary value
zero_at = np.abs(times - t_op).argmin()  #find index closes to t=0
data_trunc=np.array([ x - data[-1,0] for x in data[zero_at:,0]]) # subtract stationary value and truncate


# Fourier transform:
spectrum = np.fft.fftshift(np.fft.fft(data_trunc))
freqs = np.fft.fftshift(np.fft.fftfreq(len(data_trunc),d=dt))

# Now, plot the results:
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,2,figsize=(15,5))
ax[0].set_title("Coherence")
ax[0].set(xlabel="Time (ps)")
ax[0].plot(times, data[:,0].real, label='Re')
ax[0].plot(times, data[:,0].imag, label='Im')
ax[0].legend(loc="upper right")
ax[1].set_title("Emission spectrum")
ax[1].set(xlabel="Energy (meV)")
ax[1].set_xlim([-3,3])
ax[1].set_yscale("log")
ax[1].plot(np.multiply(2*np.pi*hbar,freqs), spectrum.real)
plt.show()


