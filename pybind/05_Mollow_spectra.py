import sys
sys.path.append('/..../ACE/pybind/') #<---plug in your directory

from ACEutils import *

# The example calculates emission spectra from strongly driven quantum dots coupled to phonons 
# by Fourier transforming the first-order coherence <sigma^+(t_op+tau) sigma^-(t_op)> 

# Parameters
ta = -2000      # start time
t_op = 0        # time at which first operator is applied
te = 2000    # end time 
dt = 0.1   # time step
t_mem = 6.4
thr = 1e-9

gamma = 0.002
hbarOmega = 0.5
hbar = 0.6582119569

plist =  [f'dt {dt}']
plist += [f'ta {ta}']
plist += [f'te {te}']
plist += [f't_mem {t_mem}']
plist += [f'threshold {thr}']
plist += [f'use_Gaussian_periodic true']
plist += [f'Boson_J_type QDPhonon']
plist += [f'Boson_omega_max 10']
plist += [f'Boson_temperature 4.2']
param =Parameters(plist)
print(param)


rho0 = np.zeros((2,2), dtype=complex)
rho0[1,1] = 1
initial = InitialState( rho0 )

sigma_plus = np.zeros((2,2), dtype=complex)
sigma_plus[1,0] = 1
sigma_minus = np.zeros((2,2), dtype=complex)
sigma_minus[0,1] = 1
sigma_x = np.array([[0,1],[1,0]], dtype=complex)


fprop = FreePropagator(param)
fprop.add_Hamiltonian(0.5*hbarOmega*sigma_x)
fprop.add_Lindblad(gamma, sigma_minus)
fprop.apply_Operator_left(t_op, sigma_minus, True)

#PT = ProcessTensors()   #<--- Use this to test phonon-free spectra
PT = ProcessTensors(param)

outfile = "05_Mollow_spectra.out"
outp  = OutputPrinter( outfile,
                      [sigma_plus])

tgrid = TimeGrid(ta, te, dt)

sim   = Simulation()

# Run:
sim.run(fprop, PT, initial, tgrid, outp)

# Read:
(times, data) = read_outfile(outfile)
print(data)


#Consider only coherences starting from t=t_op and subtract stationary value
zero_at = np.abs(times - t_op).argmin()  #find index closes to t=0
data_trunc=np.array([ x - data[-1,0] for x in data[zero_at:,0]])


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


