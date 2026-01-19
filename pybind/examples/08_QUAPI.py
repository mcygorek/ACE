import sys
sys.path.append('.../ACE/pybind/') #<---plug in your directory
from ACEutils import *
from time import perf_counter


# Parameters
te = 200    # end time 
dt = 0.5   # time step
t_mem = 4  # For fair comparison: divide-and-conquer optimal for power of 2
thr = 1e-9
T = 4.2

Omega = 0.25

sigma_x = np.array([[0,1],[1,0]], dtype=complex)

fprop = FreePropagator()
fprop.add_Hamiltonian(hbar/2*Omega*sigma_x)

# PT-MPO Simulation 
perf_start = perf_counter()

plist =  [f'dt {dt}']
plist += [f'te {te}']
plist += [f't_mem {t_mem}']
plist += [f'threshold {thr}']
plist += [f'use_Gaussian_repeat true']
plist += [f'Boson_J_type QDPhonon']
plist += [f'Boson_omega_max 10']
plist += [f'Boson_temperature {T}']
param = Parameters(plist)
print(param)
PT = ProcessTensors(param)
#PT = ProcessTensors()

outp  = OutputPrinter( KetBra(1,1,2) )
Simulation(fprop, PT, KetBra(0,0,2), TimeGrid(0, te, dt), outp)
(times_PT, data_PT) = outp.extract()

perf_stop = perf_counter()
comptime_PT = perf_stop-perf_start




# QUAPI Simulation
perf_start = perf_counter()

t_memQ = 3.5   # For fair comparison: convergence is reached
plist_QUAPI =  [f'dt {dt}']
plist_QUAPI += [f'te {te}']
plist_QUAPI += [f't_mem {t_memQ}']
plist_QUAPI += [f'use_Gaussian true']
plist_QUAPI += [f'Boson_J_type QDPhonon']
plist_QUAPI += [f'Boson_omega_max 10']
plist_QUAPI += [f'Boson_temperature 4.2']
param_QUAPI =Parameters(plist_QUAPI)

IFQ = InfluenceFunctional_QUAPI(param_QUAPI)

outpQ  = OutputPrinter( KetBra(1,1,2) )

Simulation_QUAPI(fprop, IFQ, KetBra(0,0,2), TimeGrid(0, te, dt), outpQ)
(times_QUAPI, data_QUAPI) = outpQ.extract()

perf_stop = perf_counter()
comptime_QUAPI = perf_stop-perf_start



# TEMPO simulation
perf_start = perf_counter()

thr = 1e-7
plist_TEMPO =  [f'dt {dt}']
plist_TEMPO += [f'te {te}']
plist_TEMPO += [f't_mem {t_memQ}']
plist_TEMPO += [f'threshold {thr}']
plist_TEMPO += [f'use_Gaussian true']
plist_TEMPO += [f'Boson_J_type QDPhonon']
plist_TEMPO += [f'Boson_omega_max 10']
plist_TEMPO += [f'Boson_temperature 4.2']
param_TEMPO =Parameters(plist_TEMPO)

IFQ = InfluenceFunctional_TEMPO(param_TEMPO)

outpT  = OutputPrinter(KetBra(1,1,2))
Simulation_TEMPO(fprop, IFQ, KetBra(0,0,2), TimeGrid(0, te, dt), outpT, thr)
(times_TEMPO, data_TEMPO) = outpT.extract()

perf_stop = perf_counter()
comptime_TEMPO = perf_stop-perf_start

# Now, plot the results:
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1,figsize=(12,5))
ax.set(xlabel="Time (ps)")
ax.set(ylabel="Occupations")
ax.plot(times_PT, data_PT[:,0].real, 'x', label='PT')
ax.plot(times_QUAPI, data_QUAPI[:,0].real, '+', label='QUAPI')
ax.plot(times_TEMPO, data_TEMPO[:,0].real, label='TEMPO')
ax.legend(loc="upper right")
plt.show()



print(f'Time for PT-MPO simulation: {comptime_PT}s')
print(f'Time for QUAPI simulation: {comptime_QUAPI}s')
print(f'Time for TEMPO simulation: {comptime_TEMPO}s')
