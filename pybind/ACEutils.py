import csv
import numpy as np
from ACE import *

ACE_sigma_p=np.array([[0,0],[1,0]],dtype=complex)
ACE_sigma_m=np.array([[0,1],[0,0]],dtype=complex)
ACE_sigma_x=np.array([[0,1],[1,0]],dtype=complex)
ACE_sigma_y=np.array([[0,1j],[-1j,0]],dtype=complex)
ACE_sigma_z=np.array([[-1,0],[0,1]],dtype=complex)

def Boson_create(dim):
    O = np.zeros((dim,dim),dtype=complex)
    for i in range(dim-1):
        O[i+1,i]=np.sqrt(i+1)
    return O

def Boson_destroy(dim):
    O = np.zeros((dim,dim),dtype=complex)
    for i in range(dim-1):
        O[i,i+1]=np.sqrt(i+1)
    return O

def Boson_n(dim):
    O = np.zeros((dim,dim),dtype=complex)
    for i in range(dim):
        O[i,i]=i
    return O

def Boson_vacuum(dim):
    O = np.zeros((dim,dim),dtype=complex)
    O[0,0]=1
    return O

def KetBra(i,j,dim):   # |i><j|_dim
    O = np.zeros((dim,dim),dtype=complex)
    O[i,j]=1
    return O


def run(parameter_list):
    param = Parameters(parameter_list)
    fprop = FreePropagator(param)
    PT    = ProcessTensors(param)
    initial = InitialState(param)
    outp  = OutputPrinter(param)
    tgrid = TimeGrid(param)
    sim   = Simulation(param)
    sim.run(fprop, PT, initial, tgrid, outp)

    
def read_outfile(outfile, ncol=0):
    if ncol<1:  #read number of column pairs in the outfile 
      with open(outfile,'r') as fil:
        nreal = 0
        data_iter = csv.reader(fil, delimiter = ' ')
        for line in data_iter:
            for i in range(len(line)):
                if line[i].startswith('#') or len(line[i])<1: break
                nreal += 1
            if nreal > 0: break
        ncol = int((nreal-1)/2)
    
    with open(outfile,'r') as fil:
        data_iter = csv.reader(fil, delimiter = ' ')
        times = []
        values = []
        for line in data_iter:
            if len(line)<1 or line[0].startswith('#') or len(line[0])<1: 
              continue
            times.append(float(line[0]))
            val_line = []
            for col in range(ncol):
                val_line.append(complex(float(line[1+2*col]),float(line[2+2*col])))
            values.append(val_line)
        return (np.array(times, dtype=float), np.array(values, dtype=complex))
        
        
def write_outfile(filename, data):
    if data[0].shape[0]!=data[1].shape[0]:
      print('data[0].shape[0]!=data[1].shape[0]')
      exit()
    if len(data[1].shape)==1:
      with open(filename, 'w') as fil:
        writer = csv.writer(fil, delimiter=' ')
        for i in range(data[0].shape[0]):
          string = [data[0][i], data[1][i].real, data[1][i].imag]
          writer.writerow(string)
    else:
      with open(filename, 'w') as fil:
        writer = csv.writer(fil, delimiter=' ')
        for i in range(data[0].shape[0]):
          string = []
          string.append(data[0][i])
          for j in range(data[1].shape[1]):
            string.append(data[1][i,j].real)
            string.append(data[1][i,j].imag)
          writer.writerow(string)


