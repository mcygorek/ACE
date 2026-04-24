"""
ACE utility functions and convenience helpers.

Import as:
    from ACE.utils import run, read_outfile, write_outfile
    from ACE.utils import Boson_create, Boson_destroy, Boson_n, Boson_vacuum
    from ACE.utils import ACE_sigma_p, ACE_sigma_m, ACE_sigma_x, ACE_sigma_y, ACE_sigma_z
"""

import csv
import numpy as np

# ── Operators imported from the parent package (avoids circular import) ───────
from . import _ACE as _c
Parameters     = _c.Parameters
FreePropagator = _c.FreePropagator
ProcessTensors = _c.ProcessTensors
InitialState   = _c.InitialState
TimeGrid       = _c.TimeGrid
OutputPrinter  = _c.OutputPrinter
Simulation     = _c.Simulation

# ── Common sigma matrices ─────────────────────────────────────────────
ACE_sigma_p = np.array([[0, 0], [1, 0]], dtype=complex)
ACE_sigma_m = np.array([[0, 1], [0, 0]], dtype=complex)
ACE_sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
ACE_sigma_y = np.array([[0, 1j], [-1j, 0]], dtype=complex)
ACE_sigma_z = np.array([[-1, 0], [0, 1]], dtype=complex)


def Boson_create(dim):
    O = np.zeros((dim, dim), dtype=complex)
    for i in range(dim - 1):
        O[i + 1, i] = np.sqrt(i + 1)
    return O


def Boson_destroy(dim):
    O = np.zeros((dim, dim), dtype=complex)
    for i in range(dim - 1):
        O[i, i + 1] = np.sqrt(i + 1)
    return O


def Boson_n(dim):
    O = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        O[i, i] = i
    return O


def Boson_vacuum(dim):
    O = np.zeros((dim, dim), dtype=complex)
    O[0, 0] = 1
    return O


def KetBra(i, j, dim):
    """Return |i><j| in a dim-dimensional space."""
    O = np.zeros((dim, dim), dtype=complex)
    O[i, j] = 1
    return O


# ── High-level run helper ─────────────────────────────────────────────────────
def run(parameter_list):
    """Run a simulation from a list of parameter strings and return results."""
    param   = Parameters(parameter_list)
    fprop   = FreePropagator(param)
    PT      = ProcessTensors(param)
    initial = InitialState(param)
    outp    = OutputPrinter(param)
    tgrid   = TimeGrid(param)
    sim     = Simulation(param)
    sim.run(fprop, PT, initial, tgrid, outp)


# ── File I/O ──────────────────────────────────────────────────────────────────
def read_outfile(outfile, ncol=0):
    """Read an ACE output file, returning (times, values) as numpy arrays."""
    if ncol < 1:
        with open(outfile, 'r') as f:
            nreal = 0
            for line in csv.reader(f, delimiter=' '):
                for item in line:
                    if item.startswith('#') or len(item) < 1:
                        break
                    nreal += 1
                if nreal > 0:
                    break
        ncol = int((nreal - 1) / 2)

    times, values = [], []
    with open(outfile, 'r') as f:
        for line in csv.reader(f, delimiter=' '):
            if not line or line[0].startswith('#') or not line[0]:
                continue
            times.append(float(line[0]))
            values.append([
                complex(float(line[1 + 2 * c]), float(line[2 + 2 * c]))
                for c in range(ncol)
            ])
    return np.array(times, dtype=float), np.array(values, dtype=complex)


def write_outfile(filename, data):
    """Write (times, values) to an ACE-format output file."""
    times, values = data
    if times.shape[0] != values.shape[0]:
        raise ValueError('times and values must have the same length')
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        if values.ndim == 1:
            for i in range(times.shape[0]):
                writer.writerow([times[i], values[i].real, values[i].imag])
        else:
            for i in range(times.shape[0]):
                row = [times[i]]
                for j in range(values.shape[1]):
                    row += [values[i, j].real, values[i, j].imag]
                writer.writerow(row)

def MatrixToString(mat, precision=5):
  threshold=10**(-precision)/10.*0.999
  if mat.shape[0] != mat.shape[1]:
    raise ValueError(f'mat.shape[0] != mat.shape[1]')   
  d = mat.shape[0]
  str = ''
  for i in range(d):
    for j in range(d):
      if np.abs(mat[i,j])>threshold:
        if str != '':
          str += '+'
        if mat[i,j].imag>threshold:
          str += f'({mat[i,j].real:.{precision}f}+{mat[i,j].imag:.{precision}f}*i)'
        elif mat[i,j].imag<-threshold:
          str += f'({mat[i,j].real:.{precision}f}{mat[i,j].imag:.{precision}f}*i)'
        else:
          str += f'({mat[i,j].real:.{precision}f})'
        str += f'*|{i}><{j}|_{d}'
  if str == '':
    str = f'0*|0><0|_{d}'
  return str
