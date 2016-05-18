# AGC_gain -v.1.0

import argparse
import h5py
import numpy as np

from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='create a new ouputfile gained by user defined type {constant, linear, exponential}')
parser.add_argument('filename', help='base name of output file including path')
parser.add_argument('gatelenght', help='lengh (ns) of the AGC time gate')
args = parser.parse_args()

filename = args.filename
outputfile = args.filename.replace('.out', '_AGC.out')

# Open filename and read some attributes
f = h5py.File(filename, 'r')
nrx = f.attrs['nrx']
iterations = f.attrs['Iterations']
dt = f.attrs['dt']

# Write new outputfile
fg = h5py.File(outputfile, 'w')
fg.attrs['Iterations'] = iterations
fg.attrs['dt'] = dt
fg.attrs['nrx'] = nrx
time = np.linspace(0, 1, iterations) * (iterations * dt)


