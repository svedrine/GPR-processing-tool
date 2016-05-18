# user_defined_gain -v.1.0

import argparse
import h5py
import numpy as np

from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='create a new ouputfile gained by user defined type {constant, linear, exponential}')
parser.add_argument('filename', help='base name of output file including path')
parser.add_argument('type', choices = ('constant', 'linear', 'exponential'))
parser.add_argument('-window', type=float, nargs=2)
args = parser.parse_args()

filename = args.filename
outputfile = args.filename.replace('.out', '_%s.out' % args.type)

# Open filename and read some attributes
f = h5py.File(filename, 'r')
nrx = f.attrs['nrx']
iterations = f.attrs['Iterations']
dt = f.attrs['dt']
modelruns = f.attrs['Modelruns']

# Write new outputfile
fg = h5py.File(outputfile, 'w')
fg.attrs['Iterations'] = iterations
fg.attrs['dt'] = dt
fg.attrs['nrx'] = nrx
fg.attrs['Modelruns'] = modelruns
fg.attrs['Positions'] = f.attrs['Positions']

time = np.linspace(0, 1, iterations) * (iterations * dt)

if '_merged' in filename:
	fg.attrs['nx, ny, nz'] = f.attrs['nx, ny, nz']
	fg.attrs['dx, dy, dz'] = f.attrs['dx, dy, dz']
	fg.attrs['merged_positions'] = f.attrs['merged_positions']

# Find window indexation
index = 0
while time[index] <= args.window[0]*1e-9: 
	index += 1
imin = index

while time[index] <= args.window[1]*1e-9:
	index += 1
imax = index
	
# Write the gain function
width = imax - imin

if args.type == 'constant':
	check = input('Please enter a constant value :')
	constant = float(check)
	fgain = [constant]*width
		

if args.type == 'linear':
	check = input('Please enter gradient value: ')
	gradient = float(check)
	fgain = [gradient*1e9*t for t in time[imin:imax]]

if args.type == 'exponential':
	check = input('fgain = A * exp(B*t). Please enter A, B values: ')
	A = float(check[0])
	B = float(check[2])
	fgain = [A*np.exp(B*t) for t in time[imin:imax]]

for rx in range(1, nrx + 1):

	path = '/rxs/rx%s/' % str(rx)
	availableoutputs = f[path]

	# Write properties for all availiable outputs
	grp = fg.create_group(path)
			
	if modelruns:
		for output in availableoutputs:
			for model in range(modelruns):
				# Initialize
				if model ==0: 
					grp.create_dataset(output, (iterations, modelruns), dtype=f['%s%s' %(path, output)].dtype)
				fg['%s%s' %(path, output)][:, model] = f['%s%s' %(path, output)][:,model]
				# Apply gain
				fg['%s%s' %(path, output)][imin:imax, model] *= fgain
	
	else:
		for output in availableoutputs:
				# Initialize
				grp.create_dataset(output, (f.attrs['Iterations'],), dtype=f['%s%s' %(path, output)].dtype)
				fg['%s%s' %(path, output)][:,] = f['%s%s' %(path, output)][:,]
				# Apply gain
				fg['%s%s' %(path, output)][imin:imax,] *= fgain




	
