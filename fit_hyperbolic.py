# Fit_hyperolic v.1.2

import os, argparse
import h5py
from numpy import sqrt, linspace, abs, amax, amin
import matplotlib.pyplot as plt

from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plots a B-scan image plus hyberbolic fit', usage='cd gprMax; python -m tools.fit_hyperbolic mergedfile output -p x0 t0 v')
parser.add_argument('mergedfile', help='name of merged file including path')
parser.add_argument('output', help='output to be plotted (Ex, Ey, Ez, Hx, Hy, Hz, Ix, Iy, Iz)')
parser.add_argument('-p', '--parameters', type=float, help='x0 (m), t0 (ns), 0.3 > v(m/ns) > 0.03', nargs=3)
args = parser.parse_args()

# Read user's parameters 
mergedfile = args.mergedfile
parameters = args.parameters
x0 = parameters[0]
t0 = parameters[1]
v = parameters[2]

if not 0.3 >= v >= 0.03: 
	raise CmdInputError('{} (m/ns) is out of range'.format(v))

z0 = (t0*v)/2
v *= 1e09

# Open and read some attributes
f = h5py.File(mergedfile, 'r')
positions = f.attrs['merged_positions']

nrx = f.attrs['nrx']
nx = f.attrs['nx, ny, nz'][0]
dx = f.attrs['dx, dy, dz'][0] 
x = linspace(0, 1, nx) * (nx * dx)


# Fit hyperbolic
fit = (2 / v)  * sqrt((x0-x)**2 + z0**2)

# For all the receivers
for rx in range(1, nrx + 1):
	path = '/rxs/rx%s/' % (str(rx))
	availableoutputs = list(f[path].keys())
	outputdata = f[path + args.output]
	modelruns = outputdata.shape[1]
	xcoordonates = [ positions[model][rx-1][0] for model in range(modelruns) ]
	time = outputdata.shape[0]*f.attrs['dt']
	
	if t0 > amax(time)*1e9: 
		raise CmdInputError('{} (ns) is out of range'.format(t0))

	if not amax(xcoordonates) > x0 >  amin(xcoordonates): 
		raise CmdInputError('{} (m) is out of range'.format(x0))  

	# Plot Bscan plus fit hyperbolic
	fig = plt.figure(num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(outputdata, extent=[amin(xcoordonates), amax(xcoordonates), outputdata.shape[0]*f.attrs['dt'], 0], interpolation='bicubic', aspect='auto', cmap='seismic', vmin=-amax(abs(outputdata)), vmax=amax(abs(outputdata)))
	plt.plot(x, fit, color='r', linewidth=3.0)
	plt.xlabel('Distance [m]')
	plt.ylabel('Time [s]')
	plt.xlim(amin(xcoordonates), amax(xcoordonates))
	plt.ylim(amax(time), 0)
	plt.grid()

plt.show()



