import os, argparse
import h5py
from numpy import sqrt, linspace, abs, amax, amin
import matplotlib.pyplot as plt

from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plots a B-scan image plus hyberbolic fit ', usage='cd gprMax; python -m tools.fit_hyperbolic mergedfile output -p x0 t0 v')
parser.add_argument('mergedfile', help='name of merged file including path')
parser.add_argument('output', help='output to be plotted (Ex, Ey, Ez, Hx, Hy, Hz, Ix, Iy, Iz)')
parser.add_argument('-p', '--parameters', type=float, help='x0 (m), t0 (ns), v(m/ns)', nargs=3)
args = parser.parse_args()

mergedfile = args.mergedfile
parameters = args.parameters
x0 = parameters[0]
t0 = parameters[1]
v = parameters[2]
z0 = (t0*v)/2
v *= 1e09

f = h5py.File(mergedfile, 'r')


positions = f.attrs['merged_positions']
xcoordonates = list()




nrx = f.attrs['nrx']
nx = f.attrs['nx, ny, nz'][0]
dx = f.attrs['dx, dy, dz'][0] 
x = linspace(0, 1, nx)
x *= (nx * dx)

t = (2 * sqrt((x0-x)**2 + z0**2))/v 

for rx in range(1, nrx + 1):
	path = '/rxs/rx' + str(rx) + '/'
	availableoutputs = list(f[path].keys())
	outputdata = f[path + '/' + args.output]
	modelruns = f[path + args.output].shape[1]
	for model in range(modelruns):
		xcoordonates.append(positions[model][rx-1][0])
	
	fig = plt.figure(num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(outputdata, extent=[amin(xcoordonates), amax(xcoordonates), outputdata.shape[0]*f.attrs['dt'], 0], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-amax(abs(outputdata)), vmax=amax(abs(outputdata)))
	plt.plot(x,t)
	plt.xlabel('Distance [m]')
	plt.ylabel('Time [s]')
	plt.xlim(amin(xcoordonates), amax(xcoordonates))
	plt.grid()

plt.show()



