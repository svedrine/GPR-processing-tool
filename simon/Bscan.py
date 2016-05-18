# Bscan v.2.0

import os, argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

from gprMax.exceptions import CmdInputError

"""Plots a B-scan image."""

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plots a B-scan image.')
parser.add_argument('outputfile', help='name of output file including path')
parser.add_argument('output', help='name of output component to be plotted (Ex, Ey, Ez, Hx, Hy, Hz, Ix, Iy or Iz)')
parser.add_argument('-gain', choices = ['lin', 'exp'])
parser.add_argument('-p', '-parameters', type= float, help= '(t0, coef)', nargs = 2)
args = parser.parse_args()

# Open output file and read some attributes
f = h5py.File(args.outputfile, 'r')
nrx = f.attrs['nrx']
dt = f.attrs['dt']
iterations = f.attrs['Iterations']
time = np.linspace(0, 1, iterations)
time *= (iterations * dt)

# Check there are any receivers
if nrx == 0:
    raise CmdInputError('No receivers found in {}'.format(args.outputfile))

# For all the receivers
for rx in range(1, nrx + 1):
	path = '/rxs/rx%s/' % str(rx)
	availableoutputs = list(f[path].keys())
	modelruns = f['%s%s' % (path, args.output)].shape[1]

    # Check if requested output is in file
	if args.output not in availableoutputs:
		raise CmdInputError('{} output requested to plot, but the available output for receiver 1 is {}'.format(args.output, ', '.join(availableoutputs)))
	
	# Check that there is more than one A-scan present
	if modelruns == 1:
		raise CmdInputError('{} contains only a single A-scan.'.format(args.outputfile))
	
	# If Gain is turned on then uptade outputdata
	if args.gain:
		
		# Create output file on first itteration
		if rx == 1:
			outputfile = args.outputfile.replace('.out', '_gained.out')
			fg = h5py.File(outputfile, 'w')
		
			# Init datasets' output file
			fg.attrs['Iterations'] = iterations
			fg.attrs['dt'] = dt
			fg.attrs['nrx'] = nrx
			fg.attrs['nx, ny, nz'] = f.attrs['nx, ny, nz']
			fg.attrs['dx, dy, dz'] = f.attrs['dx, dy, dz']
			fg.attrs['merged_positions'] = f.attrs['merged_positions']

			# Find initial time idexation 
			tindex = 0
			tinit = args.p[0]*1e-9
			while time[tindex] <= tinit:
				tindex += 1

			# Write the gain function to be applied
			fgain = [1]*len(time[:tindex])
			coef = args.p[1]*1e9
			if args.gain == 'lin':
				for t in time[tindex:]:
					fgain.append(coef * (t-tinit) + 1) 

		for model in range(modelruns):
			# Write properties on first iteration
			if model == 0:
				grp = fg.create_group(path)
				for output in availableoutputs:
					grp.create_dataset(output, (iterations, modelruns), dtype=f['%s%s' % (path, output)].dtype)
					print(f['%s%s' % (path, output)].dtype)
		
			# For all models
			for output in availableoutputs:
				# Initialize
				fg['%s%s' % (path, output)][:,model] = f['%s%s' % (path, output)][:,model]
				# Apply gain
				fg['%s%s' % (path, output)][tindex:,model] *= fgain[tindex:]
	
		outputdata = fg['%s%s' % (path, args.output)]
	
	else:
		outputdata = f['%s%s' % (path, args.output)]
		print(outputdata)

    # Plot B-scan image
	fig = plt.figure(num='rx' + str(rx), figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(outputdata, extent=[0, outputdata.shape[1], outputdata.shape[0]*f.attrs['dt']*1e9, 0], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(np.abs(outputdata)), vmax=np.amax(np.abs(outputdata)))
	plt.xlabel('Trace number')
	plt.ylabel('Time [ns]')
	plt.grid()
	cb = plt.colorbar()
	if 'E' in args.output:
		cb.set_label('Field strength [V/m]')
	elif 'H' in args.output:
		cb.set_label('Field strength [A/m]')
	elif 'I' in args.output:
		cb.set_label('Current [A]')

    # Save a PDF/PNG of the figure
    #fig.savefig(os.path.splitext(os.path.abspath(args.outputfile))[0] + '.pdf', dpi=None, format='pdf', bbox_inches='tight', pad_inches=0.1)
    #fig.savefig(os.path.splitext(os.path.abspath(args.outputfile))[0] + '.png', dpi=150, format='png', bbox_inches='tight', pad_inches=0.1)

plt.show()
