# -v.1.0

import argparse, glob, os
import h5py
import numpy as np
from enum import Enum

from gprMax.receivers import Rx
from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='added a gain on fields values from a base file, then return the new output file.', usage='cd gprMax; python -m tools.gain basefilename initialtime (-exp)')
parser.add_argument('basefilename', help='base name of output file including path')
parser.add_argument('initialtime', help='time (ns) at which the gain is added')
parser.add_argument('-lin','--linear', action="store_true", help='chosen an linear gain type')
args = parser.parse_args()

basefilename = args.basefilename
tinit = float(args.initialtime)* 1e-9
outputfile = basefilename + '_gained.out'


# Create output file
fout = h5py.File(outputfile, 'w')

# Open base file and read some attributes
f = h5py.File(basefilename + '.out', 'r')

iterations = f.attrs['Iterations'] 
dt = f.attrs['dt'] 
nrx = f.attrs['nrx']
time = np.linspace(0, 1, iterations)
time *= (iterations * dt)
fgain = []
modelruns= 6
# Check there are any receivers
if nrx == 0:
    raise CmdInputError('No receivers found in {}'.format(args.outputfile))	

# Init datasets' output file
fout.attrs['Iterations'] = iterations
fout.attrs['dt'] = dt
fout.attrs['nrx'] = nrx

for model in range(modelruns):
	for rx in range(1, nrx + 1):

		# Write properties on first iteration
		if model == 0:
			path = '/rxs/rx' + str(rx)
			grp = fout.create_group(path)
			availableoutputs = list(f[path].keys())
			for output in availableoutputs:
				grp.create_dataset(output, (iterations, modelruns), dtype=f[path + '/' + output].dtype)
		
		# For all receiver outputs
		for output in availableoutputs:
			fout[path + '/' + output][:,model] = f[path + '/' + output][:,model]
		

# Create gain function chosen by user
if args.linear:
	flin = []
	for t in time:
		fgain.append(6*t*1e9)
	
else:
	fgain= np.ones(np.size(time))

# Find initial time idexation 
tindex=0
while time[tindex] <= tinit:
	tindex += 1 
rxindex = 0
# For each rx, add the gain on fields component values 
#for rxindex, rx in enumerate(rxs):
for model in range(modelruns):
		#if 'Ex' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Ex'][tindex:,model] *= fgain[tindex:]
		#if 'Ey' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Ey'][tindex:,model] *= fgain[tindex:]
		#if 'Ez' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Ez'][tindex:,model] *= fgain[tindex:]
		#if 'Hx' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Hx'][tindex:,model] *= fgain[tindex:]
		#if 'Hy' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Hy'][tindex:,model] *= fgain[tindex:]
		#if 'Hz' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Hz'][tindex:,model] *= fgain[tindex:]
		#if 'Ix' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Ix'][tindex:,model] *= fgain[tindex:]
		#if 'Iy' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Iy'][tindex:,model] *= fgain[tindex:]
		#if 'Iz' in rx.outputs:
	fout['/rxs/rx' + str(rxindex + 1) + '/Iz'][tindex:,model] *= fgain[tindex:]
	
f.close()

fout.close()


	
