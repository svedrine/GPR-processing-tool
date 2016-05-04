# -v.1.0

import argparse, glob, os
import h5py
import numpy as np
from enum import Enum

# Parse command line arguments
parser = argparse.ArgumentParser(description='added a gain on fields values from a base file, then return the new output file.', usage='cd gprMax; python -m tools.gain basefile initialtime (-exp)')
parser.add_argument('basefile', help='base name of output file including path')
parser.add_argument('initialtime', help='time (ns) at which the gain is added')
parser.add_argument('-exp','--exponential', action="store_true", help='chosen an exponential gain type')
args = parser.parse_args()

basefile = args.basefile
tinit = args.initialtime * 10e-9
outputfile = basefile + '_gained.out'
print(outputfile)

# Create output file
fout = h5py.File(outputfile, 'w')

# Open base file and read some attributes
f = h5py.File(basefilename, 'r+')

itterations = fin.attrs['Iterations'] 
dt = fin.attrs['dt'] 
nrx = fin.attrs['nrx']
time = np.linspace(0, 1, iterations)
time *= (iterations * dt)

fgain = []
	
# Create gain function chosen by user
if args.exponential:

	for t in time:
		fgain.append(np.exp(t-tinit))
else:
	fgain= ones(dim(time))
	

# Find initial time idexation 
tindex=0
while time[tindex] <= tinit:
	tindex += 1 


# For each rx, add the gain on fields component values 
for rxindex, rx in enumerate(rxs):
		
	if 'Ex' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Ex'][tindex:] *= fgain[tindex:]
	if 'Ey' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Ey'][tindex:] *= fgain[tindex:]
	if 'Ez' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Ez'][tindex:] *= fgain[tindex:]
	if 'Hx' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Hx'][tindex:] *= fgain[tindex:]
	if 'Hy' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Hy'][tindex:] *= fgain[tindex:]
	if 'Hz' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Hz'][tindex:] *= fgain[tindex:]
	if 'Ix' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Ix'][tindex:] *= fgain[tindex:]
	if 'Iy' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Iy'][tindex:] *= fgain[tindex:]
	if 'Iz' in rx.outputs:
		f['/rxs/rx' + str(rxindex + 1) + '/Iz'][tindex:] *= fgain[tindex:]
	



	
