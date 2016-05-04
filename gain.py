# -v.1.0

import argparse, glob, os
import h5py
import numpy as np
from enum import Enum

# Parse command line arguments
parser = argparse.ArgumentParser(description='added a gain on fields values from a base file, then return the new output file.', usage='cd gprMax; python -m tools.gain basefile initialtime (-exp)')
parser.add_argument('basefile', help='base name of output file including path')
parser.add_argument('initialtime', help='time (ns) at which the gain is added')
parser.add_argument('-lin','--linear', action="store_true", help='chosen an linear gain type')
args = parser.parse_args()

basefile = args.basefile
tinit = float(args.initialtime)* 1e-9
print(tinit)
outputfile = basefile + '_gained.out'


# Create output file
#fout = h5py.File(outputfile, 'w')

# Open base file and read some attributes
f = h5py.File(basefile+'.out', 'r+')

iterations = f.attrs['Iterations'] 
dt = f.attrs['dt'] 
#rxs = f.attrs['rxs']
time = np.linspace(0, 1, iterations)
time *= (iterations * dt)
print(time)
fgain = []
	
# Create gain function chosen by user
if args.linear:

	for t in time:
		fgain.append(6*t*1e9)
else:
	fgain= np.ones(np.size(time))


# Find initial time idexation 
tindex=0
while time[tindex] <= tinit:
	tindex += 1 
print(tindex)
rxindex = 0
# For each rx, add the gain on fields component values 
#for rxindex, rx in enumerate(rxs):
		
	#if 'Ex' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Ex'][tindex:] *= fgain[tindex:]
	#if 'Ey' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Ey'][tindex:] *= fgain[tindex:]
	#if 'Ez' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Ez'][tindex:] *= fgain[tindex:]
	#if 'Hx' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Hx'][tindex:] *= fgain[tindex:]
	#if 'Hy' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Hy'][tindex:] *= fgain[tindex:]
	#if 'Hz' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Hz'][tindex:] *= fgain[tindex:]
	#if 'Ix' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Ix'][tindex:] *= fgain[tindex:]
	#if 'Iy' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Iy'][tindex:] *= fgain[tindex:]
	#if 'Iz' in rx.outputs:
f['/rxs/rx' + str(rxindex + 1) + '/Iz'][tindex:] *= fgain[tindex:]
	



	
