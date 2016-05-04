#-v.1.0

import argparse, glob, os
import h5py
import numpy as np

# Parse command line arguments
parser = argparse.ArgumentParser(description='Added a gain on fields values from an output file, then return the new output file.', usage='cd gprMax; python -m tools.gain basefilename')
parser.add_argument('basefilename', help='base name of output file including path')
args = parser.parse_args()

basefilename = args.basefilename
outputfile = basefilename + '_gained_out'

print(outputfile)


#Update electric field compenent by adding an exponential gain
#def apply_gain_exp(Ex, Ey, Ez, Hx, Hy, Hz):

	
