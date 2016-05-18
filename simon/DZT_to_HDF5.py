# ASCII --> HDF5 v.2.0

import os, argparse
import h5py
import numpy as np

from gprMax.exceptions import CmdInputError

# Parse command line arguments
parser = argparse.ArgumentParser(description='ASCII --> HDF5 format')
parser.add_argument('filename', help='name of ASCII file including path')
args = parser.parse_args()

filename = args.filename
outputfile = filename.replace('.DZT', '.out')



# Open the DZT file 
f = open(filename,'r', errors='replace')
print(f.read())




