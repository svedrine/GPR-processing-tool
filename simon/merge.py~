# Testmerge v.1.5

import argparse, glob, os
import h5py
import numpy as np

"""Merges traces (A-scans) from multiple output files into one new file, then removes the series of output files."""

# Parse command line arguments
parser = argparse.ArgumentParser(description='Merges traces (A-scans) from multiple output files into one new file, then removes the series of output files.', usage='cd gprMax; python -m tools.outputfiles_merge basefilename')
parser.add_argument('basefilename', help='base name of output file series including path')
args = parser.parse_args()

basefilename = args.basefilename
outputfile = '%s_merged.out' % basefilename
files = glob.glob('%s*.out' % basefilename)
outputfiles = [filename for filename in files if '_merged' not in filename]
modelruns = len(outputfiles)
print('Found {} files to merge'.format(modelruns))

# Combined output file
fout = h5py.File(outputfile, 'w')
temp = [list() for i in range(modelruns)]
# Add positional data for rxs
for model in range(modelruns):
    fin = h5py.File('%s%s.out' % (basefilename, str(model + 1)), 'r')
    nrx = fin.attrs['nrx']

    # Write properties for merged file on first iteration
    if model == 0:
        fout.attrs['Iterations'] = fin.attrs['Iterations']
        fout.attrs['dt'] = fin.attrs['dt']
        fout.attrs['nrx'] = fin.attrs['nrx']
        fout.attrs['nx, ny, nz'] = fin.attrs['nx, ny, nz']
        fout.attrs['dx, dy, dz'] = fin.attrs['dx, dy, dz']

        for rx in range(1, nrx + 1):
            path = '/rxs/rx%s/' % str(rx)
            grp = fout.create_group(path)
            availableoutputs = list(fin[path].keys())
            for output in availableoutputs:
                grp.create_dataset(output, (fout.attrs['Iterations'], modelruns), dtype=fin[path + output].dtype)

    # For all receivers
    for rx in range(1, nrx + 1):
        path = '/rxs/rx%s/' % str(rx)
        temp[model].append(fin[path].attrs['Position'])
        availableoutputs = list(fin[path].keys())
        # For all receiver outputs
        for output in availableoutputs:
            fout[path + output][:,model] = fin[path + output][:]
  
    fin.close()
fout.attrs['merged_positions'] = temp
fout.close()

check = input('Do you want to remove the multiple individual output files? [y] or n:')
if not check or check == 'y':
    for model in range(modelruns):
        file = '%s%s.out' % (basefilename, str(model + 1))
        os.remove(file)




