# -*- coding: utf-8 -*-
import pandas as pd
from netCDF4 import Dataset
import os
import sys

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input', type=str,
                    help='The input hdf5 file')
parser.add_argument('output', type=str,
                    help='The output netCDF file')

args = parser.parse_args()

try:
    hdf = pd.read_hdf(args.input).as_matrix()
except IOError:
    print('%s not found' % args.input)
    sys.exit()

if os.path.isfile(args.output):
    r = raw_input('%s already exists on disk. Overwrite? [Y/n] ' % args.output).lower()
    while r not in ['', 'y', 'n']:
        r = raw_input('"%s" is an invalid answer. Please input "y" or "n" ' % r).lower()

    if r == 'n':
        sys.exit()

nc = Dataset(args.output, 'w', format='NETCDF3_CLASSIC')

if len(hdf.shape) <= 3:
    xyz = ['x', 'y', 'z'][:len(hdf.shape)]
    for i, v in enumerate(xyz):
        nc.createDimension(v, hdf.shape[i])
else:
   print('Only 1, 2 and 3D supported')
   sys.exit()

data = nc.createVariable('data', 'f8', xyz)
data[:,:,:] = hdf.copy()

nc.close()
