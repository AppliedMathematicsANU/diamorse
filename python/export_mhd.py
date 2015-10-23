#!/usr/bin/env python

import numpy as np


def parse_options():
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [OPTIONS] INFILE [OUTPREFIX]")
    parser.add_option('-t', '--threshold', dest = 'threshold', metavar = 'X',
                      type = 'float', default = 1.0,
                      help = 'simplification threshold (default 1.0)')

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("at least one argument needed; -h shows list of options")

    return options, args


def ensure_dir(path):
    import os, os.path

    if not os.path.exists(path):
        os.makedirs(path)


def scaled(data, target, symmetric = False):
    lo = data.min()
    hi = data.max()

    if symmetric:
        hi = max(abs(lo), abs(hi))
        lo = -hi

    if hi == lo:
        return data - lo
    else:
        f = float(target) / (hi - lo)
        return (data.astype('float32') - lo) * f


def write_as_mhd(data, file_name, scale = True, symmetric = False):
    import struct

    (zs, ys, xs) = data.shape

    mhd_filename = "%s.mhd" % file_name
    raw_filename = "%s.raw" % file_name

    mhd_text = (
"""
ObjectType = Image
NDims = 3
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
CompressedSize = ???
TransformMatrix = 1 0 0 0 1 0 0 0 1
Offset = 0 0 0
CenterOfRotation = 0 0 0
ElementSpacing = 1 1 1
HeaderSize = -1
DimSize = %d %d %d
AnatomicalOrientation = ???
ElementType = MET_UCHAR
ElementDataFile = %s
""" % (xs, ys, zs, raw_filename))

    fp = open(mhd_filename, 'w')
    fp.write(mhd_text)
    fp.close()

    output = scaled(data, 254, symmetric) if scale else data

    fp = open(raw_filename, 'wb')
    fp.write('\0')
    fp.write(struct.pack('<3i', *data.shape))
    fp.write(output.astype('uint8').tostring())
    fp.close()


def getProcessedMorseData(infile, threshold):
    if infile.endswith('.nc'):
        img = VolumeImage(infile)
        morse = VectorField(img, threshold)
        critical = list(morse.criticalCells())

        data = {
            'scalars'   : img.data(),
            'basins'    : morse.basinMap(),
            'watersheds': morse.watersheds(),
            'skeleton'  : morse.skeleton(),
            'paths'     : morse.paths(),
            'critical'  : critical,
            'critval'   : map(img.scalarForCell, critical),
            'critdim'   : map(img.cellDimension, critical)
            }

        outfile = "%s.npz" % os.path.splitext(os.path.basename(infile))[0]
        np.savez(outfile, **data)
    else:
        data = np.load(infile)

    return data


def shuffle(data):
    tmp = data.astype(np.uint8)
    out = np.zeros(data.shape, np.uint8)

    for i in range(6):
        out |= ((tmp >> i) & 1) << (6 - i)

    return np.where(data < data.max(), out, 254)


if __name__ == '__main__':
    import re, os.path

    from MorseAnalysis import VolumeImage, VectorField

    inf = float('inf')

    (options, args) = parse_options()
    infile = args[0]
    if len(args) > 1:
        path = args[1]
    else:
        path = '.'

    data = getProcessedMorseData(infile, options.threshold)

    basename = re.sub('^tomo_float_', '',
                      os.path.splitext(os.path.basename(infile))[0])
    
    scalars    = data['scalars']
    basins     = shuffle(data['basins'])
    pores      = np.where(scalars > 0, 254, basins)
    skeleton   = np.where(data['skeleton'] > 0, 254, basins)
    paths      = data['paths']
    watersheds = np.where(data['watersheds'] > 0, scalars, scalars.max() * 1.1)
 
    a = scaled(scalars, 256, symmetric = True)
    b = np.where(a >= 128, a, 64)
    c = np.where(data['skeleton'] <= 0, 32, b)
    d = np.where(data['watersheds'] > 0, 96, 64)
    combined = np.where(c != 64, c, d)

    def write(data, suffix, scale = False, symmetric = False):
        fname = os.path.join(path, "%s_%s" % (basename, suffix))
        ensure_dir(os.path.dirname(fname))
        write_as_mhd(data,
                     re.sub('^\./', '', fname),
                     scale = scale,
                     symmetric = symmetric)

    write(scalars,    'scalars',  symmetric = True)
    write(basins,     'basins')
    write(pores,      'pores')
    write(skeleton,   'skeleton')
    write(paths,      'paths')
    write(watersheds, 'watersheds', symmetric = True)
    write(combined,   'combined', scale = False)
