#!/usr/bin/env python

import sys
import os

from PIL import Image
import numpy as np


if len(sys.argv) < 2:
    print("usage: %s INFILE [OUTTYPE]" % sys.argv[0])
    print("  writes raw data from image file INFILE as data type OUTTYPE " +
          "(default uint8)")
    sys.exit(1)

input  = sys.argv[1]
type   = sys.argv[2] if len(sys.argv) > 2 else 'uint8'

output = os.path.splitext(input)[0] + '.raw'
info   = os.path.splitext(input)[0] + '.info'

image  = Image.open(input)
size   = image.size
data   = np.array(image.getdata())

with open(output, 'wb') as fp:
    fp.write(data.astype(type).tostring())

with open(info, 'w') as fp:
    fp.write("size: %s x %s\n" % size)
    fp.write("type: %s\n" % type)
