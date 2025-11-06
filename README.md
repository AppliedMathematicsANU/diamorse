# diamorse

Digital image analysis using discrete Morse theory and persistent homology.

## References

Delgado-Friedrichs, O., Robins, V., & Sheppard, A. (2015). Skeletonization and partitioning of digital images using discrete Morse theory. *Pattern Analysis and Machine Intelligence, IEEE Transactions on*, 37(3), 654-666.

Robins, V., Wood, P.J., Sheppard, A.P. (2011). Theory and algorithms for constructing discrete Morse complexes from grayscale digital images. *Pattern Analysis and Machine Intelligence, IEEE Transactions on*, 33(8), 1646-1658.


# Documentation

Below is some minimal information that should help you get started. We plan to add more detail over time. Please contact us if you have any questions.

## Installation

In order to compile diamorse, you will need *git*, a GNU-compatible *make*, and a C++ compiler that supports the C++11 standard. Please contact us if any of these requirements poses a problem. The python scripts use python2.7. 

* Clone this repository to your machine: `git clone https://github.com/AppliedMathematicsANU/diamorse.git`
* Change into the newly created directory: `cd diamorse`
* Run `make main` to compile only the main analysis programs, or `make` to also compile some utility programs.
* Run `make python` to compile the python wrappers (optional, requires *cython* and *numpy*). Some of the Python scripts provided require *matplotlib*.

!!! CAUTION: Python wrappers are out of date. !!!

We have installed and run the code successfully under Linux and OS X, and expect that it will work out of the box on most Unix/Posix compatible systems. It is also possible to use diamorse under Windows via *Cygwin* (see the file [cygwin.md](https://github.com/AppliedMathematicsANU/diamorse/blob/master/cygwin.md) for details).

## File formats

All diamorse programs use the NetCDF3 file format for input and output of image data. Please note that NetCDF4 (a.k.a. HDF5) is not supported. To make the conversion into NetCDF easier, we provide a utility that converts portable graymap (.pgm) files as introduced by the Netpbm library (see http://netpbm.sourceforge.net/doc/pgm.html) into NetCDF. The .pgm format was chosen because it is already being supported by a number of libraries and tools, and is also simple enough to easily write out from within a program without having to rely on a particular library. It allows a single file to contain a sequence of (2D) images, and the conversion program uses this feature to encode a 3D dataset. In other words, if the input .pgm file contains multiple images, each one is taken as a layer in the 3D output dataset, with the first image representing the layer at z=0, the next one that at z=1, and so on. If only a single image is present, the output will effectively be a 2D dataset.

* diamorse/util/pgmtonc.C

  Converts a portable graymap (.pgm) file into NetCDF data. Both 8-bit and 16-bit images are supported.
  
  OPTION: -b (create a segmented (black and white) image with all nonzero voxels set to 1)
  
  OPTION: -t <int> (create a segmented image with all voxels larger or equal to the specified value set to 1)
  
  INPUT:  file.pgm  (a file with either a single 2d image or a stack of images of equal width and height)
  
  OUTPUT: segmentedfile.nc OR tomo_floatfile.nc (the corresponding NetCDF file)
  
  USAGE:  `diamorse $ ./bin/pgmtonc image.pgm` OR `diamorse $ ./bin/pgmtonc image.pgm -t 128`

* diamorse/main/SEDT.C

  Takes a segmented image and computes the signed Euclidean distance for each voxel using the Hirata/Meijster algorithm.

  INPUT:	file.nc (a binary image in NetCDF3 format)

  OUTPUT:	tomo_float_file_SEDT.nc  (the SEDT of the input image)  

  USAGE: `diamorse $ ./bin/SEDT segmented_sample.nc tomo_float_sample_SEDT.nc`


## Usage

Once you have a greyscale image in NetCDF3 format you can generate the persistence pairs, Morse skeleton and basins using the following. 

CAUTION: Python scripts are out of date.

* diamorse/python/persistence.py

  Python2.7 wrapper for the vector field and persistence computations
  
  INPUT:	file.nc  (greyscale NetCDF image)

  OPTION:	-t <float>  (specify the simplification threshold for the vector field computations, default is 1.0 ) 

  OPTION: -r (tells the script to write out the persistence pairs to stdout, pipe to pairs.txt) 

  check the source code for other options for input, output, and usage.  

  USAGE: `diamorse $ ./python/persistence.py -t 1.0 -r file.nc > pairs.txt`

  pairs.txt contains the persistence pairing results listed as 
  
  `<birth> <death> <dimension> <creator xyz> <destructor xyz> <weight>`

  The persistence diagram for homology in dimension k is extracted by grabbing lines with `<dimension> = k` 

  The locations of creator and destroyer critical cells are specified by the geometric center of the cell. This means that vertices in the cubical complex will have coordinates that are all integers, edges will have exactly one coordinate that is an integer plus 0.5, 2d faces (squares) will have exactly two half-integer coordinates, and 3d faces (cubes) will have three half-integer coordinates. For example, the cell represented by the coordinate pair (205.5, 169.0) is the edge connecting vertices (205,169) and (206,169).

  the `<weight>` information is an experimental feature - please ignore for now. 


* diamorse/python/plot_persistence.py

  Python2.7 scripts to provide basic plots of persistence diagrams. 


* diamorse/python/plot_basins.py

  For a 2D image this script can create figures such as Figure 4 in our 2015 IEEE TPAMI paper (reference above).  

  USAGE: `diamorse $ ./python/plot_basins -h` will display the full list of options. 


The following commandline executable programs provide lower level functionality for 3d images.  These compute the Morse vector field from a NetCDF image, simplify it to a desired threshold, compute persistence pairs, output the Morse Skeleton and pore labels as NetCDF files for visualisation. 

NOTE: 3D visualisation is not currently provided as part of diamorse.  

* `diamorse $ ./bin/VectorField tomo_float_input.nc [OUTPUT]`

  OUTPUT defaults to vector_field_input_GVF.nc

* `diamorse $ ./bin/Simplify tomo_float_input.nc vector_field_input_GVF.nc [OUTPUT]`

  OUTPUT defaults to vector_field_input_GVF_SMP.nc

  OPTION: -p <float> (default = 1.0) persistence limit for feature cancellation
  
  OPTION: -s <float> size limit for feature cancellation
  
  OPTION: -t <float> value of a level threshold to preserve

* `diamorse $ ./bin/PersistencePairs tomo_float_input.nc vector_field_input_GVF_SMP.nc [OUTPUT]`

  OUTPUT defaults to persistence_input_GVF_SMP_PP.txt

* `diamorse $ ./bin/Skeleton tomo_float_input.nc vector_field_input_GVF_SMP.nc [OUTPUT]`

  OUTPUT defaults to segmented_input_GVF_SMP_SKL.nc

  OPTION: -t <float> (default = 0.0) level set threshold 

* `diamorse $ ./bin/Pores tomo_float_input.nc vector_field_input_GVF_SMP.nc [OUTPUT]`

  OUTPUT defaults to labels_input_GVF_SMP_POR.nc

  OPTION: -t <float> (default = 0.0) level set threshold 
  


# License

The MIT License (MIT)

Copyright (c) 2015 The Australian National University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
