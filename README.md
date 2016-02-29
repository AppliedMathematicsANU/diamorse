# diamorse

Digital image analysis using discrete Morse theory and persistent homology.

## References

Delgado-Friedrichs, O., Robins, V., & Sheppard, A. (2015). Skeletonization and partitioning of digital images using discrete Morse theory. *Pattern Analysis and Machine Intelligence, IEEE Transactions on*, 37(3), 654-666.

Robins, V., Wood, P.J., Sheppard, A.P. (2011). Theory and algorithms for constructing discrete Morse complexes from grayscale digital images. *Pattern Analysis and Machine Intelligence, IEEE Transactions on*, 33(8), 1646-1658.


# Documentation

Below is some minimal information that should help you get started. We plan to add more detail over time. Please contact us if you have any questions.

## Installation

In order to compile diamorse, you will need *git*, a C++ compiler, the *Boost* library, and a GNU-compatible *make*.

* Clone this repository to your machine: `git clone https://github.com/AppliedMathematicsANU/diamorse.git`
* Run `make main` to compile only the main analysis programs, or `make` to also compile some utility programs.
* Run `make python` to compile the python wrappers (optional, requires *cython* and *numpy*).

Installation on MacBook Pro is possible with some modification to the make file.  Contact us for details if you are having problems. 

## File formats

Diamorse expects input data to be stored in the NetCDF version 3 format.

We have provided some code for converting other image formats to NetCDF3 as follows. 

* diamorse/python/img2raw.py 

  requires python package *PIL*

  INPUT:  file.img  (a 2d image file in one of the standard formats) 

  OUTPUT:	file.raw (a raw image file with default type uint8)
          file.info (a text file with info about the image size and type)

  USAGE: `diamorse $ ./python/img2raw.py sample.bmp`


* diamorse/src/util/nctomofromraw.C 

  Creates a NetCDF tomo float file from a raw 8-bit (greyscale) 2d or 3d image file.
  
  INPUT:  file.raw (a raw 8-bit 2d or 3d image file)
          xdim (OPTIONAL ydim) (the x dimension of the image, and y dimension if input is a 3d image). 
 
  OUTPUT:	file.nc (a NetCDF3 file) 

  USAGE: `diamorse $ ./bin/nctomofromraw sample.raw sample.nc XDIM [YDIM]`

* diamorse/src/util/ncfromraw.C

  Creates a NetCDF file from a raw binary (black-and-white) 2d or 3d image file. 

  INPUT:	file.raw (a raw binary 2d or 3d image file)
          xdim (OPTIONAL ydim) (the x dimension of the image, and y dimension if input is a 3d image). 

  OUTPUT:	file.nc (a NetCDF3 file)

  USAGE: `diamorse $ ./bin/ncfromraw sample.raw sample.nc XDIM [YDIM]`

* diamorse/main/SEDT.C

  Takes a segmented image and computes the signed Euclidean distance for each voxel using the Hirata/Meijster algorithm.

  INPUT:	file.nc (a binary image in NetCDF3 format)

  OUTPUT:	tomo_float_file_SEDT.nc  (the SEDT of the input image)  

  USAGE: `diamorse $ ./bin/SEDT segmented_sample.nc tomo_float_sample_SEDT.nc`


## Usage

Once you have a greyscale image in NetCDF3 format you can generate the persistence pairs, Morse skeleton and basins using the following. 

* diamorse/python/persistence.py

  Python wrapper for the vector field and persistence computations
  
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

  Python scripts to provide basic plots of persistence diagrams. 


* diamorse/python/plot_basins.py

  For a 2D image this script can create figures such as Figure 4 in our 2015 IEEE TPAMI paper (reference above).  

  USAGE: `diamorse $ ./python/plot_basins -h` will display the full list of options. 



* diamorse/bin/VectorField 

* diamorse/bin/Simplify 

* diamorse/bin/Skeleton

* diamorse/bin/Pores

  The above programs provide lower level functionality for 3d images.  These compute the Morse vector field from a NetCDF image, simplify it to a desired threshold, output the Morse Skeleton and pore labels as NetCDF files for visualisation. 

  3D visualisation is not currently provided as part of diamorse.  


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
