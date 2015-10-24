# diamorse

Digital image analysis using discrete Morse theory and persistent homology.

## Reference

Delgado-Friedrichs, O., Robins, V., & Sheppard, A. (2015). Skeletonization and partitioning of digital images using discrete Morse theory. *Pattern Analysis and Machine Intelligence, IEEE Transactions on*, 37(3), 654-666.

# Documentation

Below is some minimal information that should help you get started. We plan to add more detail over time. Please contact us if you have any questions.

## Installation

In order to compile diamorse, you will need *git*, a C++ compiler, the *Boost* library, and a GNU-compatible *make*.

* Clone this repository to your machine: `git clone https://github.com/AppliedMathematicsANU/diamorse.git`
* Run `make main` to compile only the main analysis programs, or `make` to also compile some utility programs.
* Run `make python` to compile the python wrappers.

## File formats

Diamorse expects volume data to be stored in the NetCDF version 3 format.

(in progress...)

## Usage

(in progress...)

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
