# Copyright 2014 The Australian National University
#
# MorseAnalysis.pyx
#
# Cython wrapper for volume image data and associated information such as
# Morse vector field, Morse chain complex and cell labellings.
#
# Olaf Delgado-Friedrichs feb 14

from cython.operator cimport dereference as deref
from libc.stdint cimport uint8_t, uint32_t
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

UINT8 = np.uint8
ctypedef np.uint8_t UINT8_t

FLOAT32 = np.float32
ctypedef np.float32_t FLOAT32_t


cdef extern from "ImageAnalysis.hpp" namespace "anu_am::diamorse":
    cdef cppclass ImageData[V]:
        ImageData() except +
        ImageData(string) except +
        int xdim() except +
        int ydim() except +
        int zdim() except +
        vector[float] positionForCell(size_t) except +
        size_t cellAtPosition(vector[float]) except +
        int cellDimension(size_t) except +
        V cellValue(size_t) except +

    cdef cppclass MorseData[V]:
        MorseData() except +
        MorseData(ImageData[V], float) except +
        MorseData(ImageData[V], string) except +
        unsigned char cellDirection(size_t) except +
        size_t associatedMinimum(size_t) except +
        bint isOnWatershed(size_t) except +
        float skeletonValue(size_t) except +
        bint isOnPath(size_t) except +
        size_t chainComplexSize() except +
        size_t chainComplexCell(size_t) except +
        vector[size_t] chainComplexBoundary(size_t) except +
        vector[size_t] chainComplexCoboundary(size_t) except +
        vector[pair[size_t, size_t]] birthDeathPairs() except +
        vector[pair[size_t, int]] weights() except +


cdef class VolumeImage:
    cdef ImageData[float] * _img

    def __cinit__(self, string filename):
        self._img = new ImageData[float](filename)
 
    def __dealloc__(self):
        del self._img

    def data(self):
        cdef int xdim = self._img.xdim()
        cdef int ydim = self._img.ydim()
        cdef int zdim = self._img.zdim()
        cdef np.ndarray[FLOAT32_t, ndim=3] out = np.zeros(
            [zdim, ydim, xdim], dtype=FLOAT32)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell

        for x in range(xdim):
            pos[0] = x
            for y in range(ydim):
                pos[1] = y
                for z in range(zdim):
                    pos[2] = z
                    cell = self._img.cellAtPosition(pos)
                    out[z, y, x] = self._img.cellValue(cell)

        return out

    def shape(self):
        return (self._img.xdim(), self._img.ydim(), self._img.zdim())

    def cellDimension(self, vector[float] pos):
        return self._img.cellDimension(self._img.cellAtPosition(pos))

    def scalarForCell(self, vector[float] pos):
        return self._img.cellValue(self._img.cellAtPosition(pos))


cdef class VectorField:
    cdef VolumeImage _volume
    cdef MorseData[float] * _morse

    def __cinit__(self,
                  VolumeImage volume,
                  float threshold = -1,
                  string filename = ''
                  ):
        self._volume = volume
        if filename.size() > 0:
            self._morse = new MorseData[float](deref(volume._img), filename)
        else:
            self._morse = new MorseData[float](deref(volume._img), threshold)

    def __dealloc__(self):
        del self._morse

    def data(self):
        cdef int xdim = self._volume._img.xdim()
        cdef int ydim = self._volume._img.ydim()
        cdef int zdim = self._volume._img.zdim()
        cdef np.ndarray[UINT8_t, ndim=3] out = np.zeros(
            [2 * zdim - 1, 2 * ydim - 1, 2 * xdim - 1], dtype=UINT8)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell

        for x in range(2 * xdim - 1):
            pos[0] = x / 2.0
            for y in range(2 * ydim - 1):
                pos[1] = y / 2.0
                for z in range(2 * zdim - 1):
                    pos[2] = z / 2.0
                    cell = self._volume._img.cellAtPosition(pos)
                    out[z, y, x] = self._morse.cellDirection(cell)

        return out

    def vectorForCell(self, vector[float] pos):
        cdef unsigned char d
        d = self._morse.cellDirection(self._volume._img.cellAtPosition(pos))
        return (None, (0, 0, 0),
                (0.5, 0, 0), (-0.5, 0, 0),
                (0, 0.5, 0), (0, -0.5, 0),
                (0, 0, 0.5), (0, 0, -0.5))[d]

    def criticalCells(self):
        cdef size_t i

        for i in range(self._morse.chainComplexSize()):
            yield self._volume._img.positionForCell(
                self._morse.chainComplexCell(i))

    def basinMap(self):
        cdef int xdim = self._volume._img.xdim()
        cdef int ydim = self._volume._img.ydim()
        cdef int zdim = self._volume._img.zdim()
        cdef np.ndarray[FLOAT32_t, ndim=3] out = np.zeros(
            [zdim, ydim, xdim], dtype=FLOAT32)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell
        cdef float val

        index = dict((self._morse.chainComplexCell(i), i)
                     for i in range(self._morse.chainComplexSize()))

        for x in range(xdim):
            pos[0] = x
            for y in range(ydim):
                pos[1] = y
                for z in range(zdim):
                    pos[2] = z
                    cell = self._volume._img.cellAtPosition(pos)
                    out[z, y, x] = index[self._morse.associatedMinimum(cell)]

        return out

    def watersheds(self):
        cdef int xdim = self._volume._img.xdim()
        cdef int ydim = self._volume._img.ydim()
        cdef int zdim = self._volume._img.zdim()
        cdef np.ndarray[UINT8_t, ndim=3] out = np.zeros(
            [zdim, ydim, xdim], dtype=UINT8)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell

        for x in range(xdim):
            pos[0] = x
            for y in range(ydim):
                pos[1] = y
                for z in range(zdim):
                    pos[2] = z
                    cell = self._volume._img.cellAtPosition(pos)
                    out[z, y, x] = self._morse.isOnWatershed(cell)

        return out

    def skeleton(self):
        cdef int xdim = self._volume._img.xdim()
        cdef int ydim = self._volume._img.ydim()
        cdef int zdim = self._volume._img.zdim()
        cdef np.ndarray[FLOAT32_t, ndim=3] out = np.zeros(
            [zdim, ydim, xdim], dtype=FLOAT32)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell

        for x in range(xdim):
            pos[0] = x
            for y in range(ydim):
                pos[1] = y
                for z in range(zdim):
                    pos[2] = z
                    cell = self._volume._img.cellAtPosition(pos)
                    out[z, y, x] = self._morse.skeletonValue(cell)

        return out

    def paths(self):
        cdef int xdim = self._volume._img.xdim()
        cdef int ydim = self._volume._img.ydim()
        cdef int zdim = self._volume._img.zdim()
        cdef np.ndarray[UINT8_t, ndim=3] out = np.zeros(
            [zdim, ydim, xdim], dtype=UINT8)

        cdef int x, y, z
        cdef vector[float] pos = vector[float](3)
        cdef size_t cell

        for x in range(xdim):
            pos[0] = x
            for y in range(ydim):
                pos[1] = y
                for z in range(zdim):
                    pos[2] = z
                    cell = self._volume._img.cellAtPosition(pos)
                    out[z, y, x] = self._morse.isOnPath(cell)

        return out

    def associatedMinimum(self, vector[float] pos):
        cdef size_t v = self._volume._img.cellAtPosition(pos)
        cdef size_t m = self._morse.associatedMinimum(v)
        return self._volume._img.positionForCell(m)

    def chainBoundary(self, vector[float] pos):
        cdef size_t v = self._volume._img.cellAtPosition(pos)
        cdef vector[size_t] cells = self._morse.chainComplexBoundary(v)
        cdef size_t i

        for i in range(cells.size()):
            yield self._volume._img.positionForCell(cells.at(i))

    def birthsAndDeaths(self):
        cdef ImageData[float] * img = self._volume._img
        cdef vector[pair[size_t, size_t]] pairs = self._morse.birthDeathPairs()
        cdef size_t i
        cdef pair[size_t, size_t] p

        for i in range(pairs.size()):
            p = pairs.at(i)
            a = img.positionForCell(p.first)
            b = img.positionForCell(p.second) if p.first != p.second else None
            yield (a, b)

    def weights(self):
        cdef ImageData[float] * img = self._volume._img
        cdef vector[pair[size_t, int]] pairs = self._morse.weights()
        cdef size_t i
        cdef pair[size_t, int] p

        for i in range(pairs.size()):
            p = pairs.at(i)
            a = img.positionForCell(p.first)
            v = p.second
            yield (a, v)
