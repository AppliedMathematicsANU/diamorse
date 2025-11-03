#!/usr/bin/env python

import numpy as _np


class MorseVectorField(object):
    def __init__(self, source, threshold=-1.0):
        from MorseAnalysis import read_netcdf, VolumeImage, VectorField

        if isinstance(source, str):
            data = read_netcdf(source)
            self._inputfile = source
        elif isinstance(source, _np.ndarray):
            data = source
            self._inputfile = None
        else:
            raise RuntimeError("expected a file name or Numpy array")

        self._morse = VectorField(VolumeImage(
            data.astype(_np.float32)),
            threshold=threshold
        )
        self._threshold = threshold


    def scalars(self):
        return self._morse.img_data()


    def basin_labels(self):
        return self._morse.basinMap().astype(_np.int32)


    def pore_labels(self, watermark=0.0):
        scalars = self.scalars()
        basins = self.basin_labels()
        basins[scalars > watermark] = 0

        return basins


    def births_and_deaths(self, dimension, threshold):
        if not hasattr(self, "_births_and_deaths"):
            dim = lambda v: self._morse.cellDimension(v)
            val = lambda v: self._morse.scalarForCell(v) if v else _np.inf

            self._births_and_deaths = tuple(
                (val(v), val(w), dim(v))
                for v, w in self._morse.birthsAndDeaths()
            )

        return tuple(
            (birth, death)
            for birth, death, dim in self._births_and_deaths
            if dim == dimension and death - birth > threshold
        )


    def births(self, dimension, threshold):
        return tuple(
            birth for birth, _ in self.births_and_deaths(dimension, threshold)
        )


    def deaths(self, dimension, threshold):
        return tuple(
            death for _, death in self.births_and_deaths(dimension, threshold)
        )


    def betti_numbers(self, dim, threshold):
        births = tuple((birth,  1) for birth in self.births(dim, threshold))
        deaths = tuple((death, -1) for death in self.deaths(dim, threshold))
        events = sorted(births + deaths)

        result = []

        if events:
            (x, n, n0) = (events[0][0], 0, -1)
            for (y, m) in events:
                if y != x:
                    if n != n0:
                        result.append((x, n))
                        n0 = n
                    x = y
                n += m
            if n != n0:
                result.append((x, n))

        return tuple(result)


if __name__ == "__main__":
    a = 0.5 - _np.array([[
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
    ]])

    mvf = MorseVectorField(a)

    print(f"Scalars:\n{mvf.scalars()}\n")
    print(f"Basin labels:\n{mvf.basin_labels()}\n")
    print(f"Pore labels:\n{mvf.pore_labels()}\n")

    for dim in range(4):
        print(f"Dimension {dim}")
        print(f"  Birth-death pairs: {mvf.births_and_deaths(dim, 0)}")
        print(f"  Births: {mvf.births(dim, 0)}")
        print(f"  Deaths: {mvf.deaths(dim, 0)}")
        print(f"  Betti numbers: {mvf.betti_numbers(dim, 0)}")
        print()
