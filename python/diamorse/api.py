#!/usr/bin/env python

import dataclasses as _dc
import numpy as _np


@_dc.dataclass(frozen=True)
class Cell:
    position: tuple[float, float, float]
    dimension: int
    value: float
    weight: int


class MorseVectorField(object):
    def __init__(self, source, threshold=-1.0):
        from .MorseAnalysis import read_netcdf, VolumeImage, VectorField

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


    def vector_field(self):
        return self._morse.data()


    def basin_labels(self):
        return self._morse.basinMap().astype(_np.int32)


    def pore_labels(self, watermark=0.0):
        scalars = self.scalars()
        basins = self.basin_labels()
        basins[scalars > watermark] = 0

        return basins


    def on_watershed(self):
        return self._morse.watersheds()


    def on_path(self):
        return self._morse.paths()


    def skeleton(self):
        return self._morse.skeleton()


    def _make_cell(self, pos, weights):
        if pos is None:
            return None
        else:
            return Cell(
                position=tuple(pos),
                dimension=self._morse.cellDimension(pos),
                value=self._morse.scalarForCell(pos),
                weight=weights.get(tuple(pos), 0),
            )


    def critical_cells(self):
        weights = dict((tuple(v), x) for v, x in self._morse.weights())

        return tuple(
            self._make_cell(pos, weights)
            for pos in self._morse.criticalCells()
        )


    def birth_death_pairs(self):
        weights = dict((tuple(v), x) for v, x in self._morse.weights())

        return tuple(
            (self._make_cell(v, weights), self._make_cell(w, weights))
            for v, w in self._morse.birthsAndDeaths()
        )


    def births(self, dimension, threshold=-1):
        return sorted(tuple(
            birth.value for birth, death in self.birth_death_pairs()
            if (
                birth.dimension == dimension
                and
                (death is None or death.value - birth.value > threshold)
            )
        ))


    def deaths(self, dimension, threshold=-1):
        return sorted(tuple(
            (_np.inf if death is None else death.value)
            for birth, death in self.birth_death_pairs()
            if (
                birth.dimension == dimension
                and
                (death is None or death.value - birth.value > threshold)
            )
        ))


    def betti_numbers(self, dim, threshold=-1):
        births = tuple((birth,  1) for birth in self.births(dim, threshold))
        deaths = tuple((death, -1) for death in self.deaths(dim, threshold))
        events = sorted(births + deaths)

        result = []

        if events:
            last_value = events[0][0]
            last_count = -1
            count = 0

            for value, increment in events:
                if value != last_value:
                    if count != last_count:
                        result.append((last_value, count))
                        last_count = count

                    last_value = value

                count += increment

            if count != last_count:
                result.append((last_value, count))

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
    print(f"Watersheds:\n{mvf.on_watershed()}\n")
    print(f"Paths:\n{mvf.on_path()}\n")
    print(f"Skeleton:\n{mvf.skeleton()}\n")

    print(f"Critical cells:")
    for cell in mvf.critical_cells():
        print(f"  {cell}")
    print()

    print(f"Birth-death pairs:")
    for pair in mvf.birth_death_pairs():
        print(f"  {pair}")
    print()
    print()

    for dim in range(4):
        print(f"Dimension {dim}")
        print(f"  Births: {mvf.births(dim)}")
        print(f"  Deaths: {mvf.deaths(dim)}")
        print(f"  Betti numbers: {mvf.betti_numbers(dim)}")
        print()
