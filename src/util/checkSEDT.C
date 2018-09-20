/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  checkSEDT.C
 *
 *  Sanity check for SEDT output.
 *
 *  Olaf Delgado-Friedrichs feb 15
 *
 */


#include "volume_io.hpp"
#include "collections.hpp"
#include "CubicalComplex.hpp"
#include "VertexMap.hpp"

using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef int8_t PhaseValue;
typedef float DistanceValue;
typedef VertexMap<CubicalComplex, PhaseValue> Phases;
typedef VertexMap<CubicalComplex, DistanceValue> Distances;


inline double sq(double const x) {
    return x * x;
}


int run(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage:" << argv[0] << " SEGMENTED SEDT" << std::endl;
        return 1;
    }

    char * segmentedFile = argv[1];
    char * sedtFile = argv[2];

    std::vector<size_t> const dims = readDimensions(segmentedFile);
    size_t const xdim = dims.at(0);
    size_t const ydim = dims.at(1);
    size_t const zdim = dims.at(2);

    CubicalComplex complex(xdim, ydim, zdim);

    Phases::DataPtr phaseData = readVolumeData<PhaseValue>(segmentedFile);
    Phases phases(complex, phaseData);

    Distances::DataPtr distData = readVolumeData<DistanceValue>(sedtFile);
    Distances dists(complex, distData);

    for (size_t z = 0; z < zdim; ++z)
    {
        for (size_t y = 0; y < ydim; ++y)
        {
            for (size_t x = 0; x < xdim; ++x)
            {
                Cell const v = complex.cellAt(x, y, z);
                if (dists(v) == 0)
                    std::cout << "distance is zero" << std::endl;

                if (phases(v) != (dists(v) > 0))
                    std::cout << "distance has wrong sign" << std::endl;

                std::vector<Cell> neighbors;
                if (x+1 < xdim)
                    neighbors.push_back(complex.cellAt(x+1, y, z));
                if (y+1 < ydim)
                    neighbors.push_back(complex.cellAt(x, y+1, z));
                if (z+1 < zdim)
                    neighbors.push_back(complex.cellAt(x, y, z+1));

                for (size_t i = 0; i < neighbors.size(); ++i) {
                    Cell const w = neighbors.at(i);

                    if (fabs(dists(v) - dists(w)) > 1 + 1e-6)
                    {
                        std::cout << "distances "
                                  << dists(v) << " and " << dists(w)
                                  << " at " << complex.cellPosition(v)
                                  << " and " << complex.cellPosition(w)
                                  << std::endl;
                    }

                    if (phases(v) != phases(w) && (fabs(dists(v)) != 0.5 ||
                                                   fabs(dists(w)) != 0.5))
                        std::cout << "bad value at boundary";
                }
            }
        }
    }

    for (size_t i = 0; i < 3 * xdim * ydim * zdim; ++i) {
        size_t const xv = rand() % xdim;
        size_t const yv = rand() % ydim;
        size_t const zv = rand() % zdim;
        size_t const xw = rand() % xdim;
        size_t const yw = rand() % ydim;
        size_t const zw = rand() % zdim;

        Cell const v = complex.cellAt(xv, yv, zv);
        Cell const w = complex.cellAt(xw, yw, zw);

        double const d1 = sqrt(sq(xw - xv) + sq(yw - yv) + sq(zw -zv));
        double const d2 = fabs(dists(v) - dists(w));

        if (d2 > d1 and d2 - d1 > 1e-6 * std::max(d1, d2))
            std::cout << "distances "
                      << dists(v) << " and " << dists(w)
                      << " at " << complex.cellPosition(v)
                      << " and " << complex.cellPosition(w)
                      << " (d1 = " << d1 << ", d2 = " << d2 << ")"
                      << std::endl;
    }

    return 0;
}


int main(const int argc, char* argv[])
{
  try
  {
    run(argc, argv);
  }
  catch(std::runtime_error& e)
  {
    std::clog
      << "terminate called after throwing an instance of "
      << "'std::runtime_error'\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
  catch(std::exception& e)
  {
    std::clog
      << "terminate called after throwing an exception\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
}
