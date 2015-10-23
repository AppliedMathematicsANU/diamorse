/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  ncdiffstats.C
 *
 *  Analyses the differences between two NetCDF volume datasets.
 *
 *  Olaf Delgado-Friedrichs mar 15
 *
 */

#include <cmath>

#include "netcdf.hpp"
#include "netcdfIO.hpp"

using namespace anu_am::netcdf;


Variable findVolumeVariable(NCFileInfo const info)
{
    std::vector<Variable> const vars = info.variables();

    for (size_t i = 0; i < vars.size(); ++i)
    {
        Variable const v = vars.at(i);

        if (v.dimensions().size() == 3)
            return v;
    }

    throw std::runtime_error("No appropriate variable found");
}


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage:" << argv[0] << " FILE1 FILE2" << std::endl;
        return 1;
    }

    int i = 0;
    char* path1 = argv[++i];
    char* path2 = argv[++i];

    FileBuffer data1(path1);
    FileBuffer data2(path2);

    NCFile<FileBuffer> file1(data1);
    NCFile<FileBuffer> file2(data2);

    Variable const var1 = findVolumeVariable(file1.info());
    Variable const var2 = findVolumeVariable(file2.info());

    std::vector<size_t> const dims = var1.dimensions();
    if (var2.dimensions() != dims)
        throw std::runtime_error("dimension mismatch");

    size_t const zdim = dims.at(0);
    size_t const ydim = dims.at(1);
    size_t const xdim = dims.at(2);

    double mindiff = 0, maxdiff = 0;

    for (size_t z = 0; z < zdim; ++z)
    {
        for (size_t y = 0; y < ydim; ++y)
        {
            for (size_t x = 0; x < xdim; ++x)
            {
                double const val1 = file1.getFloat(var1, x, y, z);
                double const val2 = file2.getFloat(var2, x, y, z);
                double const diff = val1 - val2;

                if (diff < mindiff) mindiff = diff;
                if (diff > maxdiff) maxdiff = diff;
            }
        }
    }

    std::cout << mindiff << " ... " << maxdiff << std::endl;
}
