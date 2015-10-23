/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  ncfrompoints.C
 *
 *  Creates a NetCDF file from an ASCII file containing point coordinates.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */


#include <fstream>
#include <sstream>
#include <string>

#include "volume_io.hpp"

using namespace anu_am::diamorse;


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT OUTPUT DIM"
                  << std::endl;
        return 1;
    }

    char* infile  = argv[1];
    char* outfile = argv[2];
    long  dim     = atol(argv[3]);

    boost::shared_ptr<std::vector<int8_t> >
        data(new std::vector<int8_t>(dim*dim*dim, 1));

    std::ifstream file(infile);
    std::string str; 
    while (std::getline(file, str))
    {
        std::istringstream ss(str);
        float xx, yy, zz;
        ss >> xx >> yy >> zz;
        int x = xx * dim, y = yy * dim, z = zz * dim;
        data->at((z * dim + y) * dim + x) = 0;
    }

    writeVolumeData(data, outfile, "segmented", dim, dim, dim);
}
