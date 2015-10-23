/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  rawfrompoints.C
 *
 *  Creates a raw data file from an ASCII file containing point coordinates.
 *
 *  Olaf Delgado-Friedrichs jul 15
 *
 */


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>


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
    long  size    = dim * dim * dim;

    std::vector<int8_t> data(size, 0);

    std::ifstream file(infile);
    std::string str; 
    while (std::getline(file, str))
    {
        std::istringstream ss(str);
        float xx, yy, zz;
        ss >> xx >> yy >> zz;
        long x = xx * dim, y = yy * dim, z = zz * dim;
        data.at((z * dim + y) * dim + x) = 1;
    }

    std::ofstream out(outfile,  std::ios::binary);
    out.write((char *) &data.at(0), size);
}
