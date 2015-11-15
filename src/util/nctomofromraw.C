/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  nctomofromraw.C
 *
 *  Creates a NetCDF tomo float file from a raw 8-bit 2d or 3d image file.
 *
 *  Olaf Delgado-Friedrichs nov 15
 *
 */


#include <fstream>

#include "volume_io.hpp"

using namespace anu_am::diamorse;


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT OUTPUT XDIM [YDIM]"
                  << std::endl;
        return 1;
    }

    char* infile = argv[1];
    char* outfile = argv[2];

    std::ifstream instream(infile, std::ifstream::binary);

    instream.seekg(0, instream.end);
    int length = instream.tellg();
    instream.seekg(0, instream.beg);

    boost::shared_ptr<std::vector<uint8_t> >
        indata(new std::vector<uint8_t>(length));

    instream.read((char*) &(indata->at(0)), length);

    boost::shared_ptr<std::vector<float_t> >
        outdata(new std::vector<float_t>(length));

    for (int i = 0; i < length; ++i)
        (*outdata)[i] = (*indata)[i];

    int xdim = atoi(argv[3]);
    int ydim = argc > 4 ? atoi(argv[4]) : length / xdim;
    int zdim = length / (xdim * ydim);

    writeVolumeData(outdata, outfile, "tomo_float", xdim, ydim, zdim);
}
