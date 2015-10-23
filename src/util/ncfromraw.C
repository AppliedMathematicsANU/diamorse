/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  ncfromraw.C
 *
 *  Creates a NetCDF file from raw 2d or 3d binary image file.
 *
 *  Olaf Delgado-Friedrichs feb 15
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

    boost::shared_ptr<std::vector<int8_t> >
        data(new std::vector<int8_t>(length));

    instream.read((char*) &(data->at(0)), length);

    int xdim = atoi(argv[3]);
    int ydim = argc > 4 ? atoi(argv[4]) : length / xdim;
    int zdim = length / (xdim * ydim);

    writeVolumeData(data, outfile, "segmented", xdim, ydim, zdim);
}
