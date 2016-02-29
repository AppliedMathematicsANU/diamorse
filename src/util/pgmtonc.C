/** -*-c++-*-
 *
 *  Copyright 2016 The Australian National University
 *
 *  pgmtonc.C
 *
 *  Convert a multi-image Netpbm .pgm file to a NetCDF volume file.
 *
 *  Olaf Delgado-Friedrichs feb 16
 *
 */


#include <fstream>
#include <vector>

#include "volume_io.hpp"

using namespace anu_am::diamorse;


typedef struct ImageDescriptor {
    size_t width;
    size_t height;
    size_t maxval;
    size_t offset;

    ImageDescriptor(
        size_t const width,
        size_t const height,
        size_t const maxval,
        size_t const offset
        )
        : width(width),
          height(height),
          maxval(maxval),
          offset(offset)
    {
    }
};


ImageDescriptor nextImage(std::ifstream& instream)
{
    size_t const pos = instream.tellg();

    checkAndSkipMagicNumber(instream);
    size_t const width  = readNumber(instream, "width");
    size_t const height = readNumber(instream, "height");
    size_t const maxval = readNumber(instream, "maximum grayscale value");

    size_t const bytesPerPixel = maxval > 255 ? 2 : 1;

    checkAndSkipWhitespaceChar(instream);
    instream.seekg(width * height * bytesPerPixel, instream.cur);

    return ImageDescriptor(width, height, maxval, pos);
}


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT OUTPUT" << std::endl;
        return 1;
    }

    char* infile = argv[1];
    char* outfile = argv[2];

    std::ifstream instream(infile, std::ifstream::binary);
}
