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


#include <cctype>
#include <fstream>
#include <string>
#include <vector>

#include "volume_io.hpp"

using namespace anu_am::diamorse;


struct ImageDescriptor {
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


void skipAndCheckMagicNumber(std::ifstream& instream)
{
    char magic[2];

    instream.read(&(magic[0]), 2);
    if (magic[0] != 'P' or magic[1] != '5')
        throw std::runtime_error("expected magic number 'P5'");
}


bool isCRorLF(int const c)
{
    return c == '\r' or c == '\n';
}


int peekNextSignificantCharacter(std::ifstream& instream)
{
    while (instream.peek() == '#')
        while (not isCRorLF(instream.get()))
            ;
    return instream.peek();
}


void skipMandatorySingleWhitespaceCharacter(std::ifstream& instream)
{
    if (not isspace(peekNextSignificantCharacter(instream)))
        throw std::runtime_error("expected whitespace");
    instream.get();
}


void skipMandatoryWhiteSpace(std::ifstream& instream)
{
    skipMandatorySingleWhitespaceCharacter(instream);

    while (isspace(peekNextSignificantCharacter(instream)))
        instream.get();
}


size_t readNumber(std::ifstream& instream, std::string const varname)
{
    size_t n = 0;

    skipMandatoryWhiteSpace(instream);

    if (not isdigit(peekNextSignificantCharacter(instream)))
        throw std::runtime_error("expected number for " + varname);

    while (isdigit(peekNextSignificantCharacter(instream)))
        n = n * 10 + instream.get() - '0';

    return n;
}


void skipData(std::ifstream& instream, ImageDescriptor const& descriptor)
{
    size_t const pixelCount = descriptor.width * descriptor.height;
    size_t const bytesPerPixel = descriptor.maxval > 255 ? 2 : 1;

    instream.seekg(pixelCount * bytesPerPixel, instream.cur);
}


ImageDescriptor nextImage(std::ifstream& instream)
{
    skipAndCheckMagicNumber(instream);

    size_t const width  = readNumber(instream, "width");
    size_t const height = readNumber(instream, "height");
    size_t const maxval = readNumber(instream, "maximum grayscale value");

    skipMandatorySingleWhitespaceCharacter(instream);

    ImageDescriptor const descriptor(width, height, maxval, instream.tellg());
    skipData(instream, descriptor);

    return descriptor;
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT OUTPUT" << std::endl;
        return 1;
    }

    char* infile = argv[1];
    char* outfile = argv[2];

    std::ifstream instream(infile, std::ifstream::binary);

    while (instream.peek() != EOF)
    {
        ImageDescriptor img = nextImage(instream);
        printf("%ld %ld %ld\n", img.width, img.height, img.maxval);
    }
}
