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

#include "json.hpp"
#include "netcdf.hpp"
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


bool isCarriageReturnOrLineFeed(int const c)
{
    return c == '\r' or c == '\n';
}


int peekNextSignificantCharacter(std::ifstream& instream)
{
    while (instream.peek() == '#')
        while (not isCarriageReturnOrLineFeed(instream.get()))
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
        throw std::runtime_error("expected a number for " + varname);

    while (isdigit(peekNextSignificantCharacter(instream)))
        n = n * 10 + instream.get() - '0';

    if (n < 1)
        throw std::runtime_error(varname + " must be at least 1");
    else if (n > 65535)
        throw std::runtime_error(varname + " must be at most 65535");

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


bool haveMatchingParameters(
    ImageDescriptor const& img1,
    ImageDescriptor const& img2)
{
    return (
        img1.width  == img2.width and
        img1.height == img2.height and
        img1.maxval == img2.maxval
        );
}


std::string stripExtension(std::string const path)
{
    size_t const pos = path.rfind(".");
    if (pos == std::string::npos)
        return path;
    else
        return path.substr(0, pos);
}


std::string makeID(std::string const path)
{
    std::string t = stripExtension(path);
    size_t const pos = t.rfind("/");
    if (pos != std::string::npos)
        t = t.substr(pos+1);

    return derivedID("tomo_float" + t, "tomo_float", "IMP");
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT [OUTPUT]" << std::endl;
        return 1;
    }

    std::string const infile = argv[1];

    std::ifstream instream(infile.c_str(), std::ifstream::binary);

    std::vector<ImageDescriptor> images;

    while (instream.peek() != EOF)
    {
        ImageDescriptor img = nextImage(instream);

        if (images.size() > 0 and not haveMatchingParameters(img, images[0]))
            throw std::runtime_error(
                "images must have matching dimensions and grayscale ranges"
                );

        images.push_back(img);
    }

    if (images.size() == 0)
        throw std::runtime_error("no images found");

    size_t const xdim = images[0].width;
    size_t const ydim = images[0].height;
    size_t const zdim = images.size();

    size_t  const m = xdim * ydim;
    size_t  const n = zdim * m;

    boost::shared_ptr<std::vector<float_t> > data(new std::vector<float_t>(n));
    size_t k = 0;

    for (size_t i = 0; i < images.size(); ++i)
    {
        instream.seekg(images[i].offset, instream.beg);

        for (size_t j = 0; j < m; ++j)
        {
            data->at(k) = (float) instream.get();
            ++k;
        }
    }

    // Generate metadata to include with the output data
    std::string const id = makeID(infile);
    std::string const outfile = argc > 2 ? argv[2] : stripTimestamp(id) + ".nc";

    js::Object const fullSpec = js::Object
        ("id"          , id)
        ("process"     , "Import Netpbm grayscale (.pgm) image(s) into NetCDF")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("predecessors", js::Array())
        ("parameters"  , js::Object("input", infile));

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting data to the output file
    writeVolumeData(
        data, outfile, "tomo_float", xdim, ydim, zdim,
        VolumeWriteOptions()
        .datasetID(id)
        .description(description));

    //writeVolumeData(data, outfile, "tomo_float", xdim, ydim, zdim);
}
