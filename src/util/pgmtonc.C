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
    if (instream.get() != 'P' or instream.get() != '5')
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


std::string makeID(std::string const path, std::string const prefix)
{
    std::string t = path;

    if (t.rfind(".") != std::string::npos)
        t = t.substr(0, t.rfind("."));
    if (t.rfind("/") != std::string::npos)
        t = t.substr(t.rfind("/") + 1);

    return derivedID(prefix + t, prefix, "IMP");
}


template<typename T>
void writeOutput(
    std::string const inpath,
    std::string const outpath,
    std::vector<ImageDescriptor> const& images,
    std::ifstream& instream,
    bool makeSegmentation = false,
    size_t threshold = 1)
{
    namespace js = anu_am::json;

    size_t const xdim = images[0].width;
    size_t const ydim = images[0].height;
    size_t const zdim = images.size();

    size_t  const m = xdim * ydim;
    size_t  const n = zdim * m;

    boost::shared_ptr<std::vector<T> > data(new std::vector<T>(n));
    size_t k = 0;

    for (size_t i = 0; i < images.size(); ++i)
    {
        instream.seekg(images[i].offset, instream.beg);

        for (size_t j = 0; j < m; ++j)
        {
            size_t const x = instream.get();
            data->at(k) = makeSegmentation ? (x >= threshold) : x;
            ++k;
        }
    }

    // Generate metadata to include with the output data
    std::string const prefix = makeSegmentation ? "segmented" : "tomo_float";
    std::string const id = makeID(inpath, prefix);
    std::string const outfile =
        outpath.size() > 0 ? outpath : (stripTimestamp(id) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , id)
        ("process"     , "Import Netpbm grayscale (.pgm) image(s) into NetCDF")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("predecessors", js::Array())
        ("parameters"  , js::Object("input", inpath));

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting data to the output file
    writeVolumeData(
        data, outfile, prefix, xdim, ydim, zdim,
        VolumeWriteOptions()
        .datasetID(id)
        .description(description));
}


void usage(char *name)
{
    std::cerr << "Usage: " << name
              << " [-b] [-t N] INPUT [OUTPUT]" << std::endl
              << "Options:" << std::endl
              << "    -b write a segmentation" << std::endl
              << "    -t threshold for segmentation" << std::endl;
}


int main(int argc, char* argv[])
{
    bool makeSegmentation = false;
    size_t threshold = 1;

    char c;
    while ((c = getopt (argc, argv, "bt:")) != -1)
    {
        switch (c)
        {
        case 'b':
            makeSegmentation = true;
            break;
        case 't':
            makeSegmentation = true;
            threshold = atoi(optarg);
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (argc < optind + 1)
    {
        usage(argv[0]);
        return 1;
    }

    std::string const infile = argv[optind];

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

    std::string const outfile = argc > optind+1 ? argv[optind+1] : "";

    if (makeSegmentation)
        writeOutput<int8_t>(infile, outfile, images, instream, true, threshold);
    else
        writeOutput<float>(infile, outfile, images, instream);
}
