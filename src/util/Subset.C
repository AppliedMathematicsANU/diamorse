/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  Subset.C
 *
 *  Extracts a subvolume from a NetCDF volume file.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>

#include "json.hpp"
#include "volume_io.hpp"
#include "stringUtils.hpp"

using namespace anu_am::diamorse;


struct Range
{
    size_t start;
    size_t stop;
    size_t stride;

    Range(size_t const start = 0, size_t const stop = 0, size_t const stride = 1)
        : start(start),
          stop(stop),
          stride(std::max((size_t)1, stride))
    {
    }

    Range(std::string const spec)
        : start(0),
          stop(0),
          stride(1)
    {
        std::vector<std::string> const parts = anu_am::stringutils::split(spec, ':');
        if (parts.size() >= 3)
            stride = std::max(1, atoi(parts.at(2).c_str()));
        if (parts.size() >= 2)
            stop = atoi(parts.at(1).c_str());
        start = atoi(parts.at(0).c_str());
    }
};


size_t stop(Range const& range, size_t const srcSize)
{
    return range.stop ? range.stop : srcSize;
}


size_t size(Range const& range, size_t const srcSize)
{
    size_t const end = stop(range, srcSize);
    return (end + range.stride - 1 - range.start) / range.stride;
}


template<typename T>
void extract(
    std::string const inpath,
    std::string const outpath,
    Range const& xrange,
    Range const& yrange,
    Range const& zrange)
{
    namespace js = anu_am::json;

    // Read the data.
    NCFileInfo const info = readFileInfo(inpath);
    Variable const var = findVolumeVariable(info);
    std::string const name = var.name();

    std::vector<size_t> dims = readDimensions(info);
    size_t const xdim = dims.at(0);
    size_t const ydim = dims.at(1);
    size_t const zdim = dims.at(2);

    size_t const xsize = size(xrange, xdim);
    size_t const ysize = size(yrange, ydim);
    size_t const zsize = size(zrange, zdim);

    size_t const xmax = stop(xrange, xdim);
    size_t const ymax = stop(yrange, ydim);
    size_t const zmax = stop(zrange, zdim);

    std::shared_ptr<std::vector<T> > const data = readVolumeData<T>(inpath);

    // Do the processing.
    std::shared_ptr<std::vector<T> > const
        output(new std::vector<T>(xsize * ysize * zsize));

    size_t k = 0;
    for (size_t z = zrange.start; z < zmax; z += zrange.stride)
    {
        for (size_t y = yrange.start; y < ymax; y += yrange.stride)
        {
            for (size_t x = xrange.start; x < xmax; x += xrange.stride)
            {
                output->at(k) = data->at((z * ydim + y) * xdim + x);
                ++k;
            }
        }
    }

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(inpath, info.attributes());
    std::string const thisID   = derivedID(parentID, name, "SS");

    std::string const outfile =
        outpath.size() > 0 ? outpath : (stripTimestamp(thisID) + ".nc");

    js::Object const parameters = js::Object
        ("start_x" , xrange.start)
        ("stop_x"  , stop(xrange, xdim))
        ("stride_x", xrange.stride)
        ("start_y" , yrange.start)
        ("stop_y"  , stop(yrange, ydim))
        ("stride_y", yrange.stride)
        ("start_z" , zrange.start)
        ("stop_z"  , stop(zrange, zdim))
        ("stride_z", zrange.stride);

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Subset")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , parameters);

    std::string const description = js::toString(fullSpec, 2);

    // Write the results.
    writeVolumeData(
        output, outfile, name, xsize, ysize, zsize,
        VolumeWriteOptions()
        .fileAttributes(inheritableAttributes(info.attributes()))
        .datasetID(thisID)
        .description(description));
}


void usage(char *name)
{
    std::cerr << "Usage: " << name 
              << " [-x RANGE] [-y RANGE] [-z RANGE] INPUT [OUTPUT]" << std::endl
              << "where" << std::endl
              << "    RANGE = [start[:stop[:stride]]]" << std::endl;
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    Range xrange, yrange, zrange;

    char c;
    while ((c = getopt (argc, argv, "x:y:z:")) != -1)
    {
        switch (c)
        {
        case 'x':
            xrange = Range(optarg);
            break;
        case 'y':
            yrange = Range(optarg);
            break;
        case 'z':
            zrange = Range(optarg);
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

    std::string const in  = argv[optind];
    std::string const out = argc > optind+1 ? argv[optind+1] : "";

    Variable const var = findVolumeVariable(in);

    switch(var.type())
    {
    case NC_BYTE  : extract<int8_t> (in, out, xrange, yrange, zrange); break;
    case NC_SHORT : extract<int16_t>(in, out, xrange, yrange, zrange); break;
    case NC_LONG  : extract<int32_t>(in, out, xrange, yrange, zrange); break;
    case NC_FLOAT : extract<float>  (in, out, xrange, yrange, zrange); break;
    case NC_DOUBLE: extract<double> (in, out, xrange, yrange, zrange); break;
    default: break;
    }
}
