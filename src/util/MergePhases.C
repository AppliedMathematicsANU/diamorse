/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  MergePhases.C
 *
 *  Remaps a range of values in a NetCDF dataset.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <stdint.h>

#include "json.hpp"
#include "volume_io.hpp"

using namespace anu_am::diamorse;


template<typename T>
void extract(
    std::string const inpath,
    std::string const outpath,
    double const from,
    double const to,
    double const into)
{
    namespace js = anu_am::json;

    // Read the data.
    NCFileInfo const info = readFileInfo(inpath);
    Variable const var = findVolumeVariable(info);
    std::string const name = var.name();

    std::vector<size_t> dims = readDimensions(info);
    size_t const n = dims.at(0) * dims.at(1) * dims.at(2);

    std::shared_ptr<std::vector<T> > const data = readVolumeData<T>(inpath);

    // Do the processing.
    for (size_t k = 0; k < n; ++k)
    {
        T const val = data->at(k);
        if (val >= from and val <= to)
            data->at(k) = into;
    }

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(inpath, info.attributes());
    std::string const thisID   = derivedID(parentID, name, "PM");

    std::string const outfile =
        outpath.size() > 0 ? outpath : (stripTimestamp(thisID) + ".nc");

    js::Object const parameters = js::Object
        ("from", from)
        ("to"  , to)
        ("into", into);

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Phase Merge")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , parameters);

    std::string const description = js::toString(fullSpec, 2);

    // Write the results.
    writeVolumeData(
        data, outfile, name, dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(inheritableAttributes(info.attributes()))
        .datasetID(thisID)
        .description(description));
}


int run(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " FROM TO INTO INPUT [OUTPUT]"
                  << std::endl;
        return 1;
    }

    double const from     = atof(argv[1]);
    double const to       = atof(argv[2]);
    double const into     = atof(argv[3]);
    std::string const in  = argv[4];
    std::string const out = argc > 5 ? argv[5] : "";

    Variable const var = findVolumeVariable(in);

    switch(var.type())
    {
    case NC_BYTE  : extract<int8_t> (in, out, from, to, into); break;
    case NC_SHORT : extract<int16_t>(in, out, from, to, into); break;
    case NC_LONG  : extract<int32_t>(in, out, from, to, into); break;
    case NC_FLOAT : extract<float>  (in, out, from, to, into); break;
    case NC_DOUBLE: extract<double> (in, out, from, to, into); break;
    default: break;
    }

    return 0;
}


int main(const int argc, char* argv[])
{
  try
  {
    run(argc, argv);
  }
  catch(std::runtime_error& e)
  {
    std::clog
      << "terminate called after throwing an instance of "
      << "'std::runtime_error'\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
  catch(std::exception& e)
  {
    std::clog
      << "terminate called after throwing an exception\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
}
