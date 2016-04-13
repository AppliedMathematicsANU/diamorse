/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  Simplify.C
 *
 *  Simplifies a gradient vector field via Morse cancellation.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <getopt.h>
#include <limits>
#include <sstream>
#include <stdlib.h>

#include "callables.hpp"
#include "collections.hpp"
#include "CubicalComplex.hpp"
#include "json.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "simplification.hpp"
#include "stringUtils.hpp"
#include "traversals.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"


using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;


void usage(char *name)
{
    std::cerr << "Usage: " << name 
              << " [-p FLOAT] [-s FLOAT] [-t FLOAT]"
              << " SCALARS FIELD [OUTPUT]"
              << std::endl
              << "Options:"
              << std::endl
              << "    -p persistence limit for feature cancellation" 
              << std::endl
              << "    -s size limit for feature cancellation"
              << std::endl
              << "    -t value level threshold to preserve"
              << std::endl;
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;
    namespace su = anu_am::stringutils;

    int c;
    float persistenceLimit = 1;
    float sizeLimit = 0;
    float threshold = -std::numeric_limits<float>::max();

    while ((c = getopt (argc, argv, "p:s:t:")) != -1)
    {
        switch (c)
        {
        case 'p':
            persistenceLimit = atof(optarg);
            break;
        case 's':
            sizeLimit = atof(optarg);
            break;
        case 't':
            threshold = atof(optarg);
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (argc - optind < 2)
    {
        usage(argv[0]);
        return 1;
    }

    char* scalarPath = argv[optind];
    char* fieldPath  = argv[optind + 1];

    // Read the data for this process.
    NCFileInfo const info = readFileInfo(fieldPath);
    Variable const var = findVolumeVariable(info);

    std::vector<size_t> dims = readDimensions(info);
    CubicalComplex complex(dims.at(0), dims.at(1), dims.at(2));
    Vertices vertices(dims.at(0), dims.at(1), dims.at(2));

    assert(dims == readDimensions(scalarPath));

    Scalars::DataPtr scalarData = readVolumeData<Value>(scalarPath);
    Scalars scalars(complex, scalarData);

    Field::DataPtr fieldData = readVolumeData<FieldItem>(fieldPath);
    Field field = Field(dims.at(0), dims.at(1), dims.at(2), fieldData);

    // Process the data.
    std::ostringstream ss;
    simplify(
        complex,
        field,
        withArgument(maxima(scalars, vertices)),
        mayCancel(complex, scalars, field,
                  persistenceLimit, sizeLimit, threshold),
        ss);
    std::string const output = ss.str();

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "vector_field", "SMP");

    std::string const outfile =
        (argc - optind > 2) ? argv[optind+2] : (stripTimestamp(thisID) + ".nc");

    js::Array const predecessors = js::Array
        (parentID)
        (guessDatasetID(scalarPath, readFileInfo(scalarPath).attributes()));

    js::Object const parameters = js::Object
        ("persistence_threshold" , persistenceLimit)
        ("feature_size_limit"    , sizeLimit)
        ("simplification_barrier", threshold);

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Gradient Vector Field Simplification")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", predecessors)
        ("parameters"  , parameters)
        ("output"      , su::split(output, '\n'));

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting gradient vector field to the output file
    writeVolumeData(
        field.data(), outfile, "vector_field",
        dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));
}
