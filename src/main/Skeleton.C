/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  Skeleton.C
 *
 *  Reads a Morse vector field and computes the Morse skeleton.
 *
 *  Olaf Delgado-Friedrichs jun 14
 *
 */

#include <getopt.h>
#include <stdint.h>

#include "CubicalComplex.hpp"
#include "json.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "traversals.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"


using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef VertexMap<CubicalComplex, int8_t> Mask;
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;


Mask skeleton(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field const& field,
    float const threshold = 0,
    int const dimension = 3)
{
    Field::Vectors V = field.V();
    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());
    Mask skel(complex, 0);


    for (Cell s = 0; s <= complex.cellIdLimit(); ++s)
    {
        if (not (complex.isCell(s) and
                 field.isCritical(s) and
                 complex.cellDimension(s) <= dimension)
            )
            continue;

        bool good = true;
        size_t const m = vertices.count(s);
        for (size_t i = 0; i < m; ++i)
            if (scalars.get(vertices(s, i)) > threshold)
                good = false;

        if (not good)
            continue;

        std::vector<std::pair<Cell, Cell> > t = flowTraversal(s, V, I);
        for (size_t j = 0; j < t.size(); ++j)
        {
            Cell const c = V(t.at(j).second);
            size_t const m = vertices.count(c);
            for (size_t i = 0; i < m; ++i)
                skel.set(vertices(c, i), 1);
        }
    }

    return skel;
}


void usage(char *name)
{
    std::cerr << "Usage: " << name 
              << " [-t FLOAT] SCALARS FIELD [OUTPUT]"
              << std::endl
              << "Options:"
              << std::endl
              << "    -t threshold (default 0)"
              << std::endl;
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    int c;
    float threshold = 0;
    int dimension = 3;

    while ((c = getopt (argc, argv, "d:t:")) != -1)
    {
        switch (c)
        {
        case 'd':
            dimension = atoi(optarg);
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
    Mask const out = skeleton(complex, scalars, field, threshold, dimension);

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "segmented", "SKL");

    std::string const outfile =
        (argc - optind > 2) ? argv[optind+2] : (stripTimestamp(thisID) + ".nc");

    js::Array const predecessors = js::Array
        (parentID)
        (guessDatasetID(scalarPath, readFileInfo(scalarPath).attributes()));

    js::Object const parameters = js::Object("threshold" , threshold);

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Skeleton")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", predecessors)
        ("parameters"  , parameters);

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting gradient vector field to the output file
    writeVolumeData(
        out.data(), outfile, "segmented", dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));
}
