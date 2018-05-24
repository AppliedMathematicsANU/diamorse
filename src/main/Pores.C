/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  Pores.C
 *
 *  Reads a Morse vector field and computes pore labels.
 *
 *  Olaf Delgado-Friedrichs jun 15
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
typedef VertexMap<CubicalComplex, int32_t> Labels;
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;


Value cellValue(Cell const v, Scalars const scalars, Vertices const vertices)
{
    Value val = scalars.get(vertices(v, 0));
    for (size_t i = 1; i < (size_t) vertices.count(v); ++i)
        val = std::max(val, scalars.get(vertices(v, i)));

    return val;
}


struct Comparator
{
    Comparator(CubicalComplex const& complex, Scalars const& scalars)
        : _complex(complex),
          _vertices(complex.xdim(), complex.ydim(), complex.zdim()),
          _scalars(scalars)
    {
    }

    bool operator()(Cell const v, Cell const w)
    {
        size_t const dv = _complex.cellDimension(v);
        size_t const dw = _complex.cellDimension(w);
        Value const sv = cellValue(v, _scalars, _vertices);
        Value const sw = cellValue(w, _scalars, _vertices);

        return dv < dw or (dv == dw and sv < sw);
    }

private:
    CubicalComplex const& _complex;
    Vertices const _vertices;
    Scalars const& _scalars;
};


std::vector<Cell> criticalCellsSorted(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field const& field,
    float const threshold,
    int const dimension)
{
    Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());
    std::vector<Cell> critical;

    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell) and
            field.isCritical(cell) and
            complex.cellDimension(cell) <= dimension and
            cellValue(cell, scalars, vertices) <= threshold
            )
            critical.push_back(cell);

    std::stable_sort(critical.begin(), critical.end(),
                     Comparator(complex, scalars));

    return critical;
}


Labels pores(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field const& field,
    float const threshold = 0)
{
    Field::Vectors coV = field.coV();
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
    Labels pores(complex, 0x7fffffff);

    std::vector<Cell> const sources =
        criticalCellsSorted(complex, scalars, field, threshold, 0);

    for (size_t i = 0; i < sources.size(); ++i)
    {
        Cell const s = sources.at(i);
        pores.set(s, i+1);
        
        std::vector<std::pair<Cell, Cell> > t = flowTraversal(s, coV, coI);
        for (size_t j = 0; j < t.size(); ++j)
        {
            Cell const b = t.at(j).second;
            Cell const c = coV(b);
            if (b != c and scalars(c) <= threshold)
                pores.set(c, i+1);
        }
    }

    return pores;
}


void usage(char *name)
{
    std::cerr << "Usage: " << name 
              << " [-t FLOAT] SCALARS FIELD [OUTPUT]"
              << std::endl
              << "Options:"
              << std::endl
              << "    -t pore inclusion threshold (default 0)"
              << std::endl;
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    int c;
    float threshold = 0;

    while ((c = getopt (argc, argv, "t:")) != -1)
    {
        switch (c)
        {
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
    Labels const out = pores(complex, scalars, field, threshold);

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "labels", "POR");

    std::string const outfile =
        (argc - optind > 2) ? argv[optind+2] : (stripTimestamp(thisID) + ".nc");

    js::Array const predecessors = js::Array
        (parentID)
        (guessDatasetID(scalarPath, readFileInfo(scalarPath).attributes()));

    js::Object const parameters = js::Object("threshold" , threshold);

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Pores")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", predecessors)
        ("parameters"  , parameters);

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting gradient vector field to the output file
    writeVolumeData(
        out.data(), outfile, "labels", dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));
}
