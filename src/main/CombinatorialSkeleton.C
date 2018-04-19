/** -*-c++-*-
 *
 *  Copyright 2018 The Australian National University
 *
 *  CombinatorialSkeleton.C
 *
 *  Reads a Morse vector field and computes the combinatorial Morse skeleton.
 *
 *  Olaf Delgado-Friedrichs apr 18
 *
 */

#include <fstream>
#include <getopt.h>
#include <stdint.h>
#include <sstream>

#include "collections.hpp"
#include "chainComplexExtraction.hpp"
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
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;


void printWithPrefix(
    std::ostream& out,
    std::string const text,
    std::string const prefix)
{
    size_t pos = 0, next = text.find('\n');
    while (next != std::string::npos)
    {
        out << prefix << text.substr(pos, (next + 1) - pos);
        pos = next + 1;
        next = text.find('\n', pos);
    }
    if (pos < text.size())
        out << prefix << text.substr(pos) << std::endl;
}


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

    typedef std::vector<std::pair<Cell, int> > Boundary;

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
    std::map<Cell, Boundary> chains = chainComplex(complex, field);

    std::vector<Cell> const critical =
        criticalCellsSorted(complex, scalars, field, threshold, dimension);

    std::stringstream ss;
    std::map<Cell, size_t> cellIndex;

    ss << "# Combinatorial skeleton for " << scalarPath
       << std::endl
       << "#   format: <index> <dimension> <value> <x> <y> <z> <boundary...>"
       << std::endl;

    for (size_t k = 0; k < critical.size(); ++k)
    {
        Cell const cell = critical.at(k);
        cellIndex[cell] = k + 1;

        std::vector<float> const p = complex.cellPosition(cell);
        int const dim = complex.cellDimension(cell);
        Value const val = cellValue(cell, scalars, vertices);
        
        ss << cellIndex.at(cell);
        ss << "   " << dim << "   " <<  val << "   " << p << "  ";

        std::vector<std::pair<Cell, int> > bnd = chains.at(cell);
        for (size_t i = 0; i < bnd.size(); ++i)
        {
            Cell const other = bnd.at(i).first;
            int const count = bnd.at(i).second;

            for (int j = 0; j < count; ++j)
                ss << " " << cellIndex.at(other);
        }
        ss << std::endl;
    }

    // Generate metadata
    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "combinatorial_skeleton",
                                           "SKC");

    std::string const outfile = (argc - optind > 2) ?
        argv[optind+2] : (stripTimestamp(thisID) + ".txt");

    js::Array const predecessors = js::Array
        (parentID)
        (guessDatasetID(scalarPath, readFileInfo(scalarPath).attributes()));

    js::Object const parameters = js::Object("threshold" , threshold);

    js::Object const description = js::Object
        ("id"          , thisID)
        ("process"     , "Combinatorial Skeleton")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("predecessors", predecessors)
        ("parameters"  , parameters);

    // Write the results to the output file
    std::ofstream ofs(outfile.c_str());
    ofs << ss.str();

    ofs << "#" << std::endl
        << "# Metadata:" << std::endl;

    Attributes const attr = inheritableAttributes(info.attributes());
    for (size_t i = 0; i < attr.size(); ++i)
    {
        std::string const key = attr.keyAt(i);
        ofs << "#" << std::endl
            << "#+ " << key << std::endl;
        if (key == "dataset_id")
            ofs << "#= " << thisID << std::endl;
        else if (key == "zdim_total")
            ofs << "#= " << dims.at(2) << std::endl;
        else if (key == "number_of_files")
            ofs << "#= 1" << std::endl;
        else if (key == "zdim_range")
            ofs << "#= 0, " << dims.at(2)-1 << std::endl;
        else
            printWithPrefix(ofs, attr(key).valuesAsString(), "#= ");
    }

    ofs << "#" << std::endl
        << "#+ " << "history_"+thisID << std::endl;
    ofs << js::toString(description, 2, "#= ") << std::endl;
}
