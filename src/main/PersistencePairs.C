/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  PersistencePairs.C
 *
 *  Reads a gradient vector field and computes persistence pairs.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <fstream>
#include <sstream>

#include "chainComplexExtraction.hpp"
#include "CubicalComplex.hpp"
#include "json.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "persistence.hpp"
#include "SimpleComplex.hpp"
#include "traversals.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"


using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;
typedef std::vector<std::pair<Cell, int> > Boundary;


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


template<typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& v)
{
    out << std::fixed << std::setprecision(1);

    for (size_t i = 0; i < v.size(); ++i)
        out << (i > 0 ? " " : "") << std::setw(6) << v.at(i);
    return out;
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
        Value const sv = cellValue(v, _scalars, _vertices);
        Value const sw = cellValue(w, _scalars, _vertices);

        return sv < sw or 
            (sv == sw and
             _complex.cellDimension(v) < _complex.cellDimension(w));
    }

private:
    CubicalComplex const& _complex;
    Vertices const _vertices;
    Scalars const& _scalars;
};


std::vector<Cell> criticalCellsSorted(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field const& field)
{
    std::vector<Cell> critical;

    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell) and field.isCritical(cell))
            critical.push_back(cell);

    std::stable_sort(critical.begin(), critical.end(),
                     Comparator(complex, scalars));

    return critical;
}


SimpleComplex simpleChainComplex(
    CubicalComplex const& complex,
    Scalars const& scalars,
    std::map<Cell, Boundary> const& chains,
    std::vector<Cell> const& sources)
{
    size_t const n = sources.size();
    Vertices const vertices(complex.xdim(), complex.ydim(), complex.zdim());

    std::map<Cell, size_t> index;
    for (size_t i = 0; i < n; ++i)
        index[sources.at(i)] = i;

    std::vector<unsigned int> dims;
    std::vector<float> values;
    std::vector<std::vector<Cell> > faceLists;

    for (size_t i = 0; i < n; ++i)
    {
        Cell const v = sources.at(i);
        dims.push_back(complex.cellDimension(v));
        values.push_back(cellValue(v, scalars, vertices));

        Boundary const flin = chains.at(v);
        std::vector<Cell> flout;
        for (size_t j = 0; j < flin.size(); ++j)
        {
            std::pair<Cell, int> p = flin.at(j);
            for (int k = 0; k < p.second; ++k)
                flout.push_back(index.at(p.first));
        }

        faceLists.push_back(flout);
    }

    return SimpleComplex(dims, values, faceLists);
}


int run(int argc, char* argv[])
{
    namespace js = anu_am::json;

    char* scalarPath = argv[1];
    char* fieldPath  = argv[2];

    if (argc < 3)
    {
        std::cerr << "Usage:" << argv[0] << " SCALARS FIELD [OUTPUT]"
                  << std::endl;
        return 1;
    }

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

    std::vector<Cell> const sources =
        criticalCellsSorted(complex, scalars, field);
    size_t const n = sources.size();

    SimpleComplex const simple = simpleChainComplex(
        complex, scalars, chains, sources);

    std::vector<Pairing<Cell> > const pairs = persistencePairing(simple);

    // Generate metadata
    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "persistence", "PP");

    std::string const outfile =
        argc > 3 ? argv[3] : (stripTimestamp(thisID) + ".txt");

    js::Array const predecessors = js::Array
        (parentID)
        (guessDatasetID(scalarPath, readFileInfo(scalarPath).attributes()));

    js::Object const description = js::Object
        ("id"          , thisID)
        ("process"     , "Critical Cell Persistence Pairs")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", predecessors)
        ("parameters"  , js::Object());

    // Write data
    std::stringstream tmp;

    tmp << "# Persistence pairs for " << scalarPath
        << std::endl
        << "#   format: <birth> <death> <dimension> <creator xyz> <destructor xyz>"
        << std::endl;

    for (size_t i = 0; i < pairs.size(); ++i)
    {
        size_t const j = pairs.at(i).partner;

        if (j > i)
        {
            Cell const v = sources.at(i);
            Cell const w = j >= n ? sources.at(i) : sources.at(j);

            tmp << std::fixed << std::setprecision(6);

            tmp << std::setw(12) << cellValue(v, scalars, vertices) << " "
                << std::setw(12);

            if (w == v)
                tmp << "inf";
            else
                tmp << cellValue(w, scalars, vertices);

            tmp << "    "
                << complex.cellDimension(v) << "    "
                << complex.cellPosition(v) << "    ";

            if (w == v)
                tmp << "   -      -      -  ";
            else
                tmp << complex.cellPosition(w);

            tmp << std::endl;
        }
    }

    std::ofstream ofs(outfile.c_str());
    ofs << tmp.str();

    // Write metadata
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
