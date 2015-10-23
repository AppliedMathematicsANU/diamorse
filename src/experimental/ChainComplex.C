/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  ChainComplex.C
 *
 *  Reads a Morse vector field and produces the associated Morse chain complex.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <fstream>
#include <sstream>

#include "boost/format.hpp"

#include "collections.hpp"
#include "chainComplexExtraction.hpp"
#include "CubicalComplex.hpp"
#include "json.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "traversals.hpp"
#include "volume_io.hpp"


using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
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


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    typedef std::vector<std::pair<Cell, int> > Boundary;

    char* infile = argv[1];

    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT [OUTPUT]" << std::endl;
        return 1;
    }

    // Read the data for this process.
    NCFileInfo const info = readFileInfo(infile);
    Variable const var = findVolumeVariable(info);

    std::vector<size_t> const dims = readDimensions(info);
    CubicalComplex const complex(dims.at(0), dims.at(1), dims.at(2));

    Field::DataPtr fieldData = readVolumeData<FieldItem>(infile);
    Field field = Field(dims.at(0), dims.at(1), dims.at(2), fieldData);

    // Process the data.
    std::map<Cell, Boundary> chains = chainComplex(complex, field);

    std::stringstream ss;
    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
    {
        if (chains.count(cell) == 0)
            continue;
        std::vector<float> const p = complex.cellPosition(cell);
        int const dim = complex.cellDimension(cell);
        ss << p << " # dim " << dim << std::endl;

        std::vector<std::pair<Cell, int> > bnd = chains.at(cell);
        for (size_t i = 0; i < bnd.size(); ++i)
        {
            Cell const other = bnd.at(i).first;
            int const count = bnd.at(i).second;

            std::vector<float> const q = complex.cellPosition(other);
            for (int j = 0; j < count; ++j)
                ss << p << " -> " << q << std::endl;
        }
    }

    // Generate metadata
    std::string const parentID = guessDatasetID(infile, info.attributes());
    std::string const thisID   = derivedID(parentID, "chain_complex", "CCX");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".txt");

    js::Object const description = js::Object
        ("id"          , thisID)
        ("process"     , "Chain Complex Extraction")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , js::Object());

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
