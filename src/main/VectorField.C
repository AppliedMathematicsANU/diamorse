/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  VectorField.C
 *
 *  Computes the Morse vector field for a volume image.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include "CubicalComplex.hpp"
#include "json.hpp"
#include "vectorFieldExtraction.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"

using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef MorseVectorField<PackedMap> Field;


int run(const int argc, char* argv[])
{
    namespace js = anu_am::json;

    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT [OUTPUT]" << std::endl;
        return 1;
    }

    char* infile = argv[1];

    // Read the data for this process.
    NCFileInfo const info = readFileInfo(infile);
    Variable const var = findVolumeVariable(info);

    std::vector<size_t> const dims = readDimensions(info);
    CubicalComplex complex(dims.at(0), dims.at(1), dims.at(2));

    Scalars::DataPtr scalarData = readVolumeData<Value>(infile);
    Scalars scalars(complex, scalarData);

    // Process the data.
    Field field = Field(complex);
    fillMorseVectorField(complex, scalars, field);

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(infile, info.attributes());
    std::string const thisID   = derivedID(parentID, "vector_field", "GVF");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Discrete Morse Gradient Vector Field")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , js::Object());

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting gradient vector field to the output file
    writeVolumeData(
        field.data(), outfile, "vector_field",
        dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));

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
