/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  Classification.C
 *
 *  Classifies cells in a cubical complex by their stable set membership.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <map>
#include <queue>
#include <sstream>
#include <stdint.h>
#include <vector>

#include "collections.hpp"
#include "CubicalComplex.hpp"
#include "json.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "Set.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"

using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef MorseVectorField<PackedMap> Field;
typedef Field::DataPtr::element_type::value_type FieldItem;

typedef uint32_t Label;
typedef VertexMap<CubicalComplex, int32_t> Labelling;


struct Legend {
    std::map<Set<Cell>, Label> label;
    std::vector<Set<Cell> > lookup;

    Set<Cell> operator()(Label const lbl) const
    {
        return lookup.at(lbl);
    }

    Label labelFor(Set<Cell> const& s) {
        if (lookup.size() == 0)
            lookup.push_back(Set<Cell>());

        if (label.count(s) == 0)
        {
            label[s] = lookup.size();
            lookup.push_back(s);
        }
        return label[s];
    }

    Label labelFor(Cell const c) {
        return labelFor(Set<Cell>(c));
    }

    size_t size() const
    {
        return lookup.size();
    }
};


std::pair<std::vector<Label>, Legend> markStable(
    CubicalComplex const& complex,
    Field const& field)
{
    Cell const n = complex.cellIdLimit();
    Field::Vectors coV = field.coV();
    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);

    std::queue<Cell> q;
    std::vector<Label> result(n);
    Legend legend;

    for (int d = 0; d <= 3; ++d)
    {
        for (Cell s = 0; s < n; ++s)
        {
            if (not complex.isCell(s)
                or not field.isCritical(s)
                or complex.cellDimension(s) != d)
            {
                continue;
            }

            result.at(s) = legend.labelFor(s);
            q.push(s);

            while (not q.empty())
            {
                Cell const cell = q.front();
                q.pop();

                int const nrCofacets = coI.count(cell);
                for (int i = 0; i < nrCofacets; ++i)
                {
                    Cell const a = coI(cell, i);

                    if (not complex.isCell(a)
                        or field.isCritical(a)
                        or result.at(a) != 0
                        or not coV.defined(a))
                    {
                        continue;
                    }

                    Cell const b = coV(a);
                    Set<Cell> sources = Set<Cell>();

                    int const nrFacets = I.count(a);
                    for (int j = 0; j < nrFacets; ++j)
                    {
                        Cell const t = I(a, j);
                        if (t == b)
                            continue;

                        Label const lbl = result.at(t);
                        if (lbl == 0) {
                            sources = Set<Cell>();
                            break;
                        }

                        Cell src = legend(lbl).elements().at(0);
                        if (complex.cellDimension(src) == d)
                            sources = sources + legend(lbl);
                    }

                    if (sources.size() > 0) {
                        result.at(a) = result.at(b) = legend.labelFor(sources);
                        q.push(a);
                        q.push(b);
                    }
                }
            }
        }
    }

    int count = 0;
    for (Cell s = 0; s < n; ++s)
        if (complex.isCell(s) and result.at(s) == 0)
            ++count;

    if (count > 0)
        std::cerr << ">>> " << count << " cells left unmarked <<<" << std::endl;

    return std::make_pair(result, legend);
}


int main(int argc, char* argv[])
{
    namespace js = anu_am::json;

    char* fieldPath  = argv[1];

    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " FIELD [OUTPUT]"
                  << std::endl;
        return 1;
    }

    // -- read the input data

    NCFileInfo const info = readFileInfo(fieldPath);
    Variable const var = findVolumeVariable(info);

    std::vector<size_t> const dims = readDimensions(info);
    CubicalComplex complex(dims.at(0), dims.at(1), dims.at(2));

    Field::DataPtr fieldData = readVolumeData<FieldItem>(fieldPath);
    Field field = Field(dims.at(0), dims.at(1), dims.at(2), fieldData);

    // -- do the computation

    std::pair<std::vector<Label>, Legend> const expansions =
        markStable(complex, field);

    std::vector<Label> const labels = expansions.first;
    Legend const legend = expansions.second;

    // -- format the legend

    std::stringstream ss;

    for (size_t i = 0; i < legend.size(); ++i)
    {
        Set<Cell> const s = legend(i);
        size_t const n = s.size();

        ss << i << " " << n;
        if (n > 0)
            ss << " " << complex.cellDimension(s.elements().at(0));
        ss << std::endl;

        for (size_t j = 0; j < n; ++j)
            ss << "   " << complex.cellPosition(s.elements().at(j))
                << std::endl;
    }

    std::string const legendAsText = ss.str();

    // -- collect the labels

    Labelling out(complex, 0);
    int const dim = (complex.zdim() == 1) ? 2 : 3;
    for (Cell c = 0; c < complex.cellIdLimit(); ++c)
    {
        if (not complex.isCell(c) or complex.cellDimension(c) != dim)
            continue;
        Cell const v = complex.cellAt(floor(complex.cellX(c)),
                                      floor(complex.cellY(c)),
                                      floor(complex.cellZ(c)));
        out.set(v, labels.at(c));
    }

    // -- generate metadata

    std::string const parentID = guessDatasetID(fieldPath, info.attributes());
    std::string const thisID   = derivedID(parentID, "labels", "CLS");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Stable Set Membership Classification")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , js::Object());

    std::string const description = js::toString(fullSpec, 2);

    VolumeWriteOptions options = VolumeWriteOptions()
        .variableAttributes(Attributes("legend", legendAsText));

    // -- write the output file
    writeVolumeData(
        out.data(), outfile, "labels", dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(inheritableAttributes(info.attributes()))
        .variableAttributes(Attributes("legend", legendAsText))
        .datasetID(thisID)
        .description(description)
        .computeHistogram(false));
}
