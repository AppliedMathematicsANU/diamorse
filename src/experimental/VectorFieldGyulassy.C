/** -*-c++-*-
 *
 *  Copyright 2016 The Australian National University
 *
 *  VectorFieldGyulassy.C
 *
 *  Gradient vector field computed via the Gyulassy-Bremer-Pascucci method.
 *
 *  Olaf Delgado-Friedrichs apr 16
 *
 */


#include <algorithm>
#include <map>
#include <queue>
#include <stdint.h>
#include <vector>
#include <utility>

#include "boost/smart_ptr.hpp"

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
typedef VertexMap<CubicalComplex, Value> Scalars;

typedef std::vector<std::pair<Cell, Value> > Distribution;


template<class S>
class SmartMorseVectorField
{
    typedef MorseVectorField<S> BaseField;

    CubicalComplex      _complex;
    Facets              _I;
    Facets              _coI;
    BaseField           _field;
    std::vector<size_t> _nrPairableFacets;

    void update(Cell const v)
    {
        _nrPairableFacets.at(v) = 0;
        for (size_t i = 0; i < _coI.count(v); ++i)
            --_nrPairableFacets.at(_coI(v, i));
    }

    void pairUp(Cell const v, Cell const w)
    {
        _field.setPartner(v, w);
        update(v);
        update(w);
    }

public:
    typedef typename BaseField::DataPtr    DataPtr;
    typedef typename BaseField::Directions Directions;
    typedef typename BaseField::Vectors    Vectors;

    SmartMorseVectorField()
    {
    }

    SmartMorseVectorField(CubicalComplex const& complex)
    : _complex(complex),
      _I(complex.xdim(), complex.ydim(), complex.zdim(), false),
      _coI(complex.xdim(), complex.ydim(), complex.zdim(), true),
      _field(complex),
      _nrPairableFacets(complex.cellIdLimit())
    {
        for (Cell cell = 0; cell < complex.cellIdLimit(); ++cell)
            if (complex.isCell(cell))
                _nrPairableFacets.at(cell) = _I.count(cell);
    }

    Vectors const V() const
    {
        return _field.V();
    }

    Vectors const coV() const
    {
        return _field.coV();
    }

    Directions getDirection(Cell const v) const
    {
        return _field.getDirection(v);
    }

    bool isCritical(Cell const n) const
    {
        return _field.isCritical(n);
    }

    Cell getPartner(Cell const n) const
    {
        return _field.getPartner(n);
    }

    void setPartner(Cell const v, Cell const w)
    {
        pairUp(v, w);

        std::queue<Cell> q;
        if (_complex.cellDimension(v) > _complex.cellDimension(w))
            q.push(v);
        else
            q.push(w);

        while (not q.empty())
        {
            Cell const cell = q.front();
            q.pop();

            for (size_t i = 0; i < _coI.count(cell); ++i)
            {
                Cell const cof = _coI(cell, i);
                if (_nrPairableFacets.at(cof) != 1)
                    continue;
                for (size_t j = 0; j < _I.count(cof); ++j)
                {
                    Cell const f = _I(cof, j);
                    if (!_field.getDirection(f))
                    {
                        pairUp(cof, f);
                        q.push(cof);
                        q.push(f);
                        break;
                    }
                }
            }
        }
    }

    bool pointsOutward(Cell const n) const
    {
        return _field.pointsOutward(n);
    }

    DataPtr const data() const
    {
        return _field.data();
    }
};


typedef SmartMorseVectorField<PackedMap> Field;


struct CellInfo
{
    size_t countDown;
    Distribution weights;
    Set<Cell> endPoints;
    bool active;

    CellInfo()
        : active(false)
    {
    }
};

typedef boost::shared_ptr<CellInfo> InfoPtr;
typedef std::vector<InfoPtr> InfoMap;


class CellQueue {
    struct Entry {
        Cell cell;
        Value value;
        size_t order;

        Entry(Cell const cell, Value const value, size_t const order)
            : cell(cell),
              value(value),
              order(order)
        {
        }
    };

    struct Comparator
    {
        bool operator()(Entry const a, Entry const b) {
            if (a.value != b.value)
                return a.value > b.value;
            else
                return a.order > b.order;
        }
    };

    std::priority_queue<Entry, std::vector<Entry>, Comparator> _pq;
    size_t _count;

public:
    CellQueue()
        : _pq(Comparator()),
          _count(0)
    {
    }

    bool empty() const
    {
        return _pq.empty();
    }

    void push(Cell const& cell, Value const& value)
    {
        _pq.push(Entry(cell, value, ++_count));
    }

    void pop()
    {
        _pq.pop();
    }

    Cell top() const
    {
        return _pq.top().cell;
    }
};


Distribution combine(Distribution const& c, InfoMap const& info)
{
    std::map<Cell, Value> tmp;

    for (size_t i = 0; i < c.size(); ++i)
    {
        Distribution const& w = info.at(c.at(i).first)->weights;
        Value const f = c.at(i).second;

        for (size_t j = 0; j < w.size(); ++j)
        {
            Cell const v = w.at(j).first;
            if (tmp.count(v) == 0)
                tmp[v] = 0;
            tmp[v] += f * w.at(j).second;
        }
    }

    Distribution result(tmp.size());

    std::map<Cell, Value>::const_iterator iter = tmp.begin();
    for (size_t i = 0; i < tmp.size(); ++i, ++iter)
        result[i] = std::make_pair(iter->first, iter->second);

    return result;
}


Value cellAverage(Cell const v, Scalars const scalars, Vertices const vertices)
{
    size_t const n = vertices.count(v);

    Value sum = 0;
    for (size_t i = 0; i < n; ++i)
        sum += scalars.get(vertices(v, i));

    return sum / n;
}


Distribution derivedDistribution(
    Cell const v,
    std::vector<std::pair<Cell, Cell> > const& nbrs,
    InfoMap const& info,
    Scalars const& values,
    Vertices const& vertices)
{
    size_t const n = nbrs.size();
    Value const val0 = cellAverage(v, values, vertices);

    std::vector<Value> weights(n);
    Value sum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        Value const val = cellAverage(nbrs.at(i).first, values, vertices);
        Value const w = std::max(0.0f, val0 - val);
        weights.at(i) = w;
        sum += w;
    }

    if (sum == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            weights.at(i) = 1;
            sum += 1;
        }
    }

    Distribution c(n);
    for (size_t i = 0; i < n; ++i)
        c.at(i) = std::make_pair(nbrs.at(i).first, weights.at(i) / sum);

    return combine(c, info);
}


bool isLocalMinimum(
    Cell const v,
    Facets const& I,
    Facets const& coI,
    Vertices const& vertices,
    Scalars const& scalars)
{
    float const h = cellAverage(v, scalars, vertices);

    size_t const m = I.count(v);
    for (size_t i = 0; i < m; ++i)
    {
        Cell const e = I(v, i);
        size_t const k = coI.count(e);
        for (size_t j = 0; j < k; ++j)
        {
            Cell const w = coI(e, j);
            if (w != v and cellAverage(w, scalars, vertices) <= h)
                return false;
        }
    }

    return true;
}


void grow(
    int const dimension,
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field& outputField)
{
    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
    Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());
    CellQueue queue;
    InfoMap info(complex.cellIdLimit());
 
    for (Cell cell = 0; cell < complex.cellIdLimit(); ++cell)
    {
        if (complex.isCell(cell) and
            complex.cellDimension(cell) == dimension and
            outputField.getDirection(cell) == 0 and
            isLocalMinimum(cell, coI, I, vertices, scalars))
        {
            info.at(cell) = InfoPtr(new CellInfo());
            queue.push(cell, cellAverage(cell, scalars, vertices) + 1);
        }
    }

    while (not queue.empty())
    {
        Cell const cell = queue.top();
        queue.pop();

        std::vector<std::pair<Cell, Cell> > lower;
        size_t nrHigher = 0;
        for (int i = 0; i < coI.count(cell); ++i)
        {
            Cell const cof = coI(cell, i);

            for (int j = 0; j < I.count(cof); ++j)
            {
                Cell const w = I(cof, j);

                if (not info.at(w))
                {
                    info.at(w) = InfoPtr(new CellInfo());
                    queue.push(w, cellAverage(cell, scalars, vertices));
                }
                else
                {
                    ++nrHigher;
                    if (info.at(w)->active)
                        lower.push_back(std::make_pair(w, cof));
                }
            }
        }

        info.at(cell)->active = true;
        info.at(cell)->countDown = nrHigher;

        if (lower.size() == 0)
        {
            info.at(cell)->weights = Distribution(1);
            info.at(cell)->weights.at(0) = std::make_pair(cell, 1);
            info.at(cell)->endPoints = Set<Cell>(cell);
            outputField.setPartner(cell, cell);
        }
        else
        {
            Distribution const dist =
                derivedDistribution(cell, lower, info, scalars, vertices);
            info.at(cell)->weights = dist;

            Set<Cell> allEndPoints;
            for (size_t i = 0; i < lower.size(); ++i) {
                Cell const w = lower.at(i).first;
                allEndPoints = allEndPoints + info.at(w)->endPoints;
            }

            //TODO determine partner for v and update endpoints here...

            // std::pair<size_t, Value> best;
            // size_t count = 0;
            // for (size_t i = 0; i < dist.size(); ++i)
            //     if (seen.count(dist.at(i).first) > 0)
            //         if (++count == 1 or dist.at(i).second > best.second)
            //             best = dist.at(i);

            // info.at(v)->basin = best.first;
        }

        for (size_t i = 0; i < lower.size(); ++i)
        {
            Cell const w = lower.at(i).first;
            if (--(info.at(w)->countDown) == 0)
                info.at(w) = InfoPtr();
        }
    }
}


void fillMorseVectorField(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field& outputField)
{
    for (int d = 0; d < complex.dimension(); ++d)
        grow(d, complex, scalars, outputField);
}


int main(const int argc, char* argv[])
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
    std::string const thisID   = derivedID(parentID, "vector_field", "VFG");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Discrete Gradient Vector Field via Gyulassy et al")
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
        .computeHistogram(false)
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));
}
