/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  Gyulassy.C
 *
 *  Computes basins using the Gyulassy-Bremer-Pascucci method.
 *
 *  Olaf Delgado-Friedrichs jun 15
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
#include "VertexMap.hpp"
#include "volume_io.hpp"


using namespace anu_am::diamorse;

typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef VertexMap<CubicalComplex, int8_t> Mask;


typedef std::vector<std::pair<Cell, Value> > Distribution;


struct VertexInfo
{
    size_t basin;
    size_t countDown;
    Distribution weights;
    bool active;

    VertexInfo()
        : active(false)
    {
    }
};

typedef boost::shared_ptr<VertexInfo> InfoPtr;
typedef VertexMap<CubicalComplex, InfoPtr> InfoMap;


class VertexQueue {
    struct Entry {
        Cell vertex;
        Value value;
        size_t order;

        Entry(Cell const vertex, Value const value, size_t const order)
            : vertex(vertex),
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
    VertexQueue()
        : _pq(Comparator()),
          _count(0)
    {
    }

    bool empty() const
    {
        return _pq.empty();
    }

    void push(Cell const& vertex, Value const& value)
    {
        _pq.push(Entry(vertex, value, ++_count));
    }

    void pop()
    {
        _pq.pop();
    }

    Cell top() const
    {
        return _pq.top().vertex;
    }
};


Distribution combine(Distribution const& c, InfoMap const& info)
{
    std::map<Cell, Value> tmp;

    for (size_t i = 0; i < c.size(); ++i)
    {
        Distribution const& w = info(c.at(i).first)->weights;
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


Distribution derivedDistribution(
    Cell const v,
    std::vector<Cell> const& nbrs,
    InfoMap const& info,
    Scalars const& values)
{
    size_t const n = nbrs.size();
    std::vector<Value> weights(n);
    Value sum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        Value const w = std::max(0.0f, values(v) - values(nbrs.at(i)));
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
        c.at(i) = std::make_pair(nbrs.at(i), weights.at(i) / sum);

    return combine(c, info);
}


Cell opposite(Cell const e, Cell const v, Facets const& I)
{
    return I(e, I(e,0) == v);
}


bool isLocalMinimum(
    Cell const v,
    Facets const& I,
    Facets const& coI,
    Scalars const& scalars)
{
    float const h = scalars(v);

    size_t const m = I.count(v);
    for (size_t i = 0; i < m; ++i)
    {
        Cell const e = I(v, i);
        size_t const k = coI.count(e);
        for (size_t j = 0; j < k; ++j)
        {
            Cell const w = coI(e, j);
            if (w != v and scalars(w) <= h)
                return false;
        }
    }

    return true;
}


Mask boundaries(
    CubicalComplex const& complex,
    Scalars const& scalars)
{
    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
    VertexQueue queue;
    InfoMap info(complex);
    Mask result(complex);

    for (Cell v = 0; v < complex.cellIdLimit(); ++v)
    {
        if (complex.isCell(v) and complex.cellDimension(v) == 0 and 
            isLocalMinimum(v, coI, I, scalars))
        {
            queue.push(v, scalars(v) + 1);
            info.set(v, InfoPtr(new VertexInfo()));
            result.set(v, 2);
        }
    }

    while (not queue.empty())
    {
        Cell const v = queue.top();
        queue.pop();

        size_t const m = coI.count(v);
        std::vector<size_t> lower;
        for (size_t i = 0; i < m; ++i)
        {
            Cell const w = opposite(coI(v, i), v, I);
            if (not info(w))
            {
                queue.push(w, scalars(w));
                info.set(w, InfoPtr(new VertexInfo()));
            }
            else if (info(w)->active)
                lower.push_back(w);
        }

        info(v)->active = true;
        info(v)->countDown = m - lower.size();

        if (lower.size() == 0)
        {
            info(v)->weights = Distribution(1);
            info(v)->weights.at(0) = std::make_pair(v, 1);
            info(v)->basin = v;
        }
        else
        {
            Distribution const dist =
                derivedDistribution(v, lower, info, scalars);
            info(v)->weights = dist;

            std::set<size_t> seen;
            for (size_t i = 0; i < lower.size(); ++i)
                seen.insert(info(lower.at(i))->basin);

            std::pair<size_t, Value> best;
            size_t count = 0;
            for (size_t i = 0; i < dist.size(); ++i)
                if (seen.count(dist.at(i).first) > 0)
                    if (++count == 1 or dist.at(i).second > best.second)
                        best = dist.at(i);

            info(v)->basin = best.first;
        }

        for (size_t i = 0; i < lower.size(); ++i)
        {
            Cell const w = lower.at(i);
            if (info(v)->basin != info(w)->basin and scalars(v) <= 0)
                result.set(std::min(v, w), 1);

            if (--(info(w)->countDown) == 0)
                info.set(w, InfoPtr());
        }
    }

    return result;
}


int main(const int argc, char* argv[])
{
    namespace js = anu_am::json;

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
    CubicalComplex complex(dims.at(0), dims.at(1), dims.at(2));

    Scalars::DataPtr scalarData = readVolumeData<Value>(infile);
    Scalars scalars(complex, scalarData);

    // Process the data.
    Mask const out = boundaries(complex, scalars);

    // Generate metadata
    std::string const parentID = guessDatasetID(infile, info.attributes());
    std::string const thisID   = derivedID(parentID, "segmented", "BBG");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Basin Boundaries via Gyulassi-Bremer-Pascucci Method")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , js::Object());

    std::string const description = js::toString(fullSpec, 2);

    // Write the results to the output file
    writeVolumeData(
        out.data(), outfile, "segmented", dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(inheritableAttributes(info.attributes()))
        .datasetID(thisID)
        .description(description)
        .computeHistogram(false));
}
