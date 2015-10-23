/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  SignedEuclideanDistanceTransform.C
 *
 *  Takes a segmented image and computes the signed Euclidean distance
 *  for each voxel using the Hirata/Meijster algorithm.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <stdint.h>
#include <vector>

#include "CubicalComplex.hpp"
#include "json.hpp"
#include "netcdf.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"

using namespace anu_am::diamorse;


typedef CubicalComplex::cell_id_type Cell;
typedef int8_t PhaseValue;
typedef float DistanceValue;
typedef VertexMap<CubicalComplex, PhaseValue> Phases;
typedef VertexMap<CubicalComplex, DistanceValue> Distances;


float const MASK_VALUE = 1.0e30f;


inline double sq(double const x)
{
    return x * x;
}


inline double sq_g(std::vector<double> const& g, size_t const i)
{
    return sq(g.at(i));
}


inline double f(std::vector<double> const& g, double const u, size_t const i)
{
    return sq(u - i) + sq_g(g, i);
}


std::vector<double> positives(std::vector<double> const& a)
{
    std::vector<double> out;
    for (size_t i = 0; i < a.size(); ++i)
        out.push_back(a.at(i) > 0 ? a.at(i) : 0);
    return out;
}


std::vector<double> negatives(std::vector<double> const& a)
{
    std::vector<double> out;
    for (size_t i = 0; i < a.size(); ++i)
        out.push_back(a.at(i) < 0 ? -a.at(i) : 0);
    return out;
}


std::vector<double> combine(
    std::vector<double> const& negatives,
    std::vector<double> const& positives)
{
    std::vector<double> out;
    for (size_t i = 0; i < negatives.size(); ++i)
        out.push_back(negatives.at(i) ? -negatives.at(i) : positives.at(i));
    return out;
}


std::vector<double> init(std::vector<double> const& a)
{
    size_t const n = a.size();
    std::vector<double> out;

    out.push_back(a.at(0) ? n : 0);
    for (size_t i = 1; i < n; ++i)
        out.push_back(a.at(i) ? 1 + out.at(i-1) : 0);

    for (size_t i = n-1; i > 0; --i)
        if (out.at(i-1) > out.at(i))
            out.at(i-1) = 1 + out.at(i);

    return out;
}


std::vector<double> propagate(std::vector<double> const& g)
{
    if (g.size() <= 1)
        return g;

    size_t const m = g.size();
    size_t q = 1;
    std::vector<size_t> s(m+1);
    std::vector<double> t(m+1);

    s.at(q) = 0;
    t.at(q) = 0;

    for (size_t u = 1; u < m; ++u)
    {
        while (q > 0 && f(g, t.at(q), s.at(q)) > f(g, t.at(q), u))
            --q;
        if (q < 1)
        {
            q = 1;
            s.at(q) = u;
            t.at(q) = 0;
        }
        else
        {
            size_t const i = s.at(q);
            double const w =
                (sq(u) - sq(i) + sq_g(g, u) - sq_g(g, i)) / (2 * (u - i));

            if (w < m) {
                ++q;
                s.at(q) = u;
                t.at(q) = w;
            }
        }
    }

    std::vector<double> result(m);
    for (size_t u = m; u > 0; --u)
    {
        while (q > 1 && u-1 < t.at(q))
            --q;
        result.at(u-1) = sqrt(f(g, u-1, s.at(q)));
    }

    return result;
}


void classify(
    CubicalComplex const& complex,
    Phases const& phases,
    Distances & distances)
{
    size_t const xmax = complex.xdim();
    size_t const ymax = complex.ydim();
    size_t const zmax = complex.zdim();

    for (size_t z = 0; z < zmax; ++z)
    {
        for (size_t y = 0; y < ymax; ++y)
        {
            for (size_t x = 0; x < xmax; ++x)
            {
                Cell v = complex.cellAt(x, y, z);
                PhaseValue const val = phases(v);
                distances.set(v, (val > 0) - (val == 0));
            }
        }
    }
}


void compute(CubicalComplex const& complex, Distances & distances)
{
    size_t const xmax = complex.xdim();
    size_t const ymax = complex.ydim();
    size_t const zmax = complex.zdim();

    for (size_t z = 0; z < zmax; ++z)
    {
        for (size_t y = 0; y < ymax; ++y)
        {
            std::vector<double> g(xmax);
            for (size_t x = 0; x < xmax; ++x)
                g.at(x) = distances.get(complex.cellAt(x, y, z));

            std::vector<double> const d = combine(
                init(negatives(g)),
                init(positives(g)));

            for (size_t x = 0; x < xmax; ++x)
                distances.set(complex.cellAt(x, y, z), d.at(x));
        }
    }

    for (size_t z = 0; z < zmax; ++z)
    {
        for (size_t x = 0; x < xmax; ++x)
        {
            std::vector<double> g(ymax);
            for (size_t y = 0; y < ymax; ++y)
                g.at(y) = distances.get(complex.cellAt(x, y, z));

            std::vector<double> const d = combine(
                propagate(negatives(g)),
                propagate(positives(g)));

            for (size_t y = 0; y < ymax; ++y)
                distances.set(complex.cellAt(x, y, z), d.at(y));
        }
    }    

    for (size_t y = 0; y < ymax; ++y)
    {
        for (size_t x = 0; x < xmax; ++x)
        {
            std::vector<double> g(zmax);
            for (size_t z = 0; z < zmax; ++z)
                g.at(z) = distances.get(complex.cellAt(x, y, z));

            std::vector<double> const d = combine(
                propagate(negatives(g)),
                propagate(positives(g)));

            for (size_t z = 0; z < zmax; ++z)
                distances.set(complex.cellAt(x, y, z), d.at(z));
        }
    }

    for (size_t z = 0; z < zmax; ++z)
    {
        for (size_t y = 0; y < ymax; ++y)
        {
            for (size_t x = 0; x < xmax; ++x)
            {
                Cell const v = complex.cellAt(x, y, z);
                if (distances(v) > 0)
                    distances.set(v, distances(v) - 0.5);
                else if (distances(v) < 0)
                    distances.set(v, distances(v) + 0.5);
                else
                    distances.set(v, MASK_VALUE);
            }
        }
    }
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

    Phases::DataPtr phaseData = readVolumeData<PhaseValue>(infile);
    Phases phases(complex, phaseData);

    // Process the data.
    Distances dist(complex);
    classify(complex, phases, dist);
    compute(complex, dist);

    // Generate metadata to include with the output data
    std::string const parentID = guessDatasetID(infile, info.attributes());
    std::string const thisID   = derivedID(parentID, "tomo_float", "SEDT");

    std::string const outfile =
        argc > 2 ? argv[2] : (stripTimestamp(thisID) + ".nc");

    js::Object const fullSpec = js::Object
        ("id"          , thisID)
        ("process"     , "Signed Euclidean Distance Transform (SEDT)")
        ("sourcefile"  , __FILE__)
        ("revision"    , js::Object("id", GIT_REVISION)("date", GIT_TIMESTAMP))
        ("parent"      , parentID)
        ("predecessors", js::Array(parentID))
        ("parameters"  , js::Object());

    std::string const description = js::toString(fullSpec, 2);

    // Write the resulting gradient vector field to the output file
    writeVolumeData(
        dist.data(), outfile, "tomo_float", dims.at(0), dims.at(1), dims.at(2),
        VolumeWriteOptions()
        .fileAttributes(info.attributes())
        .datasetID(thisID)
        .description(description));
}
