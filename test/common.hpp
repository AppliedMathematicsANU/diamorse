/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  common.hpp
 *
 *  Common code for generative testing of Morse theory algorithma
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_TEST_HPP
#define ANU_AM_DIAMORSE_TEST_HPP

#include <assert.h>
#include <iostream>
#include <map>
#include <stack>

#include "generative.hpp"

#include "chainComplexExtraction.hpp"
#include "CubicalComplex.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "persistence.hpp"
#include "SimpleComplex.hpp"
#include "vectorFieldExtraction.hpp"
#include "VertexMap.hpp"

namespace anu_am
{
namespace diamorse
{


typedef CubicalComplex::cell_id_type Cell;
typedef float Value;
typedef VertexMap<CubicalComplex, Value> Scalars;
typedef MorseVectorField<PackedMap> Field;
typedef std::vector<std::pair<Cell, int> > Boundary;


// === Infrastructure for generative testing on volume data instances.

struct VolumeData
{
    CubicalComplex complex;
    Scalars scalars;
    Field field;

    VolumeData(CubicalComplex const& complex,
               Scalars const& scalars,
               Field const& field)
        : complex(complex),
          scalars(scalars),
          field(field)
    {
    }
};

std::ostream& operator<<(std::ostream& out, VolumeData const& inst)
{
    CubicalComplex const& c = inst.complex;
    Scalars const& s = inst.scalars;

    out << c.xdim() << "x" << c.ydim() << "x" << c.zdim() 
        << " [" << *s.data() << "]";

    return out;
}


VolumeData makeVolumeData(
    CubicalComplex const& complex,
    Scalars const& scalars)
{
    Field field(complex);
    fillMorseVectorField(complex, scalars, field);

    return VolumeData(complex, scalars, field);
}


template<size_t N>
VolumeData fixedVolumeData(
    int const xd, int const yd, int const zd,
    Value const(&input)[N])
{
    assert(xd * yd * zd == N);
    CubicalComplex complex(xd, yd, zd);

    Scalars::DataPtr data(new Scalars::DataPtr::element_type(N));
    for (size_t i = 0; i < N; ++i)
        data->at(i) = input[i];

    return makeVolumeData(complex, Scalars(complex, data));
}


VolumeData randomVolumeData(int const size)
{
    using namespace anu_am::generative;

    int const d = sqrt(size);
    CubicalComplex complex(randomInt(d), randomInt(d), randomInt(d));

    size_t N = complex.xdim() * complex.ydim() * complex.zdim();
    Scalars::DataPtr data(new Scalars::DataPtr::element_type(N));
    for (size_t i = 0; i < N; ++i)
        data->at(i) = randomFloat(5.0, 0.0);

    return makeVolumeData(complex, Scalars(complex, data));
}


std::vector<VolumeData> shrinkVolumeData(VolumeData const& instance)
{
    CubicalComplex const& complex = instance.complex;
    Scalars::DataPtr data0 = instance.scalars.data();
    size_t const dim0[] = { complex.xdim(), complex.ydim(), complex.zdim() };

    std::vector<VolumeData> result;

    for (int axis = 0; axis < 3; ++axis)
    {
        for (int offset = 0; offset <= 1; ++offset)
        {
            size_t dim[] = { dim0[0], dim0[1], dim0[2] }; 
            int off[] = { 0, 0, 0 };
            --dim[axis];
            off[axis] = offset;

            CubicalComplex c(dim[0], dim[1], dim[2]);
            size_t N = dim[0] * dim[1] * dim[2];
            Scalars::DataPtr data(new Scalars::DataPtr::element_type(N));

            for (size_t x = 0; x < dim[0]; ++x)
            {
                size_t const x0 = x + off[0];
                for (size_t y = 0; y < dim[1]; ++y)
                {
                    size_t const y0 = y + off[1];
                    for (size_t z = 0; z < dim[2]; ++z)
                    {
                        size_t const z0 = z + off[2];
                        size_t const p = (z * dim[1] + y) * dim[0] + x;
                        size_t const p0 = (z0 * dim0[1] + y0) * dim0[0] + x0;
                        data->at(p) = data0->at(p0);
                    }
                }
            }

            result.push_back(makeVolumeData(c, Scalars(c, data)));
        }
    }

    return result;
}


template<typename P>
anu_am::generative::Result checkWithVolumeData(
    P const& predicate,
    int const N = 100)
{
    return checkPredicate(predicate, randomVolumeData, shrinkVolumeData, N);
}


template<typename P>
struct ForAllCells
{
    typedef anu_am::generative::Result Result;
    typedef Result result_type;
    typedef VolumeData argument_type;

    P const& predicate;

    ForAllCells(P const& predicate)
        : predicate(predicate)
    {
    }

    Result operator()(VolumeData const& candidate) const
    {
        CubicalComplex const& c = candidate.complex;

        for (Cell cell = 0; cell <= c.cellIdLimit(); ++cell)
        {
            if (not c.isCell(cell))
                continue;

            Result r = predicate(cell, candidate);
            if (!r)
            {
                std::stringstream msg;
                msg << "At cell " << c.cellPosition(cell) << ": " << r.cause();
                return anu_am::generative::failure(msg.str());
            }
        }
        return anu_am::generative::success();
    }
};

template<typename P>
ForAllCells<P> forAllCells(P const& predicate)
{
    return ForAllCells<P>(predicate);
}


// === Shared helpers.

template<typename F, typename A1>
class Bind
{
    F const& fn_;
    A1 const& arg1_;

public:
    typedef typename anu_am::generative::function_traits<F> traits;
    typedef typename traits::result_type result_type;

    Bind(F const& fn, A1 const& arg1)
        : fn_(fn),
          arg1_(arg1)
    {
    }

    template<typename A2>
    result_type operator()(A2 const& arg2) const
    {
        return fn_(arg1_, arg2);
    }

    template<typename A2, typename A3>
    result_type operator()(A2 const& arg2, A3 const& arg3) const
    {
        return fn_(arg1_, arg2, arg3);
    }
};

template<typename F, typename A>
Bind<F, A> bind(F const& fn, A const& arg)
{
    return Bind<F, A>(fn, arg);
}


std::vector<Cell> criticalCells(VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;

    std::vector<Cell> result;

    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell) and field.isCritical(cell))
            result.push_back(cell);

    return result;
}


std::vector<Cell> allCells(VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;

    std::vector<Cell> result;

    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell))
            result.push_back(cell);

    return result;
}


std::vector<Cell>
convertedBoundary(
    std::map<Cell, Boundary> chains,
    Cell const v)
{
    std::vector<std::pair<Cell, int> > tmp = chains.at(v);
    std::vector<Cell> result;

    for (size_t i = 0; i < tmp.size(); ++i)
    {
        std::pair<Cell, int> p = tmp.at(i);
        for (int j = 0; j < p.second; ++j)
            result.push_back(p.first);
    }

    return result;
}


Value cellValue(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Cell const v)
{
    Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());

    Value val = scalars.get(vertices(v, 0));
    for (int i = 1; i < vertices.count(v); ++i)
        val = std::max(val, scalars.get(vertices(v, i)));

    return val;
}


struct Comparator
{
    Comparator(CubicalComplex const& complex, Scalars const& scalars)
        : complex_(complex),
          scalars_(scalars)
    {
    }

    bool operator()(Cell const v, Cell const w)
    {
        Value const sv = cellValue(complex_, scalars_, v);
        Value const sw = cellValue(complex_, scalars_, w);

        return sv < sw
            or (sv == sw
                and complex_.cellDimension(v) < complex_.cellDimension(w));
    }

private:
    CubicalComplex const& complex_;
    Scalars const& scalars_;
};


SimpleComplex convertedChainComplex(VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;

    std::map<Cell, Boundary> chains = chainComplex(complex, field);

    std::vector<Cell> sources = criticalCells(candidate);
    std::stable_sort(sources.begin(), sources.end(),
                     Comparator(complex, candidate.scalars));

    size_t const n = sources.size();

    std::map<Cell, size_t> index;
    for (size_t i = 0; i < sources.size(); ++i)
        index[sources.at(i)] = i;

    std::vector<unsigned int> dims;
    std::vector<float> scalars;
    std::vector<std::vector<Cell> > faceLists;

    for (size_t i = 0; i < n; ++i)
    {
        Cell const v = sources.at(i);
        dims.push_back(complex.cellDimension(v));
        scalars.push_back(cellValue(complex, candidate.scalars, v));

        std::vector<Cell> const flin = convertedBoundary(chains, v);
        std::vector<Cell> flout;
        for (size_t j = 0; j < flin.size(); ++j)
            flout.push_back(index.at(flin.at(j)));

        faceLists.push_back(flout);
    }

    return SimpleComplex(dims, scalars, faceLists);
}


SimpleComplex convertedCubicalComplex(VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Facets const facets(complex.xdim(), complex.ydim(), complex.zdim(), false);

    std::vector<Cell> cells = allCells(candidate);
    std::stable_sort(cells.begin(), cells.end(),
                     Comparator(complex, candidate.scalars));

    size_t const n = cells.size();

    std::map<Cell, size_t> index;
    for (size_t i = 0; i < cells.size(); ++i)
        index[cells.at(i)] = i;

    std::vector<unsigned int> dims;
    std::vector<float> scalars;
    std::vector<std::vector<Cell> > faceLists;

    for (size_t i = 0; i < n; ++i)
    {
        Cell const v = cells.at(i);
        dims.push_back(complex.cellDimension(v));
        scalars.push_back(cellValue(complex, candidate.scalars, v));

        int const m = facets.count(v);
        std::vector<Cell> flout;
        for (int j = 0; j < m; ++j)
            flout.push_back(index.at(facets(v, j)));

        faceLists.push_back(flout);
    }

    return SimpleComplex(dims, scalars, faceLists);
}


std::ostream& operator<<(std::ostream& out, SimpleComplex const& c)
{
    out << "[";

    for (size_t i = 0; i < c.nrCells(); ++i)
    {
        if (i > 0)
            out << " ";
        out << i << "," << c.cellDimension(i) << "," << c.cellValue(i);

        out << "(";

        std::vector<Cell> fs = c.cellFaces(i);
        for (size_t j = 0; j < fs.size(); ++j)
        {
            if (j > 0)
                out << " ";
            out << fs.at(j);
        }

        out << ")";
    }

    out << "]";

    return out;
}


std::vector<float> makeVector(double const x, double const y, double const z)
{
    std::vector<float> result(3);
    result.at(0) = x;
    result.at(1) = y;
    result.at(2) = z;
    return result;
}

std::vector<float> sum(
    std::vector<float> const& a,
    std::vector<float> const& b)
{
    std::vector<float> result(std::max(a.size(), b.size()));

    for (size_t i = 0; i < result.size(); ++i)
    {
        if (i >= a.size())
            result.at(i) = b.at(i);
        else if (i >= b.size())
            result.at(i) = a.at(i);
        else
            result.at(i) = a.at(i) + b.at(i);
    }
    return result;
}

std::vector<float> directionDelta(Field::Directions const& d)
{
    assert((d & 7) > 0);

    switch (d & 7)
    {
    case Field::SELF:
        return makeVector( 0.0,  0.0,  0.0);
    case Field::XUP:
        return makeVector( 0.5,  0.0,  0.0);
    case Field::XDOWN:
        return makeVector(-0.5,  0.0,  0.0);
    case Field::YUP:
        return makeVector( 0.0,  0.5,  0.0);
    case Field::YDOWN:
        return makeVector( 0.0, -0.5,  0.0);
    case Field::ZUP:
        return makeVector( 0.0,  0.0,  0.5);
    case Field::ZDOWN:
        return makeVector( 0.0,  0.0, -0.5);
    default:
        throw "unreachable";
    }
}


bool isOutwardPointing(Cell const& cell, VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;

    std::vector<float> pos = sum(complex.cellPosition(cell),
                                 directionDelta(field.getDirection(cell)));

    return (pos.at(0) < 0 or pos.at(0) >= complex.xdim() or
            pos.at(1) < 0 or pos.at(1) >= complex.ydim() or
            pos.at(2) < 0 or pos.at(2) >= complex.zdim());
}


// === shared predicates

anu_am::generative::Result
vectorDirectionIsDefined(Cell const& cell, VolumeData const& candidate)
{
    if (candidate.field.getDirection(cell) == Field::UNDEFINED)
        return anu_am::generative::failure("direction not defined");
    else
        return anu_am::generative::success();
}


anu_am::generative::Result
vectorIsNotOutwardPointing(Cell const& cell, VolumeData const& candidate)
{
    if (candidate.field.pointsOutward(cell))
        return anu_am::generative::failure("vector marked as pointing outward");
    else
        return anu_am::generative::success();
}


anu_am::generative::Result
directionIsNotOutwardPointing(Cell const& cell, VolumeData const& candidate)
{
    if (isOutwardPointing(cell, candidate))
    {
        std::stringstream msg;
        msg << " -> " << candidate.complex.cellPosition(cell);
        return anu_am::generative::failure(msg.str());
    }
    else
        return anu_am::generative::success();
}


anu_am::generative::Result cellPartnerMatchesDirection(Cell const& cell, VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;

    std::vector<float> p0 = complex.cellPosition(cell);
    std::vector<float> d  = directionDelta(field.getDirection(cell));
    std::vector<float> p1 = complex.cellPosition(field.getPartner(cell));
    std::vector<float> p2 = sum(p0, d);

    if (p1.at(0) != p2.at(0) or
        p1.at(1) != p2.at(1) or
        p1.at(2) != p2.at(2))
    {
        std::stringstream msg;
        msg << "partner mismatch " << p1 << " <-> " << p0 << " + " << d;
        return anu_am::generative::failure(msg.str());
    }
    else
        return anu_am::generative::success();
}


anu_am::generative::Result partnerOfPartnerIsOriginalCell(Cell const& a, VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;
    Cell const b = field.getPartner(a);
    Cell const c = field.getPartner(b);

    if (c != a)
    {
        std::stringstream msg;
        msg << " -> " << complex.cellPosition(b)
            << " -> " << complex.cellPosition(c);
        return anu_am::generative::failure(msg.str());
    }
    return anu_am::generative::success();
}


anu_am::generative::Result containsNoCyclicVPaths(VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Facets const facets(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Field const& field = candidate.field;

    std::stack<std::pair<Cell, int> > path;
    std::set<Cell> seen;
    std::set<Cell> done;

    for (Cell start = 0; start < complex.cellIdLimit(); ++start)
    {
        if (not complex.isCell(start) or not field.isCritical(start))
            continue;

        assert(seen.count(start) == 0);
        path.push(std::make_pair(start, 0));
        seen.insert(start);

        while (not path.empty())
        {
            Cell const cell = path.top().first;
            int const k = path.top().second;
            path.pop();

            int const n = facets.count(cell);

            if (k >= n)
            {
                done.insert(cell);
                continue;
            }

            path.push(std::make_pair(cell, k + 1));

            Cell const f = facets(cell, k);
            Cell const next = field.getPartner(f);

            if (next == f or next == cell or done.count(next) > 0)
                continue;
            else if (seen.count(next) == 0)
            {
                path.push(std::make_pair(next, 0));
                seen.insert(next);
            }
            else
            {
                std::stringstream msg;
                msg << "cyclic V-path at " << complex.cellPosition(next);
                return anu_am::generative::failure(msg.str());
            }
        }
    }

    return anu_am::generative::success();
}


anu_am::generative::Result vImageHasCompatibleValue(
    Value const& threshold,
    Cell const& cell,
    VolumeData const& candidate)
{
    CubicalComplex const& complex = candidate.complex;
    Scalars const& scalars = candidate.scalars;
    Field const& field = candidate.field;

    if (not field.isCritical(cell))
    {
        Cell const& partner = field.getPartner(cell);
        if (complex.cellDimension(partner) > complex.cellDimension(cell))
        {
            if (cellValue(complex, scalars, partner) >
                cellValue(complex, scalars, cell) + threshold)
            {
                std::stringstream msg;
                msg << complex.cellPosition(cell)
                    << " (" << cellValue(complex, scalars, cell) << ")"
                    << " -> " << complex.cellPosition(partner)
                    << " (" << cellValue(complex, scalars, partner) << ")";
                return anu_am::generative::failure(msg.str());
            }
        }
    }

    return anu_am::generative::success();
}


std::vector<std::vector<std::pair<float, int> > >
normalizedBettiNumbers(
    SimpleComplex const& c,
    size_t const dim,
    float const threshold = 0)
{
    std::vector<std::map<float, int> > const b =
        bettiNumbers(persistencePairing(c), threshold);
    std::vector<std::vector<std::pair<float, int> > > result;

    for (size_t d = 0; d < b.size(); ++d)
    {
        std::map<float, int> const& input = b.at(d);
        std::vector<std::pair<float, int> > output;

        std::map<float, int>::const_iterator iter;
        for (iter = input.begin(); iter != input.end(); ++iter)
            output.push_back(std::make_pair(iter->first, iter->second));

        std::stable_sort(output.begin(), output.end());

        if (output.at(0).second == 0)
            output.erase(output.begin());

        result.push_back(output);
    }

    while (result.size() <= dim)
        result.push_back(std::vector<std::pair<float, int> >());

    return result;
}


anu_am::generative::Result checkPersistentHomology(
    SimpleComplex const& c1,
    SimpleComplex const& c2,
    float const threshold = 0)
{
    typedef std::vector<std::vector<std::pair<float, int> > > Betti;

    int const d = c1.dimension();

    Betti const b1 = normalizedBettiNumbers(c1, d, threshold);
    Betti const b2 = normalizedBettiNumbers(c2, d, threshold);

    if (b1 != b2)
    {
        std::stringstream msg;
        msg << "Mismatched Betti numbers " << b1 << " vs " << b2
            << " for complexes " << c1 << " and " << c2;
        return anu_am::generative::failure(msg.str());
    }
    else
        return anu_am::generative::success();
}



} // namespace diamorse
} // namespace anu_am

#endif //!ANU_AM_DIAMORSE_TEST_HPP
