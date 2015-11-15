/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  simplification.hpp
 *
 *  Simplifies a gradient vector field via repeated Morse cancellation.
 *
 *  Olaf Delgado-Friedrichs mar 14
 *
 */

#ifndef ANU_AM_DIAMORSE_SIMPLIFICATION_HPP
#define ANU_AM_DIAMORSE_SIMPLIFICATION_HPP

#include <iostream>
#include <limits>
#include <map>
#include <set>

#include "callables.hpp"
#include "chainComplexExtraction.hpp"
#include "CubicalComplex.hpp"


namespace anu_am
{
namespace diamorse
{


template<class Scalars, class Field>
class MayCancel
{
    typedef typename CubicalComplex::cell_id_type Cell;

    CubicalComplex const& complex_;
    Scalars const& scalars_;
    typename Field::Vectors const V_;
    typename Field::Vectors const coV_;
    Facets const I_;
    Facets const coI_;
    Vertices const vertices_;
    float const persistenceLimit_;
    float const sizeLimit_;
    float const threshold_;

    float value(Cell const u) const
    {
        size_t const n = vertices_.count(u);
        float val = scalars_(vertices_(u, 0));

        for (size_t i = 1; i < n; ++i)
            val = std::max(val, scalars_(vertices_(u, i)));

        return val;
    }

public:
    MayCancel(CubicalComplex const& complex,
              Scalars const& scalars,
              Field const& field,
              float const persistenceLimit = 1,
              float const sizeLimit = 0,
              float const threshold = -std::numeric_limits<float>::max())
        : complex_(complex),
          scalars_(scalars),
          V_(field.V()),
          coV_(field.coV()),
          I_(complex.xdim(), complex.ydim(), complex.zdim(), false),
          coI_(complex.xdim(), complex.ydim(), complex.zdim(), true),
          vertices_(complex.xdim(), complex.ydim(), complex.zdim()),
          persistenceLimit_(persistenceLimit),
          sizeLimit_(sizeLimit),
          threshold_(threshold)
    {
    }

    bool operator()(Cell const u, Cell const v) const
    {
        float const fu = value(u);
        float const fv = value(v);

        if (fu > threshold_ and fv < threshold_)
            return false;

        if (persistenceLimit_ >= 0 and fu - fv > persistenceLimit_)
            return false;

        if (sizeLimit_ > 0)
        {
            float size;
            if (complex_.cellDimension(v) == 0)
                size = unstableSet(v, coV_, coI_).size();
            else
                size = unstableSet(u, V_, I_).size();

            if (size > sizeLimit_)
                return false;
        }

        return true;
    }
};


template<class Scalars, class Field>
MayCancel<Scalars, Field> mayCancel(
    CubicalComplex const& complex,
    Scalars const& scalars,
    Field const& field,
    float const persistenceLimit = 1,
    float const sizeLimit = 0,
    float const threshold = -std::numeric_limits<float>::max())
{
    return MayCancel<Scalars, Field>(
        complex,
        scalars,
        field,
        persistenceLimit,
        sizeLimit,
        threshold);
}


template<typename Cell, class Boundary>
std::map<Cell, Boundary>
reverse(std::map<Cell, Boundary> const& chains)
{
    std::map<Cell, Boundary> result;

    typename std::map<Cell, Boundary>::const_iterator outer;
    for (outer = chains.begin(); outer != chains.end(); ++outer)
    {
        Cell const s = outer->first;
        Boundary const& bnd = outer->second;

        typename Boundary::const_iterator inner;
        for (inner = bnd.begin(); inner != bnd.end(); ++inner)
        {
            Cell const t = inner->first;
            int const n = inner->second;

            if (result.count(t) == 0)
                result[t] = Boundary();
            result[t].push_back(std::make_pair(s, n));
        }
    }

    return result;
}


template<typename S, typename T>
std::vector<S> firsts(std::vector<std::pair<S, T> > const& pairs)
{
    std::vector<S> result;
    for (size_t i = 0; i < pairs.size(); ++i)
        result.push_back(pairs.at(i).first);
    return result;
}


template<typename Cell, class Boundary>
std::set<Cell> dirtySet(
    std::vector<std::pair<Cell, Cell> > const& pairs,
    Boundary const& chains,
    Boundary const& cochains)
{
    std::set<Cell> dirty;

    for (size_t i = 0; i < pairs.size(); ++i)
    {
        Cell const source = pairs.at(i).first;

        if (cochains.count(source) > 0)
        {
            std::vector<Cell> const cob = firsts(cochains.at(source));
            dirty.insert(cob.begin(), cob.end());
        }

        std::vector<Cell> const bnd = firsts(chains.at(source));
        for (size_t j = 0; j < bnd.size(); ++j)
        {
            Cell const target = bnd.at(j);
            if (cochains.count(target) > 0)
            {
                std::vector<Cell> const cob = firsts(cochains.at(target));
                dirty.insert(cob.begin(), cob.end());
            }
        }
    }

    return dirty;
}


template<class Field, class Key, class Predicate>
void simplify(
    CubicalComplex const& complex,
    Field& field,
    Key const& sortingKey,
    Predicate const& allow,
    std::ostream& log = std::cerr)
{
    typedef CubicalComplex::cell_id_type Cell;
    typedef RestrictedIncidences<Facets> Incidences;
    typedef std::vector<std::pair<Cell, int> > Boundary;
    typedef std::vector<std::pair<Cell, Cell> > Pairs;

    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
    typename Field::Vectors V = field.V();
    typename Field::Vectors coV = field.coV();

    size_t const n = complex.cellIdLimit();
    boost::shared_ptr<std::vector<bool> > marked =
        connectingPaths(complex, field, V, coV, I, coI);
    Incidences rI(I, marked);
    Incidences rcoI(coI, marked);

    std::set<Cell> candidates;
    std::map<Cell, Boundary> chains;
    size_t oldsize;

    for (Cell cell = 0; cell <= n; ++cell)
    {
        if (complex.isCell(cell) and field.isCritical(cell))
        {
            candidates.insert(cell);
            chains[cell] = morseBoundary(cell, V, rI);
        }
    }
    oldsize = candidates.size();

    log << "starting with " << candidates.size() << " critical cells"
        << std::endl << std::endl;

    while (true)
    {
        std::map<Cell, Boundary> cochains = reverse(chains);
        Pairs pairs = pairsToCancel(candidates,
                                    callableMap(chains),
                                    callableMap(cochains),
                                    sortingKey,
                                    allow);
        log << pairs.size() << " cancellations" << std::endl;

        if (pairs.size() == 0)
            break;

        int count = 0;
        typename Pairs::const_iterator outer;
        for (outer = pairs.begin(); outer != pairs.end(); ++outer)
        {
            Pairs path =
                connections(outer->first, outer->second, V, coV, rI, rcoI);
            count += path.size();

            typename Pairs::const_iterator inner;
            for (inner = path.begin(); inner != path.end(); ++inner)
                field.setPartner(inner->second, inner->first);
        }
        log << "  " << count << " modifications" << std::endl;

        size_t const newsize = candidates.size() - 2 * pairs.size();
        if (newsize < 0.3 * oldsize)
        {
            marked = connectingPaths(complex, field, V, coV, rI, rcoI);
            rI = Incidences(I, marked);
            rcoI = Incidences(coI, marked);
            oldsize = newsize;
        }

        std::set<Cell> dirty = dirtySet(pairs, chains, cochains);

        std::set<Cell>::const_iterator iter;
        for (iter = dirty.begin(); iter != dirty.end(); ++iter)
            chains[*iter] = morseBoundary(*iter, V, rI);

        for (outer = pairs.begin(); outer != pairs.end(); ++outer)
        {
            candidates.erase(outer->first);
            candidates.erase(outer->second);
            chains.erase(outer->first);
            chains.erase(outer->second);
        }

        log << "  " << candidates.size() << " critical cells left" << std::endl;
    }
}



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_SIMPLIFICATION_HPP
