/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  traversals.hpp
 *
 *  Stable/unstable set traversal and derived algorithms.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#ifndef ANU_AM_DIAMORSE_TRAVERSALS_HPP
#define ANU_AM_DIAMORSE_TRAVERSALS_HPP

#include <cstddef>
#include <map>
#include <queue>
#include <set>
#include <vector>

namespace anu_am
{
namespace diamorse
{


template<class Cell, class VectorField, class Incidences>
std::vector<std::pair<Cell, Cell> > flowTraversal(
    Cell        const& s,
    VectorField const& V,
    Incidences  const& I,
    bool        const  coordinated = false)
{
    std::map<Cell, int> k;
    std::queue<Cell>    queue;
    std::vector<std::pair<Cell, Cell> > result;

    if (coordinated)
    {
        std::vector<std::pair<Cell, Cell> > t = flowTraversal(s, V, I);
        for (size_t i = 0; i < t.size(); ++i)
            ++k[V(t.at(i).second)];
    }

    queue.push(s);
    --k[s];

    while (not queue.empty())
    {
        Cell const a = queue.front();
        queue.pop();

        size_t const n = I.count(a);
        for (size_t i = 0; i < n; ++i)
        {
            Cell const b = I(a, i);

            if (V.defined(b) and V(b) != a)
            {
                result.push_back(std::make_pair(a, b));
                
                Cell const c = V(b);
                if ((k.count(c) == 0 or k[c] > 0) and c != b)
                {
                    if (coordinated and --k[c] != 0)
                        continue;
                    queue.push(c);
                    --k[c];
                }
            }
        }
    }

    return result;
}

template<class Cell, class VectorField, class Incidences>
std::set<Cell> unstableSet(
    Cell const& s,
    VectorField const& V,
    Incidences const& I)
{
    std::set<Cell> seen;
    seen.insert(s);

    std::vector<std::pair<Cell, Cell> > t = flowTraversal(s, V, I);
    for (size_t i = 0; i < t.size(); ++i)
    {
        Cell const b = t.at(i).second;
        Cell const c = V(b);
        if (b != c and seen.count(c) == 0)
            seen.insert(c);
    }

    return seen;
}

template<class Cell, class VectorField, class Incidences>
std::vector<std::pair<Cell, int> > morseBoundary(
    Cell const& s,
    VectorField const& V,
    Incidences const& I)
{
    std::map<Cell, int> counts;
    std::set<Cell> boundary;

    counts[s] = 1;

    std::vector<std::pair<Cell, Cell> > t = flowTraversal(s, V, I, true);
    for (size_t i = 0; i < t.size(); ++i)
    {
        Cell const a = t.at(i).first;
        Cell const b = t.at(i).second;
        Cell const c = V(b);

        int const n = counts[a] + (counts.count(c) > 0 ? counts[c] : 0);
        counts[c] = (n <= 3) ? n : n % 2 + 2;
        if (b == c)
            boundary.insert(c);
    }

    std::vector<std::pair<Cell, int> > result;
    typename std::set<Cell>::const_iterator iter;
    for (iter = boundary.begin(); iter != boundary.end(); ++iter)
        result.push_back(std::make_pair(*iter, counts[*iter]));

    return result;
}

template<class Cell, class VectorField, class Incidences>
std::vector<std::pair<Cell, int> > morseBoundaryFast(
    Cell const& s,
    VectorField const& V,
    Incidences const& I)
{
    std::queue<Cell> queue;
    std::map<Cell, int> counts;

    queue.push(s);

    while (not queue.empty())
    {
        Cell const a = queue.front();
        queue.pop();

        size_t const n = I.count(a);
        for (size_t i = 0; i < n; ++i)
        {
            Cell const b = I(a, i);

            if (V.defined(b) and V(b) != a )
            {
                Cell const c = V(b);

                if (c != b)
                {
                    queue.push(c);
                }
                else
                {
                    int const n = ++counts[b];
                    if (n > 3)
                        counts[b] = n % 2 + 2;
                }
            }
        }
    }

    return std::vector<std::pair<Cell, int> >(counts.begin(), counts.end());
}

template<class Cell, class VectorField, class Incidences>
std::vector<std::pair<Cell, Cell> > connections(
    Cell const& s,
    Cell const& t,
    VectorField const& V,
    VectorField const& coV,
    Incidences const& I,
    Incidences const& coI)
{
    std::set<Cell> active;
    active.insert(t);

    std::vector<std::pair<Cell, Cell> > back = flowTraversal(t, coV, coI);
    for (size_t i = 0; i < back.size(); ++i)
        active.insert(coV(back.at(i).second));

    std::vector<std::pair<Cell, Cell> > result;
    std::vector<std::pair<Cell, Cell> > forw = flowTraversal(s, V, I);
    for (size_t i = 0; i < forw.size(); ++i)
    {
        Cell const a = forw.at(i).first;
        Cell const b = forw.at(i).second;
        if (active.count(b) > 0)
            result.push_back(std::make_pair(a, b));
    }

    return result;
}

template<class Cell, class Boundary, class Key>
std::pair<Cell, int> closePartner(
    Cell const s,
    Boundary const& B,
    Boundary const& coB,
    Key const& K)
{
    typename Boundary::result_type bnd = B(s);

    if (bnd.size() == 0)
        return std::make_pair(s, 0);

    std::pair<Cell, int> best = *bnd.begin();

    typename std::vector<std::pair<Cell, int> >::const_iterator iter;
    for (iter = bnd.begin(); iter != bnd.end(); ++iter)
        if (K(iter->first) > K(best.first))
            best = *iter;

    typename Boundary::result_type cob = coB(best.first);
    for (iter = cob.begin(); iter != cob.end(); ++iter)
        if (K(iter->first) < K(s))
            return std::make_pair(s, 0);

    return best;
}

template<class Cell, class Boundary, class Key, class Predicate>
std::vector<std::pair<Cell, Cell> > pairsToCancel(
    std::set<Cell> const& S,
    Boundary const& B,
    Boundary const& coB,
    Key const& K,
    Predicate const& P)
{
    std::vector<std::pair<Cell, Cell> > result;

    typename std::set<Cell>::const_iterator outer;
    for (outer = S.begin(); outer != S.end(); ++outer)
    {
        Cell const s = *outer;
        std::pair<Cell, int> res = closePartner(s, B, coB, K);
        if (res.second == 1) {
            Cell const t = res.first;
            if (P(s, t))
            {
                result.push_back(std::make_pair(s, t));
            }
        }
    }

    return result;
}



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_TRAVERSALS_HPP
