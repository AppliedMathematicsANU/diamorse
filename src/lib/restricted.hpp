/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  restricted.hpp
 *
 *  Code based on restricted traversal algorithms.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_RESTRICTED_HPP
#define ANU_AM_DIAMORSE_RESTRICTED_HPP

#include <set>

#include <memory>


namespace anu_am
{
namespace diamorse
{


template<class Incidences>
class RestrictedIncidences
{
    Incidences I_;
    std::shared_ptr<std::vector<bool> > good_;

public:
    RestrictedIncidences(
        Incidences const& I,
        std::shared_ptr<std::vector<bool> > good)
        : I_(I),
          good_(good)
    {
    }

    RestrictedIncidences(
        Incidences const& I,
        std::vector<bool> good)
        : I_(I),
          good_(new std::vector<bool>(good))
    {
    }

    int count(size_t const id) const
    {
        int n = 0;
        for (int i = 0; i < I_.count(id); ++i)
            if (good_->at(I_(id, i)))
                ++n;
        return n;
    }

    size_t operator()(size_t const id, int const n) const
    {
        int m = -1;
        for (int i = 0; i < I_.count(id); ++i)
            if (good_->at(I_(id, i)) and ++m == n)
                return I_(id, i);
        assert(false);
    }
};


template<class VectorField, class Incidences>
void mark(
    size_t const s,
    std::shared_ptr<std::vector<bool> >& marked,
    bool const value,
    VectorField const& V,
    Incidences const& I)
{
    std::queue<size_t> queue;

    marked->at(s) = value;
    queue.push(s);

    while (not queue.empty())
    {
        size_t const a = queue.front();
        queue.pop();
        for (int i = 0; i < I.count(a); ++i)
        {
            size_t const b = I(a, i);
            if (V.defined(b) and V(b) != a)
            {
                size_t const c = V(b);
                if (c != b and marked->at(c) != value)
                {
                    queue.push(c);
                    marked->at(b) = marked->at(c) = value;
                }
            }
        }
    }
}


template<class VectorField, class Incidences>
std::shared_ptr<std::vector<bool> > connectingPaths(
    size_t const size,
    std::vector<size_t> const& downstreamSources,
    std::vector<size_t> const& upstreamSources,
    std::vector<size_t> const& downstreamTargets,
    std::vector<size_t> const& upstreamTargets,
    VectorField const& V,
    VectorField const& coV,
    Incidences const& I,
    Incidences const& coI)
{
    typedef std::vector<bool> Bits;

    std::shared_ptr<Bits> active(new Bits(size));
    std::shared_ptr<Bits> result(new Bits(size));

    for (size_t i = 0; i < downstreamSources.size(); ++i)
        mark(downstreamSources.at(i), active, true, V, I);
    for (size_t i = 0; i < upstreamSources.size(); ++i)
        mark(upstreamSources.at(i), active, true, coV, coI);

    RestrictedIncidences<Incidences> rI(I, active);
    RestrictedIncidences<Incidences> rcoI(coI, active);

    for (size_t i = 0; i < downstreamTargets.size(); ++i)
        mark(downstreamTargets.at(i), result, true, coV, rcoI);
    for (size_t i = 0; i < upstreamTargets.size(); ++i)
        mark(upstreamTargets.at(i), result, true, V, rI);

    for (size_t i = 0; i < downstreamSources.size(); ++i)
        result->at(downstreamSources.at(i)) = true;
    for (size_t i = 0; i < upstreamSources.size(); ++i)
        result->at(upstreamSources.at(i)) = true;

    return result;
}



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_RESTRICTED_HPP
