/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  Partition.hpp
 *
 *  Partitions of integer ranges of the form [0..n).
 *
 *  Olaf Delgado-Friedrichs mar 14
 *
 */

#ifndef ANU_AM_DIAMORSE_PARTITION_HPP
#define ANU_AM_DIAMORSE_PARTITION_HPP

#include <vector>


namespace anu_am
{
namespace diamorse
{


template<typename T>
class Partition
{
    mutable std::vector<T>      parent_;
    mutable std::vector<size_t> rank_;

    T hasParent(T const element) const
    {
        return parent_.at(element) > 0;
    }

    T getParent(T const element) const
    {
        return parent_.at(element) - 1;
    }

    void setParent(T const element, T const parent)
    {
        parent_.at(element) = parent + 1;
    }

    void shortcutParent(T const element, T const root) const
    {
        parent_.at(element) = root + 1;
    }

    T getRank(T const element) const
    {
        return rank_.at(element);
    }

    void setRank(T const element, size_t const value)
    {
        rank_.at(element) = value;
    }

    void connect(T const a, T const b)
    {
        setParent(a, b);
        setRank(b, getRank(a) + getRank(b) + 1);
        setRank(a, 0);
    }

public:
    Partition(T const size = 0)
        : parent_(size),
          rank_(size)
    {
    }

    T size() const
    {
        return parent_.size();
    }

    T find(T const element) const
    {
        T root = element;
        while (hasParent(root))
            root = getParent(root);

        T x = element;
        while (hasParent(x))
        {
            T y = getParent(x);
            shortcutParent(x, root);
            x = y;
        }

        return root;
    }

    void unite(T const a, T const b)
    {
        T const root_a = find(a);
        T const root_b = find(b);

        if (root_a == root_b)
            return;

        if (getRank(root_b) > getRank(root_a))
            connect(root_a, root_b);
        else
            connect(root_b, root_a);
    }
};



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_PARTITION_HPP
