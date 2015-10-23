/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  Set.hpp
 *
 *  Sets represented as sorted vectors.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#ifndef ANU_AM_DIAMORSE_SET_HPP
#define ANU_AM_DIAMORSE_SET_HPP

#include <set>
#include <vector>

namespace anu_am
{
namespace diamorse
{


template<typename T>
class Set
{
    typedef typename std::vector<T>::const_iterator Iter;

    std::vector<T> elements_;

    explicit Set(std::vector<T> const& source, bool const)
        :elements_(source)
    {
    }

public:
    // Create an empty set.
    Set()
        : elements_()
    {
    }

    // Construct a set with a single entry.
    Set(T const& t)
        : elements_()
    {
        elements_.push_back(t);
    }

    // Construct a set from a vector.
    explicit Set(std::vector<T> const& source)
        :elements_()
    {
        std::vector<T> tmp(source);
        std::sort(tmp.begin(), tmp.end());

        for (Iter it = tmp.begin(); it != tmp.end(); ++it)
            if (elements_.size() == 0 or *it != elements_.back())
                elements_.push_back(*it);
    }

    // Construct a set from an STL set.
    explicit Set(std::set<T> const& source)
        :elements_()
    {
        typename std::set<T>::const_iterator it;
        for (it = source.begin(); it != source.end(); ++it)
            elements_.push_back(*it);
    }

    // The size of the set.
    int size() const {
        return elements_.size();
    }

    // Creates a new set as the union of this set and the one given.
    Set operator+(Set const& other) const
    {
        std::vector<T> merged;

        Iter const left_end = elements_.end();
        Iter const right_end = other.elements_.end();

        Iter left = elements_.begin();
        Iter right = other.elements_.begin();
        while (left != left_end and right != right_end)
        {
            if (*left < *right)
            {
                merged.push_back(*left);
                ++left;
            }
            else if (*left > *right)
            {
                merged.push_back(*right);
                ++right;
            }
            else
            {
                merged.push_back(*right);
                ++left;
                ++right;
            }
        }
        if (left != left_end)
            merged.insert(merged.end(), left, left_end);
        else if (right != right_end)
            merged.insert(merged.end(), right, right_end);

        return Set(merged, true);
    }

    bool operator<(Set const& other) const
    {
        Iter const left_end = elements_.end();
        Iter const right_end = other.elements_.end();

        Iter left = elements_.begin();
        Iter right = other.elements_.begin();

        while (left != left_end and right != right_end)
        {
            if (*left < *right)
                return true;
            else if (*left > *right)
                return false;

            ++left;
            ++right;
        }

        return right != right_end;
    }

    std::vector<T> elements() const
    {
        return elements_;
    }
};



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_SET_HPP
