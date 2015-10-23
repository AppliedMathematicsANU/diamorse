/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  collections.hpp
 *
 *  Some utility code to make working with C++ collections a little easier.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_COLLECTIONS_HPP
#define ANU_AM_DIAMORSE_COLLECTIONS_HPP

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace anu_am
{
namespace diamorse
{


template<typename Iter>
std::ostream& print(std::ostream& out, Iter begin, Iter const end);
    
template<typename S, typename T>
std::ostream& operator<<(std::ostream& out, std::pair<S, T> const& p)
{
    out << "(" << p.first << ", " << p.second << ")";
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& v)
{
    return print(out, v.begin(), v.end());
}

template<typename T>
std::ostream& operator<<(std::ostream& out, std::set<T> const& v)
{
    return print(out, v.begin(), v.end());
}

template<typename K, typename V>
std::ostream& operator<<(std::ostream& out, std::map<K, V> const& v)
{
    return print(out, v.begin(), v.end());
}

template<typename Iter>
std::ostream& print(std::ostream& out, Iter begin, Iter const end)
{
    Iter p = begin;
    if (p != end)
    {
        out << *p++;
        while (p != end)
            out << " " << *p++;
    }
    return out;
}

    

} // namespace diamorse
} // namespace anu_am

#endif //!ANU_AM_DIAMORSE_COLLECTIONS_HPP
