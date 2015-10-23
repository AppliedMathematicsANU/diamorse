/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  OrderedMap.hpp
 *
 *  A map that preserves insertion order.
 *
 *  Olaf Delgado-Friedrichs may 15
 *
 */


#ifndef ANU_AM_DIAMORSE_ORDEREDMAP_HPP
#define ANU_AM_DIAMORSE_ORDEREDMAP_HPP

#include <map>
#include <vector>

namespace anu_am
{
namespace diamorse
{



template<typename K, typename V>
class OrderedMap
{
    std::map<K,V> _map;
    std::vector<K> _keysInOrder;

public:
    OrderedMap()
    {
    }

    OrderedMap(K const key, V const value)
    {
        this->set(key, value);
    }

    V operator()(K const key) const
    {
        return _map.at(key);
    }

    bool hasKey(K const key) const
    {
        return _map.count(key) > 0;
    }

    size_t size() const
    {
        return _keysInOrder.size();
    }

    K keyAt(size_t const i) const
    {
        return _keysInOrder.at(i);
    }

    V at(size_t const i) const
    {
        return _map.at(_keysInOrder.at(i));
    }

    void set(K const key, V const value)
    {
        if (!hasKey(key))
            _keysInOrder.push_back(key);
        _map[key] = value;
    }

    OrderedMap<K, V>& operator()(K const key, V const value)
    {
        set(key, value);
        return *this;
    }
};



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_ORDEREDMAP_HPP
