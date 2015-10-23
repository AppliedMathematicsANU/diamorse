/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  VertexMap.h
 *
 *  Stores data at the vertices of a cubical complex.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_VERTEXMAP_H
#define ANU_AM_DIAMORSE_VERTEXMAP_H

#include <vector>

#include "boost/smart_ptr.hpp"

namespace anu_am
{
namespace diamorse
{


template<class Complex, typename Value>
class VertexMap
{
    typedef typename Complex::cell_id_type Cell;

public:
    typedef Cell  argument_type;
    typedef Value value_type;
    typedef Value result_type;

    typedef boost::shared_ptr<std::vector<Value> > DataPtr;

    VertexMap()
    {
    }

    VertexMap(Complex const& complex, Value const defaultValue = Value())
        : complex_(complex),
          defaultValue_(defaultValue),
          data_(new std::vector<Value>(
                    (size_t) complex_.xdim() * complex_.ydim() * complex_.zdim(),
                    defaultValue))
    {
    }

    VertexMap(Complex const& complex,
              DataPtr data,
              Value const defaultValue = Value())
        : complex_(complex),
          defaultValue_(defaultValue),
          data_(data)
    {
    }

    void clear()
    {
        size_t const n =
            (size_t) complex_.xdim() * complex_.ydim() * complex_.zdim();
        for (size_t i = 0; i < n; ++i)
            data_->at(i) = defaultValue_;
    }

    Value get(Cell const v) const
    {
        assert((v & 7) == 0);
        return data_->at(v >> 3);
    }

    Value operator()(Cell const v) const
    {
        return get(v);
    }

    void set(Cell const v, Value const val)
    {
        assert((v & 7) == 0);
        data_->at(v >> 3) = val;
    }

    DataPtr const data() const
    {
        return data_;
    }

private:
    Complex complex_;
    Value defaultValue_;
    DataPtr data_;
};


} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_VERTEXMAP_H
