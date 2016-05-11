/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  PackedMap.h
 *
 *  A packed data structure with 4 bits for each entry.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_PACKEDMAP_H
#define ANU_AM_DIAMORSE_PACKEDMAP_H

#include <stdint.h>
#include <vector>

#include <memory>

namespace anu_am
{
namespace diamorse
{


class PackedMap
{
    typedef uint8_t Value;
    typedef uint32_t DataItem;

public:
    typedef Value value_type;
    typedef std::shared_ptr<std::vector<DataItem> > DataPtr;

    PackedMap()
    {
    }

    PackedMap(size_t const size,
              Value const defaultValue = Value())
        : defaultValue_(defaultValue * 0x11111111),
          data_(new std::vector<DataItem>((size + 7)/ 8, defaultValue_))
    {
    }

    PackedMap(DataPtr data,
              Value const defaultValue = Value())
        : defaultValue_(defaultValue),
          data_(data)
    {
    }

    void clear()
    {
        for (size_t i = 0; i < data_->size(); ++i)
            data_->at(i) = defaultValue_;
    }

    Value get(size_t const v) const
    {
        return (data_->at(v >> 3) >> shift(v)) & 15;
    }

    void set(size_t const v, Value const x)
    {
        data_->at(v >> 3) ^= ((x & 15) ^ get(v)) << shift(v);
   }

    DataPtr const data() const
    {
        return data_;
    }

private:
    DataItem defaultValue_;
    DataPtr data_;

    size_t shift(size_t const v) const
    {
        return 4 * (v & 7);
    }
};


} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_PACKEDMAP_H
