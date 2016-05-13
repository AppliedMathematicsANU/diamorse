/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  MorseVectorField.hpp
 *
 *  Gradient vector fields on cubical complexes.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_MORSEVECTORFIELD_HPP
#define ANU_AM_DIAMORSE_MORSEVECTORFIELD_HPP

#include <vector>

#include <memory>

#include "CubicalComplex.hpp"

namespace anu_am
{
namespace diamorse
{


/// This class represents a gradient vector fields on cubical complexes.
template<class S>
class MorseVectorField
{
    typedef S Storage;
    typedef CubicalComplex Complex;
    typedef Complex::cell_id_type Cell;

public:
    typedef typename Storage::DataPtr DataPtr;

    typedef enum {
        UNDEFINED = 0,
        SELF      = 1,
        XUP       = 2,
        XDOWN     = 3,
        YUP       = 4,
        YDOWN     = 5,
        ZUP       = 6,
        ZDOWN     = 7
    } Directions;

    class Vectors
    {
        Complex complex_;
        MorseVectorField field_;
        bool dual_;

    public:
        Vectors(Complex const& complex,
                MorseVectorField const& field,
                bool const dual)
            : complex_(complex),
              field_(field),
              dual_(dual)
        {
        }

        bool defined(Cell const& cell) const
        {
            if (field_.isCritical(cell))
                return true;
            else if (field_.pointsOutward(cell))
                return false;
            else if (dual_)
                return not field_.forward(cell, field_.getDirection(cell));
            else
                return field_.forward(cell, field_.getDirection(cell));
        }

        Cell operator()(Cell const& cell) const
        {
            assert(defined(cell));
            return field_.getPartner(cell);
        }
    };

    MorseVectorField()
    {
    }

    MorseVectorField(Complex const& complex)
        : xdim_(complex.xdim()),
          ydim_(complex.ydim()),
          zdim_(complex.zdim()),
          complex_(xdim_, ydim_, zdim_),
          storage_(xdim_ * ydim_ * zdim_ * 8, 0)
    {
    }

    MorseVectorField(size_t xdim, size_t ydim, size_t zdim, DataPtr data)
        : xdim_(xdim),
          ydim_(ydim),
          zdim_(zdim),
          complex_(xdim_, ydim_, zdim_),
          storage_(data, 0)
    {
    }

    Vectors const V() const
    {
        return Vectors(complex_, *this, false);
    }

    Vectors const coV() const
    {
        return Vectors(complex_, *this, true);
    }

    Directions getDirection(Cell const v) const
    {
        return (Directions) (storage_.get(v) & 7);
    }

    void setDirection(Cell const v, Directions const d)
    {
        int const out = outward(v, d) ? 8 : 0;
        storage_.set(v, (typename Storage::value_type) (d | out));
    }

    bool isCritical(Cell const n) const
    {
        return getDirection(n) == SELF;
    }

    Cell getPartner(Cell const n) const
    {
        assert(not pointsOutward(n));
        return neighbor(n, getDirection(n));
    }

    void setPartner(Cell const v, Cell const w)
    {
        for (int i = 1; i <= 7; ++i)
        {
            if (neighbor(v, i) == w)
            {
                setDirection(v, (Directions) i);
                setDirection(w, (Directions) (i == 1 ? i : i ^ 1));
                assert(getPartner(v) == w);
                assert(getPartner(w) == v);
                return;
            }
        }
        assert(false);
    }

    bool pointsOutward(Cell const v) const
    {
        return storage_.get(v) > 7;
    }

    DataPtr const data() const
    {
        return storage_.data();
    }

private:
    size_t xdim_;
    size_t ydim_;
    size_t zdim_;
    CubicalComplex complex_;
    Storage storage_;
    
    Cell neighbor(Cell const n, int const direction) const
    {
        switch (direction)
        {
        case XUP:
            return (complex_.cellDX(n) ? n + 8 : n) ^ 1;
        case XDOWN:
            return (complex_.cellDX(n) ? n : n - 8) ^ 1;
        case YUP:
            return (complex_.cellDY(n) ? n + 8 * xdim_ : n) ^ 2;
        case YDOWN:
            return (complex_.cellDY(n) ? n : n - 8 * xdim_) ^ 2;
        case ZUP:
            return (complex_.cellDZ(n) ? n + 8 * xdim_ * ydim_ : n) ^ 4;
        case ZDOWN:
            return (complex_.cellDZ(n) ? n : n - 8 * xdim_ * ydim_) ^ 4;
        default:
            return n;
        }
    }

    bool outward(Cell const v, Directions const d) const
    {
        switch (d)
        {
        case XUP:
            return complex_.cellX(v) == xdim_ - 1 and complex_.cellDX(v) == 0;
        case XDOWN:
            return complex_.cellX(v) == 0 and complex_.cellDX(v) == 0;
        case YUP:
            return complex_.cellY(v) == ydim_ - 1 and complex_.cellDY(v) == 0;
        case YDOWN:
            return complex_.cellY(v) == 0 and complex_.cellDY(v) == 0;
        case ZUP:
            return complex_.cellZ(v) == zdim_ - 1 and complex_.cellDZ(v) == 0;
        case ZDOWN:
            return complex_.cellZ(v) == 0 and complex_.cellDZ(v) == 0;
        default:
            return false;
        }
    }

    bool forward(Cell const v, Directions const d) const
    {
        switch (d)
        {
        case XUP:
            return not complex_.cellDX(v);
        case XDOWN:
            return not complex_.cellDX(v);
        case YUP:
            return not complex_.cellDY(v);
        case YDOWN:
            return not complex_.cellDY(v);
        case ZUP:
            return not complex_.cellDZ(v);
        case ZDOWN:
            return not complex_.cellDZ(v);
        default:
            return false;
        }
    }
};


} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_MORSEVECTORFIELD_HPP
