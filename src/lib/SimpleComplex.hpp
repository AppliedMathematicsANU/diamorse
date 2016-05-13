/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  SimpleComplex.hpp
 *
 *  A trivial chain complex implementation with scalar per-cell values.
 *
 *  Olaf Delgado-Friedrichs aug 13
 *
 */

#ifndef ANU_AM_DIAMORSE_SIMPLEXCOMPLEX_HPP
#define ANU_AM_DIAMORSE_SIMPLEXCOMPLEX_HPP

#include <vector>

#include <memory>

namespace anu_am
{
namespace diamorse
{


class SimpleComplex
{
public:
    typedef size_t cell_id_type;

private:
    typedef cell_id_type Cell;

    std::shared_ptr<std::vector<unsigned int> >       dims_;
    std::shared_ptr<std::vector<float> >              scalars_;
    std::shared_ptr<std::vector<std::vector<Cell> > > faceLists_;

public:
    SimpleComplex() {}

    SimpleComplex(std::vector<unsigned int> const& dims,
                  std::vector<float> const& scalars,
                  std::vector<std::vector<Cell> > const& faceLists)
        : dims_(new std::vector<unsigned int>(dims)),
          scalars_(new std::vector<float>(scalars)),
          faceLists_(new std::vector<std::vector<Cell> >(faceLists))
    {
    }

    /// The dimension of the complex is always 3.
    int dimension() const
    {
        return 3;
    }

    /// The total number of cells in the complex.
    size_t nrCells() const
    {
        return dims_->size();
    }

    /// The dimension of a given cell.
    int cellDimension(Cell const id) const
    {
        return dims_->at(id);
    }

    /// The scalar value for a given cell.
    float cellValue(Cell const id) const
    {
        return scalars_->at(id);
    }

    /// The list of (highest-dimensional proper) faces for a given cell.
    std::vector<Cell> cellFaces(Cell const id) const
    {
        return faceLists_->at(id);
    }
};



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_SIMPLEXCOMPLEX_HPP
