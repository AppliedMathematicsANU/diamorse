/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  CubicalComplex.h
 *
 *  Implicit representation of a cubical complex on a rectangular domain.
 *
 *  Olaf Delgado-Friedrichs nov 13
 *
 */

#ifndef ANU_AM_DIAMORSE_CUBICALCOMPLEX_H
#define ANU_AM_DIAMORSE_CUBICALCOMPLEX_H

#include <assert.h>
#include <cstddef>
#include <cmath>
#include <vector>


namespace anu_am
{
namespace diamorse
{


/// This class provides an implicit representation of a cubical cell complex
/// on a vertex set of the form [0,a) x [0,b) x [0,c] within N^3. In other
/// words, the class implements a mapping of subset of the range [0, 8abc)
/// onto the cells of the complex, with respect to which all further
/// topological properties and relations are computed on the fly.
///
/// A flat vector containing scalar values for each of the vertices (with x as
/// the fastest changing coordinate) can be provided on construction. The
/// value for a higher-dimensional cell is then determined as the maximum over
/// the values on its vertices.
class CubicalComplex
{
public:
    // Some typedefs to use internally and within client code.
    typedef size_t cell_id_type;
    typedef std::vector<float> cell_position_type;

private:
    typedef cell_id_type Cell;

    friend class Vertices;
    friend class Facets;

public:
    CubicalComplex()
    {
    }

    /// The constructor takes the number of vertices in x, y and z direction.
    CubicalComplex(int const xdim, int const ydim, int const zdim)
        : xdim_(xdim),
          ydim_(ydim),
          zdim_(zdim)
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
        return (2 * xdim_ - 1) * (2 * ydim_ - 1) * (2 * zdim_ - 1);
    }

    /// The highest legal cell id.
    Cell cellIdLimit() const
    {
        return 8 * xdim_ * ydim_ * zdim_;
    }

    /// Checks whether a given integer is a cell index.
    bool isCell(Cell const id) const;

    /// The dimension of a given cell.
    int cellDimension(Cell const id) const
    {
        static int dim[] = { 0, 1, 1, 2, 1, 2, 2, 3 };

        return dim[id % 8];
    }

    /// The x coordinate associated to the given cell.
    size_t cellX(Cell const id) const
    {
        return cellVertex(id) % xdim_;
    }

    /// The y coordinate associated to the given cell.
    size_t cellY(Cell const id) const
    {
        return (cellVertex(id) / xdim_) % ydim_; 
    }

    /// The z coordinate associated to the given cell.
    size_t cellZ(Cell const id) const
    {
        return cellVertex(id) / (xdim_ * ydim_);
    }

    int cellDX(Cell const id) const
    {
        return (id & 1) > 0;
    }

    int cellDY(Cell const id) const
    {
        return (id & 2) > 0;
    }

    int cellDZ(Cell const id) const
    {
        return (id & 4) > 0;
    }

    /// The position in 3d space associated to (the center of) a given cell.
    std::vector<float> cellPosition(Cell const id) const
    {
        std::vector<float> pos(3);
        pos.at(0) = cellX(id) + 0.5 * cellDX(id);
        pos.at(1) = cellY(id) + 0.5 * cellDY(id);
        pos.at(2) = cellZ(id) + 0.5 * cellDZ(id);
        return pos;
    }

    /// Retrieved the id for the cell based at the vertex at (x, y, z) with
    /// the opposite vertex indicated by the distance (dx, dy, dz).
    Cell cellAt(int const x, int const y, int const z,
               bool const dx, bool const dy, bool const dz) const
    {
        return (((z * ydim_ + y) * xdim_ + x) << 3) + dx + 2 * dy + 4 * dz;
    }

    /// Retrieved the id for the cell with its center at the given position.
    Cell cellAt(double const x, double const y, double const z) const
    {
        return cellAt((int) floor(x), (int) floor(y), (int) floor(z),
                      x > floor(x), y> floor(y), z > floor(z));
    }

    /// The largest valid x coordinate plus one.
    size_t xdim() const { return xdim_; }

    /// The largest valid y coordinate plus one.
    size_t ydim() const { return ydim_; }

    /// The largest valid z coordinate plus one.
    size_t zdim() const { return zdim_; }

private:
    // The "parent" vertex for a given cell.
    Cell cellVertex(Cell const id) const
    {
        return id >> 3;
    }

    int type(Cell const id) const
    {
        Cell   const v = id >> 3;
        size_t const x = v % xdim_;
        size_t const y = (v / xdim_) % ydim_;
        size_t const z = v / (xdim_ * ydim_);

        return ((id & 7)
                + ((x <= 0)    << 3)
                + ((x >= xdim_-1) << 4)
                + ((y <= 0)    << 5)
                + ((y >= ydim_-1) << 6)
                + ((z <= 0)    << 7)
                + ((z >= zdim_-1) << 8));
    }

    // Data members.
    size_t xdim_;
    size_t ydim_;
    size_t zdim_;
};

/// Checks whether a given integer is a cell index.
bool CubicalComplex::isCell(CubicalComplex::cell_id_type const id) const
{
    typedef CubicalComplex::cell_id_type Cell;

    if (id >= ((Cell) 8) * xdim_ * ydim_ * zdim_)
        return false;

    if ((id & 1) and cellX(id) >= xdim_ - 1)
        return false;

    if ((id & 2) and cellY(id) >= ydim_ - 1)
        return false;
    
    if ((id & 4) and cellZ(id) >= zdim_ - 1)
        return false;

    return true;
}


class Facets
{
    typedef CubicalComplex::cell_id_type Cell;

    CubicalComplex complex_;
    bool dual_;

    size_t offset_[512][8];
    int mask_[512][8];
    int count_[512];

    void initialize(size_t const xdim, size_t const ydim, bool const dual);

public:
    Facets(int const xdim, int const ydim, int const zdim, bool const dual)
        : complex_(xdim, ydim, zdim),
          dual_(dual)
    {
        initialize(xdim, ydim, dual);
    }

    int count(Cell const& id) const
    {
        return count_[complex_.type(id)];
    }

    Cell operator()(Cell const& id, int const n) const
    {
        int const t = complex_.type(id);
        assert(n < count_[t]);

        if (dual_)
            return (id | mask_[t][n]) - offset_[t][n];
        else
            return (id ^ mask_[t][n]) + offset_[t][n];
    }
};

void Facets::initialize(size_t const xdim, size_t const ydim, bool const dual)
{
    size_t increment[] = { 8, 8 * xdim, 8 * xdim * ydim };

    for (int t = 0; t < 512; ++t)
    {
        int k = 0;

        if (dual)
        {
            bool const ismin[] =
                { (t & (1 << 3)) > 0, (t & (1 << 5)) > 0, (t & (1 << 7)) > 0 };
            bool const ismax[] =
                { (t & (1 << 4)) > 0, (t & (1 << 6)) > 0, (t & (1 << 8)) > 0 };

            for (int i = 0; i < 3; ++i)
            {
                int const bit = 1 << i;
                if ((t & bit) == 0)
                {
                    if (not ismin[i])
                    {
                        offset_[t][k] = increment[i];
                        mask_[t][k] = bit;
                        ++k;
                    }
                    if (not ismax[i])
                    {
                        offset_[t][k] = 0;
                        mask_[t][k] = bit;
                        ++k;
                    }
                }
            }                
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                int const bit = 1 << i;
                if (t & bit)
                {
                    offset_[t][k] = 0;
                    mask_[t][k] = bit;
                    ++k;
                    offset_[t][k] = increment[i];
                    mask_[t][k] = bit;
                    ++k;
                }
            }
        }

        count_[t] = k;
    }
}


class Vertices
{
    typedef CubicalComplex::cell_id_type Cell;

    CubicalComplex complex_;

    size_t offset_[512][8];
    int count_[512];

    void initialize(size_t const xdim, size_t const ydim);

public:
    Vertices(int const xdim, int const ydim, int const zdim)
        : complex_(xdim, ydim, zdim)
    {
        initialize(xdim, ydim);
    }

    int count(Cell const& id) const
    {
        return count_[complex_.type(id)];
    }

    Cell operator()(Cell const& id, int const n) const
    {
        int const t = complex_.type(id);
        assert(n < count_[t]);

        return (id ^ (id & 7)) + offset_[t][n];
    }
};

void Vertices::initialize(size_t const xdim, size_t const ydim)
{
    size_t increment[] = { 8, 8 * xdim, 8 * xdim * ydim };

    for (int t = 0; t < 512; ++t)
    {
        int k = 0;

        for (int dc = 0; dc <= (t & 4); dc += 4) {
            for (int db = 0; db <= (t & 2); db += 2) {
                for (int da = 0; da <= (t & 1); da += 1) {
                    offset_[t][k] = ((da ? increment[0] : 0) +
                                     (db ? increment[1] : 0) +
                                     (dc ? increment[2] : 0));
                    ++k;
                }
            }
        }

        count_[t] = k;
    }
}


} // namespace diamorse
} // namespace anu_am


#endif // !ANU_AM_DIAMORSE_CUBICALCOMPLEX_H
