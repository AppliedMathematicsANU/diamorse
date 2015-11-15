/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  vectorFieldExtraction.hpp
 *
 *  Computes a gradient vector field for a cell complex with scalar values
 *  defined at its vertices.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#ifndef ANU_AM_DIAMORSE_VECTORFIELDEXTRACTION_HPP
#define ANU_AM_DIAMORSE_VECTORFIELDEXTRACTION_HPP

#include <exception>
#include <vector>

#include "collections.hpp"

namespace anu_am
{
namespace diamorse
{


int const MAX_STAR = 32;

class StarException: public std::exception
{
    virtual const char* what() const throw()
    {
        return "vertex star has too many cells (> 32)";
    }
} starException;


/// Helper class used in fillMorseVectorField().
template<typename Cell, class Scalars>
class StarField
{
    typedef typename Scalars::value_type Value;

    int size_;
    Cell cells_[MAX_STAR];
    Value weight_[MAX_STAR];
    int ranked_[MAX_STAR];
    int defined_[MAX_STAR];
    int incidences_[MAX_STAR][MAX_STAR];
    int incidence_counts_[MAX_STAR];

public:
    StarField(
        Cell const v,
        Facets const& cofacets,
        Vertices const& vertices,
        Scalars const& scalars)
    {
        extract(v, cofacets, vertices, scalars);
    }
    
    int size() const
    {
        return size_;
    }

    bool defined(int const cell) const
    {
        return defined_[cell];
    }

    void set(int const v, int const w)
    {
        defined_[v] = true;
        defined_[w] = true;
    }

    int indexForRank(int const i)
    {
        return ranked_[i];
    }
    
    Cell cell(int const i)
    {
        return cells_[i];
    }

    /// Returns the first face of the given cell that the vector field does
    /// not yet assign a partner to, but only if this is the case for exactly
    /// k faces of that cell. As a special case, if k = 0 and the cell has 0
    /// unused faces, the cell itself is returned.
    /// 
    /// If the cell does not have k unused faces, the function returns a
    /// negative value.
    int firstFreeFaceIfKFree(int const cell, int const k)
    {
        int partner = cell;
        int count = 0;

        int const n = incidence_counts_[cell];
        for (int j = 0; j < n; ++j)
        {
            int const f = incidences_[cell][j];
            if (not defined(f))
            {
                if (++count > k)
                    break;
                partner = f;
            }
        }

        if (count == k)
            return partner;
        else
            return -1;
    }

private:
    void extract(
        Cell const v,
        Facets const& cofacets,
        Vertices const& vertices,
        Scalars const& scalars)
    {
        Value value = scalars.get(v);

        cells_[0] = v;
        defined_[0] = false;
        weight_[0] = value;
        ranked_[0] = 0;
        incidence_counts_[0] = 0;
        size_ = 1;

        int mark = size_;
        int next = 0;
        while (next < size_)
        {
            Cell const cell = cells_[next];

            int const n = cofacets.count(cell);

            for (int i = n - 1; i >= 0; --i)
            {
                Cell const coface = cofacets(cell, i);
                int k;
                for (k = mark; k < size_; ++k)
                    if (coface == cells_[k])
                        break;
                if (k < size_)
                {
                    incidences_[k][incidence_counts_[k]] = next;
                    ++incidence_counts_[k];
                    continue;
                }

                int const m = vertices.count(coface);
                bool add = true;
                Value sum = 0.0;
                for (int j = 0; j < m; ++j)
                {
                    Cell const w = vertices(coface, j);
                    Value const d = scalars.get(w);
                    if (d > value or (d == value and w > v))
                    {
                        add = false;
                        break;
                    }
                    sum += d;
                }

                if (add)
                {
                    cells_[size_] = coface;
                    defined_[size_] = false;
                    weight_[size_] = sum;
                    ranked_[size_] = size_;
                    incidences_[size_][0] = next;
                    incidence_counts_[size_] = 1;
                    ++size_;

                    if (size_ > MAX_STAR)
                        throw starException;
                }
            }

            ++next;
            if (next == mark)
                mark = size_;
        }

        for (int i = 0; i < size_; ++i) {
            int j;
            Value x = weight_[i];
            for (j = i; j > 0 && x < weight_[j-1]; --j) {
                ranked_[j] = ranked_[j-1];
                weight_[j] = weight_[j-1];
            }
            ranked_[j] = i;
            weight_[j] = x;
        }
    }
};


/// Implements a slight variation of the ProcessLowerStar algorithm by
/// Robbins, Sheppard and Wood (TPAMI 33(8), 2011, pp. 1646-1658).
template<class Cell, class Scalars, class VectorField>
void processLowerStar(Cell const v,
                      Scalars const& scalars,
                      VectorField& outputField,
                      Facets const& cofacets,
                      Vertices const& vertices)
{
    // The lower star with the current partial field.
    StarField<Cell, Scalars> star(v, cofacets, vertices, scalars);

    // Determines whether to look for edges or singular cells next. Alternates
    // between 0 and 1.
    int k = 1;

    while (true)
    {
        int cell;
        int partner = -1;

        // Find an unused cell with k unused faces and assign the first such
        // face to partner if k > 0, or the cell itself otherwise.
        for (int i = 0; i < star.size(); ++i)
        {
            cell = star.indexForRank(i);
            if (not star.defined(cell))
                partner = star.firstFreeFaceIfKFree(cell, k);

            if (partner >= 0)
                break;
        }

        if (partner >= 0)
        {
            // New pairing found; record and resume scanning.
            star.set(cell, partner);
            outputField.setPartner(star.cell(cell), star.cell(partner));
            k = 1;
        }
        else if (k == 1) // Scan for singular cells next.
            k = 0;
        else             // Nothing found; terminate inner loop.
            break;
    }
}


/// Computes a Morse vector field for a cell complex.
template<class Complex, class Scalars, class VectorField>
void fillMorseVectorField(Complex const& source,
                          Scalars const& scalars,
                          VectorField& outputField)
{
    Facets cofacets(source.xdim(), source.ydim(), source.zdim(), true);
    Vertices vertices(source.xdim(), source.ydim(), source.zdim());

    for (size_t x = 0; x < source.xdim(); ++x)
        for (size_t y = 0; y < source.ydim(); ++y)
            for (size_t z = 0; z < source.zdim(); ++z)
                processLowerStar(source.cellAt(x, y, z),
                                 scalars, outputField,
                                 cofacets, vertices);
}


} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_VECTORFIELDEXTRACTION_HPP
