/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  persistence.hpp
 *
 *  Persistence computation for general chain complexes.
 *
 *  Olaf Delgado-Friedrichs feb 14
 *
 */

#ifndef ANU_AM_DIAMORSE_PERSISTENCE_HPP
#define ANU_AM_DIAMORSE_PERSISTENCE_HPP

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

namespace anu_am
{
namespace diamorse
{


/// This class implements a column in the incidence matrix for the purpose of
/// persistence pair computation.
template<typename T>
class Column
{
    typedef typename std::vector<T>::const_iterator Iter;

public:
    // Create an empty column.
    Column()
        : elements_()
    {
    }

    // Construct a column with a single entry.
    Column(T const& t)
        : elements_()
    {
        elements_.push_back(t);
    }

    // Construct a column from a vector. Only includes elements that occur an
    // odd number of times.
    explicit Column(std::vector<T> const& source)
        :elements_()
    {
        std::vector<T> tmp(source);
        std::sort(tmp.begin(), tmp.end());

        for (Iter it = tmp.begin(); it != tmp.end(); ++it)
            if (elements_.size() > 0 and *it == elements_.back())
                elements_.pop_back();
            else
                elements_.push_back(*it);
    }

    // Checks whether the column is empty.
    bool empty() const { return elements_.empty(); }
    
    // Returns the highest entry in the column.
    T max() const { return elements_.back(); }

    // Creates a new column with all entries that are in either this column or
    // the given one, but not in both. Equivalent to in-place addition modulo
    // 2.
    Column symmetricDifference(Column const& other) const
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
                ++left;
                ++right;
            }
        }
        if (left != left_end)
            merged.insert(merged.end(), left, left_end);
        else if (right != right_end)
            merged.insert(merged.end(), right, right_end);

        return Column(merged, true);
    }

private:
    explicit Column(std::vector<T> const& source, bool const)
        :elements_(source)
    {
    }

    std::vector<T> elements_;
};


template <typename T>
struct Pairing
{
    T partner;
    int dimension;
    float value;
    std::vector<T> chain;

    Pairing()
    {
    }

    Pairing(T const partner, int const dimension, float const value)
        : partner(partner),
          dimension(dimension),
          value(value),
          chain()
    {
    }
};


template<typename T>
int maxdim(std::vector<Pairing<T> > const& pairs)
{
    int d = -1;

    for (size_t i = 0; i < pairs.size(); ++i)
        if (pairs.at(i).dimension > d)
            d = pairs.at(i).dimension;

    return d;
}


template<typename Cell, typename T>
std::pair<Column<Cell>, std::vector<Cell> > reducedColumn(
    T const& cells,
    std::vector<Pairing<Cell> > const& pairing,
    size_t const n,
    std::map<Cell, Column<Cell> > const& M)
{
    Column<Cell> reduced(cells);
    std::vector<Cell> chains;
    Cell k;

    while (not reduced.empty() and (k = pairing.at(reduced.max()).partner) < n)
    {
        reduced = reduced.symmetricDifference(M.at(k));
        chains.push_back(k);
    }

    return std::make_pair(reduced, chains);
}


/// Computes a pairing of cells from the given complex such that the cell with
/// the lower index in a pair creates a homological d-cycle (over Z2), where d
/// is the dimension of that cell, and the one with the higher index, which
/// must be (d+1)-dimensional, destroys that cycle.
///
/// The vector returned contains the partner in this pairing for each cell,
/// with a value higher than complex.nrCells() indicating that the cycle
/// created by that cell is permanent.
///
/// The algorithm requires that cell indexes are consecutive, starting at 0,
/// and that the index of each cell in the complex is larger than the indexes
/// of all its faces.
///
/// The function uses the modified matrix-reduction method for persistent
/// homology computation as introduced by Chen and Kerber (in EuroCG 2011). In
/// order to emphasize that the algorithm as presented is restricted to
/// coefficients in Z2, we use set-theoretic rather than algebraic
/// terminology, in particular symmetric difference rather than addition
/// (modulo 2).
///
/// Columns in the incidence matrix are created on demand and only stored
/// during the lifetime of the outer loop's body, i.e. while reducing columns
/// associated to the current cell dimension. Instead of "killing" rows
/// corresponding to lower-dimensional cells, it is simply checked whether a
/// cell has already been paired when it is time for its column to be reduced.
template<typename Complex>
std::vector<Pairing<typename Complex::cell_id_type> >
persistencePairing(Complex const& complex)
{
    typedef typename Complex::cell_id_type Cell;
    typedef std::pair<Column<Cell>, std::vector<Cell> > R;

    Cell const n = complex.nrCells();
    int const dim = complex.dimension();

    // Initialise the vector that will represent the pairing by provisionally
    // assuming that all cells create permanent cycles.
    std::vector<Pairing<Cell> > pairing(n);
    for (Cell j = 0; j < n; ++j)
    {
        int const d = complex.cellDimension(j);
        float const v = complex.cellValue(j);
        pairing.at(j) = Pairing<Cell>(n+1, d, v);
    }

    // Perform the reduction by decreasing dimension.
    for (int d = dim; d > 0; --d)
    {
        // Create a container for storing incidence columns.
        std::map<Cell, Column<Cell> > M;
        
        for (Cell j = 0; j < n; ++j)
        {
            if (complex.cellDimension(j) == d and pairing.at(j).partner >= n)
            {
                std::vector<Cell> const faces = complex.cellFaces(j);
                R const r = reducedColumn(faces, pairing, n, M);
                Column<Cell> const column = r.first;

                if (not column.empty())
                {
                    // Record the new cancellation pair.
                    Cell const i = column.max();
                    pairing.at(i).partner = j;
                    pairing.at(j).partner = i;
                    pairing.at(i).chain   = r.second;

                    // Store the reduced column.
                    M[j] = column;
                }
            }
        }
    }

    // Return the resulting pairing.
    return pairing;
}

/// Computes the p-persistent Betti numbers for the level sets of the given
/// cell complex with respect to the values associated to cells by the
/// non-decreasing function f. The lifetime p is also measured in terms of
/// this function.
///
/// The output is a map in which for each relevant level the Betti numbers for
/// all relevant dimensions are listed. Only levels at which the Betti number
/// changes are included in the output.
///
/// The algorithm requires that cell indexes are consecutive, starting at 0,
/// and that the index of each cell in the complex is larger than the indexes
/// of all its faces.
///
/// The vector pairing must contain the result of running persistencePairing()
/// on the complex or a suitable equivalent.
template<typename T>
std::map<float, std::vector<int> > bettiNumbersUncollated(
    std::vector<Pairing<T> > const& pairing,
    float const p = 0)
{
    T const n = pairing.size();
    int const dim = maxdim(pairing);
    std::map<float, std::vector<int> > result;
    std::vector<int> current(dim + 1);
    std::vector<int> last(dim + 1);

    for (T cell = 0; cell < n; ++cell)
    {
        assert(cell+1 >= n or
               pairing.at(cell+1).value >= pairing.at(cell).value);

        T const partner = pairing.at(cell).partner;
        assert(partner >= n or pairing.at(partner).partner == cell);

        int const d = pairing.at(cell).dimension;

        if (partner > cell)
        {
            assert(partner >= n or pairing.at(partner).dimension == d + 1);
            if (partner >= n
                or pairing.at(partner).value > pairing.at(cell).value + p)
            {
                ++current.at(d);
            }
        }
        else if (partner < cell)
        {
            assert(pairing.at(partner).dimension == d - 1);
            if (pairing.at(partner).value < pairing.at(cell).value - p)
                --current.at(d - 1);
        }

        if ((cell+1 >= n or pairing.at(cell+1).value > pairing.at(cell).value)
            and (result.empty() or current != last))
        {
            last = current;
            result[pairing.at(cell).value] = current;
        }
    }

    return result;
}


/// Computes the p-persistent Betti numbers for the level sets of the given
/// cell complex with respect to the values associated to cells by the
/// non-decreasing function f. The lifetime p is also measured in terms of
/// this function.
///
/// The output is a vector of maps with one map for each dimension. Each map
/// contains the scalar levels at which the Betti numbers change as keys and
/// the corresponding numbers as values.
///
/// The algorithm requires that cell indexes are consecutive, starting at 0,
/// and that the index of each cell in the complex is larger than the indexes
/// of all its faces.
///
/// The vector pairing must contain the result of running persistencePairing()
/// on the complex or a suitable equivalent.
template<typename T>
std::vector<std::map<float, int> > bettiNumbers(
    std::vector<Pairing<T> > const& pairing,
    float const p = 0)
{
    typedef std::map<float, std::vector<int> > Betti;

    int const dim = maxdim(pairing);

    Betti betti = bettiNumbersUncollated(pairing, p);

    std::vector<std::map<float, int> > result;

    for (int d = 0; d <= dim; ++d)
    {
        int last = 0;
        std::map<float, int> b;

        Betti::const_iterator iter;
        for (iter = betti.begin(); iter != betti.end(); ++iter)
        {
            int current = iter->second.at(d);
            if (iter == betti.begin() or current != last)
            {
                b[iter->first] = current;
                last = current;
            }
        }

        result.push_back(b);
    }

    return result;
}



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_PERSISTENCE_HPP
