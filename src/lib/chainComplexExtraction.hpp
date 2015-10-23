/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  chainComplexExtraction.hpp
 *
 *  Computes a Morse chain complex given a gradient vector field.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_CHAINCOMPLEXEXTRACTION_HPP
#define ANU_AM_DIAMORSE_CHAINCOMPLEXEXTRACTION_HPP


#include "CubicalComplex.hpp"
#include "traversals.hpp"
#include "restricted.hpp"


namespace anu_am
{
namespace diamorse
{



template<class Field, class Vectors, class Incidences>
boost::shared_ptr<std::vector<bool> >
connectingPaths(
    CubicalComplex const& complex,
    Field const& field,
    Vectors const& V,
    Vectors const& coV,
    Incidences const& I,
    Incidences const& coI)
{
    typedef CubicalComplex::cell_id_type Cell;

    std::vector<Cell> critical;
    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell) and field.isCritical(cell))
            critical.push_back(cell);

    std::vector<Cell> downstreamSources;
    std::vector<Cell> upstreamSources;
    std::vector<Cell> downstreamTargets;
    std::vector<Cell> upstreamTargets;

    for (size_t i = 0; i < critical.size(); ++i)
    {
        Cell const cell = critical.at(i);
        switch(complex.cellDimension(cell))
        {
        case 0:
            downstreamTargets.push_back(cell);
            break;
        case 1:
            downstreamSources.push_back(cell);
            downstreamTargets.push_back(cell);
            break;
        case 2:
            downstreamSources.push_back(cell);
            upstreamSources.push_back(cell);
            break;
        case 3:
            upstreamTargets.push_back(cell);
            break;
        default:
            break;
        }
    }

    return connectingPaths(complex.cellIdLimit(),
                           downstreamSources, upstreamSources,
                           downstreamTargets, upstreamTargets,
                           V, coV, I, coI);
}


template<class Field>
std::map<CubicalComplex::cell_id_type,
         std::vector<std::pair<CubicalComplex::cell_id_type, int> > >
chainComplex(CubicalComplex const& complex, Field const& field)
{
    typedef CubicalComplex::cell_id_type Cell;
    typedef RestrictedIncidences<Facets> Incidences;
    typedef std::vector<std::pair<Cell, int> > Boundary;

    Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
    Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
    typename Field::Vectors V = field.V();
    typename Field::Vectors coV = field.coV();

    Incidences rI(I, connectingPaths(complex, field, V, coV, I, coI));

    std::map<Cell, Boundary> result;

    for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
        if (complex.isCell(cell) and field.isCritical(cell))
            result[cell] = morseBoundary(cell, V, rI);

    return result;
}



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_CHAINCOMPLEXEXTRACTION_HPP
