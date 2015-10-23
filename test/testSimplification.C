/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  testSimplification.C
 *
 *  Tests the gradient vector field simplification algorithm.
 *
 *  Olaf Delgado-Friedrichs mar 14
 *
 */

#include "generative.hpp"
#include "common.hpp"
#include "booster.hpp"

#include "simplification.hpp"

using namespace anu_am::generative;
using namespace anu_am::generative::booster;
using namespace anu_am::diamorse;


struct Simplified
{
    typedef VolumeData argument_type;
    typedef VolumeData result_type;

    Value const threshold;

    Simplified(Value const threshold)
        : threshold(threshold)
    {
    }

    VolumeData operator()(VolumeData const& original) const
    {
        CubicalComplex const& complex = original.complex;
        Scalars const& scalars = original.scalars;
        Field const& field = original.field;

        Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());

        Field::DataPtr data(new Field::DataPtr::element_type(*field.data()));
        Field clonedField(complex.xdim(), complex.ydim(), complex.zdim(), data);

        std::stringstream devNull;

        simplify(complex,
                 clonedField,
                 withArgument(maxima(scalars, vertices)),
                 mayCancel(complex, scalars, field, threshold),
                 devNull);

        return VolumeData(complex, original.scalars, clonedField);
    }
};


static float const THRESHOLD = 1;


template<typename P>
anu_am::generative::Result checkWithSimplifiedVolumeData(
    P const& predicate,
    int const N = 500)
{
    return checkWithVolumeData(composition(predicate, Simplified(THRESHOLD)),
                               N);
}


// === Tests start here.


SIMPLE_TEST_CASE(
    theVectorFieldIsComplete,
    checkWithSimplifiedVolumeData(forAllCells(vectorDirectionIsDefined)))

SIMPLE_TEST_CASE(
    noVectorsAreMarkedAsPointingOutward,
    checkWithSimplifiedVolumeData(forAllCells(vectorIsNotOutwardPointing)))

SIMPLE_TEST_CASE(
    noVectorsAreActuallyPointingOutward,
    checkWithSimplifiedVolumeData(forAllCells(directionIsNotOutwardPointing)))

SIMPLE_TEST_CASE(
    directionsAndCellPartnersMatch,
    checkWithSimplifiedVolumeData(forAllCells(cellPartnerMatchesDirection)))

SIMPLE_TEST_CASE(
    theCellPairingIsSymmetrical,
    checkWithSimplifiedVolumeData(forAllCells(partnerOfPartnerIsOriginalCell)))

SIMPLE_TEST_CASE(
    theDiscreteVectorFieldIsAGradientVectorField,
    checkWithSimplifiedVolumeData(containsNoCyclicVPaths))


Result checkOriginalVersusSimplifiedHomology(VolumeData const& candidate)
{
    return checkPersistentHomology(
        convertedChainComplex(candidate),
        convertedChainComplex(Simplified(THRESHOLD)(candidate)),
        THRESHOLD);
}

SIMPLE_TEST_CASE(
    originalAndSimplifiedComplexHaveTheSamePersistentHomoloy,
    checkWithVolumeData(checkOriginalVersusSimplifiedHomology))


SIMPLE_TEST_CASE(
    vectorFieldIsCompatibleWithScalars,
    checkWithSimplifiedVolumeData(
        forAllCells(bind(vImageHasCompatibleValue, THRESHOLD))))


Result checkAbsenceOfCancellableClosePairs(VolumeData const& candidate)
{
    typedef CubicalComplex::cell_id_type Cell;
    typedef std::vector<std::pair<Cell, int> > Boundary;

    CubicalComplex const& complex = candidate.complex;
    Field const& field = candidate.field;
    Scalars const& scalars = candidate.scalars;

    Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());
    Maxima<Scalars, Vertices> value(scalars, vertices);

    std::map<Cell, Boundary> chains = chainComplex(complex, field);
    std::map<Cell, Boundary> cochains = reverse(chains);

    for (Cell cell = 0; cell < complex.cellIdLimit(); ++cell)
    {
        if (not complex.isCell(cell) or not field.isCritical(cell))
            continue;

        std::pair<Cell, int> const res = closePartner(cell,
                                                      callableMap(chains),
                                                      callableMap(cochains),
                                                      withArgument(value));

        if (res.second == 1 and value(cell) - value(res.first) <= THRESHOLD)
        {
            std::stringstream msg;
            msg << "found cancellable pair "
                << complex.cellPosition(cell)
                << " (" << value(cell) << ")"
                << " <-> "
                << complex.cellPosition(res.first)
                << " (" << value(res.first) << ")";
            return failure(msg.str());
        }
    }

    return success();
}

SIMPLE_TEST_CASE(
    simplificationLeavesNoCancellableClosePairs,
    checkWithSimplifiedVolumeData(checkAbsenceOfCancellableClosePairs))
