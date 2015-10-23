/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  testVectorField.C
 *
 *  Tests gradient vector field computation.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#include <algorithm>
#include <math.h>
#include <sstream>
#include <vector>

#include "generative.hpp"
#include "common.hpp"
#include "booster.hpp"

using namespace anu_am::generative;
using namespace anu_am::generative::booster;
using namespace anu_am::diamorse;


// === First some sanity checks for the test harness.


BOOST_AUTO_TEST_CASE(aFixedTestInstanceCanBeMade)
{
    Value data[] = { 0, 1, 2,   3, 4, 5,   6, 7, 8,
                     0, 3, 6,   1, 4, 7,   2, 5, 8,
                     0, 1, 0,   1, 0, 1,   0, 1, 0  };

    BOOST_REQUIRE_NO_THROW(fixedVolumeData(3, 3, 3, data));
}

BOOST_AUTO_TEST_CASE(aRandomTestInstanceCanBeMade)
{
    BOOST_REQUIRE_NO_THROW(randomVolumeData(100));
}


Result alwaysTrue(VolumeData const&)
{
    return success();
}

BOOST_AUTO_TEST_CASE(aPassingPropertyProducesNoErrors)
{
    BOOST_REQUIRE(boostify(checkWithVolumeData(alwaysTrue)));
}


Result cellIsCritical(Cell const& cell, VolumeData const& candidate)
{
    if (candidate.field.isCritical(cell))
        return success();
    else
        return failure("cell is not critical");
}

BOOST_AUTO_TEST_CASE(aFailingPropertyProducesAnAppropriateErrorResult)
{
    boost::test_tools::predicate_result r =
        boostify(checkWithVolumeData(forAllCells(cellIsCritical)));
    std::string expected = "Reason: At cell ";
    BOOST_REQUIRE(!r);
    BOOST_REQUIRE_EQUAL(r.message().str().substr(1, expected.length()),
                        expected);
}


// === Tests for the gradient vector field functionality start here.

SIMPLE_TEST_CASE(
    theVectorFieldIsComplete,
    checkWithVolumeData(forAllCells(vectorDirectionIsDefined)))


SIMPLE_TEST_CASE(
    noVectorsAreMarkedAsPointingOutward,
    checkWithVolumeData(forAllCells(vectorIsNotOutwardPointing)))


SIMPLE_TEST_CASE(
    noVectorsAreActuallyPointingOutward,
    checkWithVolumeData(forAllCells(directionIsNotOutwardPointing)))


SIMPLE_TEST_CASE(
    directionsAndCellPartnersMatch,
    checkWithVolumeData(forAllCells(cellPartnerMatchesDirection)))


SIMPLE_TEST_CASE(
    theCellPairingIsSymmetrical,
    checkWithVolumeData(forAllCells(partnerOfPartnerIsOriginalCell)))


SIMPLE_TEST_CASE(
    theDiscreteVectorFieldIsAGradientVectorField,
    checkWithVolumeData(containsNoCyclicVPaths))


VolumeData withSaturatedVectorField(VolumeData const& original)
{
    CubicalComplex const& complex = original.complex;
    Facets const facets(complex.xdim(), complex.ydim(), complex.zdim(), false);

    Field::DataPtr const& oldData = original.field.data();
    Field::DataPtr newData(new Field::DataPtr::element_type(*oldData));
    Field field(complex.xdim(), complex.ydim(), complex.zdim(), newData);

    for (Cell start = 0; start < complex.cellIdLimit(); ++start)
    {
        if (not complex.isCell(start) or not field.isCritical(start))
            continue;
        int const n = facets.count(start);
        for (int i = 0; i < n; ++i)
        {
            Cell const f = facets(start, i);
            if (field.isCritical(f))
                field.setPartner(start, f);
        }
    }

    return VolumeData(complex, original.scalars, field);
}


Result containsNoCyclicVPathsAfterSaturation(VolumeData const& candidate)
{
    return containsNoCyclicVPaths(withSaturatedVectorField(candidate));
}

BOOST_AUTO_TEST_CASE(notAlwaysAGradientVectorFieldAfterSaturation)
{
    boost::test_tools::predicate_result r =
        boostify(checkWithVolumeData(containsNoCyclicVPathsAfterSaturation));
    std::string expected = "Reason: cyclic V-path at ";
    BOOST_REQUIRE(!r);
    BOOST_REQUIRE_EQUAL(r.message().str().substr(1, expected.length()),
                        expected);
}


SIMPLE_TEST_CASE(
    levelSetsAreHermetic,
    checkWithVolumeData(forAllCells(bind(vImageHasCompatibleValue, 0))))
