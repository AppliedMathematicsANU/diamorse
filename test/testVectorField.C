/** -*-c++-*-
 *
 *  Copyright 2016 The Australian National University
 *
 *  testVectorField.C
 *
 *  Tests gradient vector field computation.
 *
 *  Olaf Delgado-Friedrichs may 16
 *
 */

#include <algorithm>
#include <math.h>
#include <sstream>
#include <vector>

#include "generative.hpp"
#include "common.hpp"
#include "stringUtils.hpp"


using namespace anu_am::generative;
using namespace anu_am::diamorse;
using namespace anu_am::stringutils;


Result alwaysTrue(VolumeData const&)
{
    return success();
}


Result failsWithMessagePrefix(std::string const& prefix, Result const& r)
{
    if (r)
    {
        return failure("check should have failed");
    }
    else if (!startsWith(r.cause(), prefix))
    {
        std::stringstream msg;
        msg << "unexpected failure cause:" << std::endl;
        msg << r.cause() << std::endl;
        msg << "(was expected to start with '" << prefix << "')" << std::endl;

        return failure(msg.str());
    }
    else
    {
        return success();
    }
}


Result cellIsCritical(Cell const& cell, VolumeData const& candidate)
{
    if (candidate.field.isCritical(cell))
        return success();
    else
        return failure("cell is not critical");
}


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


int run()
{
    report("a passing property produces no errors",
           checkWithVolumeData(alwaysTrue));

    report("a failing property produces an appropriate error result",
           failsWithMessagePrefix(
               "\nReason: At cell ",
               checkWithVolumeData(forAllCells(cellIsCritical))));

    report("the vector field is complete",
           checkWithVolumeData(forAllCells(vectorDirectionIsDefined)));

    report("no vectors are marked as pointing outward",
           checkWithVolumeData(forAllCells(vectorIsNotOutwardPointing)));

    report("no vectors are actually pointing outward",
           checkWithVolumeData(forAllCells(directionIsNotOutwardPointing)));

    report("directions are consistent with cell pairings",
           checkWithVolumeData(forAllCells(cellPartnerMatchesDirection)));

    report("the call pairing is symmetrical",
           checkWithVolumeData(forAllCells(partnerOfPartnerIsOriginalCell)));

    report("the discrete vector field is a gradient vector field",
           checkWithVolumeData(containsNoCyclicVPaths));

    report("saturation does not always preserve being a gradient vector field",
           failsWithMessagePrefix(
               "\nReason: cyclic V-path at ",
               checkWithVolumeData(containsNoCyclicVPathsAfterSaturation)));

    report("level sets are hermetic",
           checkWithVolumeData(forAllCells(bind(vImageHasCompatibleValue, 0))));

    std::cerr << std::endl;

    return 0;
}


int main()
{
  try
  {
    run();
  }
  catch(std::runtime_error& e)
  {
    std::clog
      << "terminate called after throwing an instance of "
      << "'std::runtime_error'\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
  catch(std::exception& e)
  {
    std::clog
      << "terminate called after throwing an exception\n"
      << "  what():  " << e.what() << '\n';
    abort();
  }
}
