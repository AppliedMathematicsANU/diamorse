/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  common.hpp
 *
 *  Common code for generative testing of Morse theory algorithma
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */

#ifndef ANU_AM_DIAMORSE_TEST_BOOSTER_HPP
#define ANU_AM_DIAMORSE_TEST_BOOSTER_HPP

#define BOOST_TEST_MODULE vectorField
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

namespace anu_am
{
namespace generative
{
namespace booster
{

boost::test_tools::predicate_result boostify(
    anu_am::generative::Result const result)
{
    if (!result)
    {
        boost::test_tools::predicate_result r(false);
        r.message() << result.cause();
        return r;
    }
    else
        return true;
}

#define SIMPLE_TEST_CASE(name, predicate) \
    BOOST_AUTO_TEST_CASE(name) { BOOST_REQUIRE(boostify(predicate)); }


} // namespace boost
} // namespace generative
} // namespace anu_am

#endif //!ANU_AM_DIAMORSE_TEST_BOOSTER_HPP
