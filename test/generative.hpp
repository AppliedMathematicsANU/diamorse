/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  generative.hpp
 *
 *  A simple generative testing framework inspired by Haskell's QuickCheck.
 *
 *  Olaf Delgado-Friedrichs jan 14
 *
 */


#ifndef ANU_AM_GENERATIVE_HPP
#define ANU_AM_GENERATIVE_HPP

#include <iostream>
#include <sstream>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace anu_am
{
namespace generative
{

typedef boost::mt19937 RandomNumberGenerator;
RandomNumberGenerator rng(time(0));


template<class T, class D>
T random(D const& d)
{
    return boost::variate_generator<boost::mt19937&, D>(rng, d)();
}


float randomFloat(float const sigma = 5.0, float const mean = 0.0)
{
    return random<float>(boost::normal_distribution<float>(mean, sigma));
}

int randomInt(int const limit)
{
    return random<int>(boost::uniform_int<>(0, limit));
}


template<typename G>
struct function_traits
{
    typedef typename G::result_type result_type;
    typedef typename G::argument_type argument_type;
};

template<typename R>
struct function_traits<R()>
{
    typedef R result_type;
};

template<typename R, typename A>
struct function_traits<R(A)>
{
    typedef R result_type;
    typedef A argument_type;
};

template<typename R, typename A, typename B>
struct function_traits<R(A, B)>
{
    typedef R result_type;
};

template<typename R, typename A, typename B, typename C>
struct function_traits<R(A, B, C)>
{
    typedef R result_type;
};

template<class K, typename R>
struct function_traits<R(K::*)()>
{
    typedef R result_type;
};

template<class K, typename R, typename A>
struct function_traits<R(K::*)(A)>
{
    typedef R result_type;
    typedef A argument_type;
};

template<class K, typename R, typename A, typename B>
struct function_traits<R(K::*)(A, B)>
{
    typedef R result_type;
};


template<typename F, typename G>
struct Composition
{
    typedef typename function_traits<F>::result_type result_type;
    typedef typename function_traits<G>::argument_type argument_type;

    F const& f;
    G const& g;

    Composition(F const& f, G const& g)
        : f(f),
          g(g)
    {
    }

    result_type operator()(argument_type const x) const
    {
        return f(g(x));
    }
};

template<typename F, typename G>
struct Composition<F, G> composition(F const& f, G const& g)
{
    return Composition<F, G>(f, g);
}


class Result
{
    bool successful_;
    std::string cause_;

public:
    Result(bool const successful, std::string const cause = "")
        : successful_(successful),
          cause_(cause)
    {
    }

    std::string cause() const
    {
        return cause_;
    }

    operator bool() const
    {
        return successful_;
    }
};



Result failure(std::string const& cause)
{
    return Result(false, cause);
}


Result success(std::string const& cause = "")
{
    return Result(true, cause);
}


template<typename P, typename C, typename S>
std::pair<C, typename function_traits<P>::result_type>
shrink(P const& predicate, C const& candidate, S const& shrinker)
{
    C smallest = candidate;
    bool done = false;

    while (not done)
    {
        std::vector<C> const& shrunk = shrinker(smallest);
        done = true;

        for (size_t i = 0; i < shrunk.size(); ++i)
        {
            if (!predicate(shrunk.at(i)))
            {
                smallest = shrunk.at(i);
                done = false;
                break;
            }
        }
    }

    return std::make_pair(smallest, predicate(smallest));
}


template<typename P, typename G, typename S>
Result checkPredicate(
    P const& predicate,
    G const& generator,
    S const& shrinker,
    int const N = 100)
{
    typedef typename function_traits<G>::result_type Instance;

    for (int i = 0; i < N; ++i)
    {
        Instance const& candidate = generator(i);

        if (!predicate(candidate))
        {
            std::pair<Instance, Result> const shrunk =
                shrink(predicate, candidate, shrinker);

            std::stringstream msg;

            msg << std::endl
                << "Reason: " << shrunk.second.cause() << std::endl
                << "     in " << shrunk.first << std::endl
                << "  (from " << candidate << ")" << std::endl;

            return failure(msg.str());
        }
    }
    return success();
}

} // namespace generative
} // namespace anu_am

#endif //!ANU_AM_GENERATIVE_HPP
