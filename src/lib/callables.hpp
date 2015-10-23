/** -*-c++-*-
 *
 *  Copyright 2013 The Australian National University
 *
 *  callables.hpp
 *
 *  Some generic function-like objects.
 *
 *  Olaf Delgado-Friedrichs mar 14
 *
 */

#ifndef ANU_AM_DIAMORSE_CALLABLES_HPP
#define ANU_AM_DIAMORSE_CALLABLES_HPP

#include <algorithm>
#include <map>

namespace anu_am
{
namespace diamorse
{


template<typename K, typename V>
class CallableMap
{
public:
    typedef V result_type;

private:
    typedef std::map<K, V> Data;
    Data const& data_;
    V const& defaultValue_;

public:
    CallableMap(Data const& data, V const& defaultValue = V())
        : data_(data),
          defaultValue_(defaultValue)
    {
    }

    V operator()(K const& key) const
    {
        if (data_.count(key) > 0)
            return data_.at(key);
        else
            return defaultValue_;
    }
};

template<typename K, typename V>
CallableMap<K, V>
callableMap(std::map<K, V> const& data, V const& defaultValue = V())
{
    return CallableMap<K, V>(data, defaultValue);
}


template<class F>
class WithArgument
{
    F const& evaluator_;

public:
    typedef typename F::argument_type argument_type;
    typedef std::pair<typename F::result_type, argument_type> result_type;

    WithArgument(F const& evaluator)
        : evaluator_(evaluator)
    {
    }

    result_type operator()(argument_type const& arg) const
    {
        return std::make_pair(evaluator_(arg), arg);
    }
};

template<class F>
WithArgument<F> withArgument(F const& evaluator)
{
    return WithArgument<F>(evaluator);
}


template<class F, class T>
class Maxima
{
    F const& evaluator_;
    T const& sampler_;

public:
    typedef typename F::argument_type argument_type;
    typedef typename F::result_type   result_type;

    Maxima(F const& evaluator, T const& sampler)
        : evaluator_(evaluator),
          sampler_(sampler)
    {
    }

    result_type operator()(argument_type const arg) const
    {
        size_t const n = sampler_.count(arg);
        result_type val = evaluator_(sampler_(arg, 0));

        for (size_t i = 1; i < n; ++i)
            val = std::max(val, evaluator_(sampler_(arg, i)));

        return val;
    }
};

template<class F, class T>
Maxima<F, T> maxima(F const& evaluator, T const& sampler)
{
    return Maxima<F, T>(evaluator, sampler);
}



} // namespace diamorse
} // namespace anu_am

#endif //!ANU_AM_DIAMORSE_CALLABLES_HPP
