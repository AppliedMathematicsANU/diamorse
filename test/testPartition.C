/** -*-c++-*-
 *
 *  Copyright 2016 The Australian National University
 *
 *  testPartition.C
 *
 *  Tests the Partition (union-find) class.
 *
 *  Olaf Delgado-Friedrichs may 16
 *
 */

#include <stdint.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "generative.hpp"

#include "collections.hpp"
#include "Partition.hpp"


using namespace anu_am::generative;
using namespace anu_am::diamorse;


template<typename T>
class Model
{
    std::vector<std::set<T> > setFor_;

public:
    Model(T const size = 0)
        : setFor_(size)
    {
        for (T i = 0; i < size; ++i)
            setFor_[i].insert(i);
    }

    void unite(T const a, T const b)
    {
        std::set<T> sa = setFor_.at(a);
        std::set<T> sb = setFor_.at(b);
        std::set<T> setunion;
        setunion.insert(sa.begin(), sa.end());
        setunion.insert(sb.begin(), sb.end());

        typename std::set<T>::const_iterator iter;
        for (iter = setunion.begin(); iter != setunion.end(); ++iter)
            setFor_.at(*iter) = setunion;
    }

    std::vector<std::set<T> > sets() const
    {
        std::vector<std::set<T> > result;

        for (T i = 0; i < setFor_.size(); ++i)
        {
            std::set<T> s = setFor_.at(i);
            if (s.size() > 1 and *s.begin() == i)
                result.push_back(s);
        }
        return result;
    }
};


template<typename T>
class Implementation
{
    Partition<T> partition_;

public:
    Implementation(T const size = 0)
        : partition_(size)
    {
    }

    void unite(T const a, T const b)
    {
        partition_.unite(a, b);
    }

    std::vector<std::set<T> > sets() const
    {
        std::map<T, std::set<T> > tmp;

        for (T i = 0; i < partition_.size(); ++i)
            tmp[partition_.find(i)].insert(i);

        std::vector<std::set<T> > result;

        for (T i = 0; i < partition_.size(); ++i)
            if (tmp[i].size() > 1)
                result.push_back(tmp.at(i));
        std::sort(result.begin(), result.end());

        return result;
    }
};


template<typename T>
class Session
{
    Model<T> model_;
    Implementation<T> implementation_;

    T size_;
    std::vector<std::pair<T, T> > unions_;

public:
    Session(T const size, std::vector<std::pair<T, T> > const unions)
        : model_(Model<T>(size)),
          implementation_(Implementation<T>(size)),
          size_(size),
          unions_(unions)
    {
    }

    T size() const
    {
        return size_;
    }

    std::vector<std::pair<T, T> > unions() const
    {
        return unions_;
    }

    Result operator()()
    {
        for (size_t i = 0; i < unions_.size(); ++i)
        {
            T const a = unions_.at(i).first;
            T const b = unions_.at(i).second;

            model_.unite(a, b);
            implementation_.unite(a, b);
        }

        std::vector<std::set<T> > expected = model_.sets();
        std::vector<std::set<T> > found    = implementation_.sets();

        if (expected == found)
        {
            return success();
        }
        else
        {
            std::stringstream msg;
            msg << "expected " << expected << ", got " << found
                << std::endl;
            return failure(msg.str());
        }
    }
};


template<typename T>
std::ostream& operator<<(std::ostream& out, Session<T> const& session)
{
    return out << session.size() << ", " << session.unions();
}


template<typename T>
Session<T> randomSession(T const size)
{
    T const n = randomInt(size);
    size_t const m = randomInt(size);
    std::vector<std::pair<T, T> > unions;

    if (n > 0)
        for (size_t i = 0; i < m; ++i)
            unions.push_back(std::make_pair(randomInt(n-1), randomInt(n-1)));

    return Session<T>(n, unions);
}


template<typename T>
std::vector<Session<T> > shrinkSession(Session<T> const session)
{
    T const size = session.size();
    std::vector<std::pair<T, T> > const unions = session.unions();
    std::vector<Session<T> > result;

    std::vector<std::pair<T, T> > tmp;
    for (size_t i = 0; i < unions.size(); ++i)
        if (unions.at(i).first < size-1 and unions.at(i).second < size-1)
            tmp.push_back(unions.at(i));
    result.push_back(Session<T>(size-1, tmp));

    for (size_t i = 0; i < unions.size(); ++i)
    {
        std::vector<std::pair<T, T> > tmp(unions);
        tmp.erase(tmp.begin() + i);
        result.push_back(Session<T>(size, tmp));
    }

    return result;
}


template<typename F>
Result call(F fun)
{
    return fun();
}


template<typename T>
Result checkPartition()
{
    return checkPredicate(call<Session<T> >,
                          randomSession<T>,
                          shrinkSession<T>);
}


int run()
{
    report("partitions work properly", checkPartition<uint16_t>());

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
