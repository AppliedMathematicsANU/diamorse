/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  json.hpp
 *
 *  Some code for writing JSON-formatted data.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#ifndef ANU_AM_JSON_HPP
#define ANU_AM_JSON_HPP

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <memory>


namespace anu_am
{
namespace json
{


typedef enum {
    JSON_NULL   = 0,
    JSON_FALSE  = 1,
    JSON_TRUE   = 2,
    JSON_NUMBER = 3,
    JSON_STRING = 4,
    JSON_ARRAY  = 5,
    JSON_OBJECT = 6
} Type;


struct Value
{
    virtual ~Value() {};

    virtual Type type() const = 0;
    virtual void write(std::ostream&, size_t = 0, size_t = 0, std::string = "")
        const = 0;
};


typedef std::shared_ptr<Value const> ValuePointer;


struct Null : public Value
{
    Null()
    {
    }

    Type type() const
    {
        return JSON_NULL;
    }

    void write(std::ostream& out, size_t = 0, size_t = 0, std::string = "")
    const
    {
        out << "null";
    }
};


struct False : public Value
{
    False()
    {
    }

    Type type() const
    {
        return JSON_FALSE;
    }

    void write(std::ostream& out, size_t = 0, size_t = 0, std::string = "")
    const
    {
        out << "false";
    }
};


struct True : public Value
{
    True()
    {
    }

    Type type() const
    {
        return JSON_TRUE;
    }

    void write(std::ostream& out, size_t = 0, size_t = 0, std::string = "")
    const
    {
        out << "true";
    }
};


class Number : public Value
{
    double const _value;
public:
    Number(double const value) : _value(value)
    {
    }

    Type type() const
    {
        return JSON_NUMBER;
    }

    double value() const
    {
        return _value;
    }

    void write(std::ostream& out, size_t = 0, size_t = 0, std::string = "")
    const
    {
        out << _value;
    }
};


class String : public Value
{
    std::string const _value;
public:
    String(std::string const value) : _value(value)
    {
    }

    Type type() const
    {
        return JSON_STRING;
    }

    std::string value() const
    {
        return _value;
    }

    void write(std::ostream& out, size_t = 0, size_t = 0, std::string = "")
    const
    {
        out << '"';

        for (size_t i = 0; i < _value.size(); ++i)
        {
            char const c = _value.at(i);

            switch (c)
            {
            case '"' : out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\\': out << "\\\\"; break;
            default  : out << c; break;
            }
        }

        out << '"';
    }
};


void next(std::ostream& out,
          size_t const wd, size_t const left, std::string const prefix)
{
    if (wd == 0)
        out << " ";
    else
    {
        out << std::endl << prefix;
        for (size_t i = 0; i < left; ++ i)
            out << " ";
    }
}


ValuePointer makeValuePointer()
{
    return ValuePointer(new Null());
}

ValuePointer makeValuePointer(bool const val)
{
    if (val)
        return ValuePointer(new True());
    else
        return ValuePointer(new False());
}

ValuePointer makeValuePointer(int const val)
{
    return ValuePointer(new Number(val));
}

ValuePointer makeValuePointer(size_t const val)
{
    return ValuePointer(new Number(val));
}

ValuePointer makeValuePointer(double const val)
{
    return ValuePointer(new Number(val));
}

ValuePointer makeValuePointer(std::string const val)
{
    return ValuePointer(new String(val));
}

ValuePointer makeValuePointer(char const val[])
{
    return ValuePointer(new String(val));
}

class Array;
ValuePointer makeValuePointer(Array const val);

class Object;
ValuePointer makeValuePointer(Object const val);


class Array : public Value
{
    std::vector<ValuePointer> _value;

public:
    Array()
    {
    }

    template<typename T>
    Array(T const val)
    {
        _value.push_back(makeValuePointer(val));
    }

    template<typename T>
    Array(std::vector<T> const& a)
    {
        for (size_t i = 0; i < a.size(); ++i)
            _value.push_back(makeValuePointer(a.at(i)));
    }

    Type type() const
    {
        return JSON_ARRAY;
    }

    size_t size() const
    {
        return _value.size();
    }

    void write(
        std::ostream& out,
        size_t const wd = 0,
        size_t const left = 0,
        std::string const prefix = "") const
    {
        out << "[";

        for (size_t i = 0; i < size(); ++i)
        {
            if (i > 0)
                out << ",";

            next(out, wd, left + wd, prefix);
            at(i)->write(out, wd, left + wd, prefix);
        }

        next(out, wd, left, prefix);
        out << "]";
    }

    ValuePointer at(size_t i) const
    {
        return _value.at(i);
    }

    Array operator()()
    {
        _value.push_back(makeValuePointer());
        return *this;
    }

    template<typename T>
    Array operator()(T const val)
    {
        _value.push_back(makeValuePointer(val));
        return *this;
    }
};


class Object : public Value
{
    std::map<std::string, ValuePointer> _map;
    std::vector<std::string> _keysInOrder;

    void set(std::string const key, ValuePointer const value)
    {
        if (_map.count(key) == 0)
            _keysInOrder.push_back(key);
        _map[key] = value;
    }

public:
    Object()
    {
    }

    template<typename T>
    Object(std::string const key, T const val)
    {
        set(key, makeValuePointer(val));
    }

    Type type() const
    {
        return JSON_OBJECT;
    }

    size_t size() const
    {
        return _keysInOrder.size();
    }

    std::string keyAt(size_t const i) const
    {
        return _keysInOrder.at(i);
    }

    ValuePointer at(std::string const key) const
    {
        return _map.at(key);
    }

    void write(
        std::ostream& out,
        size_t const wd = 0,
        size_t const left = 0,
        std::string const prefix = "") const
    {
        out << "{";

        for (size_t i = 0; i < size(); ++i)
        {
            std::string const key = keyAt(i);

            if (i > 0)
                out << ",";

            next(out, wd, left + wd, prefix);
            out << '"' << key << "\": ";
            at(key)->write(out, wd, left + wd, prefix);
        }

        next(out, wd, left, prefix);
        out << "}";
    }

    Object operator()(std::string const key)
    {
        set(key, makeValuePointer());
        return *this;
    }

    template<typename T>
    Object operator()(std::string const key, T const val)
    {
        set(key, makeValuePointer(val));
        return *this;
    }
};


ValuePointer makeValuePointer(Array const val)
{
    return ValuePointer(new Array(val));
}

ValuePointer makeValuePointer(Object const val)
{
    return ValuePointer(new Object(val));
}


std::ostream& operator<<(std::ostream& out, Value const& val)
{
    val.write(out);
    return out;
}


std::string toString(
    Value const& val, size_t const wd = 0, std::string const prefix = "")
{
    std::ostringstream ss;
    ss << prefix;
    val.write(ss, wd, 0, prefix);
    return ss.str();
}


} // namespace json
} // namespace anu_am

#endif // ANU_AM_JSON_HPP
