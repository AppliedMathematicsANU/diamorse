/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  netcdf.hpp
 *
 *  Low-level code for reading and writing volume data from and to NetCDF files.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#ifndef ANU_AM_NETCDF_HPP
#define ANU_AM_NETCDF_HPP

#include <cassert>
#include <cstddef>
#include <map>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>

#include <memory>

#include "OrderedMap.hpp"

#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x100)


namespace anu_am
{
namespace netcdf
{

typedef enum {
  NC_BYTE      =  1,
  NC_CHAR      =  2,
  NC_SHORT     =  3,
  NC_LONG      =  4,
  NC_FLOAT     =  5,
  NC_DOUBLE    =  6,

  NC_DIMENSION = 10,
  NC_VARIABLE  = 11,
  NC_ATTRIBUTE = 12
} Tag;


template<typename T>
std::string toString(T const val)
{
    std::ostringstream ss;
    ss << val;
    return ss.str();
}


template<>
std::string toString<int8_t>(int8_t const val)
{
    std::ostringstream ss;
    ss << (int) val;
    return ss.str();
}


template<typename T>
std::string toString(std::vector<T> const& vals)
{
    std::ostringstream ss;
    for (size_t i = 0; i < vals.size(); ++i)
    {
        if (i > 0)
            ss << ", ";
        ss << vals.at(i);
    }
    return ss.str();
}


template<>
std::string toString<uint8_t>(std::vector<uint8_t> const& vals)
{
    return std::string(vals.begin(), vals.end());
}


template<Tag NCType>
struct TypeTraits
{
};

template<>
struct TypeTraits<NC_BYTE>
{
    typedef int8_t CType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        return buffer.get(offset);
    }

    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        out.write(val);
    }
};

template<>
struct TypeTraits<NC_CHAR>
{
    typedef uint8_t CType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        return buffer.get(offset);
    }


    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        out.write(val);
    }
};

template<>
struct TypeTraits<NC_SHORT>
{
    typedef int16_t CType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        return
            (buffer.get(offset)   << 8) |
            (buffer.get(offset+1) << 0);
    }

    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        out.write((val >> 8) & 0xff);
        out.write((val >> 0) & 0xff);
    }
};

template<>
struct TypeTraits<NC_LONG>
{
    typedef int32_t CType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        return
            (buffer.get(offset)   << 24) |
            (buffer.get(offset+1) << 16) |
            (buffer.get(offset+2) <<  8) |
            (buffer.get(offset+3) <<  0);
    }

    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        out.write((val >> 24) & 0xff);
        out.write((val >> 16) & 0xff);
        out.write((val >>  8) & 0xff);
        out.write((val >>  0) & 0xff);
    }
};

template<>
struct TypeTraits<NC_FLOAT>
{
    typedef float CType;
    typedef TypeTraits<NC_LONG>::CType IntType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        IntType tmp[] = { TypeTraits<NC_LONG>::get(buffer, offset) };
        CType *fp = reinterpret_cast<CType *>(tmp);

        return *fp;
    }

    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        CType tmp[] = { val };
        IntType *ip = reinterpret_cast<IntType *>(tmp);

        TypeTraits<NC_LONG>::write(out, *ip);
    }
};


template<>
struct TypeTraits<NC_DOUBLE>
{
    typedef double CType;
    typedef TypeTraits<NC_LONG>::CType IntType;

    template<class Buffer> static inline CType get(
        Buffer const& buffer,
        size_t const offset)
    {
        IntType tmp[] = { 0, 0 };

        if (IS_BIG_ENDIAN)
        {
            tmp[0] = TypeTraits<NC_LONG>::get(buffer, offset);
            tmp[1] = TypeTraits<NC_LONG>::get(buffer, offset+4);
        }
        else
        {
            tmp[0] = TypeTraits<NC_LONG>::get(buffer, offset+4);
            tmp[1] = TypeTraits<NC_LONG>::get(buffer, offset);
        }

        CType *fp = reinterpret_cast<CType *>(tmp);

        return *fp;
    }

    template<class Writer> static inline void write(
        Writer &out,
        CType const val)
    {
        CType tmp[] = { val };
        IntType *ip = reinterpret_cast<IntType *>(tmp);

        if (IS_BIG_ENDIAN)
        {
            TypeTraits<NC_LONG>::write(out, ip[0]);
            TypeTraits<NC_LONG>::write(out, ip[1]);
        }
        else
        {
            TypeTraits<NC_LONG>::write(out, ip[1]);
            TypeTraits<NC_LONG>::write(out, ip[0]);
        }
    }
};


template<typename T>
struct InverseTraits
{
};

template<>
struct InverseTraits<int8_t>
{
    static const Tag type = NC_BYTE;
};

template<>
struct InverseTraits<uint8_t>
{
    static const Tag type = NC_CHAR;
};

template<>
struct InverseTraits<int16_t>
{
    static const Tag type = NC_SHORT;
};

template<>
struct InverseTraits<int32_t>
{
    static const Tag type = NC_LONG;
};

template<>
struct InverseTraits<uint32_t>
{
    static const Tag type = NC_LONG;
};

template<>
struct InverseTraits<float>
{
    static const Tag type = NC_FLOAT;
};

template<>
struct InverseTraits<double>
{
    static const Tag type = NC_DOUBLE;
};


struct Dimension
{
    std::string name;
    size_t      size;

    Dimension() {}

    Dimension(std::string const name, size_t const size)
        : name(name),
          size(size)
    {
    }
};


struct AttrImpl
{
    virtual ~AttrImpl() {};

    virtual Tag type() const = 0;

    virtual size_t size() const = 0;

    virtual std::string valuesAsString() const = 0;

    virtual std::vector<int32_t> intValues() const = 0;

    virtual std::vector<double> floatValues() const = 0;
};


template<Tag NCType>
struct TypedAttrImpl : public AttrImpl
{
    typedef typename TypeTraits<NCType>::CType CType;

    std::vector<CType> _values;

    explicit TypedAttrImpl(std::vector<CType> const& values)
        : _values(values)
    {
    }

    Tag type() const
    {
        return NCType;
    }

    size_t size() const
    {
        return _values.size();
    }

    std::vector<CType> values() const
    {
        return _values;
    }

    std::vector<int32_t> intValues() const
    {
        std::vector<int32_t> result;
        for (size_t i = 0; i < _values.size(); ++i)
            result.push_back(_values.at(i));
        return result;
    }

    std::vector<double> floatValues() const
    {
        std::vector<double> result;
        for (size_t i = 0; i < _values.size(); ++i)
            result.push_back(_values.at(i));
        return result;
    }

    std::string valuesAsString() const
    {
        return toString(_values);
    }
};


template<typename T>
std::shared_ptr<AttrImpl> makeImpl(std::vector<T> const& values)
{
    return std::shared_ptr<AttrImpl>(
        new TypedAttrImpl<InverseTraits<T>::type>(values));
}


template<>
std::shared_ptr<AttrImpl> makeImpl<char>(std::vector<char> const& text)
{
    std::vector<uint8_t> values(text.begin(), text.end());
    return makeImpl(values);
}


template<>
std::shared_ptr<AttrImpl> makeImpl<size_t>(std::vector<size_t> const& values)
{
    std::vector<int32_t> const converted(values.begin(), values.end());
    return makeImpl(converted);
}


template<>
std::shared_ptr<AttrImpl> makeImpl<uint32_t>(std::vector<uint32_t> const& values)
{
    std::vector<int32_t> const converted(values.begin(), values.end());
    return makeImpl(converted);
}


template<typename T, size_t N>
std::shared_ptr<AttrImpl> makeImpl(T const(&a)[N])
{
    std::vector<T> values(a, a + N);

    return makeImpl(values);
}


template<typename T>
std::shared_ptr<AttrImpl> makeImpl(T const value)
{
    std::vector<T> values;
    values.push_back(value);

    return makeImpl(values);
}


template<typename T>
std::shared_ptr<AttrImpl> makeImpl(T const v1, T const v2)
{
    std::vector<T> values;
    values.push_back(v1);
    values.push_back(v2);

    return makeImpl(values);
}


template<>
std::shared_ptr<AttrImpl> makeImpl<std::string>(std::string const text)
{
    std::vector<uint8_t> values(text.begin(), text.end());
    return makeImpl(values);
}


class Attribute
{
    std::shared_ptr<AttrImpl> _impl;

public:
    Attribute()
    {
    }

    template<typename T>
    Attribute(std::vector<T> const& values)
    : _impl(makeImpl(values))
    {
    }

    template<typename T, size_t N>
    Attribute(T const(&values)[N])
    : _impl(makeImpl(values))
    {
    }

    template<typename T>
    Attribute(T const value)
    : _impl(makeImpl(value))
    {
    }

    template<typename T>
    Attribute(T const v1, T const v2)
        : _impl(makeImpl(v1, v2))
    {
    }

    Tag type() const
    {
        return _impl->type();
    }

    size_t size() const
    {
        return _impl->size();
    }

    std::vector<int32_t> intValues() const
    {
        return _impl->intValues();
    }

    std::vector<double> floatValues() const
    {
        return _impl->floatValues();
    }

    std::string valuesAsString() const
    {
        return _impl->valuesAsString();
    }
};


typedef anu_am::diamorse::OrderedMap<std::string, Attribute> Attributes;


class Variable
{
private: 
    std::string            _name;
    Tag                    _type;
    std::vector<Dimension> _dimensions;
    Attributes             _attributes;

public:
    Variable() {};

    Variable(
        std::string const name,
        Tag const type,
        std::vector<Dimension> const& dimensions,
        Attributes const& attributes = Attributes()
        )
        : _name(name),
          _type(type),
          _dimensions(dimensions),
          _attributes(attributes)
    {
    }

    Variable(
        std::string const name,
        Tag const type,
        size_t const xdim,
        size_t const ydim,
        Attributes const& attributes = Attributes()
        )
        : _name(name),
          _type(type),
          _attributes(attributes)
    {
        _dimensions.push_back(Dimension(name + "_ydim", ydim));
        _dimensions.push_back(Dimension(name + "_xdim", xdim));
    }

    Variable(
        std::string const name,
        Tag const type,
        size_t const dim,
        Attributes const& attributes = Attributes()
        )
        : _name(name),
          _type(type),
          _attributes(attributes)
    {
        _dimensions.push_back(Dimension(name + "_dim", dim));
    }

    Variable(
        std::string const name,
        Tag const type,
        size_t const xdim,
        size_t const ydim,
        size_t const zdim,
        Attributes const& attributes = Attributes()
        )
        : _name(name),
          _type(type),
          _attributes(attributes)
    {
        _dimensions.push_back(Dimension(name + "_zdim", zdim));
        _dimensions.push_back(Dimension(name + "_ydim", ydim));
        _dimensions.push_back(Dimension(name + "_xdim", xdim));
    }

    std::string name() const
    {
        return _name;
    }

    Tag type() const
    {
        return _type;
    }

    std::vector<size_t> dimensions() const
    {
        std::vector<size_t> result;

        for (size_t i = 0; i < _dimensions.size(); ++i)
            result.push_back(_dimensions.at(i).size);

        return result;
    }

    std::vector<std::string> dimensionNames() const
    {
        std::vector<std::string> result;

        for (size_t i = 0; i < _dimensions.size(); ++i)
            result.push_back(_dimensions.at(i).name);

        return result;
    }

    bool hasAttribute(std::string const name) const
    {
        return _attributes.hasKey(name);
    }

    Attribute attribute(std::string const name) const
    {
        return _attributes(name);
    }

    Attributes attributes() const
    {
        return _attributes;
    }
};


struct Header
{
    size_t nrRecords;
    std::vector<Dimension> dimensions;
    Attributes attributes;
    std::vector<Variable>  variables;
    std::vector<size_t> variableOffsets;

    Header() {};

    Header(size_t const nrRecords,
           std::vector<Dimension> const& dimensions,
           Attributes const& attributes,
           std::vector<Variable>  const& variables,
           std::vector<size_t> const& offsets)
        : nrRecords(nrRecords),
          dimensions(dimensions),
          attributes(attributes),
          variables(variables),
          variableOffsets(offsets)
    {
    }
};


template<class Buffer>
class BufferStream {
private:
    Buffer _buf;
    size_t _pos;

    void advance(size_t const n) {
        size_t const k = _pos + n;
        _pos = k + 3 - (k + 3) % 4;
    }

    template<Tag NCType>
    std::vector<typename TypeTraits<NCType>::CType> read(size_t n) {
        size_t const s = sizeof(typename TypeTraits<NCType>::CType);
        size_t i, p;

        std::vector<typename TypeTraits<NCType>::CType> result;
        for (i = 0, p = _pos; i < n; ++i, p += s)
            result.push_back(TypeTraits<NCType>::get(_buf, p));
        advance(s * n);

        return result;
    }

    template<Tag NCType>
    Attribute readAttr(size_t const size)
    {
        return Attribute(read<NCType>(size));
    }

    int32_t readInteger()
    {
        int32_t const val = TypeTraits<NC_LONG>::get(_buf, _pos);

        advance(sizeof(TypeTraits<NC_LONG>::CType));

        return val;
    }

    int32_t readNonNegative()
    {
        int32_t const val = readInteger();

        assert(val >= 0);

        return val;
    }

    std::string readString()
    {
        int32_t const n = readNonNegative();

        std::stringstream ss;
        for (int32_t i = 0; i < n; ++i)
            ss << TypeTraits<NC_CHAR>::get(_buf, _pos+i);

        advance(n * sizeof(TypeTraits<NC_CHAR>::CType));

        return ss.str();
    }

    std::vector<Dimension> readDimensions() {
        std::vector<Dimension> result;
        int32_t const tag = readInteger();
        int32_t const n = readNonNegative();

        if (tag == NC_DIMENSION)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                std::string const name = readString();
                int32_t const size = readNonNegative();
                result.push_back(Dimension(name, size));
            }
        }
        else
            assert(tag == 0 and n == 0);

        return result;
    }

    Attributes readAttributes() {
        Attributes result;
        int32_t const tag = readInteger();
        int32_t const n = readNonNegative();

        if (tag == NC_ATTRIBUTE)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                std::string const name = readString();
                Tag const type = (Tag) readInteger();
                int32_t const size = readNonNegative();

                switch (type) {
                case NC_BYTE:
                    result.set(name, readAttr<NC_BYTE>(size));
                    break;
                case NC_CHAR:
                    result.set(name, readAttr<NC_CHAR>(size));
                    break;
                case NC_SHORT:
                    result.set(name, readAttr<NC_SHORT>(size));
                    break;
                case NC_LONG:
                    result.set(name, readAttr<NC_LONG>(size));
                    break;
                case NC_FLOAT:
                    result.set(name, readAttr<NC_FLOAT>(size));
                    break;
                case NC_DOUBLE:
                    result.set(name, readAttr<NC_DOUBLE>(size));
                    break;
                default:
                    assert(false);
                }
            }
        }
        else
            assert(tag == 0 and n == 0);

        return result;
    }

    std::pair<std::vector<Variable>, std::vector<size_t> >
    readVariables(std::vector<Dimension> const dimensions)
    {
        std::vector<Variable> vars;
        std::vector<size_t> offsets;
        int32_t const tag = readInteger();
        int32_t const n = readNonNegative();

        if (tag == NC_VARIABLE)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                std::string const name = readString();

                int32_t const ndims = readNonNegative();
                std::vector<Dimension> dims;
                for (int32_t k = 0; k < ndims; ++k)
                    dims.push_back(dimensions.at(readNonNegative()));

                Attributes const attr = readAttributes();
                Tag const type = (Tag) readInteger();
                /*int32_t const size = */ readNonNegative();
                int32_t const start = readNonNegative();
            
                vars.push_back(Variable(name, type, dims, attr));
                offsets.push_back(start);
            }
        }
        else
            assert(tag == 0 and n == 0);

        return std::make_pair(vars, offsets);
    }

public:
    explicit BufferStream(Buffer const& buf)
    : _buf(buf),
      _pos(0)
    {
    }

    Header readHeader() {
        std::vector<uint8_t> const magic = read<NC_CHAR>(4);
        assert(std::string(magic.begin(), magic.end()) == "CDF\x1");

        size_t const nrRecords = readNonNegative();
        std::vector<Dimension> const dimensions = readDimensions();
        Attributes const attributes = readAttributes();
        std::pair<std::vector<Variable>, std::vector<size_t> >
            const variables = readVariables(dimensions);

        return Header(nrRecords, dimensions, attributes,
                      variables.first, variables.second);
    }
};


class NCFileInfo
{
private:
    std::vector<Dimension>           _dimensions;
    std::map<std::string, Dimension> _dimensionsByName;
    std::vector<Variable>            _variables;
    std::map<std::string, Variable>  _variablesByName;
    Attributes                       _attributes;

public:
    NCFileInfo() {}

    NCFileInfo(
        std::vector<Dimension> const& dimensions,
        std::vector<Variable>  const& variables,
        Attributes             const& attributes = Attributes()
        )
        : _attributes(attributes)
    {
        for (size_t i = 0; i < dimensions.size(); ++i)
        {
            Dimension const d = dimensions.at(i);
            _dimensions.push_back(d);
            _dimensionsByName[d.name] = d;
        }

        for (size_t i = 0; i < variables.size(); ++i)
        {
            Variable const v = variables.at(i);
            _variables.push_back(v);
            _variablesByName[v.name()] = v;
        }
    }

    explicit NCFileInfo(
        std::vector<Variable>  const& variables,
        Attributes             const& attributes = Attributes()
        )
        : _attributes(attributes)
    {
        for (size_t i = 0; i < variables.size(); ++i)
        {
            Variable const v = variables.at(i);
            _variables.push_back(v);
            _variablesByName[v.name()] = v;

            for (size_t j = 0; j < v.dimensions().size(); ++j)
            {
                std::string const name = v.dimensionNames().at(j);
                size_t const size = v.dimensions().at(j);
                Dimension const d(name, size);

                _dimensions.push_back(d);
                _dimensionsByName[d.name] = d;
            }
        }
    }

    Dimension dimension(std::string const name) const
    {
        return _dimensionsByName.at(name);
    }

    std::vector<Dimension> dimensions() const
    {
        return _dimensions;
    }

    Variable variable(std::string const name) const
    {
        return _variablesByName.at(name);
    }

    std::vector<Variable> variables() const
    {
        return _variables;
    }

    bool hasAttribute(std::string const name) const
    {
        return _attributes.hasKey(name);
    }

    Attribute attribute(std::string const name) const
    {
        return _attributes(name);
    }

    Attributes attributes() const
    {
        return _attributes;
    }
};


struct AccImpl
{
    virtual ~AccImpl() {};

    virtual int32_t getInt(
        size_t const x, size_t const y, size_t const z) const = 0;

    virtual double getFloat(
        size_t const x, size_t const y, size_t const z) const = 0;

    virtual std::string valueAsString(
        size_t const x, size_t const y, size_t const z) const = 0;
};


template<Tag NCType, class Buffer>
struct BufferAccImpl : public AccImpl
{
    typedef typename TypeTraits<NCType>::CType CType;

private:
    size_t _start;
    size_t _xdim, _ydim, _zdim;
    Buffer _buf;

    inline bool inside(size_t const x, size_t const y, size_t const z) const
    {
        return x < _xdim and y < _ydim and z < _zdim;
    }

public:
    BufferAccImpl(std::vector<size_t> const & dimensions,
                 Buffer const& buffer,
                 size_t const start)
        : _start(start),
          _buf(buffer)
    {
        size_t const n = dimensions.size();
        _xdim = dimensions.at(n-1);
        _ydim = n >= 2 ? dimensions.at(n-2) : 1;

        _zdim = 1;
        if (n >= 3)
            for (size_t i = 0; i <= n-3; ++i)
                _zdim *= dimensions.at(i);
    }

    inline CType get(size_t const x, size_t const y, size_t const z) const
    {
        if (inside(x, y, z))
        {
            size_t offset = ((z * _ydim + y) * _xdim + x) * sizeof(CType);
            return TypeTraits<NCType>::get(_buf, _start + offset);
        }
        else
            return 0;
    }

    int32_t getInt(size_t const x, size_t const y, size_t const z) const
    {
        return get(x, y, z);
    }

    double getFloat(size_t const x, size_t const y, size_t const z) const
    {
        return get(x, y, z);
    }

    std::string valueAsString(
        size_t const x, size_t const y, size_t const z) const
    {
        return toString(get(x, y, z));
    }
};


class Accessor
{
public:
    typedef std::shared_ptr<AccImpl> ImplPtr;

private:
    ImplPtr _impl;

public:
    Accessor() {}

    Accessor(ImplPtr const impl)
    : _impl(impl)
    {
    }

    int32_t getInt(size_t const x, size_t const y, size_t const z) const
    {
        return _impl->getInt(x, y, z);
    }

    double getFloat(size_t const x, size_t const y, size_t const z) const
    {
        return _impl->getFloat(x, y, z);
    }

    std::string valueAsString(
        size_t const x, size_t const y, size_t const z) const
    {
        return _impl->valueAsString(x, y, z);
    }
};


template<Tag NCType, class Buffer>
Accessor makeAccessor(
    std::vector<size_t> const& dims,
    Buffer const& buf,
    size_t const start)
{
    AccImpl *acc = new BufferAccImpl<NCType, Buffer>(dims, buf, start);
    return Accessor(std::shared_ptr<AccImpl>(acc));
}


template<class Buffer>
class NCFile
{
    Buffer _buf;
    Header _header;
    NCFileInfo _info;
    std::map<std::string, Accessor> _accessors;

    Accessor makeAcc(Variable const v, size_t const start)
    {
        std::vector<size_t> dims = v.dimensions();

        switch (v.type()) {
        case NC_BYTE:   return makeAccessor<NC_BYTE  >(dims, _buf, start);
        case NC_CHAR:   return makeAccessor<NC_CHAR  >(dims, _buf, start);
        case NC_SHORT:  return makeAccessor<NC_SHORT >(dims, _buf, start);
        case NC_LONG:   return makeAccessor<NC_LONG  >(dims, _buf, start);
        case NC_FLOAT:  return makeAccessor<NC_FLOAT >(dims, _buf, start);
        case NC_DOUBLE: return makeAccessor<NC_DOUBLE>(dims, _buf, start);
        default:
            assert(false);
        }
    }


public:
    explicit NCFile(Buffer const& buffer)
    : _buf(buffer),
      _header(BufferStream<Buffer>(_buf).readHeader()),
      _info(_header.dimensions, _header.variables, _header.attributes)
    {
        for (size_t i = 0; i < _header.variables.size(); ++i)
        {
            std::string const name = _header.variables.at(i).name();
            size_t const start = _header.variableOffsets.at(i);

            _accessors[name] = makeAcc(_info.variable(name), start);
        }
    }

    int32_t getInt(
        Variable const& var,
        size_t const x, size_t const y, size_t const z) const
    {
        return _accessors.at(var.name()).getInt(x, y, z);
    }

    double getFloat(
        Variable const& var,
        size_t const x, size_t const y, size_t const z) const
    {
        return _accessors.at(var.name()).getFloat(x, y, z);
    }

    std::string valueAsString(
        Variable const& var,
        size_t const x, size_t const y, size_t const z) const
    {
        return _accessors.at(var.name()).valueAsString(x, y, z);
    }

    Dimension dimension(std::string const name) const
    {
        return _info.dimension(name);
    }

    std::vector<Dimension> dimensions() const
    {
        return _info.dimensions();
    }

    Variable variable(std::string const name) const
    {
        return _info.variable(name);
    }

    std::vector<Variable> variables() const
    {
        return _info.variables();
    }

    bool hasAttribute(std::string const name) const
    {
        return _info.hasAttribute(name);
    }

    Attribute attribute(std::string const name) const
    {
        return _info.attribute(name);
    }

    Attributes attributes() const
    {
        return _info.attributes();
    }

    NCFileInfo info() const
    {
        return _info;
    }

    std::map<std::string, Accessor> accessors() const
    {
        return _accessors;
    }
};


template<class Writer>
class NCWriteStream
{
    Writer _out;
    size_t _pos;

public:
    explicit NCWriteStream(Writer& out)
    : _out(out),
      _pos(0)
    {
    }


private:
    void advance(size_t const n) {
        size_t const k = _pos + n;
        _pos = k + 3 - (k + 3) % 4;
        for (size_t i = 0; i < _pos - k; ++i)
            TypeTraits<NC_BYTE>::write(_out, 0);
    }


    template<Tag NCType, typename GivenType>
    void writeSingle(GivenType const val)
    {
        TypeTraits<NCType>::write(_out, val);
    }


    template<Tag NCType, typename GivenType>
    void write(std::vector<GivenType> const& values)
    {
        size_t const n = values.size();
        size_t const s = sizeof(typename TypeTraits<NCType>::CType);

        for (size_t i = 0; i < n; ++i)
            writeSingle<NCType>(values.at(i));

        advance(s * n);
    }

public:
    void writeInteger(int32_t const val)
    {
        TypeTraits<NC_LONG>::write(_out, val);
        advance(sizeof(TypeTraits<NC_LONG>::CType));
    }

    void writeNonNegative(int32_t const val)
    {
        assert(val >= 0);
        writeInteger(val);
    }

    void writeString(std::string const& val)
    {
        size_t const n = val.size();
        size_t const s = sizeof(typename TypeTraits<NC_CHAR>::CType);

        writeNonNegative(val.size());

        for (size_t i = 0; i < n; ++i)
            TypeTraits<NC_CHAR>::write(_out, val.at(i));

        advance(s * n);
    }

    void writeDimensions(std::vector<Dimension> const& dims)
    {
        size_t const n = dims.size();

        writeInteger(n > 0 ? NC_DIMENSION : 0);
        writeNonNegative(n);

        for (size_t i = 0; i < n; ++i)
        {
            writeString(dims.at(i).name);
            writeNonNegative(dims.at(i).size);
        }
    }

    void writeAttributes(Attributes const& attrs)
    {
        size_t const n = attrs.size();

        writeInteger(n > 0 ? NC_ATTRIBUTE : 0);
        writeNonNegative(n);

        for (size_t i = 0; i < n; ++i)
        {
            Attribute const& a = attrs.at(i);
            writeString(attrs.keyAt(i));
            writeInteger(a.type());

            if (a.type() == NC_CHAR)
                writeString(a.valuesAsString());
            else
            {
                writeNonNegative(a.size());

                switch (a.type())
                {
                case NC_BYTE  : write<NC_BYTE>  (a.intValues());   break;
                case NC_SHORT : write<NC_SHORT> (a.intValues());   break;
                case NC_LONG  : write<NC_LONG>  (a.intValues());   break;
                case NC_FLOAT : write<NC_FLOAT> (a.floatValues()); break;
                case NC_DOUBLE: write<NC_DOUBLE>(a.floatValues()); break;
                default:
                    assert(false);
                }
            }
        }
    }

private:
    template<Tag NCType>
    size_t elementSize() const
    {
        return sizeof(typename TypeTraits<NCType>::CType);
    }

    size_t elementSize(Tag const type) const
    {
        switch(type)
        {
        case NC_BYTE:   return elementSize<NC_BYTE>();
        case NC_CHAR:   return elementSize<NC_CHAR>();
        case NC_SHORT:  return elementSize<NC_SHORT>();
        case NC_LONG:   return elementSize<NC_LONG>();
        case NC_FLOAT:  return elementSize<NC_FLOAT>();
        case NC_DOUBLE: return elementSize<NC_DOUBLE>();
        default:
            assert(false);
        }
    }

    template<class T>
    void writeVariables(NCWriteStream<T>& out,
                        std::vector<Variable> const& vars,
                        std::vector<Dimension> const& dimensions,
                        size_t const vstart)
    {
        size_t const n = vars.size();
        size_t start = vstart;

        out.writeInteger(n > 0 ? NC_VARIABLE : 0);
        out.writeNonNegative(n);

        for (size_t i = 0; i < n; ++i)
        {
            Variable const v = vars.at(i);
            out.writeString(v.name());

            std::vector<std::string> const dims = v.dimensionNames();
            size_t const m = dims.size();
            out.writeNonNegative(m);

            size_t size = elementSize(v.type());
            for (size_t j = 0; j < m; ++j)
            {
                std::string const s = dims.at(j);
                size_t k;
                for (k = 0; k < dimensions.size(); ++k)
                    if (dimensions.at(k).name == s)
                        break;
                out.writeNonNegative(k);
                size *= dimensions.at(k).size;
            }

            out.writeAttributes(v.attributes());

            out.writeInteger(v.type());
            out.writeNonNegative(size);
            out.writeNonNegative(start);

            start += size + 3 - (size + 3) % 4;
        }
    }

    struct BufferWriter
    {
        std::vector<uint8_t> &data;

        explicit BufferWriter(std::vector<uint8_t> &data)
            : data(data)
        {
        }

        void write(uint8_t const val)
        {
            data.push_back(val);
        }
    };

public:
    void writeVariables(std::vector<Variable> const& vars,
                        std::vector<Dimension> const& dimensions)
    {
        std::vector<uint8_t> data;
        BufferWriter buf(data);
        NCWriteStream<BufferWriter> out(buf);

        writeVariables(out, vars, dimensions, 0);
        size_t const vstart = _pos + data.size();
        writeVariables(*this, vars, dimensions, vstart);
        assert(_pos == vstart);
    }

    template<class Accessor>
    void writeNCFile(NCFileInfo const& info,
                     std::map<std::string, Accessor> const& accessors)
    {
        std::string const magic = "CDF\x1";
        write<NC_CHAR>(std::vector<uint8_t>(magic.begin(), magic.end()));
        writeNonNegative(0);
        writeDimensions(info.dimensions());
        writeAttributes(info.attributes());
        writeVariables(info.variables(), info.dimensions());

        for (size_t i = 0; i < info.variables().size(); ++i)
        {
            Variable const v = info.variables().at(i);
            std::vector<size_t> const dims = v.dimensions();
            Accessor const acc = accessors.at(v.name());

            size_t xdim, ydim, zdim;
            switch (dims.size())
            {
            case 0:
                xdim = ydim = zdim = 1; break;
            case 1:
                xdim = dims.at(0); ydim = zdim = 1; break;
            case 2:
                xdim = dims.at(1); ydim = dims.at(0); zdim = 1; break;
            default:
                xdim = dims.at(2); ydim = dims.at(1); zdim = dims.at(0); break;
            }

            for (size_t z = 0; z < zdim; ++z)
                for (size_t y = 0; y < ydim; ++y)
                    for (size_t x = 0; x < xdim; ++x)
                        switch (v.type())
                        {
                        case NC_BYTE:
                            writeSingle<NC_BYTE>(acc.getInt(x, y, z));
                            break;
                        case NC_CHAR:
                            writeSingle<NC_CHAR>(acc.getInt(x, y, z));
                            break;
                        case NC_SHORT:
                            writeSingle<NC_SHORT>(acc.getInt(x, y, z));
                            break;
                        case NC_LONG:
                            writeSingle<NC_LONG>(acc.getInt(x, y, z));
                            break;
                        case NC_FLOAT:
                            writeSingle<NC_FLOAT>(acc.getFloat(x, y, z));
                            break;
                        case NC_DOUBLE:
                            writeSingle<NC_DOUBLE>(acc.getFloat(x, y, z));
                            break;
                        default:
                            assert(false);
                        }

            advance(elementSize(v.type()) * xdim * ydim * zdim);
        }
    }
};


template<class Writer, class Accessor>
void writeNCFile(Writer& writer,
                 NCFileInfo const& info,
                 std::map<std::string, Accessor> const& accessors)
{
    NCWriteStream<Writer>(writer).writeNCFile(info, accessors);
}


template<class Writer, class Buffer, class Accessor>
void writeNCFile(Writer& writer,
                 NCFile<Buffer> const& file,
                 std::map<std::string, Accessor> const& accessors)
{
    writeNCFile(writer, file.info(), accessors);
}


template<class Writer, class Accessor>
void writeNCFile(Writer& writer,
                 std::vector<Variable> const& variables,
                 Attributes const& attributes,
                 std::map<std::string, Accessor> const& accessors)
{
    writeNCFile(writer, NCFileInfo(variables, attributes), accessors);
}


template<class Writer, class Accessor>
void writeNCFile(Writer& writer,
                 std::vector<Variable> const& variables,
                 std::map<std::string, Accessor> const& accessors)
{
    writeNCFile(writer, NCFileInfo(variables), accessors);
}


} // namespace netcdf
} // namespace anu_am

#endif // ANU_AM_NETCDF_HPP
