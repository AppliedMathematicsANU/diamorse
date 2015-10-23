/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  volume_io.hpp
 *
 *  Code for reading and writing volume images.
 *
 *  Olaf Delgado-Friedrichs jun 15
 *
 */

#ifndef ANU_AM_DIAMORSE_VOLUMEIO_HPP
#define ANU_AM_DIAMORSE_VOLUMEIO_HPP

#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <boost/filesystem.hpp>

#include "netcdf.hpp"
#include "netcdfIO.hpp"
#include "stringUtils.hpp"


#ifndef GIT_REVISION
#define GIT_REVISION "unknown"
#endif

#ifndef GIT_TIMESTAMP
#define GIT_TIMESTAMP "unknown"
#endif


#ifndef __FILE__
#defined __FILE__ "unknown"
#endif


namespace anu_am
{
namespace diamorse
{


using namespace anu_am::netcdf;
namespace su = anu_am::stringutils;

template<typename T>
struct InverseTraits
{
};

template<>
struct InverseTraits<int8_t>
{
    static const Tag type = NC_BYTE;
    static int8_t mask() { return -1; }
};

template<>
struct InverseTraits<uint8_t>
{
    static const Tag type = NC_CHAR;
    static uint8_t mask() { return 0xff; }
};

template<>
struct InverseTraits<int16_t>
{
    static const Tag type = NC_SHORT;
    static int16_t mask() { return -1; }
};

template<>
struct InverseTraits<int32_t>
{
    static const Tag type = NC_LONG;
    static int32_t mask() { return 0x7fffffff; }
};

template<>
struct InverseTraits<uint32_t>
{
    static const Tag type = NC_LONG;
    static uint32_t mask() { return 0x7fffffff; }
};

template<>
struct InverseTraits<float>
{
    static const Tag type = NC_FLOAT;
    static float mask() { return 1.0e30f; }
};

template<>
struct InverseTraits<double>
{
    static const Tag type = NC_DOUBLE;
    static float mask() { return 1.0e30; }
};


std::vector<std::string> entries(std::string const filepath)
{
    using namespace boost::filesystem;

    path const p(filepath);
    std::vector<std::string> result;

    if (is_directory(p))
    {
        directory_iterator iter;
        for (iter = directory_iterator(p); iter != directory_iterator(); ++iter)
        {
            path const entry = iter->path();
            std::string const name = entry.string();
            if (name.rfind(".nc") == name.size() - 3)
                result.push_back(name);
        }
    }
    else
        result.push_back(filepath);

    return result;
}


NCFileInfo readFileInfo(std::string const path)
{
    FileBuffer const data(entries(path).at(0));
    NCFile<FileBuffer> const file(data);
    return file.info();
}


Variable findVolumeVariable(NCFileInfo const info, std::string const name = "")
{
    std::vector<Variable> const vars = info.variables();

    for (size_t i = 0; i < vars.size(); ++i)
    {
        Variable const v = vars.at(i);

        if (v.dimensions().size() == 3
            and (name.empty() or name.compare(v.name()) == 0))
        {
            return v;
        }
    }

    throw std::runtime_error("No appropriate variable found");
}


Variable findVolumeVariable(std::string const path, std::string const name = "")
{
    return findVolumeVariable(readFileInfo(path), name);
}


std::vector<size_t> readDimensions(NCFileInfo const info)
{
    std::vector<size_t> dims = findVolumeVariable(info).dimensions();
    std::reverse(dims.begin(), dims.end());

    if (info.hasAttribute("zdim_total"))
    {
        std::vector<int32_t> const zdim = info.attribute("zdim_total").intValues();
        dims.at(2) = zdim.at(0);
    }

    return dims;
}


std::vector<size_t> readDimensions(std::string const path)
{
    return readDimensions(readFileInfo(path));
}


template<typename T>
boost::shared_ptr<std::vector<T> > readVolumeData(std::string const path)
{
    typedef boost::shared_ptr<std::vector<T> > DataPtr;

    Tag const type = InverseTraits<T>::type;

    size_t xdim = 0, ydim = 0, zdim = 0;

    DataPtr result;

    std::vector<std::string> const paths = entries(path);

    for (size_t i = 0; i < paths.size(); ++i)
    {
        FileBuffer const data(paths.at(i));
        NCFile<FileBuffer> const file(data);
        Variable const v = findVolumeVariable(file.info());

        if (i == 0)
        {
            std::vector<size_t> const dims = readDimensions(file.info());
            xdim = dims.at(0);
            ydim = dims.at(1);
            zdim = dims.at(2);

            result = DataPtr(new std::vector<T>(xdim * ydim * zdim));
        }

        size_t zmin = 0;
        size_t zlen = zdim;

        if (file.hasAttribute("zdim_range"))
        {
            std::vector<int32_t> const r = file.attribute("zdim_range").intValues();
            zmin = r.at(0);
            zlen = r.at(1) - r.at(0) + 1;
        }

        size_t k = zmin * xdim * ydim;

        for (size_t z = 0; z < zlen; ++z)
        {
            for (size_t y = 0; y < ydim; ++y)
            {
                for (size_t x = 0; x < xdim; ++x)
                {
                    if (type == NC_FLOAT or type == NC_DOUBLE)
                        result->at(k) = file.getFloat(v, x, y, z);
                    else
                        result->at(k) = file.getInt(v, x, y, z);
                    ++k;
                }
            }
        }
    }

    return result;
}


template<typename T>
class VectorAccessorImpl : public AccImpl
{
public:
    typedef boost::shared_ptr<std::vector<T> > DataPtr;

private:
    DataPtr _data;
    size_t _xdim;
    size_t _ydim;
    size_t _zdim;
    size_t _xoff;
    size_t _yoff;
    size_t _zoff;
    Tag _type;

    inline bool inside(size_t const x, size_t const y, size_t const z) const
    {
        return x < _xdim and y < _ydim and z < _zdim;
    }

    inline T get(size_t const x, size_t const y, size_t const z) const
    {
        if (inside(x+_xoff, y+_yoff, z+_zoff))
            return _data->at(((z+_zoff) * _ydim + y+_yoff) * _xdim + x+_xoff);
        else
            return 0;
    }

public:
    VectorAccessorImpl() {}

    VectorAccessorImpl(
        DataPtr const data,
        size_t const xdim,
        size_t const ydim,
        size_t const zdim,
        size_t const xoff = 0,
        size_t const yoff = 0,
        size_t const zoff = 0
        )
        : _data(data),
          _xdim(xdim),
          _ydim(ydim),
          _zdim(zdim),
          _xoff(xoff),
          _yoff(yoff),
          _zoff(zoff)
    {
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


template<typename T>
Accessor makeVectorAccessor(
    boost::shared_ptr<std::vector<T> > const data,
    size_t const xdim,
    size_t const ydim,
    size_t const zdim,
    size_t const xoff = 0,
    size_t const yoff = 0,
    size_t const zoff = 0
    )
{
    AccImpl *acc =
        new VectorAccessorImpl<T>(data, xdim, ydim, zdim, xoff, yoff, zoff);

    return Accessor(boost::shared_ptr<AccImpl>(acc));
}


std::string stripPath(std::string const path)
{
    std::string t = path;
    while (t.size() > 3)
    {
        if (t.rfind(".nc") == t.size() - 3 || t.rfind("_nc") == t.size() - 3)
            t = t.substr(0, t.size() - 3);
        else
            break;
    }
    return t;
}


std::string basename(std::string const path)
{
    std::string t = stripPath(path);
    size_t const pos = t.rfind("/");
    if (pos == std::string::npos)
        return t;
    else
        return t.substr(pos+1);
}


std::string chunkPath(std::string const basePath, int const blockNr)
{
    std::ostringstream ss;
    ss << basePath << std::setw(8) << std::setfill('0') << blockNr << ".nc";
    return ss.str();
}


std::string timestamp()
{
    time_t t = time(0); // current time in seconds from epoch
    struct tm* now = gmtime(&t);

    std::ostringstream ss;
    ss << std::setfill('0')
       << now->tm_year + 1900
       << std::setw(2) << now->tm_mon + 1
       << std::setw(2) << now->tm_mday
       << "_"
       << std::setw(2) << now->tm_hour
       << std::setw(2) << now->tm_min
       << std::setw(2) << now->tm_sec;
    return ss.str();
}


bool isTimeStampChar(char const c)
{
    return std::string("0123456789_").find(c) != std::string::npos;
}


std::string stripTimestamp(std::string const str)
{
    size_t i;

    for (i = 0; i < str.size() and isTimeStampChar(str.at(i)); ++i)
        ;

    return str.substr(i);
}


std::string guessDatasetID(
    std::string const path,
    Attributes const& attrs = Attributes())
{
    if (attrs.hasKey("dataset_id"))
        return attrs("dataset_id").valuesAsString();
    else
    {
        std::string best = "";

        for (size_t i = 0; i < attrs.size(); ++i)
        {
            std::string const key = attrs.keyAt(i);
            if (su::startsWith(key, "history_")
                and not su::endsWith(key, "_output"))
            {
                std::string const id =
                    su::stripLeading(su::stripLeading(key, "history_"), "UTC_");

                if (id > best)
                    best = id;
            }
        }

        if (best.size() > 0)
            return best;
    }

    return basename(path);
}


Attributes inheritableAttributes(Attributes const& attrs)
{
    Attributes result;

    for (size_t i = 0; i < attrs.size(); ++i)
    {
        std::string const key = attrs.keyAt(i);
        if (su::startsWith(key, "history_") or
            key == "data_description" or
            key == "voxel_size_xyz" or
            key == "voxel_unit" or
            key == "coord_transform" or
            key == "total_grid_size_xyz" or
            key == "coordinate_origin_xyz")
        {
            result.set(key, attrs(key));
        }
        else if (
            key != "data_min_max" and
            key != "data_histogram_offset" and
            key != "data_histogram_binsize")
        {
            result.set(key, "");
        }
    }

    return result;
}


template<typename T, size_t N>
size_t size(T const(&)[N])
{
    return N;
}


std::string const prefixes[] = {
    "vector_field",
    "tomo_float",
    "tomo",
    "segmented",
    "proju",
    "proj",
    "pore_throat_network",
    "medial_axis",
    "labels",
    "grain_contact_network",
    "distance_map",
    "cntr_tomo"
};


std::string derivedID(
    std::string const parent,
    std::string const prefix,
    std::string const suffix)
{
    std::string base = stripTimestamp(parent);

    for (size_t i = 0; i < size(prefixes); ++i)
    {
        std::string const p = prefixes[i];

        if (su::startsWith(base, p))
        {
            base = su::stripLeading(base, p);
            break;
        }
    }

    return timestamp() + "_" + prefix + base + "_" + suffix;
}


struct Histogram
{
    double min;
    double max;
    double offset;
    double binsize;

    boost::shared_ptr<std::vector<double> > data;

    template<typename T>
    Histogram(boost::shared_ptr<std::vector<T> > const inputs,
              size_t const numberOfBins = 0x10000)
        : data(new std::vector<double>(numberOfBins))
    {
        Tag const type = InverseTraits<T>::type;
        T const mask = InverseTraits<T>::mask();

        min = max = 0;

        for (size_t i = 0; i < inputs->size(); ++i)
        {
            T const value = inputs->at(i);

            if (value != mask)
            {
                if (i == 0) min = max = value;
                if (value < min) min = value;
                if (value > max) max = value;
            }
        }

        offset = min;
        binsize = (max - min) * (1.0 + 1.0e-12) / numberOfBins;

        if (type != NC_FLOAT and type != NC_DOUBLE)
            binsize = ceil(binsize);

        for (size_t i = 0; i < inputs->size(); ++i)
        {
            T const value = inputs->at(i);
            if (value != mask)
                ++data->at((size_t) ((value - offset) / binsize));
        }
    }
};


class VolumeWriteOptions
{
    Attributes  _fileAttributes;
    Attributes  _variableAttributes;
    std::string _datasetID;
    std::string _description;
    bool        _computeHistogram;
    size_t      _fileSizeLimit;

public:
    VolumeWriteOptions()
        : _computeHistogram(true),
          _fileSizeLimit(1024 * 1024 * 1024)
    {
    }

    Attributes fileAttributes() const
    {
        return _fileAttributes;
    }

    VolumeWriteOptions& fileAttributes(Attributes const& attr)
    {
        _fileAttributes = attr;
        return *this;
    }

    Attributes variableAttributes() const
    {
        return _variableAttributes;
    }

    VolumeWriteOptions& variableAttributes(Attributes const& attr)
    {
        _variableAttributes = attr;
        return *this;
    }

    std::string datasetID() const
    {
        return _datasetID;
    }

    VolumeWriteOptions& datasetID(std::string const& text)
    {
        _datasetID = text;
        return *this;
    }

    std::string description() const
    {
        return _description;
    }

    VolumeWriteOptions& description(std::string const& text)
    {
        _description = text;
        return *this;
    }

    bool computeHistogram() const
    {
        return _computeHistogram;
    }

    VolumeWriteOptions& computeHistogram(bool const p)
    {
        _computeHistogram = p;
        return *this;
    }

    size_t fileSizeLimit() const
    {
        return _fileSizeLimit;
    }

    VolumeWriteOptions& fileSizeLimit(size_t const n)
    {
        _fileSizeLimit = n;
        return *this;
    }
};


template<typename T>
void writeVolumeData(
    boost::shared_ptr<std::vector<T> > const data,
    std::string const path,
    std::string const varname,
    size_t const xdim, size_t const ydim, size_t const zdim,
    VolumeWriteOptions const& options = VolumeWriteOptions())
{
    size_t const chunkSize =
        std::max((size_t) 1,
                 options.fileSizeLimit() / (sizeof(T) * xdim * ydim));
    size_t const nrChunks =
        std::max((size_t) 1, (chunkSize +zdim - 1) / chunkSize);

    Tag const type = InverseTraits<T>::type;

    std::string const base = stripPath(path) + (nrChunks == 1 ? ".nc" : "_nc");
    std::string dataset_id = options.datasetID();
    if (dataset_id.size() == 0)
        dataset_id = timestamp() + "_" + basename(path);

    Attributes attr = options.fileAttributes()
        ("zdim_total", zdim)
        ("number_of_files", nrChunks)
        ("dataset_id", dataset_id)
        ("history_"+dataset_id, options.description());

    if (nrChunks > 1)
        boost::filesystem::create_directory(base.c_str());

    for (size_t i = 0; i < nrChunks; ++i)
    {
        size_t const z = i * chunkSize;
        size_t const zlen = std::min(chunkSize, zdim - z);

        std::vector<Variable> vars;
        vars.push_back(Variable(varname, type, xdim, ydim, zlen,
                                options.variableAttributes()));

        std::map<std::string, Accessor> accessors;
        accessors[varname] =
            makeVectorAccessor(data, xdim, ydim, zdim, 0, 0, z);

        size_t const zdim_range[] = { z, z + zlen - 1 };
        attr.set("zdim_range", zdim_range);

        if (i == 0 and options.computeHistogram())
        {
            size_t const size = 0x10000;
            Histogram const h(data, size);
            vars.push_back(Variable("data_histogram", NC_DOUBLE, size));
            accessors["data_histogram"] =
                makeVectorAccessor(h.data, size, 1, 1, 0, 0, 0);
            attr.set("data_min_max", Attribute(h.min, h.max));
            attr.set("data_histogram_binsize", h.binsize);
            attr.set("data_histogram_offset", h.offset);
        }

        std::string const filename =
            nrChunks == 1 ? base : chunkPath(base + "/block", i);

        FileWriter writer(filename);
        writeNCFile(writer, vars, attr, accessors);
    }
}



} // namespace diamorse
} // namespace anu_am

#endif //!ANU_AM_DIAMORSE_VOLUMEIO_HPP
