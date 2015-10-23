/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  netcdfIO.hpp
 *
 *  IO classes to be used in conjection with the generic NetCDF library.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#ifndef ANU_AM_NETCDF_IO_HPP
#define ANU_AM_NETCDF_IO_HPP

#include <iostream>
#include <fstream>
#include <stdexcept>

#include "boost/smart_ptr.hpp"


namespace anu_am
{
namespace netcdf
{

class FileBuffer
{
    typedef std::vector<uint8_t> Data;

    boost::shared_ptr<Data> _data;

    Data* contents(std::string const path)
    {
        std::ifstream file(path.c_str(),
                           std::ios::in|std::ios::binary|std::ios::ate);

        if (file.is_open())
        {
            std::streampos const size = file.tellg();
            Data* result = new Data(size);

            file.seekg(0, std::ios::beg);
            file.read((char*) &result->at(0), size);
            file.close();

            return result;
        }
        else
            throw std::runtime_error("Could not open file "+path);
    }

public:
    explicit FileBuffer(std::string const path)
        : _data(contents(path))
    {
    }

    uint8_t get(size_t const offset) const
    {
        return _data->at(offset);
    }
};


struct Writer
{
    Writer() {}

    void write(uint8_t const val)
    {
        std::cout << (char) val;
    }
};


class FileWriter
{
    boost::shared_ptr<std::ofstream> _file;

public:
    FileWriter() {}

    explicit FileWriter(std::string const path)
        : _file(new std::ofstream(path.c_str(), std::ios::binary))
    {
    }

    void write(uint8_t const val)
    {
        std::vector<uint8_t> buf;
        buf.push_back(val);
        _file->write((char *) &buf.at(0), 1);
    }
};


} // namespace netcdf
} // namespace anu_am

#endif // ANU_AM_NETCDF_IO_HPP
