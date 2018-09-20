/** -*-c++-*-
 *
 *  Copyright 2015 The Australian National University
 *
 *  labyrinth.C
 *
 *  Makes binary images based on nodal approximations of minimal surfaces
 *
 *  Olaf Delgado-Friedrichs mar 15
 *
 */

#include <cmath>

#include "netcdf.hpp"
#include "netcdfIO.hpp"

using namespace anu_am::netcdf;


double const PI = 3.141592653589793;


class NodalAccessor
{
    char _type;
    float _scale;

public:
    NodalAccessor() {};

    NodalAccessor(char const type, float const scale)
        : _type(type),
          _scale(scale)
    {
    }

    int32_t getInt(size_t const x, size_t const y, size_t const z) const
    {
        double const f = 2 * PI / _scale;
        double const cpx = cos(x * f);
        double const cpy = cos(y * f);
        double const cpz = cos(z * f);
        double const spx = sin(x * f);
        double const spy = sin(y * f);
        double const spz = sin(z * f);

        switch (_type)
        {
        case 'G':
        case 'g':
            return (spy * cpz + spz * cpx + spx * cpy) < 0;
        case 'D':
        case 'd':
            return (cpx * cpy * cpz +
                    cpx * spy * spz +
                    spx * cpy * spz +
                    spx * spy * cpz) < 0;
        case 'P':
        case 'p':
            return (cpx + cpy + cpz) < 0;
        default:
            throw new std::runtime_error("unknown type");
        }
    }

    double getFloat(size_t const x, size_t const y, size_t const z) const
    {
        return getInt(x, y, z);
    }
};


int run(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " P|D|G SCALE SIZE"
                  << std::endl;
        return 1;
    }

    int i = 0;
    char const type = argv[++i][0];
    float const scale = atof(argv[++i]);
    int const size = atoi(argv[++i]);

    std::string const varname = "segmented";

    std::vector<Variable> vars;
    vars.push_back(Variable(varname, NC_BYTE, size, size, size));

    std::map<std::string, NodalAccessor> acc;
    acc[varname] = NodalAccessor(type, scale);

    Writer writer;
    writeNCFile(writer, vars, acc);

    return 0;
}


int main(const int argc, char* argv[])
{
  try
  {
    run(argc, argv);
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
