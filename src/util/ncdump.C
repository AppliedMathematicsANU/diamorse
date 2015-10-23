/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  ncdump.C
 *
 *  Prints the header of a NetCDF 3 file in a CDL-like format.
 *
 *  Olaf Delgado-Friedrichs jan 15
 *
 */

#include <libgen.h>

#include "boost/algorithm/string.hpp"

#include "collections.hpp"
#include "netcdf.hpp"
#include "netcdfIO.hpp"

using namespace anu_am::netcdf;
using namespace anu_am::diamorse;


std::string tname(Tag const type)
{
    switch (type)
    {
    case NC_BYTE  : return "byte";
    case NC_CHAR  : return "char";
    case NC_SHORT : return "short";
    case NC_LONG  : return "int";
    case NC_FLOAT : return "float";
    case NC_DOUBLE: return "double";
    default: assert(false);
    }
}


std::string stripname(std::string const path)
{
    char *s = new char[path.size()+1];
    strncpy(s, path.c_str(), path.size()+1);

    std::string const base = basename(s);
    delete s;

    return base.substr(0, base.rfind("."));
}


std::string formatString(std::string const input)
{
    std::string s(input);
    boost::replace_all(s, "\\", "\\\\");
    boost::replace_all(s, "\"", "\\\"");
    boost::replace_all(s, "\n", "\\n\",\n\t\t\t\"");

    return "\"" + s + "\"";
}


int main(int argc, char* argv[])
{
    char* infile = argv[1];

    if (argc < 2)
    {
        std::cerr << "Usage:" << argv[0] << " INPUT" << std::endl;
        return 1;
    }

    FileBuffer data(infile);
    NCFile<FileBuffer> file(data);

    std::vector<Dimension> const dims  = file.dimensions();
    std::vector<Variable>  const vars  = file.variables();
    Attributes const attrs = file.attributes();
    size_t i;

    std::cout << "netcdf " << stripname(infile) << " {" << std::endl;

    std::cout << "dimensions:" << std::endl;
    for (i = 0; i < dims.size(); ++i)
    {
        Dimension d = dims.at(i);
        std::cout << "\t" << d.name << " = " << d.size
                  << " ;" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "variables:" << std::endl;
    for (i = 0; i < vars.size(); ++i)
    {
        Variable v = vars.at(i);
        std::cout << "\t" << tname(v.type()) << " " << v.name()
                  << "(" << toString(v.dimensionNames())
                  << ") ;" << std::endl;

        Attributes const attrs = v.attributes();
        for (size_t j = 0; j < attrs.size(); ++j)
        {
            Attribute a = attrs.at(j);
            std::cout << "\t\t" << v.name() << ":" << attrs.keyAt(j) << " = "
                      << formatString(a.valuesAsString())
                      << " ;" << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "// global attributes:" << std::endl;
    for (i = 0; i < attrs.size(); ++i)
    {
        Attribute a = attrs.at(i);
        std::cout << "\t\t:" << attrs.keyAt(i) << " = "
                  << formatString(a.valuesAsString())
                  << " ;" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "data:" << std::endl;
    for (i = 0; i < vars.size(); ++i)
    {
        Variable v = vars.at(i);
        std::cout << std::endl << " " << v.name() << " =" << std::endl;
        std::cout << "  " << file.valueAsString(v, 0, 0, 0)
                  << ", " << file.valueAsString(v, 1, 0, 0)
                  << ", " << file.valueAsString(v, 2, 0, 0)
                  << ", " << file.valueAsString(v, 3, 0, 0)
                  << ", " << file.valueAsString(v, 4, 0, 0)
                  << ", ... ;"
                  << std::endl;
    }

    std::cout << "}" << std::endl;
}
