PYTHON_PREFIX  = $(shell python -c "import sys; print sys.prefix")
PYTHON_VERSION = $(shell python -c "import sys;\
print '%d.%d' % (sys.version_info.major, sys.version_info.minor)")

PYTHON_INCL = -I ${PYTHON_PREFIX}/include/python${PYTHON_VERSION}
NUMPY_INCL  = -I $(shell python -c "import numpy; print numpy.get_include()")
EXTRA_INCL  = -I ../lib
CXXOPTS     = -std=c++11 -g -O3 -Wno-cpp
CXXFLAGS    = $(CXXOPTS) $(PYTHON_INCL) $(NUMPY_INCL) $(EXTRA_INCL)
BINDIR      = ../../bin
LDFLAGS     = -lpython${PYTHON_VERSION}

MAKEDEPEND  = $(CXX) $(CXXFLAGS) -MM -MP -MT $(BINDIR)/$*.so -o $*.P $<

.PHONY:	all clean

.SUFFIXES:	.pyx .cpp .so

.pyx.cpp:
	cython --cplus $<

$(BINDIR)/%.so : %.cpp
	@$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) --shared -fPIC $< -o $@ $(LDFLAGS)


SRCS := $(wildcard *.pyx)

all:	$(BINDIR) $(SRCS:%.pyx=$(BINDIR)/%.so)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f *.cpp *.pyc


-include $(SRCS:.pyx=.P)
