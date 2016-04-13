UNAME := $(shell uname)

PYTHON_INCL = -I /usr/include/python2.7
NUMPY_INCL  = -I $(shell python -c "import numpy; print numpy.get_include()")
EXTRA_INCL  = -I ../lib
CXXOPTS     = -std=c++11 -g -O3 -Wno-cpp
CXXFLAGS    = $(CXXOPTS) $(PYTHON_INCL) $(NUMPY_INCL) $(EXTRA_INCL)
BINDIR      = ../../bin
LDFLAGS     = -lpython2.7

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
