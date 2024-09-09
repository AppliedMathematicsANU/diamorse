PYTHON_PREFIX  = $(shell python -c "import sys; print (sys.prefix)")
PYTHON_VERSION = $(shell python -c "import sys;\
print ('%d.%d' % (sys.version_info.major, sys.version_info.minor))")

PYTHON_INCL = -I ${PYTHON_PREFIX}/include/python${PYTHON_VERSION}
NUMPY_INCL  = -I $(shell python -c "import numpy; print (numpy.get_include())")
EXTRA_INCL  = -I ../lib

ifeq ($(CXX),c++)
	CXXWARN = -Wno-\#warnings
else
	CXXWARN = -Wno-cpp
endif

CXXOPTS     = -std=c++11 -g -O3 $(CXXWARN)
CXXFLAGS    = $(CXXOPTS) $(PYTHON_INCL) $(NUMPY_INCL) $(EXTRA_INCL)
BINDIR      = ../../bin
LDFLAGS     = -L${PYTHON_PREFIX}/lib -lpython${PYTHON_VERSION}

MAKEDEPEND  = $(CXX) $(CXXFLAGS) -MM -MP -MT $(BINDIR)/$*.so -o $*.P $<

.PHONY:	all clean

.SUFFIXES:	.pyx .cpp .so

.pyx.cpp:
	cython --cplus $<

$(BINDIR)/%.so : %.cpp
	@$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) --shared -fPIC $< -o $@ $(LDFLAGS)
	(cd ../../python && ln -nsf ../bin/$*.so .)


SRCS := $(wildcard *.pyx)

all:	$(BINDIR) $(SRCS:%.pyx=$(BINDIR)/%.so)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f *.cpp *.pyc


-include $(SRCS:.pyx=.P)
