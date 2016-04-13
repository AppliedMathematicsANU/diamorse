UNAME := $(shell uname)

GIT_REV  = $(shell git rev-parse HEAD || "unknown")
GIT_DATE = $(shell git show -s --format=%ci $(GIT_REV) || "unknown")
CXXDEFS  = -DGIT_REVISION="\"$(GIT_REV)\"" -DGIT_TIMESTAMP="\"$(GIT_DATE)\""

CXXSTD   = -std=c++11
CXXWARN  = -Wall -Wextra -pedantic
CXXDEBUG = -g
CXXOPTIM = -O3
CXXLIBS  = -I ../lib 
CXXFLAGS = $(CXXSTD) $(CXXDEFS) $(CXXWARN) $(CXXDEBUG) $(CXXOPTIM) $(CXXLIBS)
BINDIR   = ../../bin

MAKEDEPEND= $(CXX) $(CXXFLAGS) -MM -MP -MT $(BINDIR)/$* -o $*.P $<

.PHONY:	all clean

.SUFFIXES:	.C

$(BINDIR)/% : %.C
	@$(MAKEDEPEND);
	$(CXX) $(CXXFLAGS) $< -o $@


SRCS := $(wildcard *.C)

all:	$(BINDIR) $(SRCS:%.C=$(BINDIR)/%)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f *.P

-include $(SRCS:.C=.P)
