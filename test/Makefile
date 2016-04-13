UNAME := $(shell uname)

CXXSTD   = -std=c++11
CXXWARN  = -Wall -Wextra -pedantic
CXXOPTIM = -O3
CXXDEBUG = -g
CXXLIBS  = -I ../src/lib
CXXFLAGS = $(CXXSTD) $(CXXWARN) $(CXXDEBUG) $(CXXOPTIM) $(CXXLIBS)
BINDIR   = bin

ifeq ($(UNAME), Darwin)
LDFLAGS  = -lboost_unit_test_framework-mt
else
LDFLAGS  = -lboost_unit_test_framework
endif

MAKEDEPEND= $(CXX) $(CXXFLAGS) -MM -MP -MT $(BINDIR)/$* -o $*.P $<

.SUFFIXES:	.C .x

$(BINDIR)/% : %.C
	@$(MAKEDEPEND);
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)


SRCS := $(wildcard *.C)

PROGRAMS = $(SRCS:.C=)

all:	build
	@ \
	for file in $(PROGRAMS); do \
	  echo; \
	  echo '===' $$file '==='; \
	  $(BINDIR)/$$file; \
	done; \
	echo

build:	$(BINDIR) $(SRCS:%.C=$(BINDIR)/%)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f *.P

-include $(SRCS:.C=.P)
