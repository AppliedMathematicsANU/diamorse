.PHONY: base main all src util python extra test clean

base:	main util

all:    base python extra

main:
	$(MAKE) -C src/main

util:
	$(MAKE) -C src/util

python:
	$(MAKE) -C src/python

extra:
	$(MAKE) -C src/experimental

test:
	$(MAKE) -C test

clean:
	$(MAKE) -C src/main clean
	$(MAKE) -C src/util clean
	$(MAKE) -C src/python clean
	$(MAKE) -C src/experimental clean
	$(MAKE) -C test clean
