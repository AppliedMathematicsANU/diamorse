#!/usr/bin/env python

import sys
from MorseAnalysis import VolumeImage, VectorField


infinity = float('inf')


def fromVolumeFile(filename, options):
    img = VolumeImage(filename)
    morse = VectorField(img,
                        threshold = options.threshold,
                        filename = options.field)

    dim = lambda v: img.cellDimension(v)
    val = lambda v: img.scalarForCell(v) if v else infinity

    weights = dict((tuple(v), x) for v, x in morse.weights())

    return tuple((val(v), val(w), dim(v), v, w, weights.get(tuple(v), 0))
                 for v, w in morse.birthsAndDeaths())


def fromText(input):
    res = []
    for line in input.readlines():
        text = line.strip()
        if len(text) == 0 or text.startswith('#'):
            continue
        try:
            birth, death, dim, v0, v1, v2, w0, w1, w2, wt = text.split()
        except ValueError:
            birth, death, dim, v0, v1, v2, w0, w1, w2 = text.split()
            wt = 0

        v = map(float, [v0, v1, v2])
        w = map(float, [w0, w1, w2]) if w0 != '-' else None
        dim = int(dim)
        birth = float(birth)
        death = float(death)
        weight = int(wt)

        res.append((birth, death, dim, v, w, weight))

    return tuple(res)


def fromTextFile(filename):
    fp = open(filename)
    result = fromText(fp)
    fp.close()
    return result


def toText(output, pairs, source):
    output.write("# Persistence pairs for %s\n" % (source,))
    output.write("#   format: ")
    output.write("<birth> <death> <dimension> <creator xyz> <destructor xyz> ")
    output.write("<weight>\n")
    for birth, death, dim, v, w, wt in pairs:
        w = w or '---'
        output.write(
            "% 14.10f % 14.10f %d    % 6s % 6s % 6s    % 6s % 6s % 6s   %d\n" %
            (birth, death, dim, v[0], v[1], v[2], w[0], w[1], w[2], wt))


def toTextFile(filename, pairs, source):
    fp = open(filename, "w")
    toText(fp, pairs, source)
    fp.close()


def births(pairs, dim, threshold):
    return tuple(birth for birth, death, d, v, w, weight in pairs
                 if d == dim and death - birth > threshold)

def deaths(pairs, dim, threshold):
    return tuple(death for birth, death, d, v, w, weight in pairs
                 if d == dim and death - birth > threshold)

def weights(pairs, dim, threshold):
    return tuple(weight for birth, death, d, v, w, weight in pairs
                 if d == dim and death - birth > threshold)

def locations(pairs, dim, threshold):
    return tuple( (v, w)  for birth, death, d, v, w, weight in pairs
                 if d == dim and death - birth > threshold)


def bettiNumbers(pairs, dim, threshold):
    bs = tuple((birth,  1) for birth in births(pairs, dim, threshold))
    ds = tuple((death, -1) for death in deaths(pairs, dim, threshold))
    events = sorted(bs + ds)

    result = []

    if events:
        (x, n, n0) = (events[0][0], 0, -1)
        for (y, m) in events:
            if y != x:
                if n != n0:
                    result.append((x, n))
                    n0 = n
                x = y
            n += m
        if n != n0:
            result.append((x, n))

    return tuple(result)

def rankHomAB(pairs,dim,a,b):
    count = 0
    for birth, death, d, v, w, weight in pairs:
        if (d == dim and birth < a and death > b):
                count += 1
    return count

def rankHomFunction(pairs,dim,bvals,dvals):
    countKeys = [(b,d) for b in bvals for d in dvals if d >= b ]
    count = { ck: 0 for ck in countKeys }
    for birth, death, d, v, w, wt in pairs:
        if d == dim:
            for (b,d) in countKeys:
                if (b >= birth and d <= death):
                    count[(b,d)] += 1

    return count


def printStats(pairs):
    for dim in range(3):
        counts = defaultdict(int)
        for birth, death, d, v, w, wt in pairs:
            if d == dim:
                counts[death - birth] += 1

        cum = 0
        for x in reversed(sorted(counts.keys())):
            cum += counts[x]
            counts[x] = cum

        print ("# Numbers of cells of dimension " + str(d) +
               " by lower persistence thresholds")
        for x in sorted(counts.keys()):
            print "%10.5f %6d" % (x, counts[x])
        print


if __name__ == '__main__':
    import argparse
    from collections import defaultdict

    parser = argparse.ArgumentParser("usage: %prog [OPTIONS] INFILE")
    parser.add_argument('infile', help='file containing the field')
    parser.add_argument('-f', '--field', dest = 'field', metavar = 'FILE',
                      default = '',
                      help = 'file containing a pre-computed vector field')
    parser.add_argument('-t', '--threshold', dest = 'threshold', metavar = 'X',
                      type = float, default = 1.0,
                      help = 'simplification threshold (default 1.0)')
    parser.add_argument("-b", "--betti", dest = "betti", default = False,
                      action = "store_true", help = "output Betti numbers")
    parser.add_argument("-r", "--raw", dest = "raw", default = False,
                      action = "store_true", help = "output persistence pairs")
    parser.add_argument("-s", "--stats", dest = "stats", default = False,
                      action = "store_true", help = "output some statistics")
    options = parser.parse_args()

    infile = options.infile
    threshold = options.threshold

    if infile.endswith(".nc") or infile.endswith("_nc"):
        pairs = fromVolumeFile(infile, options)
    else:
        pairs = fromTextFile(infile)

    xth = [ "zeroth", "first", "second", "third" ]

    if options.raw:
        toText(sys.stdout, pairs, infile)

    if options.stats:
        printStats(pairs)

    if options.betti:
        for i in range(len(xth)):
            print "# The %f-persistent %s Betti numbers:" % (threshold, xth[i])
            for (val, count) in bettiNumbers(pairs, i, threshold):
                print "%10.5f %6d" % (val, count)
            print
