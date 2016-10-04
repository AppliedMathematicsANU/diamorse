#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import matplotlib.colors as col
import numpy as np

import persistence


infinity = float('inf')


def bars(pairs, dim, threshold):
    pairs = sorted(pairs)
    births = persistence.births(pairs, dim, threshold)
    deaths = persistence.deaths(pairs, dim, threshold)
    indexes = range(len(births))
    lifetimes = list(deaths[i] - births[i] for i in indexes)

    return plt.bar(indexes, lifetimes, 1.0, births)


def deathsVersusBirths(pairs, dim, threshold):
    births = persistence.births(pairs, dim, threshold)
    deaths = persistence.deaths(pairs, dim, threshold)
    return plt.plot(births, deaths, '.')


def weightsVersusPersistence(pairs, dim, threshold):
    births = persistence.births(pairs, dim, threshold)
    deaths = persistence.deaths(pairs, dim, threshold)
    spans = tuple(d - b for (b, d) in zip(births, deaths))
    weights = persistence.weights(pairs, dim, threshold)
    if max(weights) > 0:
        return plt.semilogy(spans, weights, '.')
    else:
        return plt.plot(spans, weights, '.')


def deathsVersusBirthsHistogram(pairs, dim, threshold, nbins=100):
    pairs = [ p for p in pairs if p[1] < infinity ]
    births = persistence.births(pairs, dim, threshold)
    deaths = persistence.deaths(pairs, dim, threshold)

    axmin = min([min(births), min(deaths),0])
    axmax = max([max(births), max(deaths), 0])

    return plt.hist2d(births,deaths,[nbins,nbins],[[axmin,axmax],[axmin,axmax]],False,None,1,norm=col.LogNorm())



if __name__ == '__main__':
    import sys, os.path

    import argparse
    parser = argparse.ArgumentParser(description='Process and plot.')
    parser.add_argument('infile', help='file containing the field')
    parser.add_argument('-f', '--field', metavar = 'FILE',
                        default = '',
                        help = 'file containing a pre-computed vector field')
    parser.add_argument('-t', '--threshold', metavar = 'X',
                        type = float, default = 1.0,
                        help = 'simplification threshold (default 1.0)')
    parser.add_argument('-d', '--dimensions',
                        type = int, default = 3,
                        help = 'number of dimensions (default 3)')
    parser.add_argument("-b", "--betti", dest = "betti", default = False,
                        action = "store_true", help = "output Betti numbers")
    parser.add_argument("-r", "--raw", dest = "raw", default = False,
                        action = "store_true", help = "output persistence pairs")
    parser.add_argument("-s", "--stats", dest = "stats", default = False,
                        action = "store_true", help = "output some statistics")
    args = parser.parse_args()
    infile = args.infile
    threshold = args.threshold
    dim = args.dimensions

    

    if infile.endswith(".nc") or infile.endswith("_nc"):
        pairs = persistence.fromVolumeFile(infile, args)
    else:
        pairs = fromTextFile(infile)

    plt.figure(1)
    plt.title('Cycle births and deaths, t = %.2f, d = %d' % (threshold, dim))
    bars(pairs, dim, threshold)

    plt.figure(2)
    plt.title('Deaths vs Births, t = %.2f, d = %d' % (threshold, dim))
    points = deathsVersusBirths(pairs, dim, threshold)
    plt.setp(points, color = 'black')
    plt.xlabel('Value at birth')
    plt.ylabel('Value at death')

    plt.figure(3)
    plt.title('Weight vs Persistence, t = %.2f, d = %d' % (threshold, dim))
    points = weightsVersusPersistence(pairs, dim, threshold)
    plt.setp(points, color = 'black')
    plt.xlabel('Feature persistence')
    plt.ylabel('Feature weight')
    
    plt.figure(4)
    plt.title('Deaths vs Births Histogram,  t =  %.2f, d = %d' % (threshold, dim))
    points = deathsVersusBirthsHistogram(pairs, dim, threshold)
    plt.xlabel('Level set value at birth')
    plt.ylabel('Level set value at death')
    plt.axhline(linewidth=1,color='black')
    plt.axvline(linewidth=1,color='black')

    plt.show()
