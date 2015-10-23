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

    infile = sys.argv[1]
    threshold = float(sys.argv[2])
    dim = int(sys.argv[3])

    pairs = persistence.fromTextFile(infile)

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
