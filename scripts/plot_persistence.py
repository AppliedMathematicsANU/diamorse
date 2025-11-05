#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.colors as col


infinity = float('inf')


def get_births(pairs, dim, threshold):
    return tuple(birth for birth, death, d, _ in pairs
                 if d == dim and death - birth > threshold)


def get_deaths(pairs, dim, threshold):
    return tuple(death for birth, death, d, _ in pairs
                 if d == dim and death - birth > threshold)


def get_weights(pairs, dim, threshold):
    return tuple(weight for birth, death, d, weight in pairs
                 if d == dim and death - birth > threshold)


def bars(pairs, dim, threshold):
    pairs = sorted(pairs)
    births = get_births(pairs, dim, threshold)
    deaths = get_deaths(pairs, dim, threshold)
    indexes = range(len(births))
    lifetimes = list(deaths[i] - births[i] for i in indexes)

    return plt.bar(indexes, lifetimes, 1.0, births)


def deathsVersusBirths(pairs, dim, threshold):
    births = get_births(pairs, dim, threshold)
    deaths = get_deaths(pairs, dim, threshold)

    return plt.plot(births, deaths, '.')


def weightsVersusPersistence(pairs, dim, threshold):
    births = get_births(pairs, dim, threshold)
    deaths = get_deaths(pairs, dim, threshold)
    spans = tuple(d - b for (b, d) in zip(births, deaths))
    weights = get_weights(pairs, dim, threshold)
    if max(weights) > 0:
        return plt.semilogy(spans, weights, '.')
    else:
        return plt.plot(spans, weights, '.')


def deathsVersusBirthsHistogram(pairs, dim, threshold, nbins=100):
    pairs = [ p for p in pairs if p[1] < infinity ]
    births = get_births(pairs, dim, threshold)
    deaths = get_deaths(pairs, dim, threshold)

    axmin = min([min(births), min(deaths),0])
    axmax = max([max(births), max(deaths), 0])

    return plt.hist2d(
        births,
        deaths,
        [nbins,nbins],
        [[axmin,axmax],[axmin,axmax]],
        False,
        None,
        1,
        norm=col.LogNorm()
    )



def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Process and plot.')
    parser.add_argument('infile', help='file containing the input data')
    parser.add_argument(
        '-t', '--threshold', metavar = 'X',
        type = float, default = 1.0,
        help = 'simplification threshold (default 1.0)'
    )
    parser.add_argument(
        '-d', '--dimensions',
        type = int, default = 0,
        help = 'kind of critical points to show. E.g. 0 is for minima vs. 1-saddle point'
    )

    return parser.parse_args()


if __name__ == '__main__':
    from diamorse import MorseVectorField

    args = parse_arguments()
    threshold = args.threshold
    dim = args.dimensions

    morse = MorseVectorField(args.infile, threshold=threshold)

    pairs = tuple(
        (v.value, 0 if w is None else w.value, v.dimension, v.weight)
        for v, w in morse.birth_death_pairs()
    )

    if len(get_births(pairs, dim, threshold)) == 0:
        print("Nothing to plot.")
    else:
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
