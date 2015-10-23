#!/usr/bin/env python

import sys
import pylab
import matplotlib.pyplot as plt
import matplotlib.colors as col
import numpy as np

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

DPI = 300

def set_plot_defaults():
    params = {
        'figure.figsize'       : (4.0, 3.0),
        'figure.dpi'           : DPI,
        'figure.subplot.left'  : 0.06,
        'figure.subplot.right' : 0.95,
        'figure.subplot.bottom': 0.05,
        'figure.subplot.top'   : 0.95,
        'savefig.dpi'          : DPI,
        'savefig.bbox'         : 'tight',
        'savefig.pad_inches'   : 0,
        'image.origin'         : 'lower',
        'axes.linewidth'       : 0.5,
        'xtick.labelsize'      : 3,
        'xtick.major.size'     : 2,
        'xtick.minor.size'     : 1,
        'ytick.labelsize'      : 3,
        'ytick.major.size'     : 2,
        'ytick.minor.size'     : 1
        }
    plt.rcParams.update(params)


def shuffle(data, positions = None):
    out = np.zeros(data.shape, np.uint32) + 0xff000000

    if positions is not None:
        n = len(positions)
        h = np.histogram(data, bins = n, range = (0, n))[0]

        tmp = np.zeros(data.shape, np.uint32)
        for i in range(n):
            if h[i] > 0:
                p = list(int(x) for x in positions[i])
                val = ((p[0] * 37) + p[1]) * 37 + p[2]
                tmp |= np.where(data == i, val, 0)
    else:
        tmp = data.astype(np.uint32)

    bmap = [ 7, 15, 23, 6, 14, 22, 5, 13, 21, 4, 12, 20, 3, 11, 19 ]

    for i in range(15):
        out |= ((tmp >> i) & 1) << bmap[i]

    return out


def split(data):
    out = np.zeros(data.shape + (4,), np.uint8)
    out[:,:,0] = data & 0xff
    out[:,:,1] = (data >> 8) & 0xff
    out[:,:,2] = (data >> 16) & 0xff
    out[:,:,3] = (data >> 24) & 0xff

    return out


def colorize(data, color, alpha):
    (r, g, b) = color
    mask = np.where(data, 0xff, 0)

    out = np.zeros(data.shape + (4,), np.uint8)
    out[:,:,0] = mask * r
    out[:,:,1] = mask * g
    out[:,:,2] = mask * b
    out[:,:,3] = mask * alpha

    return out


def grayscale(data, options):
    vmin = min(0, data.min())
    vmax = max(0, data.max())

    if options['split']:
        positive = np.where(data > 0, data, 0)
        negative = np.where(data < 0, data, 0)
        gray = (positive * 0xff / vmax) + ((negative - vmin) * 0xff / -vmin)
    else:
        gray = (data - vmin) * 0xff / (vmax - vmin)

    if not options['inverse']:
        gray = 0xff - gray

    out = np.zeros(data.shape + (4,), np.uint8)
    out[:,:,0] = gray
    out[:,:,1] = gray
    out[:,:,2] = gray
    out[:,:,3] = options['alpha'] * 0xff

    return out


def blend(dst, src):
    srca = src[:,:,3].astype(np.float32) / 0xff

    out = dst.copy()
    for x in (0,1,2):
        srcx = src[:,:,x]
        dstx = dst[:,:,x]
        out[:,:,x] = srcx * srca + dstx * (1 - srca)

    out[:,:,3] = 0xff

    return out


def plot_points(points, dims, dim, sz, color, alpha = 1):
    ps = list(points[i] for i in range(len(points)) if dims[i] == dim)
    if len(ps) > 0:
        plt.scatter(list(p[0] for p in ps), list(p[1] for p in ps),
                    s = sz, c = color, alpha = alpha, linewidths = 0.2)


def plot(data, options):
    z = options.z
    dim = 2 if data['scalars'].shape[2] == 1 else 3

    scalars = data['scalars'][z]
    directions = data['directions'][z]

    watersheds = colorize(data['watersheds'][z],
                          options.watersheds_color, options.watersheds_alpha)
    paths = colorize(data['paths'][z], options.paths_color, options.paths_alpha)
    skeleton = colorize(np.where(data['skeleton'][z] <= options.watermark,
                                 True, False),
                        options.skeleton_color, options.skeleton_alpha)

    basins = split(shuffle(data['basins'][z],
                           data['critical'] if options.stable_colors else None))
    basins[:,:,3] = options.basins_alpha * 0xff

    pores = basins.copy()
    pores[:,:,3] = np.where(scalars > options.watermark,
                            0, options.pores_alpha * 0xff)

    tmp = data['critical']
    idx = list(i for i in range(len(tmp)) if tmp[i][2] == z)
    critical = list(tmp[i][:2] for i in idx)
    critdim = list(data['critdim'][i] for i in idx)

    fig, ax = plt.subplots(figsize =
                           tuple(float(n) / DPI
                                 for n in (options.width, options.height)))

    img = np.zeros(basins.shape, np.uint8) + 0xff

    if options.show_scalars:
        img = blend(img, grayscale(scalars, {
                    'inverse': options.scalars_inverse,
                    'alpha'  : options.scalars_alpha,
                    'split'  : options.scalars_split }))

    if options.show_basins:
        img = blend(img, basins)
    if options.show_pores:
        img = blend(img, pores)
    if options.show_watersheds:
        img = blend(img, watersheds)
    if options.show_paths:
        img = blend(img, paths)
    if options.show_skeleton:
        img = blend(img, skeleton)

    imgplot = plt.imshow(img)
    imgplot.set_interpolation(options.interpolation)

    if options.show_contour:
        plt.contour(scalars,
                    levels = [options.watermark],
                    colors = [options.contour_color],
                    linewidths = 0.5,
                    alpha = options.contour_alpha)

    if options.show_field:
        x, y = directions.shape
        X, Y = pylab.meshgrid(pylab.arange(0, x, 1) * 0.5,
                              pylab.arange(0, y, 1) * 0.5)

        dirs = ((0, 0, 0), (0, 0, 0),
                (1, 0, 0), (-1, 0, 0),
                (0, 1, 0), (0, -1, 0),
                (0, 0, 1), (0, 0, -1))
        tmp = np.array([[dirs[d] for d in row] for row in directions])

        U = tmp[:,:,0] * 0.8
        V = tmp[:,:,1] * 0.8
        plt.quiver(X[::2, ::2], Y[::2, ::2], U[::2, ::2], V[::2, ::2],
                   angles='xy', scale_units='xy', scale=1,
                   units='dots', width=1.5,
                   headlength=3, headaxislength=3, headwidth=3)

    if options.show_minima:
        plot_points(critical, critdim, 0, options.points_size,
                    options.minima_color, options.minima_alpha)
    if options.show_saddles:
        for d in range(1, dim - 1):
            plot_points(critical, critdim, d, options.points_size,
                        options.saddles_color, options.saddles_alpha)
    if options.show_maxima:
        plot_points(critical, critdim, dim - 1, options.points_size,
                    options.maxima_color, options.maxima_alpha)

    d = options.padding + 0.5
    plt.xlim([-d, scalars.shape[0] + d - 1])
    plt.ylim([-d, scalars.shape[1] + d - 1])

    ax.get_xaxis().set_major_locator(MultipleLocator(50))
    ax.get_xaxis().set_major_formatter(FormatStrFormatter('%d'))
    ax.get_xaxis().set_minor_locator(MultipleLocator(10))
    ax.get_yaxis().set_major_locator(MultipleLocator(50))
    ax.get_yaxis().set_major_formatter(FormatStrFormatter('%d'))
    ax.get_yaxis().set_minor_locator(MultipleLocator(10))

    if options.output:
        plt.savefig(options.output)
    else:
        plt.show()


colors_by_name = {
    'white'  : (1.0, 1.0, 1.0),
    'black'  : (0.0, 0.0, 0.0),
    'grey'   : (0.5, 0.5, 0.5),
    'gray'   : (0.5, 0.5, 0.5),
    'red'    : (1.0, 0.0, 0.0),
    'green'  : (0.0, 1.0, 0.0),
    'blue'   : (0.0, 0.0, 1.0),
    'yellow' : (1.0, 1.0, 0.0),
    'magenta': (1.0, 0.0, 1.0),
    'cyan'   : (0.0, 1.0, 1.0),
    'orange' : (1.0, 0.5, 0.0),
    'brown'  : (0.4, 0.2, 0.0)
    }


def parse_color(s):
    s = s.lower()

    if colors_by_name.has_key(s):
        return colors_by_name[s]
    else:
        if s.startswith('#'):
            s = s[1:]
        elif s.startswith('0x'):
            s = s[2:]

        n = int(s, 16)

        r = (n >> 16) & 0xff
        g = (n >> 8) & 0xff
        b = n & 0xff
        
        return tuple(float(x) / 0xff for x in (r, g, b))


def add_feature_options(parser, name, default, alpha, color, description):
    parser.add_option("--show-%s" % name,
                      dest = "show_%s" % name,
                      default = default,
                      action = 'store_true',
                      help = ("show %s%s" %
                              (description, ' (default)' if default else '')))
    parser.add_option("--hide-%s" % name,
                      dest = "show_%s" % name,
                      default = default,
                      action = 'store_false',
                      help = ("do not show %s%s" %
                              (description, '' if default else ' (default)')))
    parser.add_option("--%s-alpha" % name,
                      dest = "%s_alpha" % name,
                      metavar = 'X',
                      type = 'float',
                      default = alpha,
                      help = ("opacity of %s (between 0 and 1, default %f)" %
                              (name, alpha)))
    if color is not None:
        parser.add_option("--%s-color" % name,
                          dest = "%s_color" % name,
                          metavar = 'C',
                          type = 'color',
                          default = parse_color(color),
                          help = ("color of %s (name or hexcode, default %s))" %
                                  (name, color)))


def parse_options():
    from copy import copy
    from optparse import Option, OptionParser, OptionValueError

    def check_color(option, opt, value):
        try:
            return parse_color(value)
        except ValueError:
            raise OptionValueError(
                "option %s: invalid color value: %s" % (opt, value))

    class MyOption(Option):
        TYPES = Option.TYPES + ('color',)
        TYPE_CHECKER = copy(Option.TYPE_CHECKER)
        TYPE_CHECKER['color'] = check_color


    parser = OptionParser("usage: %prog [OPTIONS] INFILE",
                          option_class = MyOption)
    parser.add_option('-f', '--field', dest = 'field', metavar = 'FILE',
                      default = '',
                      help = 'file containing a pre-computed vector field')
    parser.add_option('-t', '--threshold', dest = 'threshold', metavar = 'X',
                      type = 'float', default = 1.0,
                      help = 'simplification threshold (default 1.0)')
    parser.add_option('--watermark', dest = 'watermark', metavar = 'X',
                      type = 'float', default = 0.0,
                      help = 'upper limit for pore space (default 0.0)')
    parser.add_option('-z', dest = 'z', metavar = 'N', type = 'int', default = 0,
                      help = 'z coordinate of slice to plot (default 0)')

    parser.add_option('-o', '--output', dest = 'output', metavar = 'FILE',
                      default = None,
                      help = 'output file name (plot on screen if missing)')
    parser.add_option('--wd', '--width', dest = 'width', metavar = 'N',
                      default = 1200,
                      help = 'width of figure in pixels (default 1200)')
    parser.add_option('--ht', '--height', dest = 'height', metavar = 'N',
                      default = 900,
                      help = 'height of figure in pixels (default 900)')
    parser.add_option('--padding', dest = 'padding',
                      type = 'float', metavar = 'X', default = 5.0,
                      help = 'figure padding in source pixels (default 5.0)')
    parser.add_option('--stable-colors', dest = 'stable_colors',
                      default = False, action = 'store_true',
                      help = "color basins/pores by minima positions (slower)")
    parser.add_option('--interpolation', dest = 'interpolation',
                      metavar = 'MODE', default = 'bicubic',
                      help = "data pixel interpolation mode (default bicubic)")

    add_feature_options(parser, 'scalars',    True,  1.0, None,
                        'input image data')
    parser.add_option('--scalars-inverse', dest = 'scalars_inverse',
                      default = False, action = 'store_true',
                      help = "show as black-to-white instead of white-to-black")
    parser.add_option('--scalars-split', dest = 'scalars_split',
                      default = False, action = 'store_true',
                      help = "separate gradients for positive/negative scalars")
    add_feature_options(parser, 'basins',     False, 0.5, None,
                        'full basins')
    add_feature_options(parser, 'pores',      True,  0.7, None,
                        'pores')
    add_feature_options(parser, 'watersheds', False, 0.4, 'blue',
                        'watersheds')
    add_feature_options(parser, 'contour',    True,  0.2, 'black',
                        'outline of the pore space')
    add_feature_options(parser, 'paths',      False, 1.0, 'cyan',
                        'paths between critical cells')
    add_feature_options(parser, 'skeleton',   True,  1.0, 'white',
                        'the Morse skeleton')
    add_feature_options(parser, 'minima',     True,  1.0, 'black',
                        '0-dimensional critical cells')
    add_feature_options(parser, 'saddles',    False, 1.0, 'cyan',
                        '1-dimensional critical cells')
    add_feature_options(parser, 'maxima',     False, 1.0, 'white',
                        '2-dimensional critical cells')
    add_feature_options(parser, 'field',      False, 0.5, 'white',
                        'gradient vector field')
    parser.add_option('--points-size', dest = 'points_size', metavar = 'X',
                      type = 'float', default = 1.5,
                      help = "point size for plotting critical points")

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("at least one argument needed; -h shows list of options")

    return options, args


if __name__ == '__main__':
    import sys, os.path
    from MorseAnalysis import VolumeImage, VectorField

    (options, args) = parse_options()
    infile = args[0]
    
    if infile.endswith('.nc'):
        img = VolumeImage(infile)
        morse = VectorField(img,
                            threshold = options.threshold,
                            filename = options.field)
        critical = list(morse.criticalCells())

        data = {
            'scalars'   : img.data(),
            'directions': morse.data(),
            'basins'    : morse.basinMap(),
            'watersheds': morse.watersheds(),
            'skeleton'  : morse.skeleton(),
            'paths'     : morse.paths(),
            'critical'  : critical,
            'critval'   : map(img.scalarForCell, critical),
            'critdim'   : map(img.cellDimension, critical)
            }
    else:
        raise RuntimeException('must have an NetCDF input file')

    set_plot_defaults()
    plot(data, options)
