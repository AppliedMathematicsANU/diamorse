#!/usr/bin/env python3


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        prog='grayscales_to_pores',
        description='Computes pore labels using discrete Morse theory.',
    )
    parser.add_argument(
        "-t", "--simplification_threshold",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "-w", "--watermark",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "inputfile"
    )
    parser.add_argument(
        "outputfile",
        nargs='?'
    )

    return parser.parse_args()


def compute_pores(input, threshold, watermark):
    from diamorse import MorseVectorField

    return MorseVectorField(input, threshold).pore_labels(watermark)


def barename(path):
    import os

    return os.path.splitext(os.path.basename(path))[0]


if __name__ == "__main__":
    import numpy as _np
    from diamorse import read_netcdf, write_netcdf

    args = parse_arguments()

    input = read_netcdf(args.inputfile)

    output = compute_pores(
        input,
        threshold=args.simplification_threshold,
        watermark=args.watermark
    )

    outname = f"{barename(args.inputfile)}_pores"
    outputfile = args.outputfile or f"{outname}.nc"
    write_netcdf(outputfile, outname, output)
