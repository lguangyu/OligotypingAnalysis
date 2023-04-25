#!/usr/bin/env python3

import argparse
import io
import os
import sys

import numpy


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("oligo_output", type=str,
		help="oligotyping output directory")
	ap.add_argument("--abund-threshold", "-a", type=float, default=0.05,
		metavar="float",
		help="do not list oligos if its abundances are less than this threshold"
			" across all samples, 0 means report all [0.05]")
	ap.add_argument("--count-threshold", "-c", type=int, default=0,
		metavar="int",
		help="do not list oligos if its total count is less than this threshold"
			" across all samples, 0 measn report all [0]")
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="txt",
		help="output oligo name list [stdout]")

	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout

	return args


def get_fp(f, *ka, factory=open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument of get_fp() must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


def get_oligo_data(oligo_output_dir: str) -> (numpy.ndarray, numpy.ndarray):
	if not os.path.isdir(oligo_output_dir):
		raise IOError("'%s' is not a directory" % oligo_output_dir)
	fname = os.path.join(oligo_output_dir, "MATRIX-COUNT.txt")
	if not os.path.isfile(fname):
		raise IOError("cannot find file '%s', oligotyping output might be "
			"incomplete" % fname)

	# read data
	raw = numpy.loadtxt(fname, dtype=object, delimiter="\t")
	oligos = raw[0, 1:]
	samples = raw[1:, 0]
	counts = raw[1:, 1:].astype(int)

	return oligos, counts


def filter_and_save_oligos(f, oligos, counts, abund_thres: float = 0.0,
		count_thres: int = 0) -> None:
	oligo_abund = counts / counts.sum(axis=1, keepdims=True)
	# select oligos by max abundance across samples
	abund_mask = (oligo_abund.max(axis=0) >= abund_thres)
	count_mask = (counts.sum(axis=0) >= count_thres)
	oligo_select = oligos[abund_mask & count_mask]
	# save file
	with get_fp(f, "w") as fp:
		for i in oligo_select:
			print(i, file=fp)
	return


def main():
	args = get_args()
	oligos, counts = get_oligo_data(args.oligo_output)
	filter_and_save_oligos(args.output, oligos, counts,
		abund_thres=args.abund_threshold,
		count_thres=args.count_threshold,
	)
	return


if __name__ == "__main__":
	main()
