#!/usr/bin/env python3

import argparse
import numpy


class Fraction(float):
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if (new < 0) or (new > 1):
			raise ValueError("Fraction must be between 0 and 1 (both boundaries"
				" included), got '%f'" % new)
		return new


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("abund", type=Fraction, help="abundance (0.0-1.0) of the "
		"dominant base type")
	
	args = ap.parse_args()
	return args


def dominant_to_entropy(abund, n=5):
	abund_arr = numpy.asarray([abund] + ([(1 - abund) / (n - 1)] * (n - 1)))
	return numpy.sum(abund_arr * numpy.log(abund_arr)) * -1


def main():
	args = get_args()
	e = dominant_to_entropy(args.abund)
	print(e)
	return


if __name__ == "__main__":
	main()
