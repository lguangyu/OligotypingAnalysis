#!/usr/bin/env python3

import argparse
import collections
import itertools
import matplotlib
import matplotlib.pyplot
import numpy
import os
import re
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("dirname", type=str,
		help="directory name that contains all blastdbcmd taxonomy tables")
	ap.add_argument("-x", "--extension", type=str, default="blastdbcmd",
		metavar="str",
		help="the extension of files to scan [blastdbcmd]")
	ap.add_argument("-d", "--delimiter", type=str, default="\t",
		metavar="str",
		help="delimiter in input files [<tab>]")
	ap.add_argument("-k", "--key-field", type=int, default=2,
		metavar="int",
		help="the taxnomy value field number (0-based) in input files [0]")
	ap.add_argument("-p", "--plot", type=str, default="-",
		metavar="png",
		help="output plot image, or write to stdout by default '-'")
	# parse and refine arsg
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer
	return args


def read_oligo_tax_count(fname, key_field, *, delimiter="\t") \
		-> collections.Counter:
	with open(fname, "r") as fp:
		tax = [i.split(delimiter)[key_field] for i in fp.read().splitlines()]
	return collections.Counter(tax)


def read_all_oligo_tax_count(dirname, extension, key_field, **kw) -> dict:
	ret = dict()
	for i in os.scandir(dirname):
		if i.is_file() and i.name.endswith(extension):
			oligo = int(re.search(r"^(\d+)", i.name).group(1))
			ret[oligo] = read_oligo_tax_count(i.path, key_field, **kw)
	return ret


def setup_layout(figure, n_oligo):
	layout = dict()
	layout["figure"] = figure

	# margins
	left_margin_inch	= 0.3
	right_margin_inch	= 2.3
	top_margin_inch		= 0.1
	bottom_margin_inch	= 1.5

	# bar width
	bar_width_inch		= 0.2
	bar_height_inch		= 3.0
	axes_width_inch		= bar_width_inch * n_oligo
	axes_height_inch	= bar_height_inch

	# figure dimensions
	figure_width_inch	= left_margin_inch + axes_width_inch + right_margin_inch
	figure_height_inch	= top_margin_inch + axes_height_inch + bottom_margin_inch
	figure.set_size_inches(figure_width_inch, figure_height_inch)

	# create axes
	axes_left	= left_margin_inch / figure_width_inch
	axes_bottom	= bottom_margin_inch / figure_height_inch
	axes_width	= axes_width_inch / figure_width_inch
	axes_height	= axes_height_inch / figure_height_inch
	axes = figure.add_axes([axes_left, axes_bottom, axes_width, axes_height])
	layout["axes"] = axes
	# axes style
	#for sp in axes.spines.values():
	#	sp.set_visible(False)
	axes.set_facecolor("#E8E8F8")
	axes.tick_params(left = False, labelleft = False, bottom = False)

	return layout


def get_oligo_tax_bootstrap_matrix(oligo_tax: dict):
	# get unique oligos list
	u_oligo = sorted(oligo_tax.keys())

	# get unique taxons list
	tax_iter = itertools.chain(*[v.keys() for v in oligo_tax.values()])
	u_tax = sorted(set(tax_iter))

	# get oligo tax bootstrap matrix
	bs_mat = numpy.zeros((len(u_oligo), len(u_tax)), dtype=float)
	for i, o in enumerate(u_oligo):
		bs_total = oligo_tax[o].total()
		for j, t in enumerate(u_tax):
			bs_mat[i, j] = oligo_tax[o][t] / (bs_total or 1.0)

	return u_oligo, u_tax, bs_mat


def plot_oligo_tax_stackbar(png, oligo_tax: dict):
	u_oligo, u_tax, bs_mat = get_oligo_tax_bootstrap_matrix(oligo_tax)

	n_oligo, n_tax = len(u_oligo), len(u_tax)
	assert bs_mat.shape == (n_oligo, n_tax)

	# create layout
	figure = matplotlib.pyplot.figure()
	layout = setup_layout(figure, n_oligo)

	# plot
	colors = matplotlib.pyplot.get_cmap("Set3").colors \
		+ matplotlib.pyplot.get_cmap("Pastel2").colors
	axes = layout["axes"]
	bottom = numpy.zeros(n_oligo, dtype = float)
	handles = list()
	x = numpy.arange(n_oligo) + 0.5
	for h, t, c in zip(bs_mat.T, u_tax, colors):
		bar = axes.bar(x, h, width = 0.8, bottom = bottom,
			align = "center", edgecolor = "#404040", linewidth = 0.5,
			facecolor = c, label = t)
		bottom += h
		handles.append(bar)

	# check if 'others' needs to be added
	if len(handles) < n_tax:
		h = bs_mat[:, len(handles):].sum(axis=1)
		bar = axes.bar(x, h, width = 0.8, bottom = bottom,
			align = "center", edgecolor = "#404040", linewidth = 0.5,
			facecolor = "#ffffff", label = "others")
		bottom += h
		handles.append(bar)

	# legend
	axes.legend(handles = handles, loc = 2, fontsize = 8,
		bbox_to_anchor = (1.002, 1.02))

	# axis misc
	axes.set_xlim(0, n_oligo)
	axes.set_ylim(0.0, 1.0)
	axes.set_ylabel("Taxonomy bootstrap")
	axes.set_xticks(x)
	axes.set_xticklabels(["OLIGO_%03u" % i for i in u_oligo],
		rotation = 90, fontsize = 8,
		horizontalalignment = "center", verticalalignment = "top")

	matplotlib.pyplot.savefig(png, dpi = 300)
	matplotlib.pyplot.close()
	return


def main():
	args = get_args()
	# load data
	oligo_tax = read_all_oligo_tax_count(args.dirname, args.extension,
		key_field=args.key_field, delimiter=args.delimiter)
	# plot
	plot_oligo_tax_stackbar(args.plot, oligo_tax)
	return


if __name__ == "__main__":
	main()
