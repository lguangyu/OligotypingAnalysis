#!/usr/bin/env python3

import argparse
import io
import itertools
import matplotlib
import matplotlib.pyplot
import numpy
import os
import sys

import mpllayout


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("oligo_output", type=str,
		help="oligotyping output directory")
	ap.add_argument("-l", "--oligo-list", type=str, metavar="file",
		help="plot only listed oligos instead of all, if providied [no]")
	# ap.add_argument("-d", "--delimiter", type=str, default="\t",
	# metavar="str",
	# help="delimiter in input/output tables [<tab>]")
	ap.add_argument("-n", "--num-oligos", type=int, default=20,
		metavar="int",
		help="number of oligos to plot, only applies if -l/--oligo-list is not "
			"set [20]")
	ap.add_argument("-p", "--plot", type=str,
		metavar="png",
		help="the output plot image [stdout]")
	ap.add_argument("--dpi", type=int, default=300,
		metavar="int",
		help="plot image DPI [300]")

	# parse and refine arsg
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer

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


class OligoCountTable(object):
	@property
	def n_oligos(self) -> int:
		return len(self.oligos)

	@property
	def n_samples(self) -> int:
		return len(self.samples)

	@property
	def sample_count_sum(self) -> numpy.ndarray:
		return self.data.sum(axis=1)

	def validate_data(self) -> None:
		if self.data.ndim != 2:
			raise ValueError("self.data must be a 2-dimensional matrix")
		nrow, ncol = self.data.shape
		if len(self.samples) != nrow:
			raise ValueError("self.data number of rows must match length of"
				"self.samples")
		if len(self.oligos) != ncol:
			raise ValueError("self.data number of columns must match length of"
				"self.oligos")
		return

	def set_data(self, samples, oligos, data) -> None:
		self.samples = samples
		self.oligos = oligos
		self.data = data
		self.validate_data()
		return

	def __init__(self, samples, oligos, data, **kw):
		super().__init__(**kw)
		self.set_data(samples, oligos, data)
		return

	@classmethod
	def from_oligo_outptu_dir(cls, oligo_output_dir: str, *, delimiter="\t"):
		if not os.path.isdir(oligo_output_dir):
			raise IOError("'%s' is not a directory" % oligo_output_dir)
		fname = os.path.join(oligo_output_dir, "MATRIX-COUNT.txt")
		if not os.path.isfile(fname):
			raise IOError("cannot find file '%s', oligotyping output might be "
				"incomplete" % fname)
		return cls.from_oligo_output_matrix_count(fname, delimiter=delimiter)

	@classmethod
	def from_oligo_output_matrix_count(cls, fname: str, *, delimiter="\t"):
		# read data
		raw = numpy.loadtxt(fname, dtype=object, delimiter=delimiter)
		samples = raw[1:, 0]
		oligos = raw[0, 1:]
		data = raw[1:, 1:].astype(int)

		new = cls(samples=samples, oligos=oligos, data=data)
		return new

	def select_by_oligo_list(self, oligo_list, **kw):
		select_oligos = set(oligo_list)
		mask = [i in select_oligos for i in self.oligos]

		# create new object
		samples = self.samples.copy()
		oligos = self.oligos[mask]
		data = self.data[:, mask]

		new = type(self)(samples=samples, oligos=oligos, data=data, **kw)
		return new

	def select_by_oligo_list_file(self, oligo_list_file, **kw):
		with get_fp(oligo_list_file, "r") as fp:
			oligo_list = fp.read().splitlines()
		return self.select_by_oligo_list(oligo_list, **kw)


def setup_layout(n_samples: int) -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.3, right_margin=2.3,
		top_margin=0.1, bottom_margin=1.5,
	)

	# add axes
	ax = lc.add_frame("axes")
	ax.set_anchor("bottomleft")
	ax.set_size(0.2 * n_samples, 3.0)

	# create layout
	layout = lc.create_figure_layout()

	# apply axes style
	ax = layout["axes"]
	ax.set_facecolor("#e8e8f8")
	ax.tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=False, labelbottom=True,
		top=False, labeltop=False,
	)

	return layout


class OligoColorList(list):
	@staticmethod
	def _drop_replicates(vec: list):
		seen = set()
		ret = list()
		for i in vec:
			if i not in seen:
				ret.append(i)
				seen.add(i)
		return ret

	@classmethod
	def get_default_list(cls):
		proto = matplotlib.pyplot.get_cmap("Set3").colors \
			+ matplotlib.pyplot.get_cmap("Pastel2").colors
		new = cls(cls._drop_replicates(proto))
		return new

	def __getitem__(self, idx):
		n = super().__len__()
		if idx >= n:
			idx %= n
		return super().__getitem__(idx)

	def __iter__(self):
		return itertools.cycle(super().__iter__())


def plot_oligo_abund_stackbar(png, count_table: OligoCountTable, *,
		oligo_list_file=None, dpi=300):
	# calculate the total count from the original table
	total_counts = count_table.sample_count_sum

	# select oligos if necessary
	if oligo_list_file is not None:
		count_table = count_table.select_by_oligo_list_file(oligo_list_file)

	oligos, samples = count_table.oligos, count_table.samples
	n_oligos, n_samples = count_table.n_oligos, count_table.n_samples

	# create layout
	layout = setup_layout(n_samples)
	figure = layout["figure"]

	# plot
	colors = OligoColorList.get_default_list()

	ax = layout["axes"]
	counts = count_table.data
	fracs = counts / total_counts.reshape(-1, 1)
	bottom = numpy.zeros(n_samples, dtype=float)
	handles = list()

	x = numpy.arange(n_samples) + 0.5
	for i in range(n_oligos):
		height = fracs[:, i]
		p = ax.bar(x, height, width=0.8, bottom=bottom,
			align="center", edgecolor="#404040", linewidth=0.5,
			facecolor=colors[i], label=oligos[i])
		bottom += height
		handles.append(p)
	# add 'others'
	p = ax.bar(x, 1.0 - bottom, width=0.8, bottom=bottom,
		align="center", edgecolor="#404040", linewidth=0.5,
		facecolor="#ffffff", label="others")
	handles.append(p)

	# legend
	ax.legend(handles=handles, loc=2, fontsize=8,
		bbox_to_anchor=(1.002, 1.02))

	# axis misc
	ax.set_xlim(0, n_samples)
	ax.set_ylim(0.0, 1.0)
	ax.set_ylabel("Oligo abundances")
	ax.set_xticks(x)
	ax.set_xticklabels(samples, rotation=90, fontsize=8,
		horizontalalignment="center", verticalalignment="top")

	# save figure and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def main():
	args = get_args()
	# load data
	count_table = OligoCountTable.from_oligo_outptu_dir(args.oligo_output)
	# plot
	plot_oligo_abund_stackbar(args.plot, count_table,
		oligo_list_file=args.oligo_list, dpi=args.dpi)
	return


if __name__ == "__main__":
	main()
