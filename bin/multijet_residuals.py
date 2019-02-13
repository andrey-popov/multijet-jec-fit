#!/usr/bin/env python

"""Plots pre- and postfit residuals for multijet analysis."""

import argparse
import json
import os
import re

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import jecfit
from utils import mpl_style


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        'inputs', help='ROOT file with inputs from multijet analysis'
    )
    arg_parser.add_argument('fits', help='JSON file with results of fits')
    arg_parser.add_argument('-p', '--period', help='Period to be selected')
    arg_parser.add_argument('-m', '--method', help='Method to be selected')
    arg_parser.add_argument(
        '-o', '--output', default='fig/multijet_residuals.pdf',
        help='Name for output figure file'
    )
    args = arg_parser.parse_args()

    fig_dir = os.path.dirname(args.output)

    if fig_dir:
        try:
            os.makedirs(fig_dir)
        except FileExistsError:
            pass


    with open(args.fits) as f:
        fits = json.load(f)

    if isinstance(fits, list):
        try:
            fit = next(filter(
                lambda f: (
                    f['period'] == args.period and f['variant'] == args.method
                ), fits
            ))
        except StopIteration:
            raise RuntimeError(
                'File "{}" does not contain an entry for period "{}" and '
                'method "{}".'.format(args.fits, args.period, args.method)
            )
    else:
        # This is not a list.  Assume then that the file contains
        # results of a single fit.
        fit = fits


    poi_regex = re.compile(r'^p[0-9]+$')
    poi_values = []
    nuisances = {}

    for param in fit['parameters']:
        if poi_regex.match(param['name']):
            poi_values.append(param['value'])
        else:
            nuisances[param['name']] = param['value']


    measurement = jecfit.MultijetChi2(args.inputs, args.method, {'JER'})
    max_pt = 1.6e3
    measurement.set_pt_range(0., max_pt)

    prefit = measurement.compute_residuals(
        np.zeros_like(poi_values), np.zeros(len(nuisances))
    )
    postfit = measurement.compute_residuals(poi_values, nuisances)

    iend = np.searchsorted(prefit[0], max_pt)
    prefit = [prefit[i][:iend] for i in range(3)]
    postfit = [postfit[i][:iend] for i in range(3)]


    plt.style.use(mpl_style)

    fig = plt.figure()
    fig.patch.set_alpha(0.)
    axes = fig.add_subplot(111)

    axes.errorbar(
        prefit[0], prefit[1], yerr=prefit[2], label='Pre-fit',
        marker='o', c='black', lw=0, elinewidth=0.8
    )
    axes.errorbar(
        postfit[0], postfit[1], yerr=postfit[2], label='Post-fit',
        marker='o', c='C0', lw=0, elinewidth=0.8
    )
    axes.axhline(0., c='black', ls='dashed', lw=0.8)

    axes.legend(loc='upper right')
    axes.set_xscale('log')

    axes.xaxis.set_major_formatter(mpl.ticker.LogFormatter())
    axes.xaxis.set_minor_formatter(
        mpl.ticker.LogFormatter(minor_thresholds=(2, 0.4))
    )

    axes.set_xlabel(r'$\tau_{1}^\mathrm{L3}$ [GeV]')
    axes.set_ylabel(
        r'$\langle B^\mathrm{{L3}}_\mathrm{{{0}}}\rangle\//\/'
        r'\langle B^\mathrm{{Sim}}_\mathrm{{{0}}}\rangle - 1$'.format(
            'jet' if args.method == 'PtBal' else 'MPF'
        )
    )

    axes.text(
        1., 1.002, args.period,
        ha='right', va='bottom', transform=axes.transAxes
    )

    fig.savefig(args.output)
    plt.close(fig)

