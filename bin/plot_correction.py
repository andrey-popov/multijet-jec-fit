#!/usr/bin/env python

"""Plots fitted correction.

Uncertainties are decomposed into independent components.
"""

import argparse
import json
import os
import re

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

from config import Config
import jecfit
from utils import mpl_style


class Correction:
    """Wrapper for jet correction."""

    def __init__(self, corr_form='2p'):
        """Initialize from a C++ class."""
        self.corr = jecfit.create_correction(corr_form)


    def __call__(self, pt):
        """Evaluate correction at given pt.
        
        The argument can be a scalar or an array-like.
        """

        if np.isscalar(pt):
            return self.corr.Eval(pt)
        else:
            res = np.empty_like(pt)

            for i in range(len(pt)):
                res[i] = self.corr.Eval(pt[i])

            return res


    def set_params(self, params):
        """Set parameters of the correction."""
        self.corr.SetParams(params)


if __name__ == '__main__':

    plt.style.use(mpl_style)

    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('fits', help='JSON file with fit results')
    arg_parser.add_argument('-p', '--period', help='Period to be selected')
    arg_parser.add_argument('-m', '--method', help='Method to be selected')
    arg_parser.add_argument(
        '-o', '--output', default='fig/correction.pdf',
        help='Name for output figure file'
    )
    args = arg_parser.parse_args()

    fig_dir = os.path.dirname(args.output)

    if fig_dir:
        try:
            os.makedirs(fig_dir)
        except FileExistsError:
            pass

    if args.method == 'PtBal':
        method_label = '$p_\\mathrm{T}$ balance'
    else:
        method_label = 'MPF'

    config = Config('config/plot_config.yaml')


    with open(args.fits) as f:
        fit_infos = json.load(f)

    if isinstance(fit_infos, list):
        try:
            fit_info = next(filter(
                lambda f: (
                    f['period'] == args.period and f['variant'] == args.method
                ), fit_infos
            ))
        except StopIteration:
            raise RuntimeError(
                'File "{}" does not contain an entry for period "{}" and '
                'method "{}".'.format(args.fits, args.period, args.method)
            )
    else:
        # This is not a list.  Assume then that the file contains
        # results of a single fit.
        fit_info = fit_infos

    fit_results = jecfit.FitResults(fit_info)
    corr_form = fit_info.get('corr_form', '2p')


    # POIs always precede nuisances.  Count how many they are.
    poi_regex = re.compile(r'^p[0-9]+$')
    num_pois = sum(
        1 for p in fit_results.parameters if poi_regex.match(p.name)
    )

    corr_function = Correction(corr_form)
    nominal_params = np.asarray(
        [p.value for p in fit_results.parameters[:num_pois]]
    )
    covariance = fit_results.covariance_matrix[:num_pois, :num_pois]

    # Find shifts with respect to nominal values of POIs that correspond
    # to semiaxes of the 1 sigma ellipse.  Coordinate transformation
    # given by matrix v.T below diagonalizes the covariance matrix,
    # which becomes diag(w).  Transformation diag(w)^(-1/2) * v.T turns
    # the covariance matrix into identity one.  To find the base shifts,
    # apply the inverse of that transformation to an identity matrix.
    w, v = np.linalg.eigh(covariance)
    a = np.dot(v, np.diag(np.sqrt(w)))
    base_shifts = [a[:, i] for i in range(a.shape[1])]


    fig = plt.figure()
    fig.patch.set_alpha(0)
    axes = fig.add_subplot(111)

    pt = np.geomspace(30., 1600., num=100)

    for i, shift in enumerate(base_shifts):
        corr_function.set_params(nominal_params + shift)
        up = corr_function(pt)

        corr_function.set_params(nominal_params - shift)
        down = corr_function(pt)

        axes.fill_between(
            pt, down, up,
            color=mpl.colors.to_rgba('C{}'.format(i), 0.3), lw=0
        )

    corr_function.set_params(nominal_params)
    corr_nominal = corr_function(pt)
    axes.plot(pt, corr_nominal, color='black')

    axes.margins(x=0.)
    axes.set_xscale('log')
    axes.xaxis.set_major_formatter(mpl.ticker.LogFormatter())
    axes.xaxis.set_minor_formatter(
        mpl.ticker.LogFormatter(minor_thresholds=(2, 0.4))
    )

    axes.set_xlabel(r'$\tau_{1}^\mathrm{L2}$ [GeV]')
    axes.set_ylabel('L3Res correction')

    axes.text(
        1., 1.002, '{}, {}'.format(
            method_label, config.get_period_label(args.period)
        ),
        ha='right', va='bottom', transform=axes.transAxes
    )

    fig.savefig(args.output)
    plt.close(fig)

