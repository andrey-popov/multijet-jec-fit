#!/usr/bin/env python

"""Plots fitted values for POI and nuisances."""

import argparse
import json
import os
import re

import numpy as np

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

from config import Config
from utils import mpl_style


def plot_parameters(fits, parameter_names, config, fig_name):

    periods = list({key[0] for key in fits})
    periods.sort()
    period_labels = [config.get_period_label(p) for p in periods]

    parameter_names = config.sort_parameters(parameter_names)

    fig = plt.figure(figsize=(6., 2 * len(parameter_names) / 0.8))
    fig.patch.set_alpha(0.)
    gs = mpl.gridspec.GridSpec(len(parameter_names), 1, hspace=0.)

    for iparameter, parameter_name in enumerate(parameter_names):
        axes = fig.add_subplot(gs[iparameter, 0])

        for variant, colour, x_offset, label in [
            ('PtBal', 'C0', -0.01, r'$p_\mathrm{T}$ bal.'),
            ('MPF', 'C1', 0.01, 'MPF')
        ]:
            x = np.arange(len(periods)) + x_offset * len(periods)
            y, yerr = [], []

            for period in periods:
                parameter = next(filter(
                    lambda p: p['name'] == parameter_name,
                    fits[period, variant]['parameters'])
                )
                y.append(parameter['value'])
                yerr.append(parameter['error'])

            axes.errorbar(
                x, y, yerr, label=label,
                marker='o', lw=0, elinewidth=0.8, c=colour
            )

        axes.grid(axis='y', c='black', ls='dotted')
        axes.set_ylabel(config.get_parameter_label(parameter_name))

        # Manually set position of label for the y axis so that they are
        # aligned for all parameters
        axes.get_yaxis().set_label_coords(-0.1, 0.5)

        axes.set_xticks(range(len(periods)))

        if iparameter == len(parameter_names) - 1:
            # This is the lowest pad
            axes.set_xticklabels(period_labels)
        else:
            axes.set_xticklabels([''] * len(axes.get_xticklabels()))

        if iparameter == 0:
            axes.legend()

    fig.savefig(fig_name)
    plt.close(fig)


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('fits', help='JSON file with results of fits')
    arg_parser.add_argument(
        '--fig-dir', default='fig/postfit', help='Directory to store plots'
    )
    args = arg_parser.parse_args()

    try:
        os.makedirs(args.fig_dir)
    except FileExistsError:
        pass

    config = Config('config/plot_config.yaml')

    plt.style.use(mpl_style)


    with open(args.fits) as f:
        fits = json.load(f)

    fits = {(fit['period'], fit['variant']): fit for fit in fits}


    # Separate parameters into POI and nuisances
    poi_regex = re.compile(r'^p[0-9]+$')
    pois = set()
    nuisances = set()

    for param in next(iter(fits.values()))['parameters']:
        name = param['name']

        if poi_regex.match(name):
            pois.add(name)
        else:
            nuisances.add(name)


    plot_parameters(fits, pois, config, os.path.join(args.fig_dir, 'poi.pdf'))
    plot_parameters(
        fits, nuisances, config, os.path.join(args.fig_dir, 'nuisances.pdf')
    )

