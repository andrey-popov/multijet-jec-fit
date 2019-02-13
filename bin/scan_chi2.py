#!/usr/bin/env python

"""Plots 1D and 2D chi^2 scans around the minimum."""

import argparse
import itertools
import os

import numpy as np
from scipy.special import gammaincinv

import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import jecfit


def find_bin_edges(centres):
    """Construct array of bin edges from positions of bin centres."""
    
    edges = np.empty(len(centres) + 1)
    edges[1:-1] = (centres[:-1] + centres[1:]) / 2
    
    edges[0] = centres[0] - (centres[1] - centres[0]) / 2
    edges[-1] = centres[-1] + (centres[-1] - centres[-2]) / 2
    
    return edges


if __name__ == '__main__':
    
    mpl.rc('figure', figsize=(6.0, 4.8))
    
    mpl.rc('xtick', top=True, direction='in')
    mpl.rc('ytick', right=True, direction='in')
    mpl.rc(['xtick.minor', 'ytick.minor'], visible=True)
    
    mpl.rc('lines', linewidth=1., markersize=3.)
    mpl.rc('errorbar', capsize=1.)
    
    mpl.rc('axes.formatter', limits=[-3, 4], use_mathtext=True)
    mpl.rc('axes', labelsize='large')
    
    ROOT.gROOT.SetBatch(True)
    
    
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument(
        '--multijet', help='File with inputs from multijet analysis'
    )
    arg_parser.add_argument('-m', '--method', default='PtBal', help='Computation method')
    arg_parser.add_argument('-o', '--output', default='fig', help='Directory for produced plots')
    arg_parser.add_argument('-l', '--label', help='Label for plots')
    args = arg_parser.parse_args()
    
    if not args.multijet:
        raise RuntimeError('No inputs provided.')
    
    try:
        os.makedirs(args.output)
    except FileExistsError:
        pass
    
    if args.method == 'PtBal':
        method_label = '$p_\\mathrm{T}$ balance'
    else:
        method_label = 'MPF'
    
    
    loss_func = jecfit.MultijetChi2(args.multijet, args.method, {'JER'})
    loss_func.set_pt_range(0., 1.6e3)
    fit_results = loss_func.fit()
    
    
    # Plot 1D scans along each POI
    for ivar in range(2):
        x = np.empty((101, 2))
        centre = fit_results.parameters[ivar].value
        half_window = 2. * fit_results.parameters[ivar].error
        x[:, ivar] = np.linspace(
            centre - half_window, centre + half_window, num=len(x)
        )
        x[:, 1 - ivar] = fit_results.parameters[1 - ivar].value
        
        chi2 = np.empty(len(x))
        
        for i in range(len(x)):
            chi2[i] = loss_func(x[i])
        
        fig = plt.figure()
        fig.patch.set_alpha(0.)
        axes = fig.add_subplot(111)
        
        axes.plot(x[:, ivar], chi2)
        
        axes.margins(x=0., y=0.)
        axes.set_xlabel('$\\theta_{:d}$'.format(ivar))
        axes.set_ylabel('$\\chi^2$')
        
        axes.text(
            0.05, 0.90,
            '$\\theta_{:d} = {:.3f}$'.format(
                1 - ivar, fit_results.parameters[1 - ivar].value
            ),
            ha='left', va='bottom', transform=axes.transAxes
        )
        axes.text(0., 1.003, method_label, ha='left', va='bottom', transform=axes.transAxes)
        
        if args.label:
            axes.text(1., 1.003, args.label, ha='right', va='bottom', transform=axes.transAxes)
        
        fig.savefig(os.path.join(args.output, 'scan_p{:d}.pdf'.format(ivar)))
        plt.close(fig)
    
    
    # Plot the 2D scan
    error_sf = 2.
    v = fit_results.parameters[0]
    p0_values = np.linspace(v.value - error_sf * v.error, v.value + error_sf * v.error, num=51)
    v = fit_results.parameters[1]
    p1_values = np.linspace(v.value - error_sf * v.error, v.value + error_sf * v.error, num=51)
    chi2 = np.empty((len(p0_values), len(p1_values)))
    
    for i, j in itertools.product(range(len(p0_values)), range(len(p1_values))):
        chi2[i, j] = loss_func([p0_values[i], p1_values[j]])
    
    
    fig = plt.figure()
    fig.patch.set_alpha(0.)
    axes = fig.add_subplot(111)
    colourmap = plt.get_cmap('viridis')
    
    xx, yy = np.meshgrid(find_bin_edges(p0_values), find_bin_edges(p1_values))
    image = axes.pcolormesh(xx, yy, chi2.T, cmap=colourmap)
    fig.colorbar(image, fraction=0.05, pad=0.02, label='$\\chi^2$')
    
    axes.set_xlabel('$\\theta_0$')
    axes.set_ylabel('$\\theta_1$')
    axes.text(0., 1.003, method_label, ha='left', va='bottom', transform=axes.transAxes)
    
    if args.label:
        axes.text(1., 1.003, args.label, ha='right', va='bottom', transform=axes.transAxes)
    
    
    # Mark global minimum
    axes.plot(
        [fit_results.parameters[0].value],
        [fit_results.parameters[1].value], marker='+', color='red'
    )
    
    coord_transform = axes.transData + axes.transAxes.inverted()
    min_label_pos = list(coord_transform.transform(
        [fit_results.parameters[0].value, fit_results.parameters[1].value]
    ))
    min_label_pos[0] += 0.015
    
    axes.text(
        min_label_pos[0], min_label_pos[1],
        '{:.3f}'.format(loss_func.p_value(fit_results.min_value)),
        ha='left', va='center', color='red', transform=axes.transAxes
    )
    
    
    # Draw contours for characteristic p-values
    pvalues = [0.1, 0.01]
    chi2_levels = [2 * gammaincinv(loss_func.ndf / 2, 1 - p) for p in pvalues]
    
    xx_c, yy_c = np.meshgrid(p0_values, p1_values)
    
    for i in range(len(pvalues)):
        contour = axes.contour(
            xx_c, yy_c, chi2.T, [chi2_levels[i]], colors='white', zorder=1.5
        )
        
        if chi2_levels[i] > fit_results.min_value:
            axes.clabel(contour, fmt='{:g}'.format(pvalues[i]))
    
    
    fig.savefig(os.path.join(args.output, 'scan_2d.pdf'))
    plt.close(fig)
