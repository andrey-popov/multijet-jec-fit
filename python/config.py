import os
import re

import yaml


class Config:
    """An interface to access configuration."""

    def __init__(self, path):
        """Initialize from a YAML configuration file."""

        if not os.path.exists(path):
            # Try to resolve the path with respect to standard location
            if 'JEC_FIT_BASE' not in os.environ:
                raise RuntimeError(
                    'File "{}" not found when interpreted literally '
                    'and no standard location is provided in environment '
                    'variable "JEC_FIT_BASE".'.format(path)
                )
            else:
                try_path = os.path.join(os.environ['JEC_FIT_BASE'], path)

                if os.path.exists(try_path):
                    path = try_path
                else:
                    raise RuntimeError('File "{}" not found.'.format(path))

        with open(path) as f:
            self.config = yaml.safe_load(f)


    def get_parameter_label(self, parameter):
        """Get formatted label for given parameter.
        
        If no rules are found for given parameter, return its name
        unchanged.
        """

        # Check if this is a POI
        match = re.match(r'^p(\d+)$', parameter)

        if match:
            return r'$\theta_{}$'.format(match.group(1))

        nuisance_labels = self.config['nuisances']['labels']

        if parameter in nuisance_labels:
            return nuisance_labels[parameter]
        else:
            return parameter


    def get_period_label(self, period, full=False):
        """Get formatted label for given period."""

        period_config = self.config['periods'][period]

        if full:
            return '{} {} $\\mathrm{fb}^{-1}$ (13 TeV)'.format(
                period_config['label'], period_config['lumi']
            )
        else:
            return period_config['label']


    def sort_parameters(self, parameter_names):
        """Sort parameters as specified in the configuration.

        If some of the parameters are not mentioned in the
        configuration, append them at the end.
        """

        try:
            desired_order = self.config['nuisances']['order']
        except KeyError:
            return parameter_names

        sorted_names = []

        # First include parameters mentioned in the configuration, then
        # everything remaining
        for name in desired_order:
            if name in parameter_names:
                sorted_names.append(name)

        for name in parameter_names:
            if name not in desired_order:
                sorted_names.append(name)

        return sorted_names


