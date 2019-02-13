from collections import namedtuple
import itertools
import os

import scipy.special
import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

_location = os.path.dirname(os.path.dirname(__file__))
ROOT.gInterpreter.AddIncludePath(os.path.join(_location, 'include'))
ROOT.gInterpreter.Declare('#include <FitBase.hpp>')
ROOT.gInterpreter.Declare('#include <JetCorrDefinitions.hpp>')
ROOT.gInterpreter.Declare('#include <MultijetCrawlingBins.hpp>')
ROOT.gInterpreter.Declare('#include <PythonWrapping.hpp>')
ROOT.gSystem.Load(os.path.join(_location, 'lib', 'libjecfit.so'))
ROOT.gSystem.Load(os.path.join(
    _location, 'lib', 'libjecfit_pythonwrapping.so')
)

JetCorrStd2P = ROOT.JetCorrStd2P
JetCorrStd2P.__doc__ = """L3Res correction with two parameters."""


class MultijetChi2:
    """Python wrapper to fit corrections using multijet data.
    
    The standard correction with two parameters is used.
    """
    
    def __init__(self, file_path, method, exclude_syst=set()):
        """Initialize from input ROOT file and method label."""
        
        if method == 'PtBal':
            method_code = 0
        elif method == 'MPF':
            method_code = 1
        else:
            raise RuntimeError('Unsupported method "{}".'.format(method))

        exclude_syst_converted = ROOT.std.set('std::string')()

        for syst in exclude_syst:
            exclude_syst_converted.insert(syst)
        
        self._nuisance_defs = ROOT.NuisanceDefinitions()
        self.measurement = ROOT.MultijetCrawlingBins(
            file_path, method_code, self._nuisance_defs,
            exclude_syst_converted
        )
        self._jet_corr = ROOT.JetCorrStd2P()
        self._loss_func = ROOT.CombLossFunction(
            self._jet_corr, self._nuisance_defs
        )
        self._loss_func.AddMeasurement(self.measurement)
    
    
    def __call__(self, params, nuisances='profile'):
        """Compute chi^2 for given values of POI and nuisances.
        
        Arguments:
            params:  array_like with values of the two parameters of the
                jet correction.
            nuisances:  array_like with values of nuisances or string
                'profile'.  In the latter case nuisance parameters are
                profiled.
        
        Return value:
            Value of chi^2.
        """
        
        if (
            nuisances == 'profile' and
            self._nuisance_defs.GetNumParams() > 0
        ):
            minimizer = self._setup_minimizer()
            
            for i in range(2):
                minimizer.SetVariableValue(i, params[i])
                minimizer.FixVariable(i)
            
            minimizer.Minimize()
            return minimizer.MinValue()
        
        else:
            x = np.zeros(self._loss_func.GetNumParams())
            x[0:2] = params
            x[2:] = nuisances
            
            return self._loss_func_wrapper(x)


    def compute_residuals(self, params, nuisances):
        """Compute data-to-simulation residuals.

        Arguments:
            params:  array_like with values of POI.
            nuisances:  dict or an array_like with values of nuisances.

        Return value:
            Tuple of NumPy arrays representing a graph with residuals.
        """

        params = np.asarray(params)
        self._jet_corr.SetParams(params)

        conv_nuisances = ROOT.Nuisances(self._nuisance_defs)

        if isinstance(nuisances, dict):
            for label, value in nuisances.items():
                conv_nuisances[label] = value
        else:
            for i in range(len(nuisances)):
                conv_nuisances[i] = nuisances[i]

        graph = self.measurement.ComputeResiduals(
            self._jet_corr, conv_nuisances
        )

        n = graph.GetN()
        x, y, yerr = np.empty(n), np.empty(n), np.empty(n)

        for i in range(n):
            x_val, y_val = ROOT.Double(), ROOT.Double()
            graph.GetPoint(i, x_val, y_val)
            x[i] = x_val
            y[i] = y_val
            yerr[i] = graph.GetErrorY(i)

        return x, y, yerr
    
    
    def fit(self, print_level=3):
        """Perform the fit with all parameters floating."""
        
        minimizer = self._setup_minimizer(print_level=print_level)
        minimizer.Minimize()
        
        return FitResults(minimizer)
    
    
    @property
    def ndf(self):
        """Number of degrees of freedom."""
        
        return self._loss_func.GetNDF()


    def p_value(self, chi2):
        """Compute p-value for given chi^2."""

        return 1 - scipy.special.gammainc(self.ndf / 2, chi2 / 2)
    
    
    def set_pt_range(self, min_pt1, max_pt1):
        """Set range in pt of the leading jet used in measurement."""
        
        self.measurement.SetPtLeadRange(min_pt1, max_pt1)
    
    
    def _setup_minimizer(self, print_level=0):
        """Create and setup a minimizer.
        
        Wrapper for the loss function is stored in self, which is needed
        to prevert it from being deleted by guarbage collection.
        Because of this, only a single minimizer can be used at a time.
        """
        
        minimizer = ROOT.Minuit2.Minuit2Minimizer()
        self._loss_func_wrapper = ROOT.WrapLossFunction(self._loss_func)
        minimizer.SetFunction(self._loss_func_wrapper)
        minimizer.SetStrategy(1)
        minimizer.SetErrorDef(1.)
        minimizer.SetPrintLevel(print_level)
        
        num_params = self._loss_func.GetNumParams()
        num_poi = num_params - self._nuisance_defs.GetNumParams()
        
        for i in range(num_poi):
            minimizer.SetVariable(i, 'p{:d}'.format(i), 0., 1e-2)
            minimizer.SetVariableLimits(i, -1., 1.)
        
        for i in range(num_poi, num_params):
            minimizer.SetVariable(
                i, self._nuisance_defs.GetName(i - num_poi), 0., 1.
            )
            minimizer.SetVariableLimits(i, -5., 5.)
        
        return minimizer


class FitResults:
    """Pythonic wrapper for fit results from Minuit2Minimizer."""
    
    Variable = namedtuple('Variable', ['name', 'value', 'error'])
    
    def __init__(self, minimizer):
        
        self.status = minimizer.Status()
        self.covariance_status = minimizer.CovMatrixStatus()

        self.min_value = minimizer.MinValue()

        self.parameters = []
        
        for p in minimizer.State().MinuitParameters():
            self.parameters.append(
                FitResults.Variable(p.Name(), p.Value(), p.Error())
            )

        num_pars = len(self.parameters)
        self.covariance_matrix = np.empty((num_pars, num_pars))

        for i, j in itertools.product(range(num_pars), range(num_pars)):
            self.covariance_matrix[i, j] = minimizer.CovMatrix(i, j)


    def serialize(self):
        """Convert to a plain dictionary to store in a JSON file."""

        serialized_parameters = []

        for p in self.parameters:
            serialized_parameters.append({
                'name': p.name,
                'value': p.value,
                'error': p.error
            })

        return {
            'status': self.status,
            'covariance_status': self.covariance_status,
            'min_value': self.min_value,
            'parameters': serialized_parameters,
            'covariance_matrix': self.covariance_matrix.tolist()
        }

