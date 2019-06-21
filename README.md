# Residual jet corrections

This package implements the fit of residual (i.e. p<sub>T</sub>-dependent) jet corrections, providing an alternative to the [standard global fit](https://github.com/miquork/jecsys). It is focused on the measurement in the multijet topology, although it is generic enough to perform a simultaneous fit using inputs from different analyses. Placeholders for Z+jet and &gamma;+jet measurements are provided.

Inputs from the multijet analysis are prepared using the setup in [this repository](https://github.com/andrey-popov/multijet-jec). The computation of the &chi;<sup>2</sup> from this analysis is implemented in class [`MultijetCrawlingBins`](include/MultijetCrawlingBins.hpp). The computation is summarized in [this talk](https://indico.cern.ch/event/780845/#16-multijet-analysis-with-craw).


## Build instructions

The code base consists of C++ and Python domains. The former one implements a fast computation of the &chi;<sup>2</sup> to be minimized in the fit, while the Python code provides a convenient interface to the loss function and results of the fit.

Recommended dependencies for the C++ part are as follows:

  * Compiler with support of C++17
  * CMake 3.11
  * Boost 1.69
  * ROOT 6.16

Slightly older versions are likely to be supported. The code of this package is expected to be compliant with C++14 standard.

Requirements for the Python part: 

  * Python 3.6
  * ROOT 6.16
  * NumPy 1.14
  * Matplotlib 2.2
  * PyYAML 3.13

All the dependencies are satisfied in environment [LCG_95apython3](http://lcginfo.cern.ch/release/95apython3/), which can be set up from `/cvmfs/` with

```sh
. /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_95apython3 x86_64-slc6-gcc8-opt
```

To build the package, execute

```sh
. ./env.sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
```


## Basic fitting

Residual jet correction can be fit from the multijet analysis alone using program [`fit`](prog/fit.cpp):

```sh
inputdir="https://aapopov.web.cern.ch/aapopov/jec_inputs"
fit --multijet $inputdir/multijet.root --balance PtBal --output fit.out
```

Providing flag `--balance MPF` will run the MPF version of the measurement. The standard two-parameter functional form is used for the correction. The results, including the fitted values for the parameters of the correction, are printed in the standard output and also saved in file `fit.out`.

The same can be achieved with a Python wrapper:

```sh
./fit.py --multijet $inputdir/multijet.root --method PtBal --period 2016BCD --output fit.json
```

The data-taking period is specified for book-keeping. By default, the standard 2-parameter correction is fitted; to use a spline correction instead, provide flag `--corr spline`.

The results obtained by `fit.py` are saved in JSON format, and this is the format expected by other scripts discussed below. Program [`jq`](https://stedolan.github.io/jq/) is useful to work with such files. In particular, multiple files with fit results can be merged by running

```sh
jq -s '.' fit1.json fit2.json > fits.json
```

Usually the executable for `jq` can just be downloaded and put under `$PATH`; there is no need to build it from source.

**Important note**: Input files from the multijet analysis typically don't have any upper cut on the p<sub>T</sub> of the leading jet, but the &chi;<sup>2</sup> in bins of p<sub>T</sub> of the leading jet becomes unreliable for underpopulated bins. These should be excluded from the fit, which can be done using method `MultijetCrawlingBins::SetPtLeadRange`. The typical threshold is 1.6&nbsp;TeV (see [here](https://github.com/andrey-popov/multijet-jec/tree/Run2/analysis#inputs-for-the-fit-of-l3res-corrections)). The corresponding selection is currently hard-coded [here](https://github.com/andrey-popov/multijet-jec-fit/blob/36c35602851a514f50fb7002fbd5b0783c5ef0b4/prog/fit.cpp#L88) for C++ and [here](https://github.com/andrey-popov/multijet-jec-fit/blob/36c35602851a514f50fb7002fbd5b0783c5ef0b4/prog/fit.cpp#L88) for Python.


## Diagnostic plots

Several scripts to produce diagnostic plots are provided. Fitted jet corrections and pre- and post-fit residuals can be plotted with

```sh
plot_correction.py fits.json --period 2016BCD --method PtBal -o fig/corr.pdf
multijet_residuals.py $inputdir/multijet.root fits.json \
    --period 2016BCD --method PtBal -o fig/residuals.pdf
```

Here `fits.json` is a file with fit results produced as described above and flags `--period` and `--method` are used to identify a specific set of results within the file (and also for labels in the plots). Running

```sh
plot_parameters.py fits.json
```

generates plots with post-fit values of parameters of interest and nuisances.

Finally, 1D and 2D scans of the &chi;<sup>2</sup> for the standard 2-parameter correction can be produced with

```sh
scan_chi2.py --multijet $inputdir/multijet.root 
    --period 2016BCD --method PtBal -o fig/scans/
```

The &chi;<sup>2</sup> is minimized with respect to all nuisance parameters. Unlike all the scripts above, running the scans takes a good portion of an hour.
