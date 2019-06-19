# JEC fit

This is a prototype implementation of the global fit of residual jet corrections.
It is intended to be an alternative to the [standard global fit](https://github.com/miquork/jecsys).
The package is under development.

Dependencies:

  * CMake 2.8.12 or newer
  * Compiler with support of C++14
  * Boost 1.34 or newer
  * ROOT 6

The package is tested in environment [LCG_95apython3](http://lcginfo.cern.ch/release/95apython3/):

  * GCC 8.2
  * CMake 3.11
  * Boost 1.69
  * ROOT 6.16
  * Python 3.6
  * NumPy 1.14
  * Matplotlib 2.2

To build the package and run an example fit program, execute

```sh
source ./env.sh

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

inputdir="https://aapopov.web.cern.ch/aapopov/jec_inputs/prototype"
fit --photonjet-run1 $inputdir/photonjet_Run1.root --multijet-binnedsum $inputdir/multijet_BinnedSum.root --zjet-run1 $inputdir/Zjet_Run1.root
# An input file for the binned-sum version of the photon+jet analysis is also available and called "photonjet_BinnedSum.root"
```

There is also a Python script to run the fit, and it should be preferred. However, it only supports the multijet analysis.

Program [`jq`](https://stedolan.github.io/jq/) is useful to work with JSON file with fit results. In particular, multiple files can be merged by running

```sh
jq -s '.' fit1.json fit2.json > fits.json
```

Usually the executable for `jq` can just be downloaded and put under `$PATH`; there is no need to build it from source.
