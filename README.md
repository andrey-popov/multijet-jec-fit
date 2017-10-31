# JEC fit

This is a prototype implementation of the global fit of residual jet corrections.
It is intended to be an alternative to the [standard global fit](https://github.com/miquork/jecsys).
The package is under development.

Dependencies:
  * CMake 2.8.12 or newer
  * Compiler with support of C++14
  * Boost 1.34 or newer
  * ROOT 6

To build the package and run an example fit program, execute
```bash
cmake .
make
inputdir="https://aapopov.web.cern.ch/aapopov/jec_inputs/prototype"
bin/fit --photonjet-run1 $inputdir/photonjet_Run1.root --multijet-binnedsum $inputdir/multijet_BinnedSum.root --zjet-run1 $inputdir/Zjet_Run1.root
# An input file for the binned-sum version of the photon+jet analysis is also available and called "photonjet_BinnedSum.root"
```
