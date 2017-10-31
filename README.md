# JEC fit

This is a prototype implementation of the global fit of residual jet corrections.
It is intended to be an alternative to the [standard global fit](https://github.com/miquork/jecsys).
The package is under development.

Dependencies:
  * CMake 2.8.12 or newer
  * Compiler with support of C++14
  * ROOT 6

To build the package and run an example fit program, execute
```bash
cmake .
make
inputdir="https://aapopov.web.cern.ch/aapopov/jec_inputs/prototype"
bin/fit $inputdir/photonjet_Run1.root $inputdir/multijet_BinnedSum.root $inputdir/Zjet_Run1.root
# The input file for the binned sum version of photon+jet is called photonjet_BinnedSum.root
```
