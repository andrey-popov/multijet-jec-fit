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
bin/fit "" /gridgroup/cms/popov/Analyses/JetMET/2017.09.07_New-method-real-setup/UpdatedInputs/multijet_Run2016H.root
```
