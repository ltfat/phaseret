# PhaseReT
PhaseReT is a Matlab/GNU Octave toolbox collecting implementations of
phase-reconstruction algorithms for complex time-frequency representations
(like STFT).

## Requirements

The toolbox depends on [LTFAT - Large Time-Frequency Analysis Toolbox](http://ltfat.github.io).

Further, some functions require additional dependencies:

* [decolbfgs](http://ltfat.github.io/phaseret/mat/decolbfgs.html) requires
[minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)

## Instalation and usage

In order to use the toolbox it is necessary to run
```
phaseretstart;
```
from the [path_to_phaseret] directory. The command checks whether LTFAT has
already been started and whether MEX interfaces are already compiled.

If you are using the pre-compiled Windows package, the next step is not necessary.

The mex interfaces can be compiled by running
```
phaseretmex;
```
from within Matlab/Octave command line. Note that
the [FFTW library](http://fftw.org/) has to be installed
by e.g. `apt-get install libfftw3-dev` on Debian-based systems.
On Windows, this has to be done manually by downloading the archive from [here](http://fftw.org/install/windows.html) and placing
libfftw3-3.dll and libfftw3f-3.dll to the [path_to_phaseret]/mex subdirectory).

To start the toolbox automatically add the following lines
```
addpath(path_to_ltfat);
addpath(path_to_phaseret);
ltfatstart;
phaseretstart;
```
to your startup script [Matlab](http://de.mathworks.com/help/matlab/ref/startup.html)
[Octave](https://www.gnu.org/software/octave/doc/interpreter/Startup-Files.html).

Note that pre-compiled packages for Windows can be downloaded from
[the GitHub release page](https://github.com/ltfat/phaseret/releases).

## Documentation
Online documentation is available [here](http://ltfat.github.io/phaseret/doc).

# References

If you use this toolbox/library in your research, please cite

> Zdenek Prusa, Peter Soendergaard: TBD.

and/or relevant references found in help of the individual files.

# License
PhaseReT is distributed under terms of
[GPL3](http://www.gnu.org/licenses/gpl-3.0.en.html)
