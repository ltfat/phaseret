## PhaseReT
PhaseReT is a Matlab/GNU Octave toolbox collecting implementations of
phase-reconstruction algorithms for complex time-frequency representations
(like STFT).

## Requirements

The toolbox depends on [LTFAT - Large Time-Frequency Analysis Toolbox](http://ltfat.github.io).
LTFAT can be installed in various ways. Please follow instructions at the
[LTFAT release page](https://github.com/ltfat/ltfat/releases/latest).

Further, some functions require additional dependencies:

* [decolbfgs](http://ltfat.github.io/phaseret/mat/decolbfgs.html) requires
[minFunc](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)

## Installation and usage (Matlab and GNU Octave)

On Matlab in Windows, it is recommended to use the pre-compiled binary release package
(phaseret-[version]-win64.zip) from
[the GitHub release page](https://github.com/ltfat/phaseret/releases).

For other configurations, download the source-only package (phaseret-[version]-src.tgz).

In order to use the toolbox it is necessary to run
```
phaseretstart;
```
from the [path_to_phaseret] directory. The command checks whether LTFAT has
already been started and whether phaseret MEX interfaces are already compiled.

In order to compile the missing mex files
run
```
phaseretmex;
```
from within Matlab/Octave command line.

## Documentation
Online documentation is available [here](http://ltfat.github.io/phaseret/doc).

## References

If you use this toolbox/library in your research, please cite

> Z. Prusa: The Phase Retrieval Toolbox, AES Int. Conf. On 
> Semantic Audio, Erlangen, Germany, June, 2017 [preprint](http://ltfat.github.io/notes/ltfatnote045.pdf).

```
@inproceedings{ltfatnote045,
 year={2017},
 month={June},
 title={{The Phase Retrieval Toolbox}},
 author={Zden\v{e}k Pr\r{u}\v{s}a},
 booktitle = {{AES} International Conference On Semantic Audio},
 address = {Erlangen, Germany}}
```

and/or relevant references found in help of the individual files.

## License
PhaseReT is distributed under terms of
[GPL3](http://www.gnu.org/licenses/gpl-3.0.en.html)
