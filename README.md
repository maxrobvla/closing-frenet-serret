The code in this repository was written as part of an internship with the group for stellerator theory of the Max Planck institute for plasma physics in Greifswald mainly by me (Max-Robin Schultz). It is meant to efficiently generate curves by a method similar to the one described by G.G. Plunk and E. Rodríguez ([link to the preprint](https://doi.org/10.48550/arXiv.2508.12820)). It uses algorithms for numeric integration of ordinary differential equations provided by the [Boost C++ library (Boost.Odeint)](https://www.boost.org/library/latest/numericodeint/) and optimization algorithms provided by [SciPy](https://scipy.org/).

For more information of the structure and functionality of the code see in `DOCUMENTATION.md`

# Building

use the python build script `make` to build the solver

## CLI-arguments for the build script


`--help`: displays all implemented command line arguments

### Flags

`--enable-gradient-calculation`: enable calculation of gradient information (will possibly slow down the program, even if it is not explicitly used)

`--default-generator`: script uses `ninja` generator by default, this flag switches to the system default (if `ninja` is not found in PATH, the script switches to the system default automatically)

`--debug`: enables debug symbols

### Arguments

`--target`: choose the function of the script

- "all": (default) generates build files and builds python package
- "cmake": generates only build files without building the python package
- "rebuild": rebuilds the python package using already present build files (without (re-)generating the build files)
- "clear": deletes build directory

# License
Gnu General Public License 3.0 (or later)

see `Licence.txt`
