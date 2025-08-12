![HFB3 logo](misc/imgs/hfb3.png)

**Version:** 1.0.12  
**Authors:** N. Dubray, J.-P. Ebran, P. Carpentier, M. Frosini, A. Zdeb, N. Pillet, J. Newsome, M. Verrière, G. Accorto, D. Regnier.  
**Copyright:** CEA, DAM F-91297 Arpajon, France. All rights reserved.  
**License:** GNU GENERAL PUBLIC LICENSE v3 or higher.

This library solves the constrained Hartree-Fock-Bogoliubov equations using a cylindrical one- or two-center basis.
To see the underlying formalism, browse the generated documentation.

# How to use ?

There are several possible ways to use `HFB3`:

1. using the main command-line program (cf. `examples/cli_example.sh`)

```shell
$ cd examples; ./cli_example.sh
```

2. using a Python script (cf. `examples/16O_groundstate.py`)

```Python
import hfb3
[...]
```

3. using custom compiled `C++` code (cf. `src/main.cpp`)

```C++
#include "hfb3.h"
[...]
// compile with `-lhfb3`
```

4. using a Python `Jupyter` Notebook (cf. `examples/demo.ipynb`)

```shell
$ # dependencies
$ pip install jupyter notebook matplotlib bokeh
$ jupyter notebook examples/demo.ipynb
```

# Compilation

## Dependencies
* Main library
  * `msgpack` (included in `deps/`) [https://msgpack.org](https://msgpack.org)
  * `armadillo` (not included) [http://arma.sourceforge.net](http://arma.sourceforge.net)
* Documentation
  * `doxygen` >= 1.8.13 (not included) [http://www.stack.nl/~dimitri/doxygen](http://www.stack.nl/~dimitri/doxygen)
  * `latex` >= 3.14159265 (not included) [http://www.tug.org/texlive](http://www.tug.org/texlive)
* Unit tests
  * `GoogleTest` >= 1.10.x (included in `deps/`) [https://github.com/google/googletest](https://github.com/google/googletest)
  * `tox` (not included) [https://tox.wiki/en/latest/](https://tox.wiki/en/latest/)
* Python bindings
  * `swig` >= 4.0.0 (not included) [http://www.swig.org](http://www.swig.org)
  * `boost` (not included) [https://www.boost.org](https://www.boost.org)
* Visualization
  * `bokeh` (not included) [https://bokeh.org](https://www.bokeh.org)
  * `matplotlib` (not included) [https://matplotlib.org](https://matplotlib.org)
* Gauss quadratures generation
  * `Maxima` (not included) [http://maxima.sourceforge.net](http://maxima.sourceforge.net)

## Targets
```shell
$ make hfb3                          # generate 'bin/hfb3' (main binary)
$ make python_install                # generate and install the Python bindings
```

# Documentation (HTML)
Generate the documentation with
```shell
$ make doc
```

Then open `doc/html/index.html` with a browser.

# Python bindings

The simplest way yo use the Python bindings is to install a pre-compiled version of HFB3 with `pip`:

```shell
$ pip install hfb3
```

To re-generate the Python bindings:

```shell
$ make python_install
```

# For developers

Some useful commands if you want to tweak or contribute to HFB3.

```shell
$ make clean                                              # clean the tree
$ bin/run_tests.sh 8                                      # build and run the C++ unit tests using 8 processes
$ bin/run_tests_python.sh 8                               # build and run the Python tests using 8 processes
$ cd misc/maxima; maxima --very-quiet < tests.max         # calculate some maxima test constants
$ cd misc/maxima; maxima --very-quiet < quadratures.max   # generate quadrature constants
```
