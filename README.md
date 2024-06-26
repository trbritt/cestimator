# Cestimator
![](./cestimator.png)

![example branch parameter](https://github.com/trbritt/cestimator/actions/workflows/.github.yml/badge.svg?branch=development)
![](https://app.travis-ci.com/trbritt/cestimator.svg?token=iXZFDWLUxwEkrS9zsVhq&branch=development)

This is a first pass at my implementation of estimators of distributions. This will highlight a few different types of estimation that encompass a large class of realistic distribution, specifically elliptically distributed realisations of a random variable.

The goal is to show how these estimations perform for different regimes of data quality and covariance. 

The types of estimation that are covered are (i) non-parametric estimation (the naïve approach), (ii) shrinkage estimation, (iii) maximum likelihood estimation, and (iv) robust estimation, to discuss the impact of missing data or to cover the cases where the test distributions used in (i-iii) do not really cover the true distribution.

Details on these estimators can be found in the [docs](https://github.com/trbritt/cestimator/tree/master/docs).


## Compilation

To build the program requires only the `Eigen3` library as a prerequisite. The checks for these headers are done by the `configure` script, generated using `autoconf`. In order to compile the project, start by running `./configure --help` to see the possible build options. Notably, there is a flag `--with-visualization` that will additionally link the project against the MathGL suite, and provide visualization of the resulting analysis (in progress). Simply run:

```bash
./configure <options>
```
to generate a `Makefile`. Then run `make` to actually compile the project. All object files and the executable will be created in the `./bin` directory. 

## Usage

The programme as built is designed to take in a CSV file, where a column is the time history of the variable of interest. The CSV should have multiple columns corresponding to the quantities the user believes to be codistributed. In this case, the programme will analyse the $N$-dimensional distribution, where $N$ is the number of columns.

The programme builds a library `libcestimator.so` in the `lib` directory, that can be used to build custom scripts (don't forget to include the header `src/cestimators.hpp`!). Examples are given in the [examples](https://github.com/trbritt/cestimator/tree/master/examples) directory. These examples are also built by the make script via the target `examples`.


## Output

The output is simply formatted as text to stdout for now, but will soon include more common formats for speed.

## License

This code is distributed under the GPLv3 license.