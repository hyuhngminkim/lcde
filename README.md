# Log Concave Density Index

LCDI is a new [learned index](https://dl.acm.org/doi/pdf/10.1145/3183713.3196909) that uses nonparametric density estimation. 
At its core, LCDI uses the [CNMLCD](https://onlinelibrary.wiley.com/doi/pdf/10.1111/anzs.12232), a fast algorithm for estimating the density of an univariate random variable. 

## Requirements

1. C++ 17 or higher
2. Fortran compiler (gfortran)
3. Eigen library (during development, located in /usr/local/include/Eigen)
4. [Boost library](https://www.boost.org/) (during development, located in /usr/local/include/boost)
5. [GMP](https://gmplib.org/)
6. [sciplot](https://sciplot.github.io/) for plotting results

## How to Plot the Datasets

- We assume the necessary dataset are stored in `./data`.
- `./scripts/prepare.sh` compiles and builds the plotting script. 
- `./scripts/execute.sh` runs the plotting code. Using the `-s` flag saves the resulting plot at `./plot`. By default, the plot is saved as a "pdf" file. You can specify the file extension by setting the `-t` flag to whichever extension you like, such as "png" or "jpg".