# Log Concave Density Estimation

Briefly records what libraries are required for building and running this program. 

## Requirements

1. C++ 17 or higher
2. Fortran compiler (gfortran)
3. Eigen library (during development, located in /usr/local/include/Eigen)
4. [Boost library](https://www.boost.org/) (during development, located in /usr/local/include/boost)
5. [GMP](https://gmplib.org/)
6. [sciplot](https://sciplot.github.io/) for plotting results

## Error Rates

|Dataset|Format|Samplilng Rate|Build Time|Max Error|Min Error|
|---|---|---|---|---|---|
|books|uint64|1000|4.2317s|1,610,245|0|
|fb|uint64|1000|2.42063s|371,677|1|
|lognormal|uint64|1000|2.90579s|99,942,638|43|
|osm|uint64|1000|4.56273s|39,985,690|18|
|normal|uint64|1000|720.112s|1,766|7|
|uniform(dense)|uint64|1000|0.177492s|||
|uniform(sparse)|uint64|1000|0.872927s|6,611|0|
|wiki|uint64|1000|3.9189s|1,442,465|1|