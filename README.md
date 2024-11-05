<div style="display: flex; align-items: center;">
  <h1>circlus: Clustering and Simulation of Spherical Cauchy and PKBD Models</h1>
  <img src="logo.png" alt="circlus logo" width="100" style="margin-right: 20px;">
</div>


**circlus** is an R package for the estimation and clustering of spherical data, seamlessly integrated with the `flexmix` package. It includes the necessary M-step implementations for both Poisson Kernel-Based Distribution (PKBD) and spherical Cauchy distribution, and provides random number generators for both distributions.

## Features
- **Clustering**: Supports clustering with spherical Cauchy and PKBD models.
- **Random Number Generation**: Tools for generating random numbers from PKBD and spherical Cauchy distributions.
- **Integration with `flexmix`**: Easily integrates with the `flexmix` package to support flexible, model-based clustering workflows.

## Installation

You can install the latest release of **circlus** from CRAN:

```r
install.packages("circlus")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("lsablica/circlus")
```

## Usage

Load the **circlus** package and explore its core functionalities:

```r
library(circlus)
```

## Documentation

For a full list of functions and their usage, refer to the [circlus reference manual](https://CRAN.R-project.org/package=circlus).

## Dependencies

**circlus** requires:
- R (≥ 3.1.0)
- `Rcpp` (≥ 0.12.18)
- `Tinflex` (≥ 1.8)
- `flexmix`
- `torch`
- `methods`
  
It also links to `Rcpp` and `RcppArmadillo`.

## License

**circlus** is licensed under the GPL-3 license.

## Links

- [CRAN Package Page](https://CRAN.R-project.org/package=circlus)

