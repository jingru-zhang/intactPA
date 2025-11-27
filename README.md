## Overview

The **intactPA** package implements the **INTACT** method for harmonizing longitudinal physical activity data collected from multiple sources. It provides functions for simulation, harmonization, and evaluation of multilevel covariance structures.

## Installation

You can install the development version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install intactPA from GitHub
devtools::install_github("jingru-zhang/intactPA", build_vignettes = TRUE)
```

## Example

Here is a minimal example showing how to run **INTACT** on example data:

```r
library(intactPA)

# Load example dataset
data(example)

# Run INTACT
result <- INTACT(mydf0, Yraw, formula, needlog = 0, comparisons = 1, act = 1, wpve = 0.5, dd = NULL)
```

For more details and simulation studies, see the package vignette:
```r
vignette("intactPA", package = "intactPA")
```

## Documentation

- Reference manual and function documentation: available via `?INTACT` after installation.  
- Vignettes: introductory examples and simulation workflows.


## Acknowledgements

This package was supported by the National Natural Science Foundation of China (NSFC grant 12401388) and the Shanghai Pujiang Program (grant 23PJ1401100).

## Citation

Zhang, J., Cui, E., Li, H., & Shou, H. (2025). INTACT: A method for integration of longitudinal physical activity data from multiple sources. Biometrics.




