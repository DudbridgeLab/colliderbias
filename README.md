# colliderbias

`colliderbias` is a package for adjusting association statistics and performing Mendelian randomisation when the outcome is conditional on a selection variable, leading to collider bias.

This package is the most recent implementation, superceding the earlier [indexevent](https://github.com/DudbridgeLab/indexevent) and [ColliderBias](https://github.com/SiyangCai/ColliderBias) packages.

## Installation

```
devtools::install_github("DudbridgeLab/colliderbias")
```

## Usage
This package includes two functions to adjust for collider bias:

1. `colliderbias` Uses instrumental variables to adjust association statistics for collider bias

2. `colliderbias.mr` Uses instrumental variables to perform Mendelian randomisation in the presence of collider bias

In both functions, adjustment for weak instruments is implemented using corrected weighted least squares methods.
Alternatively, `colliderbias` can use the `mr.raps` package to adjust for weak instruments.

To see a simulated example:

```
help(testData)
```

## Citations
> Dudbridge F, Allen RJ, Sheehan NA _et al._ (2019) Adjustment for index event bias in genome-wide association studies of subsequent events. Nat Commun 10:1561. [https://doi.org/10.1038/s41467-019-09381-w](https://doi.org/10.1038/s41467-019-09381-w)
>
> Cai S, Hartley A, Mahmoud O _et al._ (2022) Adjusting for collider bias in genetic association studies using instrumental variable methods. Genetic Epidemiol 46:303-316. [https://doi.org/10.1002/gepi.22455](https://doi.org/10.1002/gepi.22455)
> 
> Cai S, Dudbridge F (2025) Estimating causal effects on a disease progression trait Using bivariate Mendelian Randomisation. Genet Epidemiol 49:e22600. [https://doi.org/10.1002/gepi.22600](https://doi.org/10.1002/gepi.22600)
