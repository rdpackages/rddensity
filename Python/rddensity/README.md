# RDDENSITY

The `rddensity` package provides tools for manipulation testing around a cutoff using local polynomial density estimators, bandwidth selection, and density plotting.

- `rddensity`: implements density discontinuity testing.
- `rdbwdensity`: implements bandwidth selectors for manipulation testing.
- `rdplotdensity`: implements density plots around the cutoff.

See Cattaneo, Jansson and Ma (2018, 2020, 2022, 2024) for references.

Website: [https://rdpackages.github.io/](https://rdpackages.github.io/).

Source code: [https://github.com/rdpackages/rddensity](https://github.com/rdpackages/rddensity).

## Authors

Matias D. Cattaneo, maintainer (<matias.d.cattaneo@gmail.com>)

Michael Jansson (<michael.jansson.berkeley@gmail.com>)

Xinwei Ma (<xinweima.pku@gmail.com>)

Ricardo Masini, Python contributor (<ricardo.masini@gmail.com>)

## Installation

To install/update use pip:
```
pip install rddensity
```

## Usage
```python
from rddensity import rddensity, rdbwdensity, rdplotdensity
```

- Replication: [rddensity illustration](https://github.com/rdpackages/rddensity/blob/main/Python/rddensity_illustration.py), [rdplotdensity illustration](https://github.com/rdpackages/rddensity/blob/main/Python/rddensity_plot_illustration.py), [senate data](https://github.com/rdpackages/rddensity/blob/main/Python/rddensity_senate.csv).

## Dependencies

- numpy
- pandas
- scipy
- sympy
- plotnine
- lpdensity

## References

For overviews and introductions, see the [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Jansson and Ma (2018): [Manipulation Testing based on Density Discontinuity](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf).<br>
_Stata Journal_ 18(1): 234-261.

- Cattaneo, Jansson and Ma (2022): [lpdensity: Local Polynomial Density Estimation and Inference](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf).<br>
_Journal of Statistical Software_ 101(2): 1-25.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf).<br>
_Journal of the American Statistical Association_ 115(531): 1449-1455.<br>
[Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).

- Cattaneo, Jansson and Ma (2024): [Local Regression Distribution Estimators](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2024_JoE.pdf).<br>
_Journal of Econometrics_ 240(2): 105074.<br>
[Supplemental appendix](https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2024_JoE--Supplement.pdf).
