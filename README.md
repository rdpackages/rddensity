# Manipulation Testing Using Local Polynomial Density Methods

The `rddensity` package implements manipulation testing procedures for regression discontinuity designs and other settings where researchers want to test for sorting around a cutoff.

- `rddensity`: density discontinuity testing with local polynomial density estimators.
- `rdbwdensity`: data-driven bandwidth selection for manipulation testing.
- `rdplotdensity`: density plots with histogram and confidence-band options.

## Python Implementation

To install/update in Python type:
```
pip install rddensity
```

- Help: [PyPI repository](https://pypi.org/project/rddensity/).

- Replication: [rddensity illustration](Python/rddensity_illustration.py), [rdplotdensity illustration](Python/rddensity_plot_illustration.py), [senate data](Python/rddensity_senate.csv).

## R Implementation

To install/update in R type:
```
install.packages('rddensity')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rddensity/rddensity.pdf), [CRAN repository](https://cran.r-project.org/package=rddensity).

- Replication: [rddensity illustration](R/rddensity_illustration.R), [rdplotdensity illustration](R/rddensity_plot_illustration.R), [senate data](R/rddensity_senate.csv).

## Stata Implementation

To install/update in Stata type:
```
net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/main/stata) replace
```

- Help: [rddensity](stata/rddensity.pdf), [rdbwdensity](stata/rdbwdensity.pdf).

- Replication: [rddensity illustration](stata/rddensity_illustration.do), [rdplotdensity illustration](stata/rddensity_plot_illustration.do), [senate data](stata/rddensity_senate.dta).

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

## Funding

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1459967](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459967), [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).
