# RDDENSITY
Density discontinuity testing (a.k.a. manipulation testing)
is commonly employed in regression discontinuity designs
and other program evaluation settings to detect perfect self-selection (manipulation)
around a cutoff where treatment/policy assignment changes.
This package implements manipulation testing procedures using
the local polynomial density estimators:
rddensity() to construct test statistics and p-values given a prespecified cutoff,
rdbwdensity() to perform data-driven bandwidth selection,
and rdplotdensity() to construct density plots.

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1459967](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459967), [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).

## Website

https://rdpackages.github.io/rddensity

## Queries and Requests

Please email: [rdpackages@googlegroups.com](mailto:rdpackages@googlegroups.com)

## Implementation

To install/update in Python type:
```
pip install rddensity
```
- Help: [PyPI](https://pypi.org/project/rddensity/2.2.0/), [Documentation](https://rdpackages.github.io/)
- Replication: [Python script](https://github.com/rdpackages/rddensity/master/Python/rddensity_illustration.py)

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Jansson and Ma (2018): [Manipulation Testing based on Density Discontinuity](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf).<br>
_Stata Journal_ 18(1): 234-261.

- Cattaneo, Jansson and Ma (2022): [lpdensity: Local Polynomial Density Estimation and Inference](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf).<br>
Working paper.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf).<br>
_Journal of the American Statistical Association_ 115(531): 1449-1455.<br>
[Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).

- Cattaneo, Jansson and Ma (2022): [Local Regression Distribution Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE.pdf).<br>
_Journal of Econometrics_, forthcoming.<br>
[Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE--Supplement.pdf).

<br><br>
