### RDDENSITY

The **rddensity** package provides Stata and R implementations of manipulation tests employing local polynomial density estimation methods. This method is useful for falsification of Regression Discontinuity Designs, as well as for testing for self-selection or sorting in other contexts. This implementation provides hypothesis tests and bandwidth selectors for manipulation testing. 

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1459967](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459967), [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).

## Major Upgrades

This package was first released in Spring 2017, and had one major upgrade in Summer 2020.

- _Summer 2020 new features include_: (i) speed improvements; (ii) improved integration with lpdensity; (iii) mass points in running variable adjustments; (iv) bandwidth selection adjustments for too few mass points in and/or overshooting of the support of the running variable; (v) density discontinuity plots with histogram and/or confidence bands; and (vi) binomial testing near cutoff as complementary discontinuity testing following results in rdlocrand methods (see references there for details).

## Stata Implementation

## R Implementation

## References

For overviews and introductions, see [rdpackages website]().

### Software and Implementation

- Cattaneo, Jansson and Ma (2018): [Manipulation Testing based on Density Discontinuity](references/Cattaneo-Jansson-Ma_2018_Stata.pdf), Stata Journal 18(1): 234-261.
- Cattaneo, Jansson and Ma (2020): [lpdensity: Local Polynomial Density Estimation and Inference.](references/Cattaneo-Jansson-Ma_2020_JSS.pdf), working paper.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](references/Cattaneo-Jansson-Ma_2020_JASA.pdf), Journal of the American Statistical Association, forthcoming. [Supplemental appendix](references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).
- Cattaneo, Jansson and Ma (2020): [Local Regression Distribution Estimators](references/Cattaneo-Jansson-Ma_2020_JoE.pdf). [Supplemental appendix](references/Cattaneo-Jansson-Ma_2020_JoE--Supplement.pdf).

<br>
<br>
