rdplotdensity
=============

Density plotting for mainpulation testing

Description
-----------
``rdplotdensity`` constructs density plots.
It is based on the local polunomial density estimator proposed in Cattaneo, Jansson and Ma (2020, 2021a).

Companion commands: ``rdbwdensity`` for bandwidth selection and ``rddensity`` for estimation.

Details
-------
Bias correction is only used for the construction of confidence intervals/bands, but not for point
estimation. The point estimates, denoted by *f_p*, are constructed using local polynomial estimates
of order *p*, while the centering of the confidence intervals/bands, denoted by *f_q*, are constructed
using local polynomial estimates of order *q*. The confidence intervals/bands take the form:
*[f_q - cv * SE(f_q) , f_q + cv * SE(f_q)]*, where *cv* denotes the appropriate critical value and *SE(f_q)*
denotes an standard error estimate for the centering of the confidence interval/band. As a result,
the confidence intervals/bands may not be centered at the point estimates because they have been bias-corrected.
Setting *q* and *p* to be equal results on centered at the point estimate confidence intervals/bands,
but requires undersmoothing for valid inference (i.e., (I)MSE-optimal bandwdith for the density point estimator
cannot be used). Hence the bandwidth would need to be specified manually when *q=p*, and the
point estimates will not be (I)MSE optimal. See Cattaneo, Jansson and Ma (2020a, 2020b) for details, and also
Calonico, Cattaneo, and Farrell (2018, 2020) for robust bias correction methods.

Sometimes the density point estimates may lie outside of the confidence intervals/bands, which can happen
if the underlying distribution exhibits high curvature at some evaluation point(s). One possible solution
in this case is to increase the polynomial order *p* or to employ a smaller bandwidth.

References
----------
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. 
`On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference <https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf>`_
*Journal of the American Statistical Association*, 113(522): 767-779.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020. 
`Coverage Error Optimal Confidence Intervals for Local Polynomial Regression <https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2020_CEopt.pdf>`_
. Working paper.

Cattaneo M. D., M. Jansson, and X. Ma. 2018.
`Manipulation Testing based on Density Discontinuity. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf>`_
*Stata Journal* 18(1): 234-261.

Cattaneo, M. D., M. Jansson, and X. Ma. 2020.
`Simple Local Polynomial Density Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf>`_.
*Journal of the American Statistical Association*, 115(531): 1449-1455.

Cattaneo, M. D., M. Jansson, and X. Ma. 2021a.
`Local Regression Distribution Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2021_JoE.pdf>`_
*Journal of Econometrics*, forthcoming.

Cattaneo M. D., M. Jansson, and X. Ma. 2022.
`lpdensity: Local Polynomial Density Estimation and Inference. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf>`_
*Journal of Statistical Software* Forthcoming.


Authors
-------
Matias D. Cattaneo, Princeton University. (cattaneo@princeton.edu).

Rajita Chandak (maintainer), Princeton University. (rchandak@princeton.edu). 

Michael Jansson, University of California Berkeley. (mjansson@econ.berkeley.edu).

Xinwei Ma (maintainer), University of California San Diego. (x1ma@ucsd.edu).

.. automodule:: rddensity.rdplotdensity
   :members:
   :undoc-members:


Example
-------
>>> import numpy as np
>>> from rddensity import rddensity
>>> from rddensity import rdplotdensity
>>> data = np.random.normal(-0.5,1,2000)
>>> rdd = rddensity(X=data, vce="jackknife")
>>> plot1 = rdplotdensity(rdd, data)
