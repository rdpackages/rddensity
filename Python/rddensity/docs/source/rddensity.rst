rddensity
=========
Manipulation testing using local polynomial density estimation


Description
-----------
``rddensity`` implements manipulation testing procedures using the
local polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020),
and implements graphical procedures with valid confidence bands using the results in
Cattaneo, Jansson and Ma (2021a,b).
In addition, the command provides complementary manipulation testing based on finite
sample exact binomial testing following the results in
Cattaneo, Frandsen and Titiunik (2015) and Cattaneo, Frandsen and Vazquez-Bare (2017).
For an introduction to manipulation testing see McCrary (2008).

Companion commands: ``rdbwdensity`` for bandwidth selection and ``rdplotdensity`` for plotting estimation results.

References
----------
Cattaneo M. D., B. Grandsen, and R. Titiunik. 2015
`Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate. <https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf>`_
*Journal of Causal Inference* 3(1): 1-24.

Cattaneo M. D., M. Jansson, and X. Ma. 2018.
`Manipulation Testing based on Density Discontinuity. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf>`_.
*Stata Journal* 18(1): 234-261.

Cattaneo, M. D., M. Jansson, and X. Ma. 2020.
`Simple Local Polynomial Density Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf>`_.
*Journal of the American Statistical Association*, 115(531): 1449-1455.

Cattaneo M. D., M. Jansson, and X. Ma. 2022.
`lpdensity: Local Polynomial Density Estimation and Inference. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf>`_.
*Journal of Statistical Software* Forthcoming.

Cattaneo M. D., R. Titiunik, and F. Vazquez-Bare. 2017.
`Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start of Child Mortality <https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf>`_
*Journal of Policy Analysis and Management*, 36(3): 643:681.

McCrary, J. 2008.
Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test.
*Journal of Econometrics* 142(2): 698-714.

Authors
-------
Matias D. Cattaneo, Princeton University. (cattaneo@princeton.edu).

Rajita Chandak (maintainer), Princeton University. (rchandak@princeton.edu).

Michael Jansson, University of California Berkeley. (mjansson@econ.berkeley.edu).

Xinwei Ma (maintainer), University of California San Diego. (x1ma@ucsd.edu).


.. automodule:: rddensity.rddensity
   :members:
   :undoc-members:

Example
-------
>>> import numpy as np
>>> from rddensity import rddensity
>>> data = np.random.normal(-0.5,1,2000)
>>> rdd = rddensity(X=data, vce="jackknife")
>>> print(repr(rdd))
