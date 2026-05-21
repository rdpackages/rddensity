rdbwdensity
===========

Bandwidth selection for manipulation testing


Description
-----------
``rdbwdensity`` implements several data-driven bandwidth selection methods useful to construct manipulation testing procedures using the local polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020).

Related ``Stata`` and ``R`` useful for inference in regression discontinuity (RD) designs are available on the `rdpackages website <https://rdpackages.github.io>`_.

Companion commands: ``rddensity`` for estimation and ``rdplotdensity`` for plotting estimation results.

References
----------
Cattaneo M. D., M. Jansson, and X. Ma. 2018.
`Manipulation Testing based on Density Discontinuity. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf>`_.
*Stata Journal* 18(1): 234-261.

Cattaneo M. D., M. Jansson, and X. Ma. 2020.
`Simple Local Polynomial Density Estimators <https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf>`_.
*Journal of the American Statistical Association*, 115(531): 1449-1455.


Authors
-------
Matias D. Cattaneo (maintainer), Princeton University. (matias.d.cattaneo@gmail.com).

Ricardo Masini, Python contributor. (ricardo.masini@gmail.com).

Rocio Titiunik, Princeton University. (rocio.titiunik@gmail.com).

Gonzalo Vazquez-Bare, UC Santa Barbara. (gvazquezbare@gmail.com).



.. automodule:: rddensity.rdbwdensity
   :members:
   :undoc-members:

Example
-------
>>> import numpy as np
>>> from rddensity import rdbwdensity
>>> data = np.random.normal(-0.5,1,2000)
>>> est = rdbwdensity(data=data, vce="jackknife")
>>> print(repr(est))
