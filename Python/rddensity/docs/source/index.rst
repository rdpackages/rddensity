.. rddensity documentation master file, created by
   sphinx-quickstart on Tue May 17 15:17:04 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

rddensity Homepage
=====================================
Density discontinuity testing (a.k.a. manipulation testing)
is commonly employed in regression discontinuity designs
and other program evaluation settings to detect perfect self-selection (manipulation)
around a cutoff where treatment/policy assignment changes.
This package implements manipulation testing procedures using
the local polynomial density estimators:
``rddensity()`` to construct test statistics and p-values given a prespecified cutoff,
``rdbwdensity()`` to perform data-driven bandwidth selection,
and ``rdplotdensity()`` to construct density plots.

Companion ``Stata`` aand ``R`` packages
and additional regression discontinuity packages
are available at `rdpackages.github.io <https://rdpackages.github.io/>`_.

Install ``rddensity`` by running
``pip install rddensity``.

Import functions by running

>>> from rddensity import rddensity

>>> from rddensity import rdbwdensity

>>> from rddensity import rdplotdensity

Source code and replication files are available in the `rddensity repository <https://github.com/rdpackages/rddensity/>`_.

References
----------
Cattaneo M. D., M. Jansson, and X. Ma. 2018.
`Manipulation Testing based on Density Discontinuity. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf>`_.
*Stata Journal* 18(1): 234-261.

Cattaneo M. D., M. Jansson, and X. Ma. 2022.
`lpdensity: Local Polynomial Density Estimation and Inference. <https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf>`_.
*Journal of Statistical Software* Forthcoming.


Authors
-------
Matias D. Cattaneo, Princeton University. (cattaneo@princeton.edu).

Rajita Chandak (maintainer), Princeton University. (rchandak@princeton.edu). 

Michael Jansson, University of California Berkeley. (mjansson@econ.berkeley.edu).

Xinwei Ma (maintainer), University of California San Diego. (x1ma@ucsd.edu).



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules


