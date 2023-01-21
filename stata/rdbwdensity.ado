********************************************************************************
* RDDENSITY STATA PACKAGE -- rdbwdensity
* Authors: Matias D. Cattaneo, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 2.4 2023-01-21

capture program drop rdbwdensity

program define rdbwdensity, eclass
syntax varlist(max=1) [if] [in] [, 	///
  C(real 0) 						///
  P(integer 2) 						///
  KERnel(string) 					///
  FITselect(string) 				///
  VCE(string)						///
  noREGularize 			    		///
  NLOCalmin (integer -1)			///
  NUNIquemin (integer -1)			///
  noMASSpoints						///
  ]
	
	marksample touse

	if ("`kernel'"=="") local kernel = "triangular"
	local kernel = lower("`kernel'")
	if ("`fitselect'"=="") local fitselect = "unrestricted"
	local fitselect = lower("`fitselect'")
	if ("`vce'"=="") local vce = "jackknife"
	local vce = lower("`vce'")

	preserve
	qui keep if `touse'

	local x "`varlist'"

	qui drop if `x'==.
	
	qui su `x'
	local x_min = r(min)
	local x_max = r(max)
	local N = r(N)

	qui su `x' if `x'<`c'
	local xl_min = r(min)
	local xl_max = r(max)
	local Nl = r(N)

	qui su `x' if `x'>=`c'
	local xr_min = r(min)
	local xr_max = r(max)
	local Nr = r(N)
	
	****************************************************************************
	*** BEGIN ERROR HANDLING *************************************************** 
	if (`c'<=`x_min' | `c'>=`x_max'){
		di "{err}{cmd:c()} should be set within the range of `x'."  
		exit 125
	}
	
	if (`Nl'<10 | `Nr'<10){
		di "{err}Not enough observations to perform calculations."  
		exit 2001
	}
	
	if (`p'!=1 & `p'!=2 & `p'!=3 & `p'!=4 & `p'!=5 & `p'!=6 & `p'!=7){
		di "{err}{cmd:p()} should be an integer value less or equal than 7."  
		exit 125
	}
		
	if ("`kernel'"!="uniform" & "`kernel'"!="triangular" & "`kernel'"!="epanechnikov"){
		di "{err}{cmd:kernel()} incorrectly specified."  
		exit 7
	}

	if ("`fitselect'"!="restricted" & "`fitselect'"!="unrestricted"){
		di "{err}{cmd:fitselect()} incorrectly specified."  
		exit 7
	}

	if ("`vce'"!="jackknife" & "`vce'"!="plugin"){ 
		di "{err}{cmd:vce()} incorrectly specified."  
		exit 7
	}
	
	if ("`regularize'" == "") {
		local regularize = 1
	}
	else {
		local regularize = 0
	}

	if ("`masspoints'" == "") {
		local masspoints = 1
	}
	else {
		local masspoints = 0
	}

	if (`nlocalmin' < 0) {
		local nlocalmin = 20 + `p' + 1
	}

	if (`nuniquemin' < 0) {
		local nuniquemin = 20 + `p' + 1
	}
	*** END ERROR HANDLING ***************************************************** 
	****************************************************************************

	qui replace `x' = `x'-`c'
	qui sort `x'

	****************************************************************************
	*** BEGIN MATA ESTIMATION ************************************************** 
	mata{
	*display("got here!")
	X = st_data(.,("`x'"), 0);
	
	XUnique   	= rddensity_unique(X)
	freqUnique  = XUnique[., 2]
	indexUnique = XUnique[., 4]
	XUnique     = XUnique[., 1]
	NUnique     = length(XUnique)
	NlUnique    = sum(XUnique :<  0)
	NrUnique    = sum(XUnique :>= 0)
	
	masspoints_flag = sum(freqUnique :!= 1) > 0 & `masspoints'
	st_numscalar("masspoints_flag", masspoints_flag)

	****************************************************************************
	** Kernel Constants
	****************************************************************************
	if ("`fitselect'"=="unrestricted") {
		if ("`kernel'"=="uniform") {
		Bsq_p=(0.24999999999999966693,0.01000000000000004878,0.00014172335600917503246,0.00000098418997230168060921,0.0000000039855627124297920874,0.000000000010481435883708594505,0.000000000000019251413808407223054,0.000000000000000026041096146069883723)
		}
		else if ("`kernel'"=="triangular") {
		Bsq_p=(0.15999999999999992006,0.0051020408163267062795,0.00006298815822632821272,0.00000039855626977185269366,0.0000000015093255191922687787,0.0000000000037733140674455142929,0.0000000000000066614606382066531783,0.00000000000000000871923521295076001)
		}
		else if ("`kernel'"=="epanechnikov") {
		Bsq_p=(0.17728531855955703689,0.0059878117913833833058,0.000076742107318398123782,0.00000049855793475223530487,0.0000000019253854002580922299,0.0000000000048868584327480008077,0.0000000000000087317551910484345913,0.000000000000000011557177676615075784)
		}
	}
	else if ("`fitselect'"=="restricted") {
		if ("`kernel'"=="uniform") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.1666666667,0.25,0.125,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667\0,0.25,0.5,0.1666666667,0.125,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545\0,0.125,0.1666666667,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846\0,0.1,0.125,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571\0,0.08333333333,0.1,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333\0,0.07142857143,0.08333333333,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125\0,0.0625,0.07142857143,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471\0,0.05555555556,0.0625,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778\0,0.05,0.05555555556,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947\0,0.04545454545,0.05,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947,0.025\0,0.04166666667,0.04545454545,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947,0.025,0.02380952381)
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.03333333333,0.05208333333,0.02430555556,0.01904761905,0.015625,0.01322751323,0.01145833333,0.0101010101,0.009027777778,0.008158508159,0.00744047619\0,0.05208333333,0.08333333333,0.0375,0.02916666667,0.02380952381,0.02008928571,0.01736111111,0.01527777778,0.01363636364,0.01231060606,0.01121794872\0,0.02430555556,0.0375,0.01785714286,0.0140625,0.01157407407,0.009821428571,0.008522727273,0.007523148148,0.006730769231,0.006087662338,0.005555555556\0,0.01904761905,0.02916666667,0.0140625,0.01111111111,0.009166666667,0.007792207792,0.006770833333,0.005982905983,0.005357142857,0.004848484848,0.004427083333\0,0.015625,0.02380952381,0.01157407407,0.009166666667,0.007575757576,0.006448412698,0.005608974359,0.00496031746,0.004444444444,0.004024621212,0.003676470588\0,0.01322751323,0.02008928571,0.009821428571,0.007792207792,0.006448412698,0.005494505495,0.004783163265,0.004232804233,0.003794642857,0.003437738732,0.003141534392\0,0.01145833333,0.01736111111,0.008522727273,0.006770833333,0.005608974359,0.004783163265,0.004166666667,0.003689236111,0.003308823529,0.002998737374,0.00274122807\0,0.0101010101,0.01527777778,0.007523148148,0.005982905983,0.00496031746,0.004232804233,0.003689236111,0.003267973856,0.002932098765,0.002658160553,0.002430555556\0,0.009027777778,0.01363636364,0.006730769231,0.005357142857,0.004444444444,0.003794642857,0.003308823529,0.002932098765,0.002631578947,0.002386363636,0.002182539683\0,0.008158508159,0.01231060606,0.006087662338,0.004848484848,0.004024621212,0.003437738732,0.002998737374,0.002658160553,0.002386363636,0.002164502165,0.001980027548\0,0.00744047619,0.01121794872,0.005555555556,0.004427083333,0.003676470588,0.003141534392,0.00274122807,0.002430555556,0.002182539683,0.001980027548,0.001811594203)
		}
		else if ("`kernel'"=="triangular") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.08333333333,0.1666666667,0.05,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641\0,0.1666666667,0.5,0.08333333333,0.05,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576\0,0.05,0.08333333333,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495\0,0.03333333333,0.05,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762\0,0.02380952381,0.03333333333,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667\0,0.01785714286,0.02380952381,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588\0,0.01388888889,0.01785714286,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856\0,0.01111111111,0.01388888889,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608\0,0.009090909091,0.01111111111,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947\0,0.007575757576,0.009090909091,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947,0.002380952381\0,0.00641025641,0.007575757576,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947,0.002380952381,0.002164502165)
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.01031746032,0.02222222222,0.005853174603,0.003736772487,0.002579365079,0.001881914382,0.001430976431,0.001123413623,0.0009046509047,0.0007437007437,0.0006219474969\0,0.02222222222,0.05,0.0123015873,0.007738095238,0.005291005291,0.003835978836,0.002904040404,0.002272727273,0.001825951826,0.001498501499,0.001251526252\0,0.005853174603,0.0123015873,0.003373015873,0.002175925926,0.001512746513,0.001109307359,0.0008466070966,0.0006664631665,0.0005377955378,0.0004428210678,0.0003707893414\0,0.003736772487,0.007738095238,0.002175925926,0.001414141414,0.0009884559885,0.0007277444777,0.0005570818071,0.0004395604396,0.0003553391053,0.0002930035651,0.000245621753\0,0.002579365079,0.005291005291,0.001512746513,0.0009884559885,0.0006937506938,0.000512384441,0.0003931914646,0.0003108465608,0.0002516764281,0.0002077851343,0.0001743612425\0,0.001881914382,0.003835978836,0.001109307359,0.0007277444777,0.000512384441,0.0003793825222,0.0002917139078,0.0002309951758,0.0001872718784,0.000154780147,0.0001299991432\0,0.001430976431,0.002904040404,0.0008466070966,0.0005570818071,0.0003931914646,0.0002917139078,0.0002246732026,0.0001781499637,0.0001445917726,0.0001196172249,0.0001005451663\0,0.001123413623,0.002272727273,0.0006664631665,0.0004395604396,0.0003108465608,0.0002309951758,0.0001781499637,0.0001414210909,0.0001148916061,9.512417407e-05,8.001258001e-05\0,0.0009046509047,0.001825951826,0.0005377955378,0.0003553391053,0.0002516764281,0.0001872718784,0.0001445917726,0.0001148916061,9.341535657e-05,7.739735012e-05,6.514127067e-05\0,0.0007437007437,0.001498501499,0.0004428210678,0.0002930035651,0.0002077851343,0.000154780147,0.0001196172249,9.512417407e-05,7.739735012e-05,6.416508393e-05,5.403303328e-05\0,0.0006219474969,0.001251526252,0.0003707893414,0.000245621753,0.0001743612425,0.0001299991432,0.0001005451663,8.001258001e-05,6.514127067e-05,5.403303328e-05,4.552211074e-05)
		}
		else if ("`kernel'"=="epanechnikov") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.1,0.1875,0.0625,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429\0,0.1875,0.5,0.1,0.0625,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049\0,0.0625,0.1,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692\0,0.04285714286,0.0625,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571\0,0.03125,0.04285714286,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941\0,0.02380952381,0.03125,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333\0,0.01875,0.02380952381,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848\0,0.01515151515,0.01875,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667\0,0.0125,0.01515151515,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496\0,0.01048951049,0.0125,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496,0.003409090909\0,0.008928571429,0.01048951049,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496,0.003409090909,0.003105590062)
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.01428571429,0.028515625,0.008515625,0.005627705628,0.003984375,0.002963702964,0.002287946429,0.001818181818,0.001478794643,0.001225832991,0.001032366071\0,0.028515625,0.05892857143,0.01666666667,0.01088169643,0.007643398268,0.005654761905,0.004348776224,0.00344629329,0.002797202797,0.002315067745,0.001947317388\0,0.008515625,0.01666666667,0.005140692641,0.003426339286,0.002440268065,0.001822916667,0.001411713287,0.001124526515,0.0009162895928,0.0007606325966,0.0006413091552\0,0.005627705628,0.01088169643,0.003426339286,0.002297702298,0.001643813776,0.001232101232,0.0009566326531,0.0007635501753,0.000623139881,0.0005179340783,0.0004371279762\0,0.003984375,0.007643398268,0.002440268065,0.001643813776,0.00118006993,0.0008868781888,0.0006900452489,0.0005516943994,0.0004508513932,0.0003751456876,0.000316903077\0,0.002963702964,0.005654761905,0.001822916667,0.001232101232,0.0008868781888,0.0006679594915,0.0005206118906,0.0004168174447,0.0003410218254,0.0002840296958,0.0002401244589\0,0.002287946429,0.004348776224,0.001411713287,0.0009566326531,0.0006900452489,0.0005206118906,0.0004063467492,0.0003257181187,0.0002667514374,0.0002223557692,0.0001881158642\0,0.001818181818,0.00344629329,0.001124526515,0.0007635501753,0.0005516943994,0.0004168174447,0.0003257181187,0.0002613485586,0.0002142160239,0.0001786923984,0.0001512691854\0,0.001478794643,0.002797202797,0.0009162895928,0.000623139881,0.0004508513932,0.0003410218254,0.0002667514374,0.0002142160239,0.0001757110167,0.0001466644151,0.0001242236025\0,0.001225832991,0.002315067745,0.0007606325966,0.0005179340783,0.0003751456876,0.0002840296958,0.0002223557692,0.0001786923984,0.0001466644151,0.0001224862094,0.0001037942608\0,0.001032366071,0.001947317388,0.0006413091552,0.0004371279762,0.000316903077,0.0002401244589,0.0001881158642,0.0001512691854,0.0001242236025,0.0001037942608,8.799171843e-05)
		}
	}
	Psi=(0,-1,0,0,0,0,0,0,0,0,0,0\-1,0,0,0,0,0,0,0,0,0,0,0\0,0,1,0,0,0,0,0,0,0,0,0\0,0,0,1,0,0,0,0,0,0,0,0\0,0,0,0,-1,0,0,0,0,0,0,0\0,0,0,0,0,1,0,0,0,0,0,0\0,0,0,0,0,0,-1,0,0,0,0,0\0,0,0,0,0,0,0,1,0,0,0,0\0,0,0,0,0,0,0,0,-1,0,0,0\0,0,0,0,0,0,0,0,0,1,0,0\0,0,0,0,0,0,0,0,0,0,-1,0\0,0,0,0,0,0,0,0,0,0,0,1)

	****************************************************************************
	** Select preliminary bandwidths.
	****************************************************************************
	mu = mean(X); sd = (variance(X))^(1/2)

	fhatb = sd^(2*`p'+5) * normalden(-mu/sd) / (rddensity_h(-mu/sd,`p'+2) * normalden(-mu/sd))^2
	C_b = (25884.444444494150957,3430865.4551236177795,845007948.04262602329,330631733667.03808594,187774809656037.3125,145729502641999264,146013502974449876992)
	b = ((2*`p'+1)/4 * fhatb * C_b[`p']/`N')^(1/(2*`p'+5))
	
	fhatc = sd^(2*`p'+1) * normalden(-mu/sd) / (rddensity_h(-mu/sd,`p') * normalden(-mu/sd))^2
	C_c = (4.8000000000000246914,548.57142857155463389,100800.00000020420703,29558225.458100609481,12896196859.612621307,7890871468221.609375,6467911284037581)
	c = (1/(2*`p') * fhatc * C_c[`p']/`N')^(1/(2*`p'+1))
	
	// b is for higher-order derivative estimation
    // c is for density estimation
	
	if (`regularize') {
		
		// bandwidth should not exceed the range of data
		b = min( (b, max(abs(XUnique))) )
		c = min( (c, max(abs(XUnique))) )

		// nlocalmin check
		
		if (`nlocalmin' > 0) {
			b = max((b, sort(abs(X[selectindex(X :< 0)]), 1)[min((20+`p'+2+1, `Nl'))], (X[selectindex(X :>= 0)])[min((20+`p'+2+1, `Nr'))]))
			c = max((c, sort(abs(X[selectindex(X :< 0)]), 1)[min((20+`p'+  1, `Nl'))], (X[selectindex(X :>= 0)])[min((20+`p'  +1, `Nr'))]))
		}

		// nuniquemin check
		if (`nuniquemin' > 0) {
			b = max((b, sort(abs(XUnique[selectindex(XUnique :< 0)]), 1)[min((20+`p'+2+1, NlUnique))], (XUnique[selectindex(XUnique :>= 0)])[min((20+`p'+2+1, NrUnique))]))
			c = max((c, sort(abs(XUnique[selectindex(XUnique :< 0)]), 1)[min((20+`p'  +1, NlUnique))], (XUnique[selectindex(XUnique :>= 0)])[min((20+`p'  +1, NrUnique))]))
		}
	}
	
	st_numscalar("BW_b", b)
	st_numscalar("BW_c", c)

	****************************************************************************
	** Estimate main bandwidths.
	****************************************************************************
	Xb = select(X, -b:<=X :& X:<=b)
	Nlb = sum(-b:<=X :& X:<0)
	Nrb = rows(Xb) - Nlb
	
	Xc = select(X, -c:<=X :& X:<=c)
	Nlc = sum(-c:<=X :& X:<0)
	Nrc = rows(Xc) - Nlc
	
	Ytemp = (0..(`N'-1))' :/ (`N'-1)
	if (`masspoints') {
		Ytemp = rddensity_rep(Ytemp[indexUnique], freqUnique)
	}
	Yb = select(Ytemp, -b:<=X :& X:<=b)
	Yc = select(Ytemp, -c:<=X :& X:<=c)

	h = J(4,3,0)

	fV_b = rddensity_fv(Yb, Xb, `Nl', `Nr', Nlb, Nrb, b, b, `p'+2 , `p'+1, "`kernel'", "`fitselect'", "`vce'", `masspoints')
	fV_c = rddensity_fv(Yc, Xc, `Nl', `Nr', Nlc, Nrc, c, c, `p'   , 1    , "`kernel'", "`fitselect'", "`vce'", `masspoints')
	

	h[.,2] = `N'*c*fV_c[.,2]

	if ("`fitselect'"=="unrestricted") {
		h[1,3] = fV_b[1,3] * Bsq_p[`p']^(1/2) * (-1)^`p' * factorial(`p'+1)
		h[2,3] = fV_b[2,3] * Bsq_p[`p']^(1/2) * factorial(`p'+1)
	}
	else if ("`fitselect'"=="restricted") {
		Psi = Psi[1..`p'+2,1..`p'+2];
		Gplus = Gplus[1..`p'+2,1..`p'+2]; Gminus = Psi*Gplus*Psi;
		vplus = Splus[1..`p'+2,`p'+3]; vminus = Psi*vplus;
		Splus = Splus[1..`p'+2,1..`p'+2]; Sminus = Psi*Splus*Psi;
		S = invsym(fV_c[2,1] * Splus + fV_c[1,1] * Sminus);
		B = fV_b[1,3] * S[1..2,] * (fV_c[1,1] * (-1)^(`p'+1) * vminus + fV_c[2,1] * vplus) 
		h[1,3] = B[1,1]
		h[2,3] = B[2,1]
	}

 	h[3,3] = h[2,3] - h[1,3]; h[4,3] = h[2,3] + h[1,3]; h[.,3] = h[.,3]:^2;
	h[.,1] = ((1/(2*`p')) * (h[.,2]:/h[.,3]) * (1/`N')):^(1/(2*`p'+1));
	
	if (`regularize') {
		
		for (i=1; i<=4; i++) {
			if (h[i, 2] < 0) { 
				h[i, 1] = 0
				h[i, 2] = . 
			}
			if (h[i, 1] == .) { 
				h[i, 1] = 0 
			}
		}
	
		// bandwidth should not exceed the range of data
		h[1,1] = min((h[1,1], abs(XUnique[1])))
		h[2,1] = min((h[2,1], XUnique[NUnique]))
		h[3,1] = min((h[3,1], max((abs(XUnique[1]), XUnique[NUnique]))))
		h[4,1] = min((h[4,1], max((abs(XUnique[1]), XUnique[NUnique]))))

		// nlocalmin check
		if (`nlocalmin' > 0) {
			hlMin = sort(abs(X[selectindex(X :< 0)]), 1)[min((`Nl', `nlocalmin'))]
			hrMin = (X[selectindex(X :>= 0)])[min((`Nr', `nlocalmin'))]
			h[1,1] = max((h[1,1], hlMin))
			h[2,1] = max((h[2,1], hrMin))
			h[3,1] = max((h[3,1], hlMin, hrMin))
			h[4,1] = max((h[4,1], hlMin, hrMin))
		}

		// nuniquemin check
		if (`nuniquemin' > 0) {
			hlMin = sort(abs(XUnique[selectindex(XUnique :< 0)]),1)[min((NlUnique, `nuniquemin'))]
			hrMin = (XUnique[selectindex(XUnique :>= 0)])[min((NrUnique, `nuniquemin'))]
			h[1,1] = max((h[1,1], hlMin))
			h[2,1] = max((h[2,1], hrMin))
			h[3,1] = max((h[3,1], hlMin, hrMin))
			h[4,1] = max((h[4,1], hlMin, hrMin))
		}
	}
	
	st_matrix("h", h);

	*display("Estimation completed.");
	}
	*** END MATA ESTIMATION **************************************************** 
	****************************************************************************

	****************************************************************************
	*** BEGIN OUTPUT TABLE *****************************************************
	if (masspoints_flag == 1) {
		disp ""
		disp "Point estimates and standard errors have been adjusted for repeated observations."
		disp "(Use option {it:nomasspoints} to suppress this adjustment.)"
	}
	
	disp ""
	disp "Bandwidth selection for manipulation testing." 

	disp ""
	disp in smcl in gr "Cutoff " in ye "c = " %10.3f `c'  _col(22) " {c |} " _col(23) in gr "Left of " in ye "c" _col(36) in gr "Right of " in y "c"  _col(58) in gr "Number of obs = "  in ye %12.0f `N'
	disp in smcl in gr "{hline 22}{c +}{hline 22}"                                                                                                    _col(58) in gr "Model         = "  in ye "{ralign 12:`fitselect'}"
	disp in smcl in gr "{ralign 21:Number of obs}"        _col(22) " {c |} " _col(23) as result %9.0f `Nl'       _col(37) as result %9.0f  `Nr'       _col(58) in gr "Kernel        = "  in ye "{ralign 12:`kernel'}"
	disp in smcl in gr "{ralign 21:Min Running var.}"     _col(22) " {c |} " _col(23) as result %9.3f `xl_min'   _col(37) as result %9.3f  `xr_min'   _col(58) in gr "VCE method    = "  in ye "{ralign 12:`vce'}"
	disp in smcl in gr "{ralign 21:Max Running var.}"     _col(22) " {c |} " _col(23) as result %9.3f `xl_max'   _col(37) as result %9.3f  `xr_max'
	disp in smcl in gr "{ralign 21:Order loc. poly. (p)}" _col(22) " {c |} " _col(23) as result %9.0f `p'        _col(37) as result %9.0f  `p'

	disp ""
	disp "Running variable: `x'."
	disp in smcl in gr "{hline 22}{c TT}{hline 34}"
	disp in smcl in gr "{ralign 21:Target}"              _col(22) " {c |} " _col(23) in gr "Bandwidth"          _col(37) " Variance"                 _col(49) "   Bias^2" 
	disp in smcl in gr "{hline 22}{c +}{hline 34}"
	disp in smcl in gr "{ralign 21:left density}"        _col(22) " {c |} " _col(23) as result %9.3f h[1,1]     _col(37) as result %9.3f h[1,2]      _col(49) as result %9.3f h[1,3]
	disp in smcl in gr "{ralign 21:right density}"       _col(22) " {c |} " _col(23) as result %9.3f h[2,1]     _col(37) as result %9.3f h[2,2]      _col(49) as result %9.3f h[2,3]
	disp in smcl in gr "{ralign 21:difference densities}" _col(22) " {c |} " _col(23) as result %9.3f h[3,1]     _col(37) as result %9.3f h[3,2]      _col(49) as result %9.3f h[3,3]
	disp in smcl in gr "{ralign 21:sum densities}"      _col(22) " {c |} " _col(23) as result %9.3f h[4,1]     _col(37) as result %9.3f h[4,2]      _col(49) as result %9.3f h[4,3]
	disp in smcl in gr "{hline 22}{c BT}{hline 34}"
	disp ""
	*** END OUTPUT TABLE ******************************************************* 
	****************************************************************************

	restore

	ereturn clear
	ereturn scalar c = `c'
	ereturn scalar p = `p'
	ereturn scalar N_l = `Nl'
	ereturn scalar N_r = `Nr'
	mat rown h = f_left f_right f_diff f_sum
	mat coln h = bandwidth var bias2
	ereturn matrix h = h
	ereturn scalar BW_b = BW_b
	ereturn scalar BW_c = BW_c
	
	ereturn local runningvar "`x'"
	ereturn local kernel = "`kernel'"
	ereturn local fitmethod = "`fitselect'"
	ereturn local vce = "`vce'"

	mata: mata clear
	
end
	
