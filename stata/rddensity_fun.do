********************************************************************************
* RDDENSITY STATA PACKAGE -- rddensity -- Mata functions
* Authors: Matias D. Cattaneo, Michael Jansson, Xinwei Ma
********************************************************************************
*!version 2.4 2023-01-21

** NOTE: DATA IS ASSUMED TO BE IN ASCENDING ORDER

** do rddensity_fun.ado 

********************************************************************************
* Extracting unique elements
********************************************************************************

capture mata mata drop rddensity_unique()

mata
real matrix rddensity_unique(real colvector x){
// Note: x should be in ascending order

n = length(x)

// x has one or no element
if (n == 0) {
	out = J(0, 4, .)
	return(out)
}
if (n == 1) {
	out = (x, 1, 1, 1)
	return(out)
}

// else
numIndexTemp = selectindex(x[2..n] :!= x[1..(n-1)])
if (length(numIndexTemp) == 0) {
	numIndexLast  = J(1, 1, 1)
	numIndexFirst = J(1, 1, 1)
} else {
	numIndexLast  = (numIndexTemp \ n)
	numIndexFirst = (1 \ numIndexTemp:+1)
}

freq = numIndexLast - numIndexFirst :+ 1
unique = x[numIndexLast]

out = (unique, freq, numIndexFirst, numIndexLast)
return(out)
}

mata mosave rddensity_unique(), replace
end

********************************************************************************
* Replicating each element in x y times (elementwise)
********************************************************************************
capture mata mata drop rddensity_rep()

mata
real matrix rddensity_rep(real colvector x, real colvector y){
// Warning: x and y should have the same length
//          y should contain strictly positive integers
//          y cannot be empty

nout = sum(y)
if (length(y) == 1) {
	out = J(nout, 1, x)
}
else {
	if (all(y :== 1)) {
		out = x
	}
	else {
		out = J(nout, 1, .)
	    out[1..y[1], 1] = J(y[1], 1, x[1])
	    indextemp = y[1]
	    for (i=2; i<=length(y); i++) {
			out[(indextemp+1)..(indextemp+y[i])] = J(y[i], 1, x[i])
			indextemp = indextemp + y[i]

		}
	}
}

return(out)
}
/* // testing
mata: lpdensity_rep(5, 1)
mata: lpdensity_rep(5, 3)
mata: lpdensity_rep((5\ 6), (1\ 1))
mata: lpdensity_rep((5\ 6), (1\ 2))
mata: lpdensity_rep((5\ 6), (2\ 1))
mata: lpdensity_rep((5\ 6), (2\ 3))
mata: lpdensity_rep((5\ 6\ 7\ 8), (1\ 1\ 1\ 1))
mata: lpdensity_rep((5\ 6\ 7\ 8), (1\ 2\ 2\ 3))
mata: lpdensity_rep((5\ 6\ 7\ 8), (3\ 1\ 2\ 1))
mata: lpdensity_rep((5\ 6\ 7\ 8), (3\ 2\ 2\ 3))
*/ 
mata mosave rddensity_rep(), replace
end

********************************************************************************
* Main estimation function
********************************************************************************

capture mata mata drop rddensity_fv()

mata
real matrix rddensity_fv(real colvector Y, real colvector X,
                         real scalar Nl, real scalar Nr, real scalar Nlh, real scalar Nrh,
                         real scalar hl, real scalar hr,
                         real scalar p, real scalar s,
						 string scalar kernel, string fitselect, string scalar vce, 
						 real scalar masspoints){
	N = Nl + Nr; Nh = Nlh + Nrh;
	W = J(Nh, 1, 1);
	if (kernel=="triangular") {W = W:-(X[1..Nlh]/(-hl)\X[(Nlh+1)..Nh]/hr);}
	else if (kernel=="epanechnikov") {W = W:-(X[1..Nlh]/hl\X[(Nlh+1)..Nh]/hr):^2;}

	if (fitselect=="restricted"){
		Xp = J(Nh,1+p+1,0); Xp[1..Nlh,1] = X[1..Nlh]; Xp[(Nlh+1)..Nh,2] = X[(Nlh+1)..Nh];
		Xp[.,3] = J(Nh, 1, 1);
		for (j=2; j<=p; j++){Xp[.,j+2] = X:^j;}
	}
	else if (fitselect=="unrestricted"){
		Xp = J(Nh,2*(p+1),0); Xp[1..Nlh,1] = X[1..Nlh]; Xp[(Nlh+1)..Nh,2] = X[(Nlh+1)..Nh];
		Xp[1..Nlh,3] = J(Nlh, 1, 1); Xp[(Nlh+1)..Nh,4] = J(Nrh, 1, 1);
		for (j=2; j<=p; j++){
			Xp[1..Nlh     ,2*j+1] = X[1..Nlh]:^j;
			Xp[(Nlh+1)..Nh,2*j+2] = X[(Nlh+1)..Nh]:^j;
		}
	}

	S = cross(Xp,W,Xp); Sinv = invsym(S);
	
	XpWY = cross(Xp,W,Y); b = Sinv * XpWY;
	
	out = J(4,3,0); out[.,1] = (b[1,1]\b[2,1]\b[2,1]-b[1,1]\b[2,1]+b[1,1]);
	
	if (fitselect=="restricted"){out[.,3] = (b[s+2,1]\b[s+2,1]\0\2*b[s+2,1]);}
	else if (fitselect=="unrestricted"){out[.,3] = (b[2*s+1,1]\b[2*s+2,1]\b[2*s+2,1]-b[2*s+1,1]\b[2*s+2,1]+b[2*s+1,1]);}

	if (vce=="jackknife"){
		XpW = Xp:*W; L = J(Nh,cols(Xp),0)
		
		if (masspoints) {
			XUnique     = rddensity_unique(X)
			freqUnique  = XUnique[., 2]
			indexUnique = XUnique[., 3]
			
			for (jj=1; jj<=cols(L); jj++) {
				L[., jj] = rddensity_rep(((runningsum((0\ XpW[Nh..1, jj])) :/ (N - 1))[Nh..1])[indexUnique], freqUnique)
			}
		}
		else {
			for (i=1; i<=Nh; i++) { 
				L[1,.] = L[1,.]   + XpW[i,.]   / (N-1) 
			}
			for (i=2; i<=Nh; i++) { 
				L[i,.] = L[i-1,.] - XpW[i-1,.] / (N-1) 
			}
		}

		V = Sinv[1..2,.]*cross(L,L)*Sinv[.,1..2];
	}
	else if (vce=="plugin" & fitselect=="unrestricted"){
		if (kernel=="uniform") {Varv0=(1.2000000000000046185,5.4857142857155736237,14.28571428575378377,29.090909094894868758,51.398601335826242575,82.707663360750302672,124.51643660507397726,178.21357816224917769);}
		else if (kernel=="triangular") {Varv0=(1.3714285714285674445,5.7142857142855163488,14.025974025962113956,27.41258741335866489,46.993007009769826254,73.889594174983358243,109.22603222282486968,154.12574072671122849);}
		else if (kernel=="epanechnikov") {Varv0=(1.345468935496637819,5.7720057720059116946,14.467758368262224167,28.72478523773543202,49.848500126756334794,79.148187721824797336,117.93483411672059447,167.52591790934093297);}
		V = diag((out[1,1] * Varv0[p] / (N*hl), out[2,1] * Varv0[p] / (N*hr)))
	}
	else if (vce=="plugin" & fitselect=="restricted"){
		if (kernel=="uniform") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.1666666667,0.25,0.125,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667\0,0.25,0.5,0.1666666667,0.125,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545\0,0.125,0.1666666667,0.1,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846\0,0.1,0.125,0.08333333333,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571\0,0.08333333333,0.1,0.07142857143,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333\0,0.07142857143,0.08333333333,0.0625,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125\0,0.0625,0.07142857143,0.05555555556,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471\0,0.05555555556,0.0625,0.05,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778\0,0.05,0.05555555556,0.04545454545,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947\0,0.04545454545,0.05,0.04166666667,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947,0.025\0,0.04166666667,0.04545454545,0.03846153846,0.03571428571,0.03333333333,0.03125,0.02941176471,0.02777777778,0.02631578947,0.025,0.02380952381)[1..p+2,1..p+2]
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.03333333333,0.05208333333,0.02430555556,0.01904761905,0.015625,0.01322751323,0.01145833333,0.0101010101,0.009027777778,0.008158508159,0.00744047619\0,0.05208333333,0.08333333333,0.0375,0.02916666667,0.02380952381,0.02008928571,0.01736111111,0.01527777778,0.01363636364,0.01231060606,0.01121794872\0,0.02430555556,0.0375,0.01785714286,0.0140625,0.01157407407,0.009821428571,0.008522727273,0.007523148148,0.006730769231,0.006087662338,0.005555555556\0,0.01904761905,0.02916666667,0.0140625,0.01111111111,0.009166666667,0.007792207792,0.006770833333,0.005982905983,0.005357142857,0.004848484848,0.004427083333\0,0.015625,0.02380952381,0.01157407407,0.009166666667,0.007575757576,0.006448412698,0.005608974359,0.00496031746,0.004444444444,0.004024621212,0.003676470588\0,0.01322751323,0.02008928571,0.009821428571,0.007792207792,0.006448412698,0.005494505495,0.004783163265,0.004232804233,0.003794642857,0.003437738732,0.003141534392\0,0.01145833333,0.01736111111,0.008522727273,0.006770833333,0.005608974359,0.004783163265,0.004166666667,0.003689236111,0.003308823529,0.002998737374,0.00274122807\0,0.0101010101,0.01527777778,0.007523148148,0.005982905983,0.00496031746,0.004232804233,0.003689236111,0.003267973856,0.002932098765,0.002658160553,0.002430555556\0,0.009027777778,0.01363636364,0.006730769231,0.005357142857,0.004444444444,0.003794642857,0.003308823529,0.002932098765,0.002631578947,0.002386363636,0.002182539683\0,0.008158508159,0.01231060606,0.006087662338,0.004848484848,0.004024621212,0.003437738732,0.002998737374,0.002658160553,0.002386363636,0.002164502165,0.001980027548\0,0.00744047619,0.01121794872,0.005555555556,0.004427083333,0.003676470588,0.003141534392,0.00274122807,0.002430555556,0.002182539683,0.001980027548,0.001811594203)[1..p+2,1..p+2]
		}
		else if (kernel=="triangular") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.08333333333,0.1666666667,0.05,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641\0,0.1666666667,0.5,0.08333333333,0.05,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576\0,0.05,0.08333333333,0.03333333333,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495\0,0.03333333333,0.05,0.02380952381,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762\0,0.02380952381,0.03333333333,0.01785714286,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667\0,0.01785714286,0.02380952381,0.01388888889,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588\0,0.01388888889,0.01785714286,0.01111111111,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856\0,0.01111111111,0.01388888889,0.009090909091,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608\0,0.009090909091,0.01111111111,0.007575757576,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947\0,0.007575757576,0.009090909091,0.00641025641,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947,0.002380952381\0,0.00641025641,0.007575757576,0.005494505495,0.004761904762,0.004166666667,0.003676470588,0.003267973856,0.002923976608,0.002631578947,0.002380952381,0.002164502165)[1..p+2,1..p+2]
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.01031746032,0.02222222222,0.005853174603,0.003736772487,0.002579365079,0.001881914382,0.001430976431,0.001123413623,0.0009046509047,0.0007437007437,0.0006219474969\0,0.02222222222,0.05,0.0123015873,0.007738095238,0.005291005291,0.003835978836,0.002904040404,0.002272727273,0.001825951826,0.001498501499,0.001251526252\0,0.005853174603,0.0123015873,0.003373015873,0.002175925926,0.001512746513,0.001109307359,0.0008466070966,0.0006664631665,0.0005377955378,0.0004428210678,0.0003707893414\0,0.003736772487,0.007738095238,0.002175925926,0.001414141414,0.0009884559885,0.0007277444777,0.0005570818071,0.0004395604396,0.0003553391053,0.0002930035651,0.000245621753\0,0.002579365079,0.005291005291,0.001512746513,0.0009884559885,0.0006937506938,0.000512384441,0.0003931914646,0.0003108465608,0.0002516764281,0.0002077851343,0.0001743612425\0,0.001881914382,0.003835978836,0.001109307359,0.0007277444777,0.000512384441,0.0003793825222,0.0002917139078,0.0002309951758,0.0001872718784,0.000154780147,0.0001299991432\0,0.001430976431,0.002904040404,0.0008466070966,0.0005570818071,0.0003931914646,0.0002917139078,0.0002246732026,0.0001781499637,0.0001445917726,0.0001196172249,0.0001005451663\0,0.001123413623,0.002272727273,0.0006664631665,0.0004395604396,0.0003108465608,0.0002309951758,0.0001781499637,0.0001414210909,0.0001148916061,9.512417407e-05,8.001258001e-05\0,0.0009046509047,0.001825951826,0.0005377955378,0.0003553391053,0.0002516764281,0.0001872718784,0.0001445917726,0.0001148916061,9.341535657e-05,7.739735012e-05,6.514127067e-05\0,0.0007437007437,0.001498501499,0.0004428210678,0.0002930035651,0.0002077851343,0.000154780147,0.0001196172249,9.512417407e-05,7.739735012e-05,6.416508393e-05,5.403303328e-05\0,0.0006219474969,0.001251526252,0.0003707893414,0.000245621753,0.0001743612425,0.0001299991432,0.0001005451663,8.001258001e-05,6.514127067e-05,5.403303328e-05,4.552211074e-05)[1..p+2,1..p+2]
		}
		else if (kernel=="epanechnikov") {
		Splus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.1,0.1875,0.0625,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429\0,0.1875,0.5,0.1,0.0625,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049\0,0.0625,0.1,0.04285714286,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692\0,0.04285714286,0.0625,0.03125,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571\0,0.03125,0.04285714286,0.02380952381,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941\0,0.02380952381,0.03125,0.01875,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333\0,0.01875,0.02380952381,0.01515151515,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848\0,0.01515151515,0.01875,0.0125,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667\0,0.0125,0.01515151515,0.01048951049,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496\0,0.01048951049,0.0125,0.008928571429,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496,0.003409090909\0,0.008928571429,0.01048951049,0.007692307692,0.006696428571,0.005882352941,0.005208333333,0.004643962848,0.004166666667,0.003759398496,0.003409090909,0.003105590062)[1..p+2,1..p+2]
		Gplus=(0,0,0,0,0,0,0,0,0,0,0,0\0,0.01428571429,0.028515625,0.008515625,0.005627705628,0.003984375,0.002963702964,0.002287946429,0.001818181818,0.001478794643,0.001225832991,0.001032366071\0,0.028515625,0.05892857143,0.01666666667,0.01088169643,0.007643398268,0.005654761905,0.004348776224,0.00344629329,0.002797202797,0.002315067745,0.001947317388\0,0.008515625,0.01666666667,0.005140692641,0.003426339286,0.002440268065,0.001822916667,0.001411713287,0.001124526515,0.0009162895928,0.0007606325966,0.0006413091552\0,0.005627705628,0.01088169643,0.003426339286,0.002297702298,0.001643813776,0.001232101232,0.0009566326531,0.0007635501753,0.000623139881,0.0005179340783,0.0004371279762\0,0.003984375,0.007643398268,0.002440268065,0.001643813776,0.00118006993,0.0008868781888,0.0006900452489,0.0005516943994,0.0004508513932,0.0003751456876,0.000316903077\0,0.002963702964,0.005654761905,0.001822916667,0.001232101232,0.0008868781888,0.0006679594915,0.0005206118906,0.0004168174447,0.0003410218254,0.0002840296958,0.0002401244589\0,0.002287946429,0.004348776224,0.001411713287,0.0009566326531,0.0006900452489,0.0005206118906,0.0004063467492,0.0003257181187,0.0002667514374,0.0002223557692,0.0001881158642\0,0.001818181818,0.00344629329,0.001124526515,0.0007635501753,0.0005516943994,0.0004168174447,0.0003257181187,0.0002613485586,0.0002142160239,0.0001786923984,0.0001512691854\0,0.001478794643,0.002797202797,0.0009162895928,0.000623139881,0.0004508513932,0.0003410218254,0.0002667514374,0.0002142160239,0.0001757110167,0.0001466644151,0.0001242236025\0,0.001225832991,0.002315067745,0.0007606325966,0.0005179340783,0.0003751456876,0.0002840296958,0.0002223557692,0.0001786923984,0.0001466644151,0.0001224862094,0.0001037942608\0,0.001032366071,0.001947317388,0.0006413091552,0.0004371279762,0.000316903077,0.0002401244589,0.0001881158642,0.0001512691854,0.0001242236025,0.0001037942608,8.799171843e-05)[1..p+2,1..p+2]
		}
		Psi=(0,-1,0,0,0,0,0,0,0,0,0,0\-1,0,0,0,0,0,0,0,0,0,0,0\0,0,1,0,0,0,0,0,0,0,0,0\0,0,0,1,0,0,0,0,0,0,0,0\0,0,0,0,-1,0,0,0,0,0,0,0\0,0,0,0,0,1,0,0,0,0,0,0\0,0,0,0,0,0,-1,0,0,0,0,0\0,0,0,0,0,0,0,1,0,0,0,0\0,0,0,0,0,0,0,0,-1,0,0,0\0,0,0,0,0,0,0,0,0,1,0,0\0,0,0,0,0,0,0,0,0,0,-1,0\0,0,0,0,0,0,0,0,0,0,0,1)[1..p+2,1..p+2]
		Sminus = Psi*Splus*Psi; Gminus = Psi*Gplus*Psi
		S = invsym(out[2,1] * Splus + out[1,1] * Sminus)
		V = S[1..2,]*(out[2,1]^3 * Gplus + out[1,1]^3 * Gminus)*S[,1..2] / (N*hl)
	}
	out[.,2] = (V[1,1]\V[2,2]\(-1,1)*V*(-1\1)\(1,1)*V*(1\1))
	
	return(out)
}
mata mosave rddensity_fv(), replace
end

********************************************************************************
* Preliminary bandwidth selection
********************************************************************************

capture mata mata drop rddensity_h()

mata
real scalar rddensity_h(real scalar x, real scalar p){
	if (p==0)  out = 1
	if (p==1)  out = x
	if (p==2)  out = x^2 - 1
	if (p==3)  out = x^3 - 3*x
	if (p==4)  out = x^4 - 6*x^2 + 3
	if (p==5)  out = x^5 - 10*x^3 + 15*x
	if (p==6)  out = x^6 - 15*x^4 + 45*x^2 - 15
	if (p==7)  out = x^7 - 21*x^5 + 105*x^3 - 105*x
	if (p==8)  out = x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
	if (p==9)  out = x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
	if (p==10) out = x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
	return(out)
}
mata mosave rddensity_h(), replace
end

********************************************************************************
* Empirical quantile
********************************************************************************
capture mata mata drop rddensity_quantile()

mata
real scalar rddensity_quantile(real colvector x, real scalar p){

x = sort(x, 1)
n = length(x)
q = ceil(p * n)
if (q < 1) q = 1
if (q > n) q = n
out = x[q]

return(out)

}
mata mosave rddensity_quantile(), replace
end
