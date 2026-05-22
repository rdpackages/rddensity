version 16
clear all
set more off
set graphics off

args repo_root
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
adopath ++ "`repo_root'/stata"

capture mkdir "`repo_root'/stata/tests"
capture mkdir "`repo_root'/stata/tests/fixtures"

tempfile results
tempname handle
postfile `handle' str16 language str16 precision str64 case str32 metric double value using "`results'", replace

capture program drop _rdd_post_density
program define _rdd_post_density
    syntax, CASE(string) HANDLE(name) PRECISION(string) [ALL]

    local vce_suffix "jk"
    if "`e(vce)'" == "plugin" local vce_suffix "asy"

    post `handle' ("stata") ("`precision'") ("`case'") ("hat_left") (e(f_ql))
    post `handle' ("stata") ("`precision'") ("`case'") ("hat_right") (e(f_qr))
    post `handle' ("stata") ("`precision'") ("`case'") ("hat_diff") (e(f_qr) - e(f_ql))
    post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_left") (e(se_ql))
    post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_right") (e(se_qr))
    post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_diff") (e(se_q))
    post `handle' ("stata") ("`precision'") ("`case'") ("test_t_`vce_suffix'") (e(T_q))
    post `handle' ("stata") ("`precision'") ("`case'") ("test_p_`vce_suffix'") (e(pv_q))
    post `handle' ("stata") ("`precision'") ("`case'") ("h_left") (e(h_l))
    post `handle' ("stata") ("`precision'") ("`case'") ("h_right") (e(h_r))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_full") (e(N_l) + e(N_r))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_left") (e(N_l))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_right") (e(N_r))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_eff_left") (e(N_h_l))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_eff_right") (e(N_h_r))

    if "`all'" != "" {
        post `handle' ("stata") ("`precision'") ("`case'") ("hat_p_left") (e(f_pl))
        post `handle' ("stata") ("`precision'") ("`case'") ("hat_p_right") (e(f_pr))
        post `handle' ("stata") ("`precision'") ("`case'") ("hat_p_diff") (e(f_pr) - e(f_pl))
        post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_p_left") (e(se_pl))
        post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_p_right") (e(se_pr))
        post `handle' ("stata") ("`precision'") ("`case'") ("sd_`vce_suffix'_p_diff") (e(se_p))
        post `handle' ("stata") ("`precision'") ("`case'") ("test_p_t_`vce_suffix'") (e(T_p))
        post `handle' ("stata") ("`precision'") ("`case'") ("test_p_p_`vce_suffix'") (e(pv_p))
    }
end

capture program drop _rdd_post_bandwidth
program define _rdd_post_bandwidth
    syntax, CASE(string) HANDLE(name) PRECISION(string)

    matrix H = e(h)
    local rownames : rownames H
    local colnames : colnames H
    forvalues i = 1/`=rowsof(H)' {
        local row : word `i' of `rownames'
        if "`row'" == "f_left" local rowkey "l"
        if "`row'" == "f_right" local rowkey "r"
        if "`row'" == "f_diff" local rowkey "diff"
        if "`row'" == "f_sum" local rowkey "sum"
        forvalues j = 1/`=colsof(H)' {
            local col : word `j' of `colnames'
            if "`col'" == "bandwidth" local colkey "bw"
            if "`col'" == "var" local colkey "variance"
            if "`col'" == "bias2" local colkey "biassq"
            post `handle' ("stata") ("`precision'") ("`case'") ("h_`rowkey'_`colkey'") (H[`i', `j'])
        }
    }
    post `handle' ("stata") ("`precision'") ("`case'") ("n_full") (e(N_l) + e(N_r))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_left") (e(N_l))
    post `handle' ("stata") ("`precision'") ("`case'") ("n_right") (e(N_r))
end

capture program drop _rdd_make_data
program define _rdd_make_data
    syntax, STORAGE(string)
    clear
    set obs 200
    gen `storage' x = -1.5 + 3 * (_n - 1) / 199 + 0.05 * sin(_n - 1)
end

capture program drop _rdd_make_masspoint_data
program define _rdd_make_masspoint_data
    syntax, STORAGE(string)
    clear
    set obs 240
    gen `storage' x = round(-1.4 + 3 * (_n - 1) / 239 + 0.08 * sin((_n - 1) / 2), .1)
end

foreach precision in single double {
    local storage_type "double"
    if "`precision'" == "single" local storage_type "float"

    _rdd_make_data, storage(`storage_type')
    rddensity x, h(0.6 0.6) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_fixed") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rddensity x, nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_estimated") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rddensity x, h(0.7 0.7) kernel(uniform) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_uniform") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rddensity x, h(0.8 0.8) fitselect(restricted) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_restricted") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rddensity x, h(0.6 0.6) all nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_all") handle(`handle') precision(`precision') all

    _rdd_make_data, storage(`storage_type')
    rddensity x, h(0.75 0.65) kernel(epanechnikov) vce(plugin) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_epanechnikov_plugin") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rddensity x, c(0.15) h(0.5 0.7) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_nonzero_cutoff") handle(`handle') precision(`precision')

    _rdd_make_masspoint_data, storage(`storage_type')
    rddensity x, h(0.8 0.8) nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_masspoints_adjusted") handle(`handle') precision(`precision')

    _rdd_make_masspoint_data, storage(`storage_type')
    rddensity x, h(0.8 0.8) nomasspoints nobinomial precision(`precision')
    _rdd_post_density, case("rddensity_masspoints_unadjusted") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rdbwdensity x, precision(`precision')
    _rdd_post_bandwidth, case("rdbwdensity_default") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rdbwdensity x, kernel(uniform) vce(plugin) precision(`precision')
    _rdd_post_bandwidth, case("rdbwdensity_uniform_plugin") handle(`handle') precision(`precision')

    _rdd_make_data, storage(`storage_type')
    rdbwdensity x, c(0.15) p(3) kernel(epanechnikov) precision(`precision')
    _rdd_post_bandwidth, case("rdbwdensity_nonzero_cutoff") handle(`handle') precision(`precision')

    _rdd_make_masspoint_data, storage(`storage_type')
    rdbwdensity x, precision(`precision')
    _rdd_post_bandwidth, case("rdbwdensity_masspoints_adjusted") handle(`handle') precision(`precision')

    _rdd_make_masspoint_data, storage(`storage_type')
    rdbwdensity x, nomasspoints precision(`precision')
    _rdd_post_bandwidth, case("rdbwdensity_masspoints_unadjusted") handle(`handle') precision(`precision')
}

_rdd_make_data, storage(double)
rddensity x, h(0.6 0.6) nobinomial
if "`e(precision)'" != "double" {
    display as error "rddensity default precision should be double."
    exit 9
}
rdbwdensity x
if "`e(precision)'" != "double" {
    display as error "rdbwdensity default precision should be double."
    exit 9
}

capture which lpdensity
if _rc == 0 {
    _rdd_make_data, storage(double)
    quietly rddensity x, h(0.6 0.6) plot plot_range(-0.8 0.8) hist_range(-0.8 0.8) genvars(rdd_double) precision(double) nobinomial
    confirm double variable rdd_double_grid
    confirm double variable rdd_double_bw
    confirm double variable rdd_double_f
    confirm double variable rdd_double_hist_center
    capture drop rdd_double_*

    _rdd_make_data, storage(float)
    quietly rddensity x, h(0.6 0.6) plot plot_range(-0.8 0.8) hist_range(-0.8 0.8) genvars(rdd_single) precision(single) nobinomial
    confirm float variable rdd_single_grid
    confirm float variable rdd_single_bw
    confirm float variable rdd_single_f
    confirm float variable rdd_single_hist_center
    capture drop rdd_single_*
}

postclose `handle'

use "`results'", clear
format value %21.16g
export delimited using "`repo_root'/stata/tests/fixtures/numerical-baseline.csv", replace
