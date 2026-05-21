version 16
clear all
set more off

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
postfile `handle' str32 case str32 metric double value using "`results'", replace

capture program drop _rdd_post_density
program define _rdd_post_density
    syntax, CASE(string) HANDLE(name) [ALL]

    post `handle' ("`case'") ("hat_left") (e(f_ql))
    post `handle' ("`case'") ("hat_right") (e(f_qr))
    post `handle' ("`case'") ("hat_diff") (e(f_qr) - e(f_ql))
    post `handle' ("`case'") ("sd_jk_left") (e(se_ql))
    post `handle' ("`case'") ("sd_jk_right") (e(se_qr))
    post `handle' ("`case'") ("sd_jk_diff") (e(se_q))
    post `handle' ("`case'") ("test_t_jk") (e(T_q))
    post `handle' ("`case'") ("test_p_jk") (e(pv_q))
    post `handle' ("`case'") ("h_left") (e(h_l))
    post `handle' ("`case'") ("h_right") (e(h_r))
    post `handle' ("`case'") ("n_full") (e(N_l) + e(N_r))
    post `handle' ("`case'") ("n_left") (e(N_l))
    post `handle' ("`case'") ("n_right") (e(N_r))
    post `handle' ("`case'") ("n_eff_left") (e(N_h_l))
    post `handle' ("`case'") ("n_eff_right") (e(N_h_r))

    if "`all'" != "" {
        post `handle' ("`case'") ("hat_p_left") (e(f_pl))
        post `handle' ("`case'") ("hat_p_right") (e(f_pr))
        post `handle' ("`case'") ("hat_p_diff") (e(f_pr) - e(f_pl))
        post `handle' ("`case'") ("sd_jk_p_left") (e(se_pl))
        post `handle' ("`case'") ("sd_jk_p_right") (e(se_pr))
        post `handle' ("`case'") ("sd_jk_p_diff") (e(se_p))
        post `handle' ("`case'") ("test_p_t_jk") (e(T_p))
        post `handle' ("`case'") ("test_p_p_jk") (e(pv_p))
    }
end

capture program drop _rdd_post_bandwidth
program define _rdd_post_bandwidth
    syntax, CASE(string) HANDLE(name)

    matrix H = e(h)
    local rownames : rownames H
    local colnames : colnames H
    forvalues i = 1/`=rowsof(H)' {
        local row : word `i' of `rownames'
        forvalues j = 1/`=colsof(H)' {
            local col : word `j' of `colnames'
            post `handle' ("`case'") ("h_`row'_`col'") (H[`i', `j'])
        }
    }
    post `handle' ("`case'") ("n_full") (e(N_l) + e(N_r))
    post `handle' ("`case'") ("n_left") (e(N_l))
    post `handle' ("`case'") ("n_right") (e(N_r))
end

set obs 200
gen double x = -1.5 + 3 * (_n - 1) / 199 + 0.05 * sin(_n - 1)

rddensity x, h(0.6 0.6) nobinomial
_rdd_post_density, case("rddensity_fixed") handle(`handle')

rddensity x, h(0.7 0.7) kernel(uniform) nobinomial
_rdd_post_density, case("rddensity_uniform") handle(`handle')

rddensity x, h(0.8 0.8) fitselect(restricted) nobinomial
_rdd_post_density, case("rddensity_restricted") handle(`handle')

rddensity x, h(0.6 0.6) all nobinomial
_rdd_post_density, case("rddensity_all") handle(`handle') all

rdbwdensity x
_rdd_post_bandwidth, case("rdbwdensity_default") handle(`handle')

rdbwdensity x, kernel(uniform) vce(plugin)
_rdd_post_bandwidth, case("rdbwdensity_uniform_plugin") handle(`handle')

postclose `handle'

use "`results'", clear
format value %21.16g
export delimited using "`repo_root'/stata/tests/fixtures/numerical-baseline.csv", replace
