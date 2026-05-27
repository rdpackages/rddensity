version 16
clear all
set more off
set graphics off

args repo_root install_dir
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all

if `"`install_dir'"' == "" {
    local install_dir "`c(tmpdir)'rddensity-plus"
}
local install_dir : subinstr local install_dir "\" "/", all
capture mkdir "`install_dir'"

sysdir set PLUS "`install_dir'"
net install rddensity, from("`repo_root'/stata") replace

capture confirm file "`install_dir'/l/lrddensity.mlib"
if _rc {
    capture confirm file "`install_dir'/r/lrddensity.mlib"
}
if _rc {
    display as error "Installed package did not include lrddensity.mlib."
    exit 601
}

local mo_files rddensity_mlib_loaded.mo rddensity_unique.mo rddensity_rep.mo rddensity_fv.mo rddensity_h.mo rddensity_quantile.mo
foreach mo_file of local mo_files {
    capture confirm file "`install_dir'/r/`mo_file'"
    if !_rc {
        display as error "Installed package should not include `mo_file'."
        exit 9
    }
    capture confirm file "`install_dir'/l/`mo_file'"
    if !_rc {
        display as error "Installed package should not include `mo_file'."
        exit 9
    }
}

which rddensity
which rdbwdensity

clear
set obs 200
gen double x = -1.5 + 3 * (_n - 1) / 199 + 0.05 * sin(_n - 1)

quietly rddensity x, h(0.6 0.6) nobinomial
assert e(N_l) + e(N_r) == 200

quietly rdbwdensity x
assert e(N_l) + e(N_r) == 200

mata: st_numscalar("rddensity_mlib_check", rddensity_mlib_loaded())

display as text "Stata package install check passed."
