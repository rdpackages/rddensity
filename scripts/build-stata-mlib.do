version 16
clear all
set more off

args repo_root
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
local stata_dir "`repo_root'/stata"

cd "`stata_dir'"

local functions rddensity_mlib_loaded rddensity_unique rddensity_rep rddensity_fv rddensity_h rddensity_quantile

capture erase "lrddensity.mlib"
do "rddensity_fun.do"

mata: mata mlib create lrddensity, replace
foreach function of local functions {
    mata: mata mlib add lrddensity `function'()
}
mata: mata mlib index

foreach function of local functions {
    capture erase "`function'.mo"
}

confirm file "lrddensity.mlib"
display as text "Built stata/lrddensity.mlib with Stata `c(stata_version)'."
