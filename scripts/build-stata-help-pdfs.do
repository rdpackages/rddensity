version 13
clear all
set more off

args repo_root
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
local stata_dir "`repo_root'/stata"

cd "`stata_dir'"

local helpfiles rddensity rdbwdensity
foreach helpfile of local helpfiles {
    display as text "Building `helpfile'.pdf"
    translate "`helpfile'.sthlp" "`helpfile'.pdf", replace translator(smcl2pdf)
}

display as text "Stata help PDF generation complete."
