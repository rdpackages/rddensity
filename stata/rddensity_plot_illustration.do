******************************************************************************** 
** RDDENSITY Stata Package
** Do-file for RDDENSITY Plot Illustration
** Authors: Matias D. Cattaneo, Michael Jansson and Xinwei Ma
********************************************************************************
** net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
** net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace
********************************************************************************
clear all
set more off
 
********************************************************************************
** Load data 
********************************************************************************
use "rddensity_senate.dta", clear

********************************************************************************
** Default RDDENSITY Plot
********************************************************************************
preserve

capture drop temp_*
qui rddensity margin, plot plot_range(-50 50) hist_range(-50 50) genvars(temp)

local ci_plot_region_l = `"(rarea temp_cil temp_cir temp_grid if temp_group == 0, sort lcolor(white%0) color(red%30))"'
local ci_plot_region_r = `"(rarea temp_cil temp_cir temp_grid if temp_group == 1, sort lcolor(white%0) color(blue%30))"'

local es_plot_line_l = `"(line temp_f temp_grid if temp_group == 0, sort lcolor(red) lwidth("medthin") lpattern(solid))"'
local es_plot_line_r = `"(line temp_f temp_grid if temp_group == 1, sort lcolor(blue) lwidth("medthin") lpattern(solid))"'

qui su temp_hist_width if temp_hist_group == 0
local hist_width_l = r(mean)
qui su temp_hist_width if temp_hist_group == 1
local hist_width_r = r(mean)

local plot_histogram_l = `"(bar temp_hist_height temp_hist_center if temp_hist_group == 0, barwidth(`hist_width_l') color(red%20))"'
local plot_histogram_r = `"(bar temp_hist_height temp_hist_center if temp_hist_group == 1, barwidth(`hist_width_r') color(blue%20))"'

local graph_opt = `"xline(0, lcolor(black) lwidth(medthin) lpattern(solid)) legend(off) title("Manipulation Testing Plot", color(gs0)) xtitle("margin") ytitle("")"'

twoway 	`plot_histogram_l' 	/// 
				`plot_histogram_r' 	/// 
				`ci_plot_region_l' 	///
				`ci_plot_line_l'   	///
				`ci_plot_ebar_l'   	///
				`ci_plot_region_r' 	///
				`ci_plot_line_r'   	///
				`ci_plot_ebar_r'   	///
				`es_plot_line_l' 	///
				`es_plot_point_l' 	///
				`es_plot_line_r' 	///
				`es_plot_point_r' 	///
				,					///
				`graph_opt'
restore
