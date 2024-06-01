*------------------------------------------------------------------------------*
**************************** Microeconometrics 20295 ***************************
******************************** PROBLEM  SET 3 ********************************
********************************************************************************
******************* Authors: Barbas, Frey, Ghelfi, Lucchini ********************
*------------------------------------------------------------------------------*
clear all 
global path "/mnt/data/Uni/2024/Micrometrics/problem_set/ps4"
*global path "C:\Users\nikib\OneDrive\Desktop\ESS\Year 1 semester 2\Microeconometrics\micrometrics_ps\ps3"
*global path "/Users/marcellalucchini/Desktop/MIRCOECONOMETRICS/PROBLEM SETS/micrometrics_ps/ps4"

global data "$path/data"
global temp "$path/output/temp"
global output "$path/output"
global tables "$path/output/tables"
global graphs "$path/output/graphs"

*uncomment the following line to install the packages

*ssc install diff, replace 
*ssc install twowayfeweights, replace
*ssc install gtools, replace
*ssc install bacondecomp, replace
*ssc install reghdfe, replace
*ssc install ftools, replace
*ssc install coefplot, replace
*ssc install avar
*ssc install eventstudyinteract


*------------------------------------------------------------------------------*
*********************************  EXCERCISE 1  ********************************
*------------------------------------------------------------------------------*

*********************************  QUESTION 1  *********************************

use "$data/pset_4", replace

*(a)

/* In our dataset, the variable stpop reports the US state population. This variable is key to weighting descriptive output and analysis by state population. For instance, Wolfers (2006) consistently report estimates using state population weights. Hence, we must do the same in reporting the evolution of divorce rates and in regressing divorce rates on unilateral divorce laws.

In principle, STATA allows four different alternative weighting methods. These are frequency weights, analytic weights, sampling weights and importance weights (Dupraz, 2013). In our case, since the divorce rates are state level averages, the analytic weights (aweight) are the relevant ones.  In fact, these should be used to weight observations as if each observation is a mean computed from a sample size n, where n is the weight variable. In brief, analytic weights should be used in samples where each observation is a group mean. 

Hence, the model will be estimated in state means. As explained in (Dupraz, 2013), this is favorable for two reasons. First, for computational efficiency. Second, because the random unit is the state. Finally, note that it can be shown that analytic weights (aweight) and frequency weights (fweight) yield the same point estimates, while they are differentiated in terms of standard errors.
*/

*(b)
* Generate the reform dummy =1 if state is reformed between 1968 and 1988	
gen reform=0
replace reform=1 if lfdivlaw>=1968 & lfdivlaw<=1988

preserve	
gen div_rate_r=div_rate if reform==1
gen div_rate_c=div_rate if reform==0

collapse div_rate_r div_rate_c [aw=stpop], by (year)	
gen div_rate_diff=div_rate_r-div_rate_c

generate upper = 8
local barcall upper year if inrange(year, 1969, 1977), bcolor(gs14) 
twoway bar  `barcall' || line  div_rate_r year, lwidth(thick) lcolor(black)|| line  div_rate_c year, lcolor(gs12) lwidth(thick) || line  div_rate_diff year, lpattern(dash) lcolor(black) xtitle("Year") ytitle("Divorce rate" "Divorces per 1,000 persons per year") xlab(, nogrid) xlabel(1956 (2) 1998, angle(45)) legend(order(1 "Reform Period (1969-1977)" 2 "Reformed States" 3 "Control States" 4 "Difference in Divorce Rates") pos(1) ring(0) col(1)) title("Figure 1:  Divorce Rates in Reformed and Control States", size(medium)) 

graph export "$output/graphs/plot_b1.png", replace

restore 

* Generate the reform dummy =1 if state is reformed between 1969 and 1973
preserve

gen reform2=. 
replace reform2=1 if lfdivlaw>=1969 & lfdivlaw<=1973
replace reform2=0 if lfdivlaw==2000

gen div_rate_r2=div_rate if reform2==1
gen div_rate_c2=div_rate if reform2==0

drop if year>1978
collapse div_rate_r2 div_rate_c2 [aw=stpop], by (year)	

line  div_rate_r2 year, xline(1969) lwidth(thick) lcolor(black) || line  div_rate_c2 year, lwidth(thick) lcolor(gs12) xtitle("Year") ytitle("Divorce rate" "Divorces per 1,000 persons per year") xlab(, nogrid) xlabel(1956 (2) 1978, angle(45)) legend(order(1 "States Treated: reformed in 1969-1977" 2 "Control States: reformed in 2000") pos(4) ring(0) col(1)) title("Figure 2: Average Divorce Rates in Reformed and Control States", size(medium)) 
graph export "$output/graphs/plot_b2.png", replace

restore 

/* The parallel pre-trends assumption appears to be valid. In fact, the graph depicting the disparity in the average outcome variable between the two groups before the introduction of unilateral divorce laws in 1968 remains relatively consistent.
*/

*(c)
preserve 
keep if year == 1968 | year == 1978
keep if (lfdivlaw >=1968 & lfdivlaw<=1973) | (lfdivlaw==2000)

gen unilateral = 0
replace unilateral = 1 if lfdivlaw >= 1969 & lfdivlaw <= 1973

gen post = 0
replace post = 1 if year == 1978
gen post_unilateral = post*unilateral

*regression (i)
reg div_rate post_unilateral post [aweight = stpop], vce(robust)
outreg2 using "$output/tables/table_c.xls", title("Regression Table Question C")  label excel replace

*regression (ii)
reg div_rate post_unilateral post unilateral [aweight = stpop], vce(robust)
outreg2 using "$output/tables/table_c.xls", title("Regression Table Question C")  label excel append

*for roboustness we also run the regression with the interaction term
reg div_rate i.post##i.unilateral [aweight = stpop], vce(robust)
outreg2 using "$output/tables/table_c.xls", title("Regression Table Question C") label excel append

*and similarly we use the diff command
diff div_rate,  t(unilateral) p(post)
outreg2 using "$output/tables/table_c.xls", title("Regression Table Question C")  label excel append

/* In the absence of state fixed effects in the model, the pooled Ordinary Least Squares (OLS) estimates reduce to solely the disparity in average outcomes between the treatment and control groups in 1978. This oversight neglects group-specific effects that persist consistently over time, potentially introducing selection bias. Drawing from the graph presented in point a), we anticipate the coefficient linked with the interaction term to be elevated in the first regression due to the presence of positive selection bias within the treatment group, as evidenced by higher mean divorce rates.

Consequently, the estimation of the causal effect stemming from the implementation of unilateral divorce laws, while maintaining the assumption of parallel trends, consistently hinges on the coefficient of the interaction term in the second regression, which is approximately -0.0050148. However, this coefficient does not significantly deviate from zero (p-value = 0.993). Furthermore, it is worth noting that the coefficient of the interaction term in regression 1 conflates the causal impact of the law with group fixed effects, which seem to be positive, yielding a positive and significant point estimate (coefficient = 1.700743 and p-value = 0.000).
 */

*(d)
matrix table_1 = J(3, 3, .)

* unilateral = 1, post = 1
qui sum  div_rate if unilateral==1 & post==1 [aw=stpop]
matrix table_1[1,1]=round(r(mean),.001)

* unilateral = 1, post = 0
qui sum  div_rate if unilateral==1 & post==0 [aw=stpop]
matrix table_1[1,2]=round(r(mean),.001)

* unilateral = 0, post = 1
qui sum  div_rate if unilateral==0 & post==1 [aw=stpop]
matrix table_1[2,1]=round(r(mean),.001)

*unilateral = 0, post = 0
qui sum  div_rate if unilateral==0 & post==0 [aw=stpop]
matrix table_1[2,2]=round(r(mean),.001)

matrix table_1[1,3]=table_1[1,1]-table_1[1,2]
matrix table_1[2,3]=table_1[2,1]-table_1[2,2]
matrix table_1[3,1]=table_1[1,1]-table_1[2,1]
matrix table_1[3,2]=table_1[1,2]-table_1[2,2]
matrix table_1[3,3]=table_1[3,1]-table_1[3,2]

matrix colnames table_1= "UNILATERAL=1" "UNILATERAL=0" "Difference 2"
matrix rownames table_1= "POST=1" "POST=0" "Difference 1"
	
* Exporting the table
putexcel set "$output/tables/table_1.xlsx", replace
putexcel A1=matrix(table_1) , names nformat(number_d2)
putexcel (A1:A4), overwr bold border(right thick) 
putexcel (A1:D1), overwr bold border(bottom thick) 
putexcel (C1:C4), border(right thick) 
putexcel (A4:D4), border(top thick)	
	
restore

*(e)
clear
use "$data/pset_4", replace
encode st, gen(state)
order state, after(st)
drop st

xtset state year
keep if year>=1956 &  year<=1988	

gen imp_unilateral = 0
replace imp_unilateral = 1 if year >= lfdivlaw
label variable imp_unilateral "State already introduced unilateral divorce laws"

qui tab year , gen(year_)
local year_fe year_*
qui tab state , gen(state_)
local state_fe state_*

forval i=1/51{
	bysort state (year): gen timetrend_lin_`i'=_n if state==`i' 
	replace timetrend_lin_`i'=0 if timetrend_lin_`i'==.
	bysort state (year): gen timetrend_sq_`i'=_n^2 if state==`i' 
	replace timetrend_sq_`i'=0 if timetrend_sq_`i'==.
}
local state_timetrend timetrend_lin_*
local state_timetrend_sq timetrend_sq_*

*regression (i)
reg div_rate imp_unilateral `year_fe' `state_fe'  [aweight = stpop], vce(cluster state)
outreg2 using "$output/tables/table_e.xls", title("Regression Table Question E") keep(imp_unilateral) label excel replace

*regression (ii) linear time trends
reg div_rate imp_unilateral `year_fe' `state_fe' `state_timetrend' [aweight = stpop], vce(cluster state)

outreg2 using "$output/tables/table_e.xls", title("Regression Table Question E") keep(imp_unilateral) label excel append

*regression (iii) quadratic time trends
gen year2 = year^2
reg div_rate imp_unilateral  `year_fe' `state_fe' `state_timetrend' `state_timetrend_sq' [aweight = stpop], vce(cluster state)

outreg2 using "$output/tables/table_e.xls", title("Regression Table Question E") keep(imp_unilateral) label excel append

/* The initial specification's results are obtained without accounting for state-specific time trends in divorce rates, thereby relying solely on a stringent assumption of parallel trends for accurately identifying the policy's causal effect. This assumption posits that, on average, the hypothetical time trend of the treatment group aligns with the observed time trend in the control group. In subsequent regression, adjustments are introduced for state-specific linear and quadratic time trends in the second and third specifications, respectively. Discrepancies among these specifications may stem from either a breach of the parallel trend assumption or a misalignment in the polynomial fitting of state-specific time trends, potentially elucidating disparities between regression 2 and 3.Should the coefficients across the three specifications remain identical, indicating no significant discrepancies, it implies that the average time trends between the groups are insubstantial, thereby bolstering the parallel trend assumption. In our specific scenario, the coefficients for the three specifications manifest as follows: (1) coefficient = -0.0511, lacking statistical significance; (2) coefficient = 0.488, significant at the 0.01% level; (3) coefficient = 0.323, also significant at the 0.01% level. The observed inconsistency among the three specifications suggests a potential violation of the parallel trend assumption, hinting that the consistent coefficient may reside in the POST_UNILATERAL interaction term in the third regression. This supposition operates under the presumption that a quadratic fitting accurately captures state-specific time trends. However, these assertions will be subject to scrutiny in point (k), where we will contend that estimating the counterfactual time trends in the outcome variable might be confounded by the dynamic effects of the treatment over time.
*/ 

*(f)
/* Creates simulated observations */
preserve 

clear
set obs 6
gen obs = _n
gen state = floor(.9 + obs/3)
bysort state: gen year = _n
gen D = state == 1 & year == 3
replace D = 1 if state == 2 & (year == 2 | year == 3)

/* Creates simulated outcomes */
set seed 1000
gen Y = 0.1 + 0.02 * (year == 2) + 0.05 * (D == 1) + uniform() / 100
set seed 1001
gen Y2 = 0.1 + 0.02 * (year == 2) + 0.05 * (D == 1) + 0.3 * (state == 2 & year == 3) + uniform() / 100
set seed 1002
gen Y3 = 0.1 + 0.02 * (year == 2) + 0.05 * (D == 1) + 0.4 * (state == 2 & year == 3) + uniform() / 100
set seed 1003
gen Y4 = 0.1 + 0.02 * (year == 2) + 0.05 * (D == 1) + 0.5 * (state == 2 & year == 3) + uniform() / 100

reg Y i.D i.state i.year 
outreg2 using "$output/tables/table_f.xls", title("Regression Table Question F")  label excel replace

reg Y2 i.D i.state i.year
outreg2 using "$output/tables/table_f.xls", title("Regression Table Question F")  label excel append

reg Y3 i.D i.state i.year
outreg2 using "$output/tables/table_f.xls", title("Regression Table Question F")  label excel append

reg Y4 i.D i.state i.year
outreg2 using "$output/tables/table_f.xls", title("Regression Table Question F")  label excel append

/* We will explore four distinct outcome variables within the initial model specification outlined in point (e). Our specific focus will be on the first specification, which exclusively integrates the UNILATERAL dummy variable. It becomes apparent that achieving consistent estimation of treatment effects proves challenging in the presence of state-specific time trends. This limitation extends to the subsequent three simulated dependent variables, where any correlation between these varied time trends and treatment assignment remains unaccounted for.

Consistent with our expectations, the coefficients derived from the regression analysis on the initial simulated dependent variable closely mirror the true treatment effect. This correspondence is attributable to the alteration in treatment assignment status for state 1 during the transition from year 2 to year 3, simultaneous with a corresponding alteration in the outcome variable for the other state, theoretically serving as its counterfactual. Consequently, given the positive treatment effect and the upward trend in state 2 during year 3, a downward bias in regression estimates ensues, a phenomenon exacerbated by larger disparities in absolute values between the trend differentials.
*/

*(g)
log using "$output//tables/table_g.log", replace
twowayfeweights Y state year D, type(feTR)
twowayfeweights Y2 state year D, type(feTR)
twowayfeweights Y3 state year D, type(feTR)
twowayfeweights Y4 state year D, type(feTR)
log close

restore

*There is no better way to our knowledge to save the results of the twowayfeweights command, so we save the log file and then extract the results from there. We opened a github issue to ask for a better solution.

*TODO add interpretation

*(h)
*all states observation start on 1956, so we can use the stpop of 1956 as the initial population
bysort state: gen init_stpop=stpop[1] if year == 1956
bysort state: replace init_stpop = init_stpop[_n-1] if missing(init_stpop)
order init_stpop, after(stpop)

*check if the initial population is correct for AK and AL
list state year init_stpop in 1/100

reg div_rate i.state i.year imp_unilateral [aweight = init_stpop], vce(cluster state)
outreg2 using "$output/tables/table_h.xls", title("Regression Table Question H")  label excel replace

bacondecomp div_rate imp_unilateral [aweight = init_stpop]
graph export "$output/graphs/plot_h.png", replace

/* This section aims at revisiting our analysis through the use of the Goodman-Bacon decomposition. The canonical DiD estimator contains two time periods and two groups. Yet, it is often the case that units receive treatment at different times. In this case, the two-way fixed effects estimator is a weighted average of all possible two-group and two-period DiD estimators in the data (DD decomposition). Weights are proportional to timing group sizes and the variance of the treatment dummy in each pair. 

Goodman and Bacon (2021) use the DD decomposition to show that the two-way fixed effects estimator estimates a variance-weighted average of the treatment effect parameters, sometimes using negative weights. Such weights are positive if treatment effects do not change over time. Instead, if average treatment effects are time varying, the weights could be negative. This is explained by the DD decomposition: when already-treated units act as controls, the changes in their outcomes are subtracted and these changes can include treatment effects that vary in time. The authors underline that while this does not imply a failure in the design in terms of parallel trends, it suggests to be cautious in the way of aggregating treatment effects when using this type of estimators. In the paper, the authors employ the DD decomposition to define common trends when using two-way fixed effects estimators to identify the variance-weighted treatment effect parameter. Finally, they develop tools to describe the estimator and discuss how estimates can change across specifications (methods are implemented using command bacondecomp). 

Running bacondecomp, we analyze the decomposition of the treatment effect. The graph shows the relationship between the 2x2 DD treatment effect estimates (on the vertical axes) and the assigned weights (horizontal axes). The graph distinguishes between three types of 2x2 comparisons: circles are terms in which one timing group acts as the treatment group and the pre-1964 reform states act as the control group. The triangles are terms in which one timing group acts as the treatment group and the non-reform states act as the control group. The x's are the timing-only terms. The dashed line indicates the two-way fixed effects estimate, that is the average of the y-axis values weighted by their x-axis values. This value is the same as the coefficient from the regression in (ii). 

We observe that the most influential 2x2 DDs are the ones determined by the "never treated vs timing". We observe that instead, the majority of the weights are concentrated around zero, and are a combination of all three types explained before. From the graph, we observe that there could be a proportion of negative weights that are however negligible in magnitude. In absolute value, the positive weights seem to be of larger magnitude than the negative ones. 

In general, the discussion of negative weights must consider (i) the presence of negative weights and (ii) the magnitude of the associated 2x2 DD estimate. A negative weight is problematic if it is associated to a large positive 2x2 DD estimate, since it would end up creating a large negative summand that could bias the final weighted average. However, this does not seem to be the case. In fact, the largest 2x2 DD estimates (above 0.5) are associated with positive weights. If anything, the presence of negative weights would be associated with the 2x2 DD estimates of lower magnitude. Hence, negative weights may not be a major problem in the present analysis. 

Technical note from Goodman and Bacon (2021) : why can these weights be negative? The probability limit of the two-way fixed effect DD estimator is decomposable into: (i) the variance weighted average treatment effect on the treated (ii) the variance weighted common trends assumed to be zero and (iii) the weighted sum of the change in treatment effects within each timing group's before and after a later treatment time. The latter can be different from zero with time-varying treatment effects, and this is a source of bias inducing negative weights.
*/

*(i)
xtset state year

*Creates variable for relative time since union entry
gen relative_divlaw = year - lfdivlaw

*Creates dummy for people who never entered the union
gen never_divlaw = 0
replace never_divlaw = 1 if lfdivlaw > 1988

tab relative_divlaw

*Creates dummies for leads and lags
forvalues k = 10(-1)2 {
gen g_`k' = relative_divlaw == -`k'
}
forvalues k = 0/15 {
gen g`k' = relative_divlaw == `k'
}

reghdfe div_rate g_* g0-g15 [aweight=stpop], absorb(`year_fe' `state_fe') vce(cluster state)
estimates store reg1
outreg2 using "$output/tables/table_i.xls",title("Regression Table Question I")  label excel replace

reghdfe div_rate g_* g0-g15  `state_timetrend' [aweight=stpop], absorb(`year_fe' `state_fe') vce(cluster state)
estimates store reg2
outreg2 using "$output/tables/table_i.xls",title("Regression Table Question I")  label excel append

reghdfe div_rate g_* g0-g15  `state_timetrend' `state_timetrend_sq' [aweight=stpop], absorb(`year_fe' `state_fe') vce(cluster state)
estimates store reg3
outreg2 using "$output/tables/table_i.xls",title("Regression Table Question I")  label excel append

*TODO add interpretation

*(j)
coefplot reg1, keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, angle(45) alternate)
graph export "$output/graphs/plot_j1.png", replace

STOP

coefplot reg2, keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, angle(45) alternate)
graph export "$output/graphs/plot_j2.png", replace

coefplot reg3, keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, angle(45) alternate)
graph export "$output/graphs/plot_j3.png", replace


*(k)

*TODO do question

*(l)
eventstudyinteract div_rate g_* g0-g15  [aweight=stpop], cohort(lfdivlaw) control_cohort(never_divlaw) absorb(`year_fe' `state_fe') vce(cluster state)

matrix C = e(b_iw)
mata st_matrix("A",sqrt(diagonal(st_matrix("e(V_iw)"))))
matrix C = C \ A'
matrix list C
coefplot matrix(C[1]), se(C[2]) keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, alternate)

graph export "$output/graphs/plot_l1.png", replace

eventstudyinteract div_rate g_* g0-g15  [aweight=stpop], cohort(lfdivlaw) covariates(`state_timetrend') control_cohort(never_divlaw) absorb(`year_fe' `state_fe') vce(cluster state)

matrix C = e(b_iw)
mata st_matrix("A",sqrt(diagonal(st_matrix("e(V_iw)"))))
matrix C = C \ A'
matrix list C
coefplot matrix(C[1]), se(C[2]) keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, alternate)

graph export "$output/graphs/plot_l2.png", replace

eventstudyinteract div_rate g_* g0-g15  [aweight=stpop], absorb(`year_fe' `state_fe') cohort(lfdivlaw) covariates(`state_timetrend' `state_timetrend_sq') control_cohort(never_divlaw)  vce(cluster state)

matrix C = e(b_iw)
mata st_matrix("A",sqrt(diagonal(st_matrix("e(V_iw)"))))
matrix C = C \ A'
matrix list C
coefplot matrix(C[1]), se(C[2]) keep(g*) vertical yline(0) xtitle("Years after law") ytitle("Estimated effect") ///
				title("Dependent variable") xlabel(, alternate)

graph export "$output/graphs/plot_l3.png", replace


*TODO add interpretation


