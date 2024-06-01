*****************************************
******** Microeconometrics 20295*********
************PROBLEM SET 1****************
*****************************************
** Authors: Barbas, Frey, Ghelfi, Lucchini

clear all 


    global path "/Users/marcellalucchini/Desktop/MIRCOECONOMETRICS/micrometrics_ps/ps1"
	*global path "D:/Uni/2024/Micrometrics/problem_set/ps1"
	global data "$path/DATA"
	global temp "$path/output/temp"
	global output "$path/output"
	
	*uncomment to install
	*uncomment to install
	*ssc install estout, replace
	*ssc install randomizr, replace
	*ssc install ritest, replace
	*ssc install outreg2, replace
	
*---------------------------------------*
**************  QUESTION 1  *************
*---------------------------------------*
use "$data/jtrain2.dta", replace

* (a)*
quietly {
matrix ex_1a =(.,.,.,.,.,.)

local i=1

foreach var of varlist age educ black hisp re74 re75 nodegree {

	qui ttest `var', by(train) unequal 
	
	matrix ex_1a[`i',1]=r(mu_2)
	matrix ex_1a[`i',2]=r(sd_2) 
	matrix ex_1a[`i',3]=r(mu_1)
	matrix ex_1a[`i',4]=r(sd_1) 
	matrix ex_1a[`i',5]= r(mu_1)-r(mu_2)
	matrix ex_1a[`i',6]=r(se) 
	matrix list ex_1a
		local i=`i'+1
	if `i'<=7 matrix ex_1a=(ex_1a \ .,.,.,.,.,.) 
	
}


matrix rownames ex_1a= "Age" "Education" "Black" "Hispanic" "Real Earn 74" "Real Earn 75" "No Degree"
matrix colnames ex_1a= "Mean Trt (1)" "StDev Trt" "Mean Ctrl (2)" "StDev Ctrl" "(1)-(2)" "StDev (1)-(2)"

local rows = rowsof(ex_1a)
local cols = colsof(ex_1a)

*Loop through each element of the matrix and round it to two decimal places
forval i = 1/`rows' {
    forval j = 1/`cols' {
        matrix ex_1a[`i',`j'] = round(ex_1a[`i',`j'], 0.001)
    }
}

	
	matrix list ex_1a, f(%9.3f) title("Balance check")
	
	putexcel set "$output/Table_1.xlsx", replace
	
	putexcel A1=matrix(ex_1a), names nformat(number_d2)
	putexcel (A2:A8), overwr bold border(right thick) 
	putexcel (B1:G1), overwr bold border(bottom thick) 

}

/*The treatment and control groups do not exhibit perfect balance across all seven variables. In particular, five of seven variables are balanced (Age, Education, Black, Real Earnings in 1974, Real Earnings in 1975), while two (Hispanic and Nodegree) are problematic. That is, 71.42% of the variables are balanced in the sample, which is slightly below the 80% expected and reccomended threshold. This suggests that treatment and control groups are mostly comparable. 
      
By inspecting the difference in means its standard errors, we witness that the difference in percentage of people with No Degree in the two groups is statistically significant at the 1% level. Similarly, the difference in percentage of Hispanic people in the two groups is significant at the 10% level. Instead, the remaining variables do not present statistically significant differences at standard significance levels (of 1%, 5%, 10%). 

The implication of the such differences is that there might be a correlation between being selected into the traineeship treatment and respectively, having a Degree or being Hispanic. This could be especially problematic for the Degree variable since it exhibits statistical differences at 1% significance level. This is intuitive, since having a degree is an important determinant on the decision to participate in a training program aiming at seizing regular job opportunities, and this mechanisms may yield selection. */

*(b)*
reg re78 train, vce(robust)
outreg2 using "$temp/reg_1_table", replace dta 
use "$temp/reg_1_table_dta"
export excel using "$output/Reg1_table", replace

/*The coefficient of the regression indicates that being in the treatment group (train=1), relative to the control group (train=0), increases earnings by an average of $1,794, indexed at the 1982 US price level. 

Note that the coefficient is significant at the 5% level, yielding a relatively large confidence interval of [550.57; 3038.11]. Yet, in general, the subsidised work experience seems to have a positive impact on real earnings, given the experimental framework. */


* (c) * 
use "$data/jtrain2.dta", replace
quietly {
local x_1 "train"
qui reg re78 `x_1', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (1) replace dta


local x_1 "train"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (2) append dta

local x_1 "train"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3' , vce(robust) 
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (3) append dta

use "$temp/Table_2_dta"
export excel using "$output/Table_2", replace 

}

/*The training coefficient exhibits minimal sensitivity to the introduction of covariates. In fact, the coefficient remains almost unchanged (1.794, 1.686, 1.680, respectively) and maintains significance at the 1% level. Note that also the standard errors of the coefficient are constant over time. This stability suggests that the relationship between the dependent variable and the treatment variable is robust. 

Concerning the coefficients relating to the remaining regressors, we find that they are not substantially changed in the different specifications. We observe a link between the significance of coefficients and the balance table results obtained in (1a). For example, the coefficients of Real Earnings in 1974 and 1975 are not significant; This may seem peculiar, in that the use of past outcomes is usually a strong predictor of the future outcomes. Yet, coefficients are insignificant because the treatment and control groups are balanced across those two Real Earnings variables. Similar reasonings apply to the other regressors. */

*(d)* 
use "$data/jtrain2.dta", replace

reg re78 train age educ black hisp re74 re75
dfbeta, stub(dfbeta1)

rename dfbeta11 influence_train

sort influence_train
hist influence_train 

preserve
keep if _n < _N - 3  & _n > 2
reg re78 train age educ black hisp re74 re75, vce(robust)  
outreg2 using "$temp/Table_3", ctitle (Removing highest and lowest 3) replace dta
restore 

preserve
keep if _n < _N - 5  & _n > 4
reg re78 train age educ black hisp re74 re75, vce(robust)   
outreg2 using "$temp/Table_3", ctitle (Removing highest and lowest 5) append dta
restore

preserve
keep if _n < _N - 10  & _n > 9
reg re78 train age educ black hisp re74 re75, vce(robust) 
outreg2 using "$temp/Table_3", ctitle (Removing highest and lowest 10) append dta
restore

use "$temp/Table_3_dta", replace
export excel using "$output/Table_3", replace

/* In this case, DFbeta measures the extent to which the average treatment effect estimate changes when we remove the most extreme individual treatment effects. It captures the distributional drivers of the ATE estimate. The results are slightly sensitive to influential observations. This can be seen by plotting a histogram of DFbeta. We see that the histogram is centred at zero standard deviations, suggesting that on average, the individual treatment effects coincide. We also observe that the most isolated individual treatment effect realisations can be captured by the 5 lowest and largest values. 

This is reflected in the table coefficients: the coefficients of train are relatively constant with the removal of the 3 and 5 lowest and largest individual treatment effects. Instead, the coefficient slightly decreases with the removal of the 10 lowest and largest effects, since excluding a  progressively larger number of extreme values yields a convergence of ATE towards zero. Still, the absolute difference is not too substantial. 
      
Finally, the only significant coefficients are for Black (across removals of 3,5,10 lowest and largest observations) and Educ (for the removal of 3 lowest and largest observations). Note that by removing influential observations we are shifting the regression line, therefore altering the constant. 
*/

*---------------------------------------*
**************  QUESTION 2  *************
*---------------------------------------*
use "$data/jtrain3.dta", clear 

*(a)*
quietly {
matrix ex_2a=(.,.,.,.,.,.)

local i=1

foreach var of varlist age educ black hisp re74 re75 {

	qui ttest `var', by(train) unequal
	
	matrix ex_2a[`i',1]=r(mu_2)
	matrix ex_2a[`i',2]=r(sd_2) 
	matrix ex_2a[`i',3]=r(mu_1)
	matrix ex_2a[`i',4]=r(sd_1) 
	
	matrix ex_2a[`i',5]= r(mu_1)-r(mu_2)
	matrix ex_2a[`i',6]=r(se) 
	matrix list ex_2a
		local i=`i'+1
	if `i'<=7 matrix ex_2a=(ex_2a \ .,.,.,.,.,.) 
	
}

matrix rownames ex_2a= "Age" "Education" "Black" "Hispanic" "Real Earn 74" "Real Earn 75" "No Degree"
matrix colnames ex_2a= "Mean Trt (1) 2" "StDev Trt 2" "Mean Ctrl (2) 2" "StDev Ctrl 2" "(1)-(2) 2" "StDev (1)-(2) 2"

local rows = rowsof(ex_2a)
local cols = colsof(ex_2a)

*Loop through each element of the matrix and round it to two decimal places
forval i = 1/`rows' {
    forval j = 1/`cols' {
        matrix ex_2a[`i',`j'] = round(ex_2a[`i',`j'], 0.001)
    }
}


matrix list ex_1a ex_2a, f(%9.3f) title("Balance check")
	
	putexcel set "$output/Table_1.xlsx", replace
	
		putexcel A2=matrix(ex_1a), names nformat(number_d2)
		putexcel (A3:A9), overwr bold border(right thick) 
		putexcel (B2:G2), overwr bold 
		
		putexcel I2=matrix(ex_2a), colnames nformat(number_d2)
		putexcel (I2:N2), overwr bold
		 
		putexcel B1 = "Table 1 - jtrain2.dta"
		putexcel I1 = "Table 1 - jtrain3.dta"
		putexcel B1, overwr bold 
		putexcel I1, overwr bold 
	
}
	

*(b)*
use "$data/jtrain3.dta", replace
gen treated=.
    set seed 63285
    gen random=uniform()
    sort random
	egen random_order=rank(random)
	qui sum random 
	gen N =r(N)
	replace treated=0 if random_order<=(N/2)
	replace treated=1 if random_order>(N/2) & random_order<=N

*(c)*
ssc install randtreat
randtreat, generate(treated_2) misfits(global)
pwcorr treated treated_2, sig star(.05)
tab treated_2, missing

randtreat, generate(treated_3) se(63285) misfits(global)
pwcorr treated treated_3, sig star(.05)
tab treated_2, missing

/* The correlation between treated and treated_2 is not statistically significant at the 95% level. This is indeed expected, as in treated_2 we are not assigning any value to the seed variable. If we assign the same seed, as we did for treated_3, we still do not find any significant correlation between treated and treated_3. 

The seed is in fact used to initialize the pseudorandom generator, that is used in both procedures to assign the treatment and control. The function runiform() uses a KISS alghorithm, while the function randtreat() uses a Linear Congruential Generator (LCG). Despite KISS being based on LGC, it will handle differently the same seed. In conclusion, this will produce different pseudorandom numbers. */

 
*(d)*
quietly {
matrix ex_2d=(.,.,.,.,.,.)

local i=1

foreach var of varlist age educ black hisp re74 re75 {

	qui ttest `var', by(treated) unequal
	
	matrix ex_2d[`i',1]=r(mu_2)
	matrix ex_2d[`i',2]=r(sd_2) 
	
	matrix ex_2d[`i',3]=r(mu_1)
	matrix ex_2d[`i',4]=r(sd_1) 
	
	matrix ex_2d[`i',5]= r(mu_1)-r(mu_2)
	matrix ex_2d[`i',6]=r(se) 
	matrix list ex_2d
		local i=`i'+1
	if `i'<=7 matrix ex_2d=(ex_2d \ .,.,.,.,.,.) 
	
}

matrix rownames ex_2d= "Age" "Education" "Black" "Hispanic" "Real Earn 74" "Real Earn 75" "No Degree"
matrix colnames ex_2d= "Mean Trt (1) 3" "StDev Trt 3" "Mean Ctrl (2) 3" "StDev Ctrl 3" "(1)-(2) 3" "StDev (1)-(2) 3"

local rows = rowsof(ex_2d)
local cols = colsof(ex_2d)

*Loop through each element of the matrix and round it to two decimal places
forval i = 1/`rows' {
    forval j = 1/`cols' {
        matrix ex_2d[`i',`j'] = round(ex_2d[`i',`j'], 0.001)
    }
}


matrix list ex_1a ex_2a ex_2d, f(%9.3f) title("Balance check")
	
	putexcel set "$output/Table_1.xlsx", replace
	
		putexcel A2=matrix(ex_1a), names nformat(number_d2)
		putexcel (A3:A9), overwr bold border(right thick) 
		putexcel (B2:G2), overwr bold 
		
		putexcel I2=matrix(ex_2a), colnames nformat(number_d2)
		putexcel (I2:N2), overwr bold
		
		putexcel P2=matrix(ex_2d), colnames nformat(number_d2)
		putexcel (P2:U2), overwr bold
		 
		putexcel B1 = "Table 1 - jtrain2.dta"
		putexcel I1 = "Table 1 - jtrain3.dta"
		putexcel P1 = "Table 1 - jtrain3.dta randomized"
		putexcel B1, overwr bold 
		putexcel I1, overwr bold 
		putexcel P1, overwr bold 

}

/* The findings are in line with our expectations. This fully randomized framework produces a perfectly balanced sample: across all variables, the treatment and control variables do not present statistically significant differences in means. This contrasts our findings of table (2.a) of observational data, where the sample was imbalanced across all the considered variables a part from Hispanic at the 5% level. This is intuitive for the observational versus randomic nature of the assignment. */

*(e)*
quietly {
local x_1 "treated"
qui reg re78 `x_1', vce(robust)
count if e(sample) & treated == 0
local n_ctrl= r(N)
count if e(sample) & treated == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Randomised Treatment 1) append dta

local x_1 "treated"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(robust)
count if e(sample) & treated == 0
local n_ctrl= r(N)
count if e(sample) & treated == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Randomised Treatment 2) append dta

local x_1 "treated"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3', vce(robust)
count if e(sample) & treated == 0
local n_ctrl= r(N)
count if e(sample) & treated == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Randomised Treatment 3) append dta

}

use "$temp/Table_2_dta", clear 
export excel using "$output/Table_2", replace

/* We find that the coefficient for "treated" is not significant even after the introduction of the other control variables. This is expected because the treatment was randomly assigned. As we add the observational covariates, some of them become significant, in that they become gradually more useful in explaining the results - compensating the evidenced lack of explanatory power of the treatment dummy. */

*(f)*

use "$data/jtrain3.dta", clear 

quietly {
local x_1 "train"
qui reg re78 `x_1', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Observational Control 1) append dta 

local x_1 "train"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Observational Control 2)append dta 

local x_1 "train"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_2", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Observational Control 3) append dta 

}

use "$temp/Table_2_dta", clear 
export excel using "$output/Table_2", replace

/* By contrasting columns (1)-(3) to columns (7)-(9), we observe substantial differences in all coefficients. In (1), the treatment coefficient is 1.794 (and significant at 1% level), while it changes sign and magnitude reaching -15.20 and maintaining its significance in (7). Concerning the sensitivity of coefficients to the introduction of covariates, we have that the first model features relatively stable coefficients; instead, the second model exhibits a gradual convergence of the treatment effect to a weakly positive value, together with a loss of significance of coefficients (mimiking the dynamics shown in point (2.e)). The sensitivity of coefficients to the introduction of covariates is also evident across the other variables. These differences are in line with our expectations given the fact that jtrain2 and jtrain3 are datesets of experimental versus non-experimental data. */
 
*---------------------------------------*
**************  QUESTION 3  *************
*---------------------------------------*
* (a)*
/*In the context of Neyman inference in randomized controlled experiments, the relevant estimator is the difference in sample average outcomes for the treatment and control group. Allowing for heterogeneous treatment effects, this estimator is unbiased for the Average Treatment Effect (ATE) if pure randomisation holds. This is reflected in the independence of treatment assignment and potential outcomes. Note that this is true both (i) if we directly estimate ATE of a sample, and (ii) if we have a larger population from which the sample was randomly drawn, and we are interested in such population's ATE. Instead, for the Neyman standard errors to be unbiased, we need the sample at hand to be envisioned as being randomly sampled from an infinite population. Then, inference must be carried out on this population's average treatment effect. */
	
*(b)*
use "$data/jtrain2.dta", replace

*** Hess code ***
ritest train _b[train]: ///
	reg re78 train

* running same test with 10,000 permutations
ritest train _b[train], reps(10000): ///
	reg re78 train

*** Our code ***
mat M = (.,.,.)

forvalue i=1(1)10000 {
	quietly{
	*replacign previous iterations variable
	cap drop select sample uni_sel
	*assign pseudorandom float according to uniform distribution
	gen uni_sel = runiform(0,1)
	*sort based on uni_sel variable
	gsort uni_sel, generate(select) 
	gen sample = .
	*split the sample in treatement and control
	replace sample = 1 if select <=185
	replace sample = 0 if select >185
	regres re78 sample, vce(robust)
	*create resault matrix
	mat M = M \ (scalar(`i'), scalar(_b[sample]), scalar(_se[sample]))
	display in smcl as text "`value_label'"
	}
	_dots `i' 0
}

cap drop rep coef std c

mat A=M[2..10000,.]
svmat A
rename A1 rep
rename A2 coef 
rename A3 std	

reg re78 train
gen c = .
replace c = 0 if abs(_b[train])>abs(coef)
replace c = 1 if abs(_b[train])<abs(coef)
egen n = total(c)
gen p_value = n/10000
display "p_value: " p_value "  n (number of c):" n

/* Fisher Inference provides an alternative way of computing the p-values and testing significance of treatment effects in Randomised Controlled Trials (RCTs). This method consists in re-randomizing the data by repeatedly assigning a pseudo treatment status. All the values of pseudo-randomized treatment effect will produce a distribution on which it is possible to test the null hypothesis that the observed realization is in line with the distribution of the pseudo-randomized treatment effect.
The strength of this approach comes from it not needing assumptions on a super population distribution as the total number of permutations of the data can, theoretically,be achieved.

Conversely, as the number of permutations is quite large, the approximation of Fisher inference is to settle on a number of permutations that we will see produce robust results (In our case: 10,000).
We compute the p-values according to this methodology by comparing randomized treatment effects with the observed one. Specifically, we count the number of times that the randomized treatment effects are larger than the observed treatment in absolute value. We obtain the p-value by dividing the former result by the total number of permutations. Using this approach, we obtain a p-value of 0.004, meanwhile both Athey and Imbens and Heß package obtain a value of 0.0044. This small discrepancy is indeed to be expected because of the random nature of the permutations drawings as we are far from reaching the total number of possible permutations. */
 
*(c)*
/* Athey and Imbens (2017) employ Fisher inference through repeated treatment re-assignment across the entire sample. However, this method may not accurately represent the randomization process in the NSW experiment as outlined in LaLonde (1986), where random assignment likely occurred at the site level due to the program being conducted by multiple sites at different times.
Stratification at that level reduces variance in the sampling distribution of the statistic Delta Y under the null hypothesis. Preventing extreme observations, such as when all treated units come from a site with a high Y value, is possible with random assignment at the site level. These less exterme values lead to lower Fisher p-values, making rejection of the sharp null hypothesis even less difficult compared to Athey and Imbens (2017). */

*(d)*
*1*
/* The standard STATA method of calculation for heteroscedasticity-consistent standard errors (HC1) engages in the estimation of the asymptotic variance of OLS coefficients through the substitution of the error term's variance-covariance matrix sigma^2, with the diagonal matrix of squared residuals S^2, multiplying the outcome by a correction factor equal to the degrees of freedom of the residual vector, n/n – k. HC1 ensures consistency by using the consistent estimator S^2, for the error term's variance sigma^2. 
Alternatively, defining X as the matrix of covariates; H as the projection matrix H= X(X'X)^(-1)X' with ith diagonal entry h_ii. HC3 adjusts each squared residual by a factor dependent on the projection matrix H: S^2 (e_i^2) by 1/(1-h_ii)^2. This tuning aims to reduce the impact of outliers and enhance standard error estimation, particularly in sample sizes with observations lower or equal than 250, and when heteroskedasticity is prominent.*/

*2*
* Standard regression
use "$data/jtrain2.dta", replace
quietly {
local x_1 "train"
qui reg re78 `x_1', vce(robust)  
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Robust SE1) replace dta


local x_1 "train"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Robust SE 2) append dta

local x_1 "train"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3', vce(robust)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Robust SE 3) append dta
}

use "$temp/Table_4_dta", clear 
export excel using "$output/Table_4", replace

use "$data/jtrain2.dta", replace
reg re78 train age educ black hisp re74 re75
dfbeta, stub(dfbeta1)

rename dfbeta11 influence_train

sort influence_train
hist influence_train 

preserve
keep if _n < _N - 3  & _n > 2
reg re78 train age educ black hisp re74 re75, vce(robust)  
outreg2 using "$temp/Table_5", ctitle (Robust SE 1) replace dta
restore 

preserve
keep if _n < _N - 5  & _n > 4
reg re78 train age educ black hisp re74 re75, vce(robust)   
outreg2 using "$temp/Table_5", ctitle (Robust SE 2) append dta
restore

preserve
keep if _n < _N - 10  & _n > 9
reg re78 train age educ black hisp re74 re75, vce(robust) 
outreg2 using "$temp/Table_5", ctitle (Robust SE 3) append dta
restore

use "$temp/Table_5_dta", replace
export excel using "$output/Table_5", replace


* Using HC3
use "$data/jtrain2.dta", replace
quietly {
local x_1 "train"
qui reg re78 `x_1', vce(hc3)   
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (HC3 SE 1) append dta


local x_1 "train"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(hc3) 
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (HC3 SE 2) append dta

local x_1 "train"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3', vce(hc3)
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (HC3 SE 3) append dta
}

use "$temp/Table_4_dta", clear 
export excel using "$output/Table_4", replace

use "$data/jtrain2.dta", replace
reg re78 train age educ black hisp re74 re75
dfbeta, stub(dfbeta1)

rename dfbeta11 influence_train

sort influence_train
hist influence_train 

preserve
keep if _n < _N - 3  & _n > 2
reg re78 train age educ black hisp re74 re75, vce(hc3)  
outreg2 using "$temp/Table_5", ctitle (HC3 SE 1) append dta
restore 

preserve
keep if _n < _N - 5  & _n > 4
reg re78 train age educ black hisp re74 re75, vce(hc3)   
outreg2 using "$temp/Table_5", ctitle (HC3 SE 2) append dta
restore

preserve
keep if _n < _N - 10  & _n > 9
reg re78 train age educ black hisp re74 re75, vce(hc3) 
outreg2 using "$temp/Table_5", ctitle (HC3 SE 3) append dta
restore

use "$temp/Table_5_dta", replace
export excel using "$output/Table_5", replace

*3*
*Bootstrap
use "$data/jtrain2.dta", replace
quietly {
local x_1 "train"
qui reg re78 `x_1', vce(bootstrap, reps(1000))   
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Bootstrap SE 1) append dta


local x_1 "train"
local x_2 "age educ black hisp"
reg re78 `x_1' `x_2', vce(bootstrap, reps(1000)) 
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Bootstrap SE 2) append dta

local x_1 "train"
local x_2 "age educ black hisp"
local x_3 "re74 re75"
reg re78 `x_1' `x_2' `x_3', vce(bootstrap, reps(1000))
count if e(sample) & train == 0
local n_ctrl= r(N)
count if e(sample) & train == 1
local n_trt= r(N)
outreg2 using "$temp/Table_4", addstat("Number Treated",`n_trt', "Number Control",`n_ctrl') ctitle (Bootstrap SE 3) append dta
}

use "$temp/Table_4_dta", clear 
export excel using "$output/Table_4", replace

use "$data/jtrain2.dta", replace
reg re78 train age educ black hisp re74 re75
dfbeta, stub(dfbeta1)

rename dfbeta11 influence_train

sort influence_train
hist influence_train 

preserve
keep if _n < _N - 3  & _n > 2
reg re78 train age educ black hisp re74 re75, vce(bootstrap, reps(1000))  
outreg2 using "$temp/Table_5", ctitle (Bootstrap SE 1) append dta
restore 

preserve
keep if _n < _N - 5  & _n > 4
reg re78 train age educ black hisp re74 re75, vce(bootstrap, reps(1000))   
outreg2 using "$temp/Table_5", ctitle (Bootstrap SE 2) append dta
restore

preserve
keep if _n < _N - 10  & _n > 9
reg re78 train age educ black hisp re74 re75, vce(bootstrap, reps(1000)) 
outreg2 using "$temp/Table_5", ctitle (Bootstrap SE 3) append dta
restore

use "$temp/Table_5_dta", replace
export excel using "$output/Table_5", replace

/*Bootstrapping is a non-parametric method used to calculate standard errors for a statistic by resampling from a sample. The central assumption is that the sample represents the population well in terms of distribution of the variable of interest. By repeatedly drawing and computing statistics from resampled samples, a valid estimate of the sampling distribution and moments can be obtained, with its  variance and standard deviation. If the sample accurately reflects the population, the standard deviation of the bootstrapping distribution can be used as the standard error for the statistic.*/

*4* 
/* The regression performed with alternative methods of HC3 and HC1 correcting for potential heteroskedasticity yield the same results as before. In fact, we observe the same magnitude in the coefficents and the statistical significance for the regression adding the controls decreases from 1% to 5%. The standard errors change as expected, and the most significant change is witnessed between the standard regression and the HC1. Note that HC1 and HC3 yield similar standard errors. These results are consistent with the Data Colada post, since the sample size in this exercise is above 250, meaning that HC1 and HC3 should not be significantly different. HC3 is designed to mitigate the excessive impact of observations with extremely high variance (Long, Ervin 2000). This implies that observations heavily shaping the regression coefficient may have their influence reduced by HC3 standard errors. Given that the effect of outliers on standard errors was already minimal even with baseline robust standard errors, it's not surprising that the introduction of HC3 has little impact.

The bootstrapping procedure also yielded minimal deviations in point estimates of standard errors and significance levels. This bolsters the robustness of the baseline estimates as they are not dependent on any Central Limit Theorem (CLT) or parametric assumptions unlike the baseline ones. Since there was hardly any alteration in the standard errors and the p-values, our inferential assessment regarding the impact of NSW on earnings remains consistent.     
*/






