*------------------------------------------------------------------------------*
**************************** Microeconometrics 20295 ***************************
******************************** PROBLEM  SET 2 ********************************
********************************************************************************
******************* Authors: Barbas, Frey, Ghelfi, Lucchini ********************
*------------------------------------------------------------------------------*
clear all 

global path "/Users/marcellalucchini/Desktop/MIRCOECONOMETRICS/PROBLEM SETS/micrometrics_ps/ps2"
global data "$path/data"
global temp "$path/output/temp"
global output "$path/output"
adopath + "$path/data"
	
*uncomment to install
* ssc install ivreg2
* ssc install ranktest

*------------------------------------------------------------------------------*
*********************************  EXCERCISE 1  ********************************
*------------------------------------------------------------------------------*

*********************************  QUESTION 1  *********************************
* (a) *
use "$data/pset_2_q_1", replace
sort birthdate

preserve
keep if birthdate >= 1930 & birthdate < 1940
collapse (mean) Education_mean = Education (mean) birthqtr_1 = birthqtr, by(birthdate)
label var Education_mean "Years of Completed Education" 
label var birthdate "Year of Birth" 
twoway connected Education birthdate, xtitle("Year of Birth") ytitle("Education") mlabel(birthqtr_1) msymbol(S) lcolor(gs0) mcolor(gs0) mlabcolor(gs0) mlabposition(6) plotregion(lcolor(black) lwidth(thin)) ylab(, nogrid) xlab(, nogrid) 
graph export "$output/fig1.png", replace
restore

preserve
keep if birthdate >= 1940 & birthdate < 1950
collapse (mean) Education_mean = Education (mean) birthqtr_1 = birthqtr, by(birthdate)
label var Education_mean "Years of Completed Education"
label var birthdate "Year of Birth" 
twoway connected Education birthdate, xtitle("Year of Birth") ytitle("Education") mlabel(birthqtr_1) msymbol(S) lcolor(gs0) mcolor(gs0) mlabcolor(gs0) mlabposition(6) plotregion(lcolor(black) lwidth(thin)) ylab(, nogrid) xlab(, nogrid) 
graph export "$output/fig2.png", replace
restore

preserve
keep if birthdate >= 1950 & birthdate <= 1960
collapse (mean) Education_mean = Education (mean) birthqtr_1 = birthqtr, by(birthdate)
label var Education_mean "Years of Completed Education"
label var birthdate "Year of Birth" 
twoway connected Education birthdate, xtitle("Year of Birth") ytitle("Education") mlabel(birthqtr_1) msymbol(S) lcolor(gs0) mcolor(gs0) mlabcolor(gs0) mlabposition(6) plotregion(lcolor(black) lwidth(thin)) ylab(, nogrid) xlab(, nogrid) 
graph export "$output/fig3.png", replace
restore

/* Figures I, II, III document the relationship between education and season of birth for for men born between 1930 and 1959. The figures plot average years of completed schooling on the vertical axis, and the year and quarter of birth on the horizontal axis. Angrist and Kruger used this as a first graphical representation of the relevance assumption, to identify a pattern between the instrument of interest (quarter of birth) and the endogenous variable (years of completed education). The idea is that individuals born in the beginning of the year start school at an older age, and can therefore drop out after completing less schooling than individuals born near the end of the year.
We observe three main trends. First, average schooling is higher for individuals born near the end of the year (Q3, Q4) than for those born earlier in the same year (Q1, Q2). Secondly, a part from a unique outlier in 1945/Q4-1946/Q1, men born in the last quarter have more education than men born in the beginning of the following year. In most cases, this holds also for third quarter versus next first quarter comparisons. Finally, the seasonal education patter is exhibited irrespective of the secular trend in education (increasing in first two decades and decreasing in the third). Note that in the original paper, the authors produce a second detrended version of the figures by removing the trend in years of education across cohorts. 
Hence, we observe a first quarter effect of lower education relative to subsequent quarters within the same year, which provides graphical evidence in favor of the relevance assumption. */ 

* (b) *
/* In a model that estimates the effect of education on health, the former variable may be endogeneous. For instance, genetic predisposition may be a latent unobservable variable in the error term, affecting simultaneously health outcomes and education. Hence, an IV strategy may be considered to instrument education through quarter of birth. Relevance of the instrument was assessed in point (1a), while exogeneity is left to be discussed. Exogeneity is composed of random assignment and exclusion restriction. 
Concerning exclusion, there is some evidence of differences in physical and mental health of individuals born at different times of the year. This would suggest that the quarter of birth is not an exogenous instrument. For instance, individuals born in earlier quarters are more likely to suffer from schizophrenia (O'Callaghan et al., 1991). Similarly, there are quarter of year effects on mental retardation (Knoblock and Pasamanick, 1958), autism (Gillberg 1990), dyslexia (Livingston, Adam, Bracha, 1993), multiple sclerosis (Templer et al., 1991), manic depression (Hare,1975) and IQ (Whorton and Karnes, 1981). Hence, the quarter of birth instrument does not exhibit a unique educational causal channel towards health outcomes, but rather multiple ones, compromising the validity of exclusion and exogeneity assumptions. 
Concerning as good as random assignment, some minor issues may arise from manipulation of quarter of birth. In principle, exact manipulation is not expected. However, manipulation could occur for several reasons related to parental leave being more convenient in specific quarters; or, being aware of literature related to seasonality of health disorders in relation with pregnancy - as described above. 

In the Angrist and Kruger setting with earning returns to education, the discussion is similar. Random assignment could still be subject to the threat of manipulation for the mechanisms described above. An example could be that high-income households could spend more resources in family planning, resulting in non-random assignment in quarter of birth. 
Instead, the issue of exclusion, in addition to the health channels illustrated above, could be threatened by three alternative channels directly linking quarter of birth to earning outcomes, as argued by Bound and Jaeger (1995).  First, there is evidence of quarter of birth affecting students' school performance through a variety of channels as attendance rates or behavioral diffulties (Carroll, 1992). Second, they suggest the presence of regional patterns in birth seasonality (Lam and Miron,1991). Finally, there is evidence on high income families being less likely to be born in winter months (Kestenbaum, 1987). This could be a threat, since being born in the first quarter may simultaneously be related with lower education and lower family income, both leading to lower income returns through the intergenerational and educational channels, respectively. Bound and Jaeger's analysis suggest that differences in family income at time of birth seem to account for virtually all the association between quarter of birth and wages. 
Hence, exogeneity and exclusion are likely violated. */

* (c) *
/* The previous discussion established imperfect exogeneity in OLS and IV settings, introducing bias in both cases. Inconsistency in OLS and IV estimates may arise from, respectively, non-zero corr(epsilon,D) and non-zero corr(epsilon,Z). The sign and magnitude of such correlations determins the overestimation or underestimation of the coefficients. 
In a simple regression model of health outcomes and education, the OLS estimates will likely be positively biased due to the discussed endogeneity issue. In fact, education may be positively correlated with essential unobserved characteristics of the individual (ex. latent family income). This positive correlation, mediated by a standard deviation ratio, will induce an overestimation of coefficients. 

The 2SLS should eliminate this bias to the extent that exogeneity and relevance assumptions hold, inducing 2SLS coefficients to be of lower magnitude compared to the biased OLS. However, as we argued before, the exogeneity and relevance assumptions could be problematic if we instrument education with quarter of birth. Hence, our framework does not entail a perfectly exogenous instrument, introducing a bias also in the 2SLS estimates. As the issues related to the health are prominent for those born in the first quarters of the year, this will lead 2SLS coefficients to also be upwardly biased. Note that comparing the relative biases of OLS and 2SLS will reveal the best estimator in this case of imperfect exogeneity, and using the instrument is not necessarily a dominating strategy. In this framework, compliers are those that change their years of schooling (D) as their quarter of birth deviates from the first (Z). Specifically, an individual is a complier if, upon reaching the age of statutory compulsory education would drop out if born in Q1 while would be forced to stay in school if born later in the year (Q2, Q3, Q4). */

*********************************  QUESTION 2  *********************************
use "$data/pset_2_q_2_and_3", replace
label var Education "Years of Education"
label var Healthy "Health"

* (a) *
sum Healthy, d
return list
scalar mu_y=r(mean)

sum Education, d
return list
scalar mu_x=r(mean)

* (b) *
tabulate birthqtr, generate(Quarter)
tabulate region, generate(region_)
local Controls "Central Married region_*"
tabulate birthyear, generate(birthyear_)
local Birth_Year_FEs "birthyear_*"

* (c) and (d) *
qui reg Healthy Education, vce(robust)
outreg2 using "$temp/Table_Q_2", label cttop(OLS) replace dta ///
keep(Education) nocons addtext(Controls, NO, Year of Birth FEs, NO) addstat("Mean y", mu_y, "Mean x", mu_x) 

reg Healthy Education `Controls', vce(robust)
outreg2 using "$temp/Table_Q_2", label cttop(OLS) append dta ///
keep(Education) nocons addtext(Controls, YES, Year of Birth FEs, NO) addstat("Mean y", mu_y, "Mean x", mu_x)

reg Healthy Education `Controls' `Birth_Year_FEs', vce(robust) 
outreg2 using "$temp/Table_Q_2", label cttop(OLS) append dta ///
keep(Education) nocons addtext(Controls, YES, Year of Birth FEs, YES) addstat("Mean y", mu_y, "Mean x", mu_x)

use "$temp/Table_Q_2_dta", replace
export excel using "$output/Table_Q_2", replace 

* (e) and (f) *
use "$data/pset_2_q_2_and_3", replace
label var Education "Years of Education"
label var Healthy "Health" 

sum Healthy, d
return list
scalar mu_y=r(mean)

sum Education, d
return list
scalar mu_x=r(mean)

tabulate birthqtr, generate(Quarter)
tabulate region, generate(region_)
local Controls "Central Married region_*"
tabulate birthyear, generate(birthyear_)
local Birth_Year_FEs "birthyear_*"

ivreg2 Healthy (Education = Quarter1 Quarter2 Quarter3), robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_2", label cttop(IV) nor2 append dta ///
keep(Education) nocons addtext(Controls, NO, Year of Birth FEs, NO) addstat("Mean y", mu_y, "Mean x", mu_x, "F-statistic IVs", F_weak)

ivreg2 Healthy (Education = Quarter1 Quarter2 Quarter3) `Controls', robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_2", label cttop(IV) nor2 append dta ///
keep(Education) nocons addtext(Controls, YES, Year of Birth FEs, NO) addstat("Mean y", mu_y, "Mean x", mu_x, "F-statistic IVs", F_weak)

ivreg2 Healthy (Education = Quarter1 Quarter2 Quarter3) `Controls' `Birth_Year_FEs', robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_2", label cttop(IV) nor2 append dta ///
keep(Education) nocons addtext(Controls, YES, Year of Birth FEs, YES) addstat("Mean y", mu_y, "Mean x", mu_x, "F-statistic IVs", F_weak)

use "$temp/Table_Q_2_dta", replace
export excel using "$output/Table_Q_2", replace  

*********************************  QUESTION 3  *********************************
use "$data/pset_2_q_2_and_3", replace

tabulate birthqtr, generate(Quarter)
tabulate region, generate(region_)
local Controls "Central Married region_*"
tabulate birthyear, generate(birthyear_)
local Birth_Year_FEs "birthyear_*"

label var Education "Years of Education"
label var Healthy "Health"
label var Quarter1 "Quarter 1"
label var Quarter2 "Quarter 2"
label var Quarter3 "Quarter 3"

sum Healthy, d
return list
scalar mu_y=r(mean)

sum Education, d
return list
scalar mu_x=r(mean)

* (a) *
reg Healthy Education `Controls' `Birth_Year_FEs', vce(robust) 
outreg2 using "$temp/Table_Q_3", label cttop(OLS) replace dta ///
keep(Education) nocons addstat("Mean Y",mu_y, "Mean X", mu_x) addtext(Controls, YES, Year of Birth FEs, YES)

* (b) *
ivreg2 Healthy (Education= Quarter1 Quarter2 Quarter3) `Controls' `Birth_Year_FEs', robust first savefirst
scalar F_weak = e(widstat)
est restore _ivreg2_Education
outreg2 using "$temp/Table_Q_3", label cttop(First Stage) append dta ///
keep(Quarter*) addtext(Controls, YES, Year of Birth FEs, YES) addstat("Mean X", mu_x, "F-Statistic IVs", F_weak) 

* (c) *
reg Healthy  Quarter1 Quarter2 Quarter3 `Controls' `Birth_Year_FEs', vce(robust) 
outreg2 using "$temp/Table_Q_3", label cttop(Reduced Form) append dta ///
keep(Quarter1 Quarter2 Quarter3) nocons addtext(Controls, YES, Year of Birth FEs, YES) addstat("Mean Y", mu_y, "Mean X", mu_x)

/* Based on the IV regression results (2e), we can see that education has a positive impact on health, when instrumenting education through quarter of birth. Given this result, we expect the first quarter to have a negative influence on health (through a channel of reduced education), with a progressive phasing out of the effect for subsequent quarters (Q2, Q3). This is verified in the data. Moreover, we observe that as quarters of birth increase, the significance levels fall, until reaching insignificance at Q3. Note that in this regression, the fourth quarter is taken as the reference, which motivates the negative signs. Note that we entered the realm of multiple instruments for education. In this case, we have positive IV coefficients and negative reduced form coefficients, inducing us to expect that reduced form coefficients will be negative. */

* (d) *
ivreg2 Healthy (Education = Quarter1 Quarter2 Quarter3) `Controls' `Birth_Year_FEs', robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_3", label cttop(IV) append dta ///
keep(Education) nocons addstat("Mean Y", mu_y, "Mean X", mu_x, "F-Statistic IVs", F_weak) addtext(Controls, YES, Year of Birth FEs, YES)

use "$temp/Table_Q_3_dta", replace
export excel using "$output/Table_Q_3", replace 

* (e) * 
/* Bound et al. (1995) discuss the problems of weak instruments and finite sample bias. Instruments are weak when they are weakly correlated with the endogenous variable. When instruments are weak, even a slight endogeneity of the IV can lead to large inconsistencies of IV estimates. Moreover, in finite samples, IV estimates are biased in the same direction as OLS, and the bias magnitude converges to that of OLS as the first stage R2 approaches zero (i.e., the weaker the instrument). They suggest that large sample sizes do not necessarily solve such issues, and that the partial R2 and F statistics can be used to assess the quality of instruments. 
In our case, the birth quarter instrument has a low correlation with the years of education, despite being significant. In fact, the weak explanatory power of the quarter of birth for education is reflected in the first stage R2, equal to 0.0004. Hence, we seem to be in the realm of weak instruments, and this can exacerbate any problems associated with a correlation between the instrument and the error. 
Moreover, as argued before, we have an imperfectly exogenous instrument, as there is a correlation between quarter of birth and health outcomes. If such correlation is sufficiently large, this may yield largely biases in IV estimates. Instead, the finite sample bias arises as we need to estimate the first stage coefficients. Such bias is exacerbated by small samples (which is not our case) and weak instruments. 

In practice, the presence of weak instruments can be detected through an F test of joint significance of the instrumental variables in the first stage. According to Staiger and Stock (1997), an instrument is strong if the statistic is larger than 10. This criterion is met in 2SLS regression, with an F statistic of 61.05. Hence, we reject the null of joint insignificance of the instruments, meaning that the IVs jointly explain the endogenous variable. This would suggest the presence of a strong instrument in this overidentified model, which contradicts the assessment made using the R2. However, one must note that the F-test is a test of significance whose power is significantly affected by the sample size (larger size increases the likelihood of rejecting the null). Hence, in this sense, we argue that the pure F-test is a noisy indicator of the presence of weak instruments. For this reason, Stock and Yogo (2005) propose alternative critical values accounting for different sample sizes, number of instruments and covariates. Such critical values slightly attenuate the weakness of the instruments. 

The F statistic contains information on the magnitude of the finite sample bias. Being the statistic far from unity, such bias should not be a concern, according to Bound and Jaeger (1995). Note that the F statistic is useful to estimate the ratio of the concentration parameter to the number of parameters. The latter captures finite sample bias, hence such statistic suggests a low likelihood of such bias. */

* (f) *
use "$data/pset_2_q_2_and_3", replace
tabulate birthqtr, generate(Quarter)
tabulate region, generate(region_)
local Controls "Central Married region_*"
tabulate birthyear, generate(birthyear_)
local Birth_Year_FEs "birthyear_*"

label var Education "Years of Education"
label var Healthy "Health"
label var Quarter1 "Quarter 1"
label var Quarter2 "Quarter 2"
label var Quarter3 "Quarter 3"

sum Healthy, d
return list
scalar mu_y=r(mean)

sum Education, d
return list
scalar mu_x=r(mean)

sum bpl, d
return list 
tab bpl if bpl != r(max), gen (State_FE_)
local State_FEs "State_FE_*"

sum birthdate, d
return list 
tab birthdate if birthdate != r(max), gen (Year_Quarter_)
local Year_Quarter_FEs "Year_Quarter_*"

egen State_Quarter = group(birthqtr bpl)
sum State_Quarter, d
return list
tab State_Quarter if State_Quarter != r(max), gen (State_Quarter_)
local State_Quarter_FEs "State_Quarter_*"

* (g) and (h) *
ivreg2 Healthy (Education = `Year_Quarter_FEs') `Controls' `Birth_Year_FEs', robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_3_g", label cttop(IV) nor2 replace dta ///
keep(Education) nocons addstat("Mean Y", mu_y, "Mean X", mu_x, "F-Statistic IVs", F_weak) addtext(Controls, YES, Year of Birth FEs, YES, State FEs, NO, Instrument, Year Quarter FEs)

ivreg2 Healthy (Education = `State_Quarter_FEs') `Controls' `Birth_Year_FEs' `State_FEs', robust first savefirst
scalar F_weak = e(widstat)
outreg2 using "$temp/Table_Q_3_g", label cttop(IV) nor2 append dta ///
keep(Education) nocons addstat("Mean Y", mu_y, "Mean X", mu_x, "F-Statistic IVs", F_weak) addtext(Controls, YES, Year of Birth FEs, YES, State FEs, YES, Instrument, State Quarter FEs)

use "$temp/Table_Q_3_g_dta", replace
export excel using "$output/Table_Q_3_g", replace 

/* The F statistic for the excluded instruments from the previous regression are, respectively, 7.672 and 3.205. Compared to the same statistic computed in (3e), the estimate falls much closer to the cutoff for finite sample bias, indicating a higher risk for it. To account for the State and Year trends in quarter of birth, interaction dummies are included as instruments. This poses a trade off between improving the first stage and decreasing the degrees of freedom that could likely exacerbate finite sample bias. The fall in the F statistic is due to such finite sample bias. */ 

*------------------------------------------------------------------------------*
*********************************  EXCERCISE 2 *********************************
*------------------------------------------------------------------------------*

*********************************  QUESTION 1  *********************************
* (a) *
/* Autor et al. (2013) seek to identify the causal effect of Chinese import exposure on US local labor markets. They leverage on the nature of the causes of the increase in Chinese imports, that are supply-driven (for instance, the rising competitiveness of Chinese manufacturers and the lowering of China's trade barriers). The employed measure of local labor market exposure to import competition is the change in Chinese import exposure per worker in a region. This is computed as the sum of the changes in US imports from China in the various industries, using as weights the region's share of national industry employment in industry j. In order to solve the likely endogeneity related to the positive correlation between US imports from China and industry import demand shock, they use an IV approach. Specifically, they instrument the growth in Chinese imports to US using the change in Chinese import exposure per worker in a region computed this time using contemporaneous industry-level growth of Chinese imports in eight other developed countries and employment levels from the prior decade. The latter aims at mitigating the potential simultaneity bias coming from potential employment response in anticipation of China trade shocks. 

Using this methodology requires that import demand shocks in high-income countries are not the primary cause of China's export surge. Such main causes must be China's productivity growth triggered by the transition to a market economy and a reduction in China's trade costs. In this context, the IV strategy identifies the components of US import growth to the extent that the common within-industry component of rising Chinese imports to the United States and other high-income countries stems from China's rising comparative advantage and (or) falling trade costs in these sectors. Hence, identification requires the inspection of three aspects. First, product demand shocks should not be highly correlated across high-income countries. Second, Chinese import growth must be driven by productivity growth in China, and not by negative productive shock in US. Finally, Chinese import rise should not be driven by the technology shocks of high-income countries negatively affecting their labor-intensive industries. 

Now, Pinkham et al. (2020) shows that the endogenous variable and the instrument are of Bartik form. A Bartik instrument is an IV that uses the inner product structure of the endogenous variable to construct an instrument. This instrument or similar ones are constructed as a sum across industries of products between shares (the shares of county i's employment in industry j) and shocks (US imports from industry j in China). The source of exogeneity can come either from the shares or from the shocks. Symmetrically, identification and consistency can come either from the shares - if they are randomly assigned across counties (as considered in Pinkham et al. (2020)), or from the shocks - if they are random across industries (as considered in Borusyak et al. (2020)). Pinkham et al. (2020) maintain that the paper can be assessed following the shares approach. In fact, the paper does not emphasize having a large number of independent shocks; the focus of the paper is on particular industries - and this is  inconsistent with the random shock explanation of the consistency of the estimators; finally, they suggest that in a setting with a fixed number of time periods and a large number of locations and industries it is appropriate to establish consistency in terms of shares. 

Consider the difference between the TSLS estimator and the parameter of interest, and let the conditions for consistency be stated in terms of the shares. In this context, the relevance assumption maintains that the denominator of the ratio converges to a nonzero term. More precisely, there must be an industry and a time period when the industry share in national employment has predictive power for the import exposure and the growth of US imports from China cannot weight the covariances between industries and time periods. Moreover, strict exogeneity must hold, meaning that the numerator of the abovementioned ratio must converge to zero. This assumption is equivalent to an exclusion restriction, in that the industry share must be uncorrelated with the structural error term after having controlled for those industries that have nonzero growth rates in imports from China. That is, the effect on employment of exposure to Chinese imports should only arise through the effect on import competition. */

* (b) *
/* In answering this question, we will be following the analysis of Pinkham et al. (2020). The authors impose a restricted form of linear heteroegenity, with constant effects within each county. In this scenario, the model must be expanded through the inclusion of county i specific coefficients. In addition, the model must specify a linear relationship between the import exposure in a county and the industry-location employment shares for all industries in that county. The additional assumptions necessary for consistency are: (i) for each industry k, the coefficients of the shares of that idustry in all the regressions for each county, have all the same sign (weakly); (ii) the conditional expectation of the product between the industry-location shares, the coefficients of interest and the location-industry specific errors, is equal to zero.

Under these assumptions, the estimator of the industry specific coefficient (beta_k hat) converges in probability to the expected value of a convex combination of the beta l, being written as the product of the beta_l times a term omega_kl, which is proportional to the quadratic mean deviation of z_lk. This points to the existence of an hetereogenous treatment effect, since the coefficients of the beta_k converge to a different combination, which varies in omega_kl. */

*********************************  QUESTION 2  *********************************
clear all
set matsize 2000

* (a) * 
/*** AKM ADH Data **/
insheet using "$data/ADHdata_AKM.csv", clear
gen year = 1990 + (t2=="TRUE")*10
drop t2

/*** BHJ SHARES **/
merge 1:m czone year using "$data/Lshares.dta", gen(merge_shares)
/*** BHJ SHOCKS **/
merge m:1 sic87dd year using "$data/shocks.dta", gen(merge_shocks)

rename ind_share share_emp_ind_bhj_
gen z_ = share_emp_ind_bhj_ * g
rename g g_
drop g_emp_ind-g_importsUSA
reshape wide share_emp_ind_bhj_ g z_, i(czone year) j(sic87dd)
egen z = rowtotal(z_*)

local controls reg_* l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource l_shind_manuf_cbp t2
local weight weight

local y d_sh_empl_mfg 
local x shock
local z z

local ind_stub share_emp_ind_bhj_
local growth_stub g_

local time_var year
local cluster_var czone

levelsof `time_var', local(years)

/** g_2141 and g_3761 = 0 for all years **/
drop g_2141 `ind_stub'2141
drop g_3761 `ind_stub'3761

forvalues t = 1990(10)2000 {
	foreach var of varlist `ind_stub'* {
		gen t`t'_`var' = (year == `t') * `var'
		}
	foreach var of varlist `growth_stub'* {
		gen t`t'_`var'b = `var' if year == `t'
		egen t`t'_`var' = max(t`t'_`var'b), by(czone)
		drop t`t'_`var'b
		}
	}

tab division, gen(reg_)
drop reg_1
tab year, gen(t)
drop t1

drop if czone == .

foreach var of varlist `ind_stub'* {
	if regexm("`var'", "`ind_stub'(.*)") {
		local ind = regexs(1) 
		}
	tempvar temp
	qui gen `temp' = `var' * `growth_stub'`ind'
	qui regress `x' `temp' `controls' [aweight=`weight'], cluster(czone)
	local pi_`ind' = _b[`temp']
	qui test `temp'
	local F_`ind' = r(F)
	qui regress `y' `temp' `controls' [aweight=`weight'], cluster(czone)
	local gamma_`ind' = _b[`temp']
	drop `temp'
	}

foreach var of varlist `ind_stub'3571 `ind_stub'3944 `ind_stub'3651 `ind_stub'3661 `ind_stub'3577 {
	if regexm("`var'", "`ind_stub'(.*)") {
		local ind = regexs(1) 
		}
	tempvar temp
	qui gen `temp' = `var' * `growth_stub'`ind'
	ch_weak, p(.05) beta_range(-10(.1)10)   y(`y') x(`x') z(`temp') weight(`weight') controls(`controls') cluster(czone)
	disp r(beta_min) ,  r(beta_max)
	local ci_min_`ind' =string( r(beta_min), "%9.2f")
	local ci_max_`ind' = string( r(beta_max), "%9.2f")
	disp "`ind', `beta_`ind'', `t_`ind'', [`ci_min_`ind'', `ci_max_`ind'']"
	drop `temp'
	}

preserve
keep `ind_stub'* czone year `weight'
reshape long `ind_stub', i(czone year) j(ind)
gen `ind_stub'pop = `ind_stub'*`weight'
collapse (sd) `ind_stub'sd = `ind_stub' (rawsum) `ind_stub'pop `weight' [aweight = `weight'], by(ind year)
tempfile tmp
save `tmp'
restore

bartik_weight, z(t*_`ind_stub'*)    weightstub(t*_`growth_stub'*) x(`x') y(`y') controls(`controls'  ) weight_var(`weight')

mat beta = r(beta)
mat alpha = r(alpha)
mat gamma = r(gam)
mat pi = r(pi)
mat G = r(G)
qui desc t*_`ind_stub'*, varlist
local varlist = r(varlist)

clear
svmat beta
svmat alpha
svmat gamma
svmat pi
svmat G

gen ind = ""
gen year = ""
local t = 1
foreach var in `varlist' {
	if regexm("`var'", "t(.*)_`ind_stub'(.*)") {
		qui replace year = regexs(1) if _n == `t'
		qui replace ind = regexs(2) if _n == `t'
		}
	local t = `t' + 1
	}

* new code * 
kdensity alpha1, normal addplot(histogram alpha1, below bin(50) color(lime)) title("Distributions of Rotemberg weights") graphregion (color(white)) legend(order(3 1 2)) 
graph export "$output/pset_2_exercise_2_question_2_a.pdf", replace

* (b) and (c) *
************** ADH code
/** Calculate Panel C: Variation across years in alpha **/
quietly {
total alpha1 if year == "1990"
mat b = e(b)
local sum_1990_alpha = string(b[1,1], "%9.3f")
total alpha1 if year == "2000"
mat b = e(b)
local sum_2000_alpha = string(b[1,1], "%9.3f")

sum alpha1 if year == "1990"
local mean_1990_alpha = string(r(mean), "%9.3f")
sum alpha1 if year == "2000"
local mean_2000_alpha = string(r(mean), "%9.3f")

destring ind, replace
destring year, replace
merge 1:1 ind year using `tmp'
gen beta2 = alpha1 * beta1
gen indshare2 = alpha1 * (`ind_stub'pop/`weight')
gen indshare_sd2 = alpha1 * `ind_stub'sd
gen G2 = alpha1 * G1
collapse (sum) alpha1 beta2 indshare2 indshare_sd2 G2 (mean) G1 , by(ind)
gen agg_beta = beta2 / alpha1
gen agg_indshare = indshare2 / alpha1
gen agg_indshare_sd = indshare_sd2 / alpha1
gen agg_g = G2 / alpha1
rename ind sic
merge 1:1 sic using "sic_code_desc"
rename sic ind
keep if _merge == 3
gen ind_name = subinstr(description, "Not Elsewhere Classified", "NEC", .)
replace ind_name = subinstr(ind_name, ", Except Dolls and Bicycles", "", .)

gsort -alpha1

/** Panel A: Negative and Positive Weights **/
total alpha1 if alpha1 > 0
mat b = e(b)
local sum_pos_alpha = string(b[1,1], "%9.3f")
total alpha1 if alpha1 < 0
mat b = e(b)
local sum_neg_alpha = string(b[1,1], "%9.3f")

sum alpha1 if alpha1 > 0
local mean_pos_alpha = string(r(mean), "%9.3f")
sum alpha1 if alpha1 < 0
local mean_neg_alpha = string(r(mean), "%9.3f")

local share_pos_alpha = string(abs(`sum_pos_alpha')/(abs(`sum_pos_alpha') + abs(`sum_neg_alpha')), "%9.3f")
local share_neg_alpha = string(abs(`sum_neg_alpha')/(abs(`sum_pos_alpha') + abs(`sum_neg_alpha')), "%9.3f")

/** Panel B: Correlations of Industry Aggregates **/
gen F = .
gen agg_pi = .
gen agg_gamma = .
levelsof ind, local(industries)
foreach ind in `industries' {
	capture replace F = `F_`ind'' if ind == `ind'
	capture replace agg_pi = `pi_`ind'' if ind == `ind'
	capture replace agg_gamma = `gamma_`ind'' if ind == `ind'		
	}
corr alpha1 agg_g agg_beta F agg_indshare_sd
mat corr = r(C)
forvalues i =1/5 {
	forvalues j = `i'/5 {
		local c_`i'_`j' = string(corr[`i',`j'], "%9.3f")
		}
	}

/** Panel  D: Top 5 Rotemberg Weight Inudstries ** (OUR PANEL B)**/
foreach ind in 3571 3944 3651 3661 3577 {
	qui sum alpha1 if ind == `ind'
   local alpha_`ind' = string(r(mean), "%9.3f")
	qui sum agg_g if ind == `ind'	
	local g_`ind' = string(r(mean), "%9.3f")
	qui sum agg_beta if ind == `ind'	
	local beta_`ind' = string(r(mean), "%9.3f")
	qui sum agg_indshare if ind == `ind'	
	local share_`ind' = string(r(mean)*100, "%9.3f")
	tempvar temp
	qui gen `temp' = ind == `ind'
	gsort -`temp'
	local ind_name_`ind' = ind_name[1]
	drop `temp'
	}


/** Over ID Figures **/
gen omega = alpha1*agg_beta
total omega
mat b = e(b)
local b = b[1,1]

gen label_var = ind 
gen beta_lab = string(agg_beta, "%9.3f")

gen abs_alpha = abs(alpha1) 
gen positive_weight = alpha1 > 0
gen agg_beta_pos = agg_beta if positive_weight == 1
gen agg_beta_neg = agg_beta if positive_weight == 0
twoway (scatter agg_beta_pos agg_beta_neg F if F >= 5 [aweight=abs_alpha ], msymbol(Oh Dh) ), legend(label(1 "Positive Weights") label( 2 "Negative Weights")) yline(`b', lcolor(black) lpattern(dash)) xtitle("First stage F-statistic")  ytitle("{&beta}{subscript:k} estimate")
graph export "$output/ps2_ex2_q2c_figureA2.pdf", replace

gsort -alpha1
twoway (scatter F alpha1 if _n <= 5, mcolor(dblue) mlabel(ind_name  ) msize(0.5) mlabsize(2) ) (scatter F alpha1 if _n > 5, mcolor(dblue) msize(0.5) ), name(a, replace) xtitle("Rotemberg Weight") ytitle("First stage F-statistic") yline(10, lcolor(black) lpattern(dash)) legend(off)
graph export "$output/ps2_ex2_q2c_figureA3.pdf", replace


/** Panel E: Weighted Betas by alpha weights ** (OUR PANEL A)**/
preserve
gen agg_beta_weight = agg_beta * alpha1

collapse (sum) agg_beta_weight alpha1 (mean)  agg_beta, by(positive_weight)
egen total_agg_beta = total(agg_beta_weight)
gen share = agg_beta_weight / total_agg_beta
gsort -positive_weight
local agg_beta_pos = string(agg_beta_weight[1], "%9.3f")
local agg_beta_neg = string(agg_beta_weight[2], "%9.3f")
local agg_beta_pos2 = string(agg_beta[1], "%9.3f")
local agg_beta_neg2 = string(agg_beta[2], "%9.3f")
local agg_beta_pos_share = string(share[1], "%9.3f")
local agg_beta_neg_share = string(share[2], "%9.3f")
restore
}

* new code
gen agg_beta_weight = agg_beta * alpha1
egen total_agg_beta = total(agg_beta_weight)
gen share = agg_beta_weight/total_agg_beta

foreach ind in 3571 3944 3651 3661 3577 {
	qui sum alpha1 if ind == `ind'
	local alpha_`ind' = string(r(mean), "%9.3f")
	qui sum agg_g if ind == `ind'
	local g_`ind' = string(r(mean), "%9.3f")
	qui sum agg_beta if ind == `ind'
	local beta_`ind' = string(r(mean), "%9.3f")
	qui sum agg_indshare if ind == `ind'
	local share_`ind' = string(r(mean)*100, "%9.3f")
	qui sum share if ind == `ind'
	local sharebeta_`ind' = string(r(mean), "%9.3f")
	tempvar temp
	qui gen `temp' = ind == `ind'
	gsort -`temp'
	local ind_name_`ind' = ind_name[1]
	drop `temp'
	}


/*** Write final table **/
capture file close fh
file open fh  using "$output/rotemberg_summary_adh.tex", write replace
file write fh "\begin{landscape}"
file write fh "\begin{table}"
file write fh "\begin{tabular}{l cccccc}"
file write fh "\toprule" _n

/** Panel A**/
file write fh "\multicolumn{7}{l}{\textbf{Panel A: Estimates of $\beta_{k}$ for positive and negative weights} }\\" _n
file write fh  " & $\alpha$-weighted Sum & Share of overall $\beta$ & Mean & & & \\ \cmidrule(lr){2-3}" _n
file write fh  " Negative & `agg_beta_neg' & `agg_beta_neg_share' &`agg_beta_neg2' & & &\\" _n
file write fh  " Positive & `agg_beta_pos' & `agg_beta_pos_share' & `agg_beta_pos2' & & &\\" _n

/** Panel B **/
file write fh "\multicolumn{7}{l}{\textbf{Panel B: Top 5 Rotemberg weight industries} }\\" _n
file write fh  " & $\hat{\alpha}_{k}$ & \$g_{k}$ & $\hat{\beta}_{k}$ & 95 \% CI & Ind Share & Share of overall $\beta$ \\ \cmidrule(lr){2-6}" _n
foreach ind in 3571 3944 3651 3661 3577 {
	if `ci_min_`ind'' != -10 & `ci_max_`ind'' != 10 {
		file write fh  "`ind_name_`ind'' & `alpha_`ind'' & `g_`ind'' & `beta_`ind'' & (`ci_min_`ind'',`ci_max_`ind'')  & `share_`ind'' & `sharebeta_`ind'' \\ " _n
		}
	else  {
		file write fh  "`ind_name_`ind'' & `alpha_`ind'' & `g_`ind'' & `beta_`ind'' & \multicolumn{1}{c}{N/A}  & `share_`ind'' & `sharebeta_`ind'' \\ " _n
		}
	}
	file write fh  "\bottomrule" _n
	file write fh "\end{tabular}"
	file write fh "\end{table}"
	file write fh "\end{landscape}"
file close fh

/* 
Under homogeneity, we can use the beta_k associated to each sector as individual instruments in the regression instead of the Bartik instrument. In fact, under such assumption, all the beta_k should probabilistically converge to a single beta, which is also the one related to the Bartik instrument. If we allow for heterogeneous treatment effects, then each beta_k will converge in probability to a convex combination of the location specific treatment effects. The Bartik IV is a weighted average of such coefficients, so it converges to another convex combination of the treatment effects, provided that the Rotemberg weights are all non negative.

The dispersion observed in figure A2 can be interpreted as either heterogeneity or endogeneity. In the case of endogenous shares, aymptotic bias in the sector specific coefficients would generate heterogeity. If we accept the argument presented in the paper that heterogeneity drives the dispersion because the model is well-specified (rouling out ednogenenity), we identify two types of treatment effect heterogeneity. The first one is in terms of magnitude, relating to deviations from the Bartik estimator. The second arising, from the sign of the coefficients, relating to the presence of both positive and negative Rotenberg weights. 

In this case, we can deduce heterogeneity in the location-specific parameters. In fact, the beta_k, which converge under heterogeneity to different convex combinations of the beta_l, exhibit significant dispersion, both in sign and magnitude. Indeed, the greater the disparity among the beta_l, the more diverse the convex combinations will be, resulting in dispersed beta_k. Additionally, negative weights to the estimates are also present, indicating heterogeneity at the sign level. In this instance, we must check that 2 conditions hold in order to avoid inconsistency in the Bartik estimator. First, we need to have minimal vertical dispersion. Secondly, negative weights should not exceed a certain threshold. In fact, if the dispersion is not excessive, the convex weights will eventually become identical in the limit. Consequently, the signs associated to each location parameter beta_l will be uniquely determined by the sum of the location specific Rotenberg weights, which can be still be positive if the negative weights are not too frequent and large in magnitude.

In our scenario, there are some negative Rotenberg weights associated with certain point estimates. This raises the possibility (but does not necessarily imply) of nonconvex weights on the beta_l. This would suggest that the combination (Bartik estimate) is not convex, and the LATE-like interpretation is invalidated. However, this is not necessarily true, as the negative weights might be small and offset by positive ones in the sum, maintaining a convex combination and preserving the LATE interpretation. This is likely the case here, as the negative weights are outliers and carry minimal weight.

However, it is worth noticing that the estimates with the largest weights are concentrated around a weak F-statistic. This can be understood though figure A3. In particular, we are in a situation in which the Bartik estimate overly depends on some specific industries, equivalent to the ones associated with high Rotemberg weights. Such overweighted industries should be coherent with their relationship with trade shocks. For instance, it would be suspicious to observe that one sector unrelated to trade is largely contributing to our Bartik estimate. Instead, we want that the sectors with the largest shocks are the ones driving the estimates. Such instruments are electronic computers, games and toys, household audio and video, telephone apparatus and computer equipment. As suggested in Pinkham Appendix (2020), these sectors are higher-skill technologically innovative industries, which could be problematic if changes in technology are the main driver of industry shocks. Note that this is different from what Autor et al. (2013) intended to achieve, since the optimal setting would have been characterised by higher weights associated with relatively low-skill technologically stagnant industries, where chinese trade shocks are the most impactful. 

Finally, note that such limitation could be linked to the weakness of the instruments highlighted both in figure A2 and A3.  This could pose a limitation to the interpretation of the Bartik coefficient (both under homogeneity and heterogeneity of treatment effect). In fact, even if the threats associated to negative weights are not too pronounced, the coefficient is largely determined by shocks that are not related to instrumented variable. 
*/

*********************************  QUESTION 3  *********************************
/*  Considering IV estimators, under the assumption of homogeneous treatment effects, the identification assumptions that must hold for the estimates to be consistent are:
(1) Relevance, the instrument should capture a part of the variation in the endogenous variable. (For the interpetation, refer to question 2.1a and 2.1b)
(2) Exogeneity which consists of:
2.1 Random assignment, the instrument should be as good as randomly allocated. (For the interpetation, refer to question 2.1a and 2.1b)
2.2 Exclusion restriction, the instrument should affect the outcome only through the effect on the exogenous variable. In the case presented, this amounts to saying that the effect of Chinese imports exposure of the "control" countries on vote shares should only arise through the effect on import competition in the US.

However, the context of the IV estimates presented in Section V of Autor et al. (2020) indicates the existence of heterogeneous treatment effects. Thus, this setting requires a further assumption:
(3) Monotonicity: Given the exposure to the instrument, that is, the exposure to Chinese imports of high income countries other than the US, there can not be a treatment effect that goes in the opposite direction as the one that you would have without the exposure. This amounts to maintaing a "no defiers" assumption. In the case presented, the growth of Chinese import penetration in the US for an industry k in a period t cannot a different effect in sign compared to same shocks at same time t for other industries.

The main difference between the models presented in Autor et al. (2013) and Autor et al. (2020) is limited to the outcome variable (the first being labour market outcomes the second electoral outcomes). The relevance, randomness and monotonicity assumptions present no difference as the instrument and instrumented variables are the same. Thus, the different outcome variable between the two papers, impacts the assumption relating to the exogeneity of the instrument. As for the issues arising in both papers, as explained in Autor et al. (2013), several possible threats to the identification strategy presented in the paper exist:
(a) Product demand shock correlated across high income countries.
(b) Productivity shocks in the US that could be driving growth in imports from China. 
(c) Technological shocks in China may cause an increase in demand for Chinese imports by high income countries in their labor-intensive industries of US and other Western countries directly impacting the employment of workers in these industries and indirectly their voting preferences. This last one being relevant only for Autor et al. (2020).
While (a) and (b) still constitute threats for both models, the relative impact of (c) cannot be determined a priori in Autor et al. (2020). Because we are studying a different outcome, that is mainly affected indirectly through labour market outcomes, one argument is that the effect of layoffs in labor-intensive industries may be less relevant for the election outcome estimation. On the other hand, one could argue that this mechanism is related to, or majorly affects political preferences.
 
Furthermore, additional threats originate from the particular outcome variable present in the Autor et al. (2020). The US election outcome may be correlated  with high income countries exposure to Chinese imports. This could happen through two opposing channels. On the one hand, US elections might have an effect on worldwide trade policy stance and indeed cause shifts to high income countries imports. On the other hand, high income countries exposure to Chinese imports might influence the US elections as observing a dire economic situation in other countries could provoke fear of replication in the US. If one or both of those scenarios were to happen, the exclusion restriction of the model would not be satisfied.
 */














