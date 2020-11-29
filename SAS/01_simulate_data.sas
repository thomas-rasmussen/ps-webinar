/* Simulation of data. It is way out of scope of this presentation to try
and explain how this works. Just run the syntax if you need to generate 
the data used in the example. */

/* Specify the number of observations used to estimate parameters and specify
the convergence criteria. Using a dataset with 1 million observations and a 
convergence griteria of 0.00001 is sufficient to estimate parameters inducing 
the necessary precision. */
%let EstParm_NObs = 10**6; 
%let EstParm_ConvCri = 0.00001; 

/*******************************************************************************
ESTIMATE PARAMETERS
*******************************************************************************/

/*** Simulate covariate data ***/

/* Number of observations in simulated data. */
%let NObs = 10000; 
/* Number of covariates to simulate from multivariate normal distribution. */
%let NVars = 5;	
/* Number of simulated covariates that are dichotomized. */
%let ND = 4;

%SimVarsMultiNormal(
	OutDS = par_dat1,
	NObs = &EstParm_NObs,
	NSim = 1,
	NVars = &NVars,
	ND = &ND,
	VarCorr = 0,
	BlockSize = 1,
	seed = 1,
	DEL = Y
);

/*** Estimate intercept parameter in PS-model ***/

/* Transform variables */
data par_dat2;
  set par_dat1;
  age = max(0, 50 + sqrt(100) * x5);
  age_trans = age**1.1;
  risk_score = x1;
  comor_1 = x2;
  comor_2 = x3;
  male = x4;
  drop x:;
run;

/* Induced proportion of treated patients. */
%let E_TargetProp = 0.5;
/* Linear predictor minus intercept, used in true PS model */
%let E_LP_NoInt = 
  log(1.3) * comor_1
  + log(1.3) * comor_2
  + log(1.5) * male  
  + log(1.02) * age_trans
  + log(1.02) * male * age_trans
  ;


%SimBinVar(
	InDS = par_dat2,
  OutParms = par1,
	TargetProp = &E_TargetProp,
	LP_NoInt = &E_LP_NoInt,
	ConvCri = &EstParm_ConvCri,
	MaxIte = 30,
	seed = 2,
	DEL = Y
);

/*** Induce desired marginal exposure effect ***/

/* Parameters specifying the Weibull distributed baseline hazard function. */
%let lambda = 0.00002;
%let eta = 1.5;
/* Marginal HR exposure effect that is to be induced. */
%let MarEffect = 0.5; 
/* Linear predictor minus exposure term, used in the true outcome Cox model. */
%let Cox_LP_NoE = 
  log(1.3) * risk_score 
  + log(1.3) * comor_1
  + log(1.3) * comor_2
  + log(1.5) * male  
  + log(1.02) * age_trans
  + log(1.02) * male * age_trans
  ;


ods listing;
ods results off;
%SimTimeToEvent(
	InDS = par_dat2,
	OutParms = par2,
	TargetMar = &MarEffect,
	LP_NoInd = &Cox_LP_NoE,
	ConvCri = &EstParm_ConvCri,
  MaxIte = 30,
	lambda = &lambda,
	eta = &eta,
	seed = 3,
	DEL = Y
);

/* Combine parameter datasets */
data par3;
	merge par1 par2;
run;

/*******************************************************************************
SIMULATE DATA
*******************************************************************************/

proc sql noprint;
	select Alpha0Est, Beta0Est 
		into :Alpha0Est, :Beta0Est
		from par3;
quit;

/* Simulate covariate data */
%SimVarsMultiNormal(
	OutDS = cov_dat1,
	NObs = &NObs,
	NSim = 1,
	NVars = &NVars,
	ND = &ND,
	VarCorr = 0,
	BlockSize = 1,
	seed = 4,
	DEL = Y
	);

/* Simulate exposure, outcome and censoring. */
data data.population;
  retain id treatment male age risk_score comorbidity_score time death;
	set cov_dat1;
  format comorbidity_score $3.;
	call streaminit(5);
  /* transform variables */
  age = max(0, 50 + sqrt(100) * x5);
  age_trans = age**1.1;
  risk_score = x1;
  comor_1 = x2;
  comor_2 = x3;
  if comor_1 + comor_2 = 0 then comorbidity_score = "0";
  if comor_1 + comor_2 = 1 then comorbidity_score = "1-2";
  if comor_1 + comor_2 = 2 then comorbidity_score = "3+";
  male = x4;
	/* Calculation of true PS */
	ps_true = 1 / (1 + exp(-&Alpha0Est - (&E_LP_NoInt)));
	/* Simulation of exposure based on the true ps */
	treatment = rand("bernoulli", ps_true);
	/* Simulation of time-to-event */
	u = rand("uniform");
	time = (-log(u) / (0.00002 * exp(&Beta0Est * treatment + &Cox_LP_NoE))) ** (1 / 1.5);
  death = 1;
  time = round(time);
  /* Censor all time-to-events above 365. */
  if time > 365 then do;
    time = 365;
    death = 0;
  end;
  id = _n_;
  keep id treatment male risk_score age comorbidity_score time death;
run;

/*******************************************************************************
ASSESS INDUCED EFFECTS
*******************************************************************************/

%SimVarsMultiNormal(
	OutDS = test_dat1,
	NObs = 10**6,
	NSim = 1,
	NVars = &NVars,
	ND = &ND,
	VarCorr = 0,
	BlockSize = 1,
	seed = 6,
	DEL = Y
	);

data test_dat2;
  retain treatment male age risk_score comorbidity_score time death;
	set test_dat1;
	call streaminit(7);
  /* transform variables */
  age = max(0, 50 + sqrt(100) * x5);
  age_trans = age**1.1;
  risk_score = x1;
   comor_1 = x2;
  comor_2 = x3;
  if comor_1 + comor_2 = 0 then comorbidity_score = "0";
  if comor_1 + comor_2 = 1 then comorbidity_score = "1-2";
  if comor_1 + comor_2 = 2 then comorbidity_score = "3+";
  male = x4;
	/* Calculation of true PS */
	ps_true = 1 / (1 + exp(-&Alpha0Est - (&E_LP_NoInt)));
	/* Simulation of exposure based on the true ps */
	treatment = rand("bernoulli", ps_true);
	/* Simulation of time-to-event */
	u = rand("uniform");
	time = (-log(u) / (0.00002 * exp(&Beta0Est * treatment + &Cox_LP_NoE))) ** (1 / 1.5);
  death = 1;
  time = round(time);
  /* Censor all time-to-events above 365. */
  if time > 365 then do;
    time = 365;
    death = 0;
  end;
  weight_ate = treatment / ps_true + (1 - treatment) / (1 - ps_true);
  weight_att = treatment + (1 - treatment) * ps_true / (1 - ps_true);
run;

/* Check induced conditional relative effects. */
proc phreg data = test_dat2 plots = none;
  class comorbidity_score(ref = "0") / param = ref;
  model time*death(0) = treatment risk_score comorbidity_score male|age_trans / risklimits;
  ods output ParameterEstimates = est_rel_cond;
run;

/* Check induced crude relative treatment effect. */
proc phreg data = test_dat2 plots = none;
  model time*death(0) = treatment / risklimits;
  ods output ParameterEstimates = est_rel_crude;
run;

/* Check induced relative ATE and ATT */
proc phreg data = test_dat2 plots = none;
  model time * death(0) = treatment / risklimits;
  weight weight_ate;
  ods output ParameterEstimates = est_rel_ate;
run;

proc phreg data = test_dat2 plots = none;
  model time*death(0) = treatment / risklimits;
  weight weight_att;
  ods output ParameterEstimates = est_rel_att;
run;

data data.induced_relative_effects;
  format check E_LP_NoInt Cox_LP_NoE parameter $200.;
  set
    est_rel_cond(in = q1)
    est_rel_crude(in = q2)
    est_rel_ate(in = q3)
    est_rel_att(in = q4)
  ;
  if q1 then check = "Induced relative conditional effects";
  if q2 then check = "Induced relative crude treatment effect";
  if q3 then check = "Induced relative ATE";
  if q4 then check = "Induced relative ATT";

  MarEffect = &MarEffect;
  E_LP_NoInt = "&E_LP_NoInt"; 
  Cox_LP_NoE = "&Cox_LP_NoE";

  keep check MarEffect E_LP_NoInt Cox_LP_NoE 
       parameter classval0 hazardratio hrlowercl hruppercl;
run;

/* Check induced risk 1-year risk differences */
proc lifetest data = test_dat2
    outsurv = risk_diff_crude1(where = (time = 365 and survival ne .))
    notable noprint plots = none method = km;
  strata treatment;
  time time * death(0);
run;

proc lifetest data = test_dat2
    outsurv = risk_diff_att1(where = (time = 365 and survival ne .))
    notable noprint plots = none method = km;
  strata treatment;
  time time * death(0);
  weight weight_att;
run;

proc lifetest data = test_dat2
    outsurv = risk_diff_ate1(where = (time = 365 and survival ne .))
    notable noprint plots = none method = km;
  strata treatment;
  time time * death(0);
  weight weight_ate;
run;

data risk_diff2;
  set risk_diff_crude1(in = q1)
      risk_diff_ate1(in = q2)
      risk_diff_att1(in = q3);
  if q1 then effect = "crude";
  if q2 then effect = "ate";
  if q3 then effect = "att"; 
run;

proc sort data = risk_diff2;
  by effect descending treatment;
run;

data data.induced_absolute_effects;
  retain effect risk_diff;
  set risk_diff2;
  by effect;
  if first.effect then risk_diff = survival;
  else risk_diff = survival - risk_diff;
  if last.effect;
  keep effect risk_diff time;
run;



proc datasets library = work kill nolist;
quit;


