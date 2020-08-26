/*******************************************************************************
Assess proportional hazards assumption
*******************************************************************************/

/* Estimate survival curves for treated and untreated patients. */
proc lifetest data = data.analysis_dat outsurv = ph1
     notable noprint method = km;
  strata pop treatment;
  time time * death(0);
  weight w;
run;

/* Calculate log(time) and log(-log(survival)) */
data ph2;
  set ph1;
  by pop treatment;
  retain surv;
  if first.treatment then surv = 1;
  if survival ne . then surv = survival;
  if time > 0 then log_time = log(time);
  else log_time = .;
  if surv < 1 then lml_surv = log(-log(surv));
  else lml_surv = .;
  keep pop treatment log_time lml_surv;
run;

proc sort data =ph2 nodupkeys;
  by _all_;
run;

/* Save data for slides */
data data.assess_ph;
  set ph2;
run;

proc format;
  value treatment
    0 = "Untreated"
    1 = "Treated"
  ;
run;

ods pdf file = "&path_output\assess_ph.pdf" pdftoc = 1 style = statistical;
  proc sgpanel data = ph2;
    panelby pop / onepanel;
    step x = log_time y = lml_surv / group = treatment;
    colaxis label = "log(time)";
    rowaxis label = "log(-log(survival))";
    format treatment treatment.;
    keylegend / title = "";
  run;
ods pdf close;


/*******************************************************************************
Bootstrap
*******************************************************************************/

/* Since we are using bootstrap sampling to estimate confidence intervals, we
start by making the bootstrap samples and add them to the original population. 
Then we estimate of ps's, make matched and weighted populations etc. like 
before, now just done in each of the bootstrap samples.*/

/* Make bootstrap samples */
proc surveyselect data = data.population out = bs1(drop = numberhits) 
  method = urs seed = 1 noprint reps = 200 samprate = 1 outhits;
run;

/* Add to original data */
data bs2;
  set data.population(in = q1) bs1;
  if q1 then do;
    replicate = 0;
  end;
run;

/* Estimate ps's */
proc logistic data = bs2 noprint;
  by replicate;
  class agegroup(ref = "0-18") / param = ref;
  model treatment(ref = "0") = male risk_score agegroup;
  output out = ps(drop = _level_) prob = ps;
run;

/* Make ps-matched populations. Takes 5-10 mins on SAS server. */
%ps_match(
  in_ds = ps, 
  out_pf = ps, 
  group_var = treatment, 
  ps_var = ps,
  by = replicate,
  seed = 1
);

data data.bootstrap_matched_dat;
  set ps_matches;
run;

data ps_matches;
  set data.bootstrap_matched_dat;
run;

/* Calculate weights */
data ps_weights;
  set ps;
  weight_ate = treatment / ps + (1 - treatment) / (1 - ps);
  weight_att = treatment + (1 - treatment) * ps / (1 - ps);
run;

/* Restructure and combine data */
data dat1;
  format pop $20.;
  set ps(in = q1) ps_matches(in = q2) ps_weights(in = q3);
  if q1 then do;
    pop = "Original";
    w   = 1;
    output;
  end;
  else if q2 then do;
    pop = "Matched";
    w = 1;
    output;
  end;
  else if q3 then do;
    pop = "ATT-weighted";
    w = weight_att;
    output;
    pop = "ATE-weighted";
    w = weight_ate;
    output;
  end;
  drop weight_ate weight_att __match;
run;

proc sort data = dat1;
  by pop replicate treatment;
run;


/*******************************************************************************
Relative conditional effect
*******************************************************************************/

/* Estimate conditional effect of treatment on 1 year mortality. */
proc phreg data = dat1(where = (pop = "Original" and replicate = 0)) 
    plots = none;
  by pop;
  class agegroup(ref = "0-18") / param = ref;
  model time * death(0) = treatment male risk_score agegroup / risklimits;
  ods output ParameterEstimates = rel_cond1(where = (parameter = "treatment"));
run;

data rel_cond2;
  set rel_cond1;
  where parameter = "treatment";
  keep pop hazardratio hrlowercl hruppercl;
run;


/*******************************************************************************
Relative marginal effects
*******************************************************************************/

proc phreg data = dat1 plots = none;
  by pop replicate;
  model time * death(0) = treatment / risklimits;
  weight w;
  ods output ParameterEstimates = rel_mar1;
run;

/* Use 2.5 and 97.5 percentile of bootstrap estimates as CI limits */
proc univariate data = rel_mar1(where = (replicate ne 0));
  by pop;
  var hazardratio;
  output out = rel_mar_ci1 pctlpts = 2.5 97.5 pctlpre = p;
run;

data rel_mar2;
  merge rel_mar1(where = (replicate = 0)) rel_mar_ci1(in = q2);
  by pop;
  if q2 then do;
    hrlowercl = p2_5;
    hruppercl = p97_5;
  end;
  keep pop hazardratio hrlowercl hruppercl;
run;

/* Combine effect estimates and save for slides. */
data all_rel_est;
  format effect $50.;
  set rel_cond2(in = q1) rel_mar2(in = q2);
  if q1 then effect = "Conditional";
  else if q2 then effect = "Marginal";
run;

proc sort data = all_rel_est out = data.relative_effect_estimates;
  by pop effect;
run;


/*******************************************************************************
Cumulative incidence curves
*******************************************************************************/

/* Estimate cumulative incidence curves for treated and untreated patients. */
proc lifetest data = dat1 outsurv = cum_inc1
    notable noprint plots = none method = km;
  by pop replicate treatment;
  time time * death(0);
  weight w;
run;

/* Expand the data so we have estimates on all days. */
data cum_inc2;
  set cum_inc1;
  by pop replicate treatment;
  retain cum_inc;
  if first.treatment then cum_inc = .;
  if survival ne . then cum_inc = 1 - survival;
  keep pop replicate treatment time cum_inc;
run;

data cum_inc3;
  set cum_inc2(rename = (time = tmp));
  by pop replicate treatment;
  time_last = lag(tmp);
  if first.treatment then time_last = 0;
  time = tmp;
  do i = time_last + 1 to tmp;
    time = i;
    output;
  end;
  drop tmp time_last i;
run;

/* Estimate CI's from bootstrap samples */
proc sort data = cum_inc3(where = (replicate ne 0)) out = cum_inc_ci1;
  by pop treatment time;
run;

proc univariate data = cum_inc_ci1;
  by pop treatment time;
  var cum_inc;
  output out = cum_inc_ci2 pctlpts = 2.5 97.5 pctlpre = p;
run;

/* Combine and save data for slides. */
data data.cumulative_incidences;
  merge 
    cum_inc3(where = (replicate = 0)) 
    cum_inc_ci2(rename = (p2_5 = lcl p97_5 = ucl))
  ;
  by pop treatment time;
  drop replicate;
run;


/*******************************************************************************
Estimate absolute ATE and ATT
*******************************************************************************/

/* Estimate risk difference */
proc sort data = cum_inc3 out = cum_inc_diff1;
  by pop replicate time descending treatment;
run;

data cum_inc_diff2;
  set cum_inc_diff1;
  by pop replicate time;
  retain cum_inc_diff;
  if first.time then cum_inc_diff = cum_inc;
  else cum_inc_diff = cum_inc_diff - cum_Inc;
  if last.time;
  drop cum_inc treatment;
run;

/* Estimate CI's */  
proc sort data = cum_inc_diff2(where = (replicate ne 0)) out = cum_inc_diff_ci1;
  by pop time;
run;

proc univariate data = cum_inc_diff_ci1;
  by pop time;
  var cum_inc_diff;
  output out = cum_inc_diff_ci2 pctlpts = 2.5 97.5 pctlpre = p;
run;

/* Combine and save data for slides */
data data.absolute_effect_estimates;
  merge 
    cum_inc_diff2(where = (replicate = 0)) 
    cum_inc_diff_ci2(rename = (p2_5 = lcl p97_5 = ucl))
  ;
  by pop time;
  drop replicate;
run;

/* 1 year effects estimates */
data abs_eff1;
  set data.absolute_effect_estimates;
  where time = 365;
run;


proc datasets library = work kill nolist;
quit;
