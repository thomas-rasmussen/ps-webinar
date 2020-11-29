/*******************************************************************************
Descriptive summary table
*******************************************************************************/

/* Make a descriptive summary table */
%pt_char(
  in_ds     = data.analysis_dat,
  out_ds    = sum_tbl1,
  var_list  = male risk_score comorbidity_score age,
  var_types = d cat cat cont,
  by        = pop,
  strata    = treatment,
  weight    = w
);

/* Remove overall treatment group made by the macro and sort table data. 
Note that we are including deaths in this example because we need it for
the slides, but we should not include the outcome in a real analysis, where 
we should blind ourselves to the outcomes, so as not to let them influence 
what analysis we choose. */
data sum_tbl2;
  set sum_tbl1;
  where treatment ne .;
  if pop = "Original" then sort = 1;
  else if pop = "Matched" then sort = 2;
  else if pop = "ATT-weighted" then sort = 3;
  else if pop = "ATE-weighted" then sort = 4;
  sort2 = _n_;
run;

proc sort data = sum_tbl2;
  by sort sort2;
run;

/* Save table data for use in slides. */
data data.summary_tbl;
  set sum_tbl2;
    keep pop treatment __label __stat_char;
run;

/* Note that we have two apparantly identical formatted values in $labelk, which is
not possible. This was achieved by using the invisible character Alt-255 in the 
second "  0" instead of spaces. */
proc format;
  value treatment
    0 = "Untreated"
    1 = "Treated"
  ;

  value $ label
    "__n" = "^S = {font_weight = bold}Number of patients, N (%)"
    "male" = "^S = {font_weight = bold}Male, N (%)"
    "age" = "^S = {font_weight = bold}Age, median (Q1-Q3)"
    "risk_score: title" = "^S = {font_weight = bold}Risk score, N (%)"
    "risk_score: 0" = "  0"
    "risk_score: 1" = "  1"
    "comorbidity_score: title" = "^S = {font_weight = bold}Comorbidity score, N (%)"
    "comorbidity_score: 0"    = "  0"
    "comorbidity_score: 1-2"  = "  1-2"
    "comorbidity_score: 3+"   = "  3+"
  ;
run;

ods escapechar = "^";
ods pdf file = "&path_output\summary_table.pdf" notoc;
  proc report data = sum_tbl2
      style(report) = {font_size = 7pt}
      style(header) = {font_size = 7pt font_weight = bold}
      style(column) = {font_size = 7pt just = c}
    ;
    columns __label pop, treatment, (__stat_char __report_dummy);
    define __label / "" group order = data format = $label.
      style(column) = {just = l asis = on};
    define pop / "Population" across order = data;
    define treatment / "" across order = data format = treatment.;
    define __stat_char / "" display;
    define __report_dummy / noprint;
  run;
ods pdf close;


/*******************************************************************************
PS distribution
*******************************************************************************/

/* Restructure analysis data */
data ps_dist;
  set data.analysis_dat;
  if treatment = 0 then do;
    ps0 = ps;
    w0 = w;
  end;
  if treatment = 1 then do;
    ps1 = ps;
    w1 = w;
  end;
run;

proc sort data = ps_dist;
  by pop;
run;


/* Using the group option in sgplot/sgpanel does not work as 
we want it to. Work around by using sgplot and overlay separate plots for 
treated and untreated patients. */
ods pdf file="&path_output/ps_distribution.pdf"
	style = statistical startpage = yes pdftoc = 2;
  proc sgplot data = ps_dist noautolegend;
    by pop;
    density ps0 / weight = w0 type = kernel(c = 0.8) 
      name = "untreated" legendlabel = "Untreated";
    density ps1 / weight = w1 type = kernel(c = 0.8) 
      name = "treated" legendlabel = "Treated";
    keylegend "untreated" "treated"; 
    xaxis label = "ps";
  run;
ods pdf close;


/*******************************************************************************
Standardized differences
*******************************************************************************/

/* Calculate SD's in each analysis */
%calculate_sd(
  in_ds     = data.analysis_dat,
  out_ds    = sd_overall,
  group_var = treatment,
  var       = male risk_score age comorbidity_score,
  weight    = w,
  by        = pop
);

/* Repeat, but now in stratas of the male variable */
%calculate_sd(
  in_ds     = data.analysis_dat,
  out_ds    = sd_strata,
  group_var = treatment,
  var       = risk_score age comorbidity_score,
  weight    = w,
  by        = pop male
);

/* Combine and save data for slides */
data data.standardized_differences;
  format strata $10.;
  set sd_overall(in = q1) sd_strata(in = q2);
  if q1 then strata = "overall";
  if q2 then strata = "male";
run;

proc sort data = sd_strata;
  by male;
run;

ods pdf file="&path_output/standardized_differences.pdf" 
  style = statistical pdftoc = 2;
  ods proclabel = "Overall";
  proc sgplot data = sd_overall;
    scatter x = __sd y = __var / group = pop;
    refline 0.1 / axis = x; 
		yaxis type = discrete display = (nolabel);
		xaxis label = "Absolute standardized difference";
	run;
  ods proclabel = "In male strata";
  proc sgplot data = sd_strata;
    by male;
    scatter x = __sd y = __var / group = pop;
    refline 0.1 / axis = x; 
		yaxis type = discrete display = (nolabel);
		xaxis label = "Absolute standardized difference";
	run;
ods pdf close;


/*******************************************************************************
Weight statistics
*******************************************************************************/

ods pdf file = "&path_output/weight_stats.pdf" style = statistical notoc;
	proc means data = data.analysis_dat min p1 p99 max mean stddev nway;
    class pop treatment;
		var w;
    output out = data.weight_stats(drop = _type_ _freq_)
      n = n min = min p1 = p1 p99 = p99 max = max mean = mean stddev = stddev
      / noinherit;
	run;
ods pdf close;


/*******************************************************************************
Empirical CDF
*******************************************************************************/

/* Calculate empirical CDF of age for treated and untreated patients 
in each analysis. */
%empirical_cdf(
  in_ds   = data.analysis_dat,
  out_ds  = cdf_overall,
  var     = age,
  strata  = pop treatment, 
  weight  = w
);

/* Repeat, now in stratas of the variable male. */
%empirical_cdf(
  in_ds   = data.analysis_dat,
  out_ds  = cdf_strata,
  var     = age,
  strata  = pop male treatment, 
  weight  = w
);

/* Combine and save data for slides. */
data data.empirical_cdf;
  format strata $10.;
  set cdf_overall(in = q1) cdf_strata(in = q2);
  if q1 then strata = "overall";
  if q2 then strata = "male";
run;

ods pdf file="&path_output/empirical_cdf.pdf" style = statistical pdftoc = 1;
  ods proclabel = "Overall";
  proc sgpanel data = cdf_overall;
    panelby pop / onepanel;
    step x = __x y = __cdf / group = treatment;
    colaxis label = "Age";
    rowaxis label = "CDF";
    keylegend / title = "";
    format treatment treatment.;
  run;
 ods proclabel = "Stratify by agegroup";
  proc sgpanel data = cdf_strata;
    panelby pop male / onepanel layout = lattice;
    step x = __x y = __cdf / group = treatment;
    colaxis label = "Age";
    rowaxis label = "CDF";
    keylegend / title = "";
    format treatment treatment.;
  run;
ods pdf close;



proc datasets library = work kill nolist;
quit;
