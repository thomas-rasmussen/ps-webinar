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

/* Make plots to grahically assess PH assumption */
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


proc datasets library = work kill nolist;
quit;
