
/* Estimate ps's using logistic regression. */
proc logistic data = data.population noprint;
  class agegroup(ref = "0-18") / param = ref;
  model treatment(ref = "0") = male risk_score agegroup / risklimits;
  output out = ps(drop = _level_) prob = ps;
run;

/* Make ps-matched population. Matching is done on logit(ps) using a 
caliper of 0.2 times the standard deviation of logit(ps). Matching is
done with replacement. */
%ps_match(
  in_ds = ps, 
  out_pf = ps, 
  group_var = treatment, 
  ps_var = ps, 
  seed = 1
);

/* Calculate ATE and ATT weights. */
data ps_weights;
  set ps;
  weight_ate = treatment / ps + (1 - treatment) / (1 - ps);
  weight_att = treatment + (1 - treatment) * ps / (1 - ps);
run;

/* Restructure and combine the original population with the matched and 
weighted populations. */
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
    pop = "ATE-weighted";
    w = weight_ate;
    output;
    pop = "ATT-weighted";
    w = weight_att;
    output;
  end;
  drop weight_ate weight_att;
run;

proc sort data = dat1 out = data.analysis_dat;
  by pop id;
run;
  


proc datasets library = work kill nolist;
quit;
