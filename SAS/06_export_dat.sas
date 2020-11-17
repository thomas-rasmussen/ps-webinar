/* Make .csv versions of the data to be version-controlled and used in the slides. */

%macro export_csv(dat);
  %local i i_dat;
  %do i = 1 %to %sysfunc(countw(&dat, %str( )));
    %let i_dat = %scan(&dat, &i, %str( ));
    proc export data = data.&i_dat outfile = "&path_data\&i_dat..csv" 
      dbms = csv replace;
    run;
  %end;
%mend export_csv;


%export_csv(
  absolute_effect_estimates
  analysis_dat
  assess_ph
  cumulative_incidences
  empirical_cdf
  induced_absolute_effects
  induced_relative_effects
  population
  relative_effect_estimates
  standardized_differences
  summary_tbl
  weight_stats
  bootstrap_est
);
