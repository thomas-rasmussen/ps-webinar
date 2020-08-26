/* To reduce file-sizes and avoid binary files in git, we make .csv versions
of the data we use in the slides. */

%macro export_csv(data);
  proc export data = data.&data outfile = "&path_data\&data..csv" 
    dbms = csv replace;
  run;
%mend export_csv;


%export_csv(absolute_effect_estimates);
%export_csv(analysis_dat);
%export_csv(assess_ph);
%export_csv(cumulative_incidences)
%export_csv(empirical_cdf);
%export_csv(induced_absolute_effects);
%export_csv(induced_relative_effects);
%export_csv(population);
%export_csv(relative_effect_estimates);
%export_csv(standardized_differences);
%export_csv(summary_tbl);
%export_csv(weight_stats);
