
/* Path to webinar folder */
%let path_main = ;
/*%let path_main = S:\Thomas Rasmussen\github_repos\ps-webinar;*/

/* Paths to subfolders */
%let path_syntax = &path_main\SAS;
%let path_macros = &path_main\SAS-macros;
%let path_output = &path_main\output;
%let path_data   = &path_main\data;

libname data   "&path_data";
libname output "&path_output";

/* Load macros */
%include "&path_macros\calculate_sd.sas";
%include "&path_macros\empirical_cdf.sas";
%include "&path_macros\ps_match.sas";
%include "&path_macros\pt_char.sas";
%include "&path_macros\SimBinVar.sas";
%include "&path_macros\SimTimeToEvent.sas";
%include "&path_macros\SimVarsMultiNormal.sas";

/* Simulate data */
%include "&path_syntax\01_simulate_data.sas";

/* Estimate ps and make matched and weighted populations */
%include "&path_syntax\02_make_analysis_dat.sas";

/* Assess balance */
%include "&path_syntax\03_assess_balance.sas";

/* Assess proportional hazards assumption */
%include "path_syntax\04_assess_ph.sas";

/* Estimate treatment effect */
%include "path_syntax\05_estimate_treatment_effect.sas";

/* Export data */
%include "path_syntax\06_export_dat.sas";
