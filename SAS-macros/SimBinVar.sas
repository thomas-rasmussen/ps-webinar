/*******************************************************************************
**	DESCRIPTION:
**	The macro takes an input dataset with variable information for a number of
**	observations and predicts a new binary variable using a logistic regression
**	model. Using an iterative bisection approach the intercept of the model 
**	is chosen to induce a speficied proportion of the observations having the
**	binary vairable equal to one. The estimated parameters is saved in an output 
**	dataset.
**
**	In the first step of the macro the logistic regression model
**		logit(p_i)=alpha0+alpha1*x1_i+...+alphak*xk_i, i=1,..,n
**	if fit. Here alpha_0 (=0) is the intercept, x1_i,...,xk_i is the k 
**	independent variables and alpha1,...,alphak is the specicified log(OR) 
**	associations. Using this model p_i is calculated for each observation,
**	and the binary variable is simulated from a Bernoulli distribution with 
**	parameter p_i. Finally we calculate the proportion of observations with 
**	the binary variable equal to one using this initial model.
**
**	In the second step the intercept alpha0 is iteratively increased(decreased) 
**	when the initial proportion is too low(high) until the induced proportion 
**	is too high(low). Thus an interval [MinAlpha0,MaxAlpha0] can been found, 
**	in which an alpha0-value exists that induces the desired proportion of 
**	observations with the simulated binary variable equal to one. 
**
**	In the third step an iterative bisection approach is used to determine the 
**	alpha0-value that induces the desired proportion. Using this approach the 
**	simulated binary variable can be simulated so that the proportion of
**	observations with the binary variable equal to one is arbitrarily close to 
**	the target value, or at least as close as the data allows, eg if there are 
**	only ten observations and the target proportion is 51%, 50% is the closest 
**	we can get.
**
** PARAMETERS:
**	InDS:		Name of input dataset
**	OutParms:	Name of output dataset with estimated parameter values. 
**	TargetProp:	Target proportion of observations with the simulated binary 
**				variable equal to one
**	LP_NoInt:	Specify the linear predictor 				
**					alpha1*x1_i+...+alphak*xk_i
**				of the logistic regression model (no intercept)	
**	ConvCri:	Convergence crieteria, calculated as the absolut difference
				between the induced and the target proportion of observations
				with the binary variable equal to one
**	MaxIte:		Maximum number of iterations to use in step 2 and step 3
**				before the macro is terminated.
**	Seed:		Specify the seed to be used. Default is 0. Choose a positive 
**				integer as seed to get reproducible simulations.
**	DEL:		Specify if the temporary datasets used in the macro is to be
**				deleted at the end of the macro. Y(Yes)/N(No). Default is Y.
**
**	NOTES:
**	The implemented approach is based on a approach described and implemented in
**	multiple studies by Austin, PC See for example
**	1. "Variance estimation when using inverse probability of treatment
**	   weighting (IPTW) with survival analysis" (2016) Austin. PC
**	2. "The performance of different porpensity score methods for 
**	   estimating marginal hazard ratios" (2013) Austin, PC
**
**	EXAMPLE OF USE:
**	%SimBinVar(
**		InDS=InDS
**		,OutDS=OutDS
**		,OutParms=OutDS_parms
**		,TargetProp=0.5
**		,LP_NoInt=log(2)*x1+log(3)*x2+log(0.7)*x3
**		,ConvCri=0.00001
**		,MaxIte=50
**		,seed=4545
**		,DEL=Y
**		);
**
********************************************************************************
**	AUTHOR:			Thomas Bøjer Rasmussen (TBR)
**	CREATION DATE:	2018/01/17
********************************************************************************
*	HISTORY:
*	1.0 MODIFIED:	2018/01/17  BY:TBR
*	- Created.
*******************************************************************************/	
%macro SimBinVar(
	InDS=
	,OutParms=
	,TargetProp=
	,LP_NoInt=
	,ConvCri=
	,MaxIte=50
	,seed=0
	,DEL=Y
	);

%put;
%put Macro &sysmacroname is executing;
%put;
%put;
%put Parameter list:;
%put _local_;
%put;

%local StartTime;
%let StartTime = %sysfunc(datetime(),datetime19.);
%put WARNING- Start Date: %sysfunc(datepart("&StartTime."dt),date9.);
%put WARNING- Start Time: %sysfunc(timepart("&StartTime."dt),time.);

/*******************************************************************************
STEP 1: INITIALIZE ALPHA_0
*******************************************************************************/
/* Setting the initial intercept in the model to zero, calculate the
predicted values, and use these to simulate an initial binary variable
for all observations. */
data _SB_data;
	call streaminit(&seed.);
	set &InDS.;
	length BinVar 3;
	format BinVar 1.;
	/* Calculate the linear predictor using a zero intercept value */ 
	LP_NoInt=&LP_NoInt.;
	/* Calculate predicted value */
	p=1/(1+exp(-0 - (LP_NoInt)));
	/* Simulation of binary variable */
	BinVar=rand("Bernoulli",p);
	keep BinVar LP_NoInt;
run;

/* Calculate the proportion of obervations with the binary variable */
proc means data=_SB_data mean noprint;
	var Binvar;
	output out=_SB_Prop(keep=Prop) 
		mean(BinVar)=Prop / noinherit;
run;

/* Create a dataset with initial parameter values */
data _SB_InitParms;
	set _SB_Prop;
	format	InitAlpha0Est InitProp 15.10;
	/* Initial alpha0 estimate */
	InitAlpha0Est=0; 
	/* Initial induced proportion of observations with the binary variable */
	InitProp=Prop; 
	drop Prop;
run;



/*******************************************************************************
STEP 2: FIND MIN AND MAX ALPHA_0 VALUES
*******************************************************************************/
/* Consider the total proportion of observations with the binary variable equal
to one. If the initial proportion is too low we add one to the value of alpha0, 
simulate BinVar again for all observations using the new alpha0-value 
(using the same seed), and calculate the new induced proportion. This is 
repeated until we get an alpha0 value that induces a proportion that is too high. 
Assume this happens in iteration n. We then know that the alpha0-value inducing 
the desired proportion is smaller than alpha0Max=alpha0(n) and larger than 
alpha0Min=alpha0(n-1), where alpha0(k) denotes the alpha0-value in iteration k. 
Similarly, if the initial proportion is too large, we can subtract one from 
alpha0 in each step and find alpha0Max=alpha0(n-1) and alpha0Min=alpha0(n). */

/* Create a dataset with intial parameter values */
%local Stop;
data _SB_MinMaxParms;
	set _SB_InitParms;
	format	Ite 3. CurrentAlpha0Est CurrentProp TargetProp Alpha0Min 
			Alpha0Max 15.10 Stop 1.;
	/* Iteration number */
	Ite=0;
	/* Current alpha0 estimate */
	CurrentAlpha0Est=InitAlpha0Est; 
	/* Current induced proportion of observations with the 
	binary variable equal to one */
	CurrentProp=InitProp; 
	/* Target proportion */
	TargetProp=&TargetProp.;
	/* Dummy variable used to indicate when a minimum and maximum alpha0
	value has been found */
	Stop=0;

	Alpha0Min=.;
	Alpha0Max=.;

	/* If the initial induced proportion is exactly the desired proportion
	we can set the min and max values equal to the initial alpha0
	estimate */
	if InitProp=TargetProp then do;
		Alpha0Min=InitAlpha0Est;
		Alpha0Max=InitAlpha0Est;
		Stop=1;
	end;

	/* Save value of Stop in a macro variable */
	call symput("Stop",Stop);
run;

/* Determine a minimum and maximum alpha0 value */

%local CurrentAlpha0Est PreviousAlpha0Est CurrentProp k OptionNotes;

/*Save current value of the Notes option in a macro variable */
%let OptionNotes=%sysfunc(getoption(notes));

options nonotes;
%let k=0;
%put WARNING- Determining alpha0 minimum and maximum value:;
%put WARNING- Iteration 0; 
%do %while (&Stop.=0 and &k.<=&MaxIte.);
	%let k=%sysevalf(&k.+1);

	/* Load the parameters from the previous iteration and
	calculate a new alpha0 estimate */
	data _SB_MinMaxParms_ite;
		set _SB_MinMaxParms end=eof;
		if eof;
		ite=&k.;
		/* Save the alpha0 estimate from last iteration in a 
		macro variable */
		call symput("PreviousAlpha0Est",CurrentAlpha0Est);
		/* Update the alpha0 estimate */
		if InitProp<TargetProp then CurrentAlpha0Est=CurrentAlpha0Est+1;
		if InitProp>TargetProp then CurrentAlpha0Est=CurrentAlpha0Est-1;
		/* Save the updated alpha0 estimate in a macro variable */
		call symput("CurrentAlpha0Est",CurrentAlpha0Est);
	run;

	/* Resimulate BinVar using the new alpha0-value */
	data _SB_data;
		call streaminit(&seed.);
		set _SB_data;
		p=1/(1+exp(-&CurrentAlpha0Est. - (LP_NoInt)));
		BinVar=rand("Bernoulli",p);
		drop p;
	run;

	/* Calculate the proportion of observations with the binary 
	variable equal to one */
	proc means data=_SB_data mean noprint;
		var Binvar;
		output out=_SB_Prop(keep=Prop) 
			mean(BinVar)=Prop / noinherit;
	run;

	/* Save the proportion in a macro variable */
	proc sql noprint;
		select Prop into :CurrentProp
			from _SB_Prop;
	quit;

	/* Evaluate if a minimum and maximum alpha0-value has been
	found.  */
	data _SB_MinMaxParms_ite;
		set _SB_MinMaxParms_ite;
		CurrentProp=&CurrentProp.;
		if	(InitProp<TargetProp and CurrentProp>=TargetProp) 
		then do;
			Alpha0Min=&PreviousAlpha0Est.;
			Alpha0Max=&CurrentAlpha0Est.;
			Stop=1;
		end;
		if	(InitProp>TargetProp and CurrentProp<=TargetProp) 
		then do;
			Alpha0Min=&CurrentAlpha0Est.;
			Alpha0Max=&PreviousAlpha0Est.;
			Stop=1;
		end;
		/* Save value of Stop in a macro variable */
		call symput("Stop",Stop);
	run;

	/* Add iteration parameters to parameter dataset */
	data _SB_MinMaxParms;
		set _SB_MinMaxParms _SB_MinMaxParms_ite;
	run;

	%put WARNING- Iteration &k.; 
%end; /*End of while-loop*/
option &OptionNotes.;


/*******************************************************************************
STEP 3: ESTIMATE ALPHA_0
*******************************************************************************/
/* A bisection approach can now be used to estimate the alpha0-value 
that induces the desired proportion arbitrarily close 
to the target proportion by using the following procedure:
1.	if CurrentProp<TargetProp then set
		CurrentAlpha0Est=CurrentAlpha0Est+abs(alpha0Max-alpha0Min)/2
	if CurrentProp>TargetProp then set
		CurrentAlpha0Est=CurrentAlpha0Est-abs(alpha0Max-alpha0Min)/2
2.	Calculate p, simulate BinVar , and calculate the proportion, CurrentProp.
3.	If CurrentProp>TargetProp then set Alpha0Max=CurrentAlpha0Est.
	If CurrentProp<TargetProp then set Alpha0Min=CurentAlpha0Est.
4.	Repeat 1.-3. until alpha0Est induces a proportion with an
	absolute difference from the target proportioon less than &ConvCri.: 
	abs(CurrentProp-TargetProp)<&ConvCri. */	


/* Create a dataset with initial parameter values */
data _SB_EstAlpha0Parms;
	set _SB_MinMAxParms end=eof;
	if eof;
	ite=0;
	Stop=0;

	/* If the current induced proportion is exactly the desired proportion
	we can set the min and max values equal to the initial alpha0
	estimate */
	if CurrentProp=TargetProp then do;
		Alpha0Min=CurrentAlpha0Est;
		Alpha0Max=CurrentAlpha0Est;
		Stop=1;
	end;

	/* Save value of Stop in a macro variable */
	call symput("Stop",Stop);

	drop InitAlpha0Est InitProp;
run;

%let k=0;
%put WARNING- Estimating alpha0 using bisection method:;
%put WARNING- Iteration 0; 

/* Iteratively estimate the alpha0-value that induces the desired 
proportion of observations with BinVar=1  */
options nonotes;
%local Diff;
%do %while (&Stop.=0 and &k.<=&MaxIte.);
	%let k=%sysevalf(&k.+1);
	/* Loading the parameters from the previous iteration and
	calculate a new alpha0 estimate from the min and max values
	of the previous iteration */
	data _SB_EstAlpha0Parms_ite;
		set _SB_EstAlpha0Parms end=eof;
		if eof;
		ite=&k.;
		/* Update the alpha0 estimate */
		if CurrentProp<TargetProp then do;
			CurrentAlpha0Est=CurrentAlpha0Est+abs(Alpha0Max-Alpha0Min)/2;
		end;
		if CurrentProp>TargetProp then do;
			CurrentAlpha0Est=CurrentAlpha0Est-abs(Alpha0Max-Alpha0Min)/2;
		end;		
		/* Save updated alpha0 estimate in macro variable */
		call symput("CurrentAlpha0Est",CurrentAlpha0Est);

	run;

	/* Resimulate BinVar using the new alpha0-value */
	data _SB_data;
		call streaminit(&seed.);
		set _SB_data;
		p=1/(1+exp(-&CurrentAlpha0Est. - (LP_NoInt)));
		BinVar=rand("Bernoulli",p);
		drop p;
	run;

	/* Calculate the proportion of observations with the binary 
	variable equal to one */
	proc means data=_SB_data mean noprint;
		var Binvar;
		output out=_SB_Prop(keep=Prop) 
			mean(BinVar)=Prop / noinherit;
	run;

	/* Save the proportion in a macro variable */
	proc sql noprint;
		select Prop into :CurrentProp
			from _SB_Prop;
	quit;

	/* Evaluate if the alpha0 estimate induces a proportion that is
	close enought to the target */
	data _SB_EstAlpha0Parms_ite;
		set _SB_EstAlpha0Parms_ite;
		format Diff 15.10;
		CurrentProp=&CurrentProp.;
		if CurrentProp<TargetProp then Alpha0Min=CurrentAlpha0Est;
		if CurrentProp>TargetProp then Alpha0Max=CurrentAlpha0Est;

		Diff=abs(CurrentProp-TargetProp);
		if Diff<&ConvCri. then Stop=1;

		/* Save value of Stop in a macro variable */
		call symput("Stop",Stop);
		/* Save value of Diff in a macro variable */
		call symput("Diff",put(Diff,15.10));
	run;

	/* Add iteration parameters to parameter dataset */
	data _SB_EstAlpha0Parms;
		set _SB_EstAlpha0Parms _SB_EstAlpha0Parms_ite;
	run;

	%put WARNING- Iteration: &k.;
	%put WARNING- Diff: &Diff.;

%end; /*End of until-loop*/
options &OptionNotes.;


/* Create dataset with final parameters */
data &OutParms.(rename=(CurrentAlpha0Est=Alpha0Est CurrentProp=PropEst));
	set _SB_EstAlpha0Parms end=eof;
	if eof;
	keep CurrentAlpha0Est CurrentProp TargetProp;
run;

%if &DEL=Y %then %do;
	proc datasets nodetails nolist;
		delete _SB_:;
	run;
	quit;
%end;

%local EndTime;
%let EndTime = %sysfunc(datetime(),datetime19.);
%put WARNING- End Date: %sysfunc(datepart("&EndTime."dt),date9.);
%put WARNING- End Time: %sysfunc(timepart("&EndTime."dt),time.);

%put;
%put Macro &sysmacroname has ended.;
%put;

%mend SimBinVar;











	





