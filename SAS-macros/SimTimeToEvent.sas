/*******************************************************************************
**	DESCRIPTION:
**	The macro takes an input dataset with variable information on a set of
**	observations and simulates a binary indicator variable and a 
**	time-to-event variable with a specified marginal effect of the binary
**	variable on the time to event. Estimated parameters are saved in an output
**	dataset.
**
**	Given the input data, a Time-to-event variable and a binary variable (BinVar)
**	is simulated so that the following Cox proportional hazard models hold:
**		h(t|x)=h0(t)*exp(beta0*BinVar+beta1*x1+...+betak*xk)
**	Here h0(t) is a specified Weibull distribution, beta0,...,betak is 
**	regression parameters, and BinVar, x1,...,xk is the associated independent
**	variables. Beta0 is the CONDITIONAL log(HR) effect of Binvar on the 
**	time-to-event variable. The macro uses a iterative bisection approach to 
**	estimate the conditonal effect that will induce the desired MARGINAL
**	effect. The macro achives this in three steps:
**	1. An initial model is specified with beta0 chosen as the logarithm of
**	the specified marginal HR to be induced, ie the conditional effect is set
**	to the marginal effect of interest. 
**	2. An interval [beta0Min,beta0Max] of conditonal effects are found in which
**	we know there is a contional effect thast will induce the desried marginal 
**	effect
**	3. The beta0 parameter that induces the desired marginal effect is estimated
**	using a iterative bisection approach.
**
** PARAMETERS:
**	InDS:		Name of input dataset
**	OutParms:	Name of output dataset with (estimated) parameter values. 
**	TargetMar:	Target marginal effect of BinVar on the time-to-event 
**				variable
**	LP_NoInd:	Specify the linear predictor (without the indicator variable)				
**					beta1*xk+...+betak*xk
**	ConvCri:	Convergence crieteria, calculated as the absolut difference
**				between the induced and the target marginal effect of BinVar
**	MaxIte:		Maximum number of iterations to use in step 2 and step 3
**				before the macro is terminated.
**	lambda:		Specify the lambda paramter from the Weibull distribution
**	eta:		Specify the eta paramter from the Weibul distribtuion
**	Seed:		Specify the seed to be used. Default is 0. Choose a positive 
**				integer as seed to get reproducible simulations.
**	DEL:		Specify if the temporary datasets used in the macro is to be
**				deleted at the end of the macro. Y(Yes)/N(No). Default is Y.
**
**	NOTES:
**	The implemented method is based on the methods described in
**	1.	"The performance of different porpensity score methods for 
**		estimating marginal hazard ratios" (2013) Austin, PC
**	2.	"Generating survival times to simulate Cox proportional hazard 
**		models" (2005) Bender, R et al
**
**	EXAMPLE OF USE:
**	%SimTimeToEvent(
**		InDS=InDS
**		,OutParms=OutDS
**		,TargetMar=3
**		,LP_NoInd=log(1.5)*x1+log(2)*x2
**		,ConvCri=0.000001
**		,MaxIte=30
**		,lambda=0.00002
**		,eta=2
**		,seed=543435
**		,DEL=Y
**		);
********************************************************************************
**	AUTHOR:			Thomas Bøjer Rasmussen (TBR)
**	CREATION DATE:	2018/01/19
********************************************************************************
*	HISTORY:
*	1.0 MODIFIED:	2018/01/19  BY:TBR
*	- Created.
*******************************************************************************/
%macro SimTimeToEvent(
	InDS=
	,OutParms=
	,TargetMar=
	,LP_NoInd=
	,ConvCri=
	,MaxIte=
	,lambda=0.00002
	,eta=2
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
STEP 1: INITIALIZE BETA0
*******************************************************************************/
/* Load data and simulate the potential time-to-event for both
values of a simulated binary variable */
%local i;
data _ST_data;
	call streaminit(&seed.);
	set &InDS.;
	LP_NoInd=&LP_NoInd.;
	%do i=0 %to 1;
		IndVar=&i.;
		u=rand("Uniform");
		/* Simulation of a time-to-event outcome from a Cox PH model with
		a Weibull baseline hazard function,
		and with the specified conditional effects. See Bender (2005) */
		TimeToEvent=(-log(u)/(&lambda.*exp(log(&TargetMar.)*IndVar+LP_NoInd)))**(1/&eta.);
		output;
	%end;
	keep IndVar LP_NoInd u TimeToEvent;
run;

/* Estimate the marginal effect of IndVar that the conditional effect induces */
ods listing exclude all;
proc phreg data=_ST_data;
	class IndVar(ref="0") / param=ref;
	model TimeToEvent=IndVar;
	ods output ParameterEstimates=_ST_InducedMarEffect(keep=Estimate);
run;
ods listing select all;

/* make a dataset with initial parameter aestimates */
data _ST_InitParms;
	format  TargetMarEffect TargetLogMarEffect 15.10;
	set _ST_InducedMarEffect(rename=(Estimate=InitLogMarEffect));
	label InitLogMarEffect=" ";
	InitBeta0Est=log(&TargetMar.);
	TargetMarEffect=&TargetMar.;
	TargetLogMarEffect=log(&TargetMar.);
run;
	

/*******************************************************************************
STEP 2: FIND MIN AND MAX BETA0 VALUES
*******************************************************************************/
/* Consider the induced marginal effect of Binvar. The marginal effect is 
always closer to the the null-association than the conditional effect is
(this is not strictly true, but it is not relevant here).
If the initially induced marginal effect is too small we can increase the
beta0 parameter and recaulculate hte induced marginal effect. This can be
repeated until we have found a beta-value that induces a marginal effect that is
higher than desired. We then know that the beta0-value inducing 
the desired marginal effect is smaller than beta0Max=beta0(n) and larger than 
beta0Min=beta0(n-1), where beta0(k) denotes the beta0-value in iteration k. 
Similarly, if the initially induced marginal effect is too large, we can 
decrease  beta0 in each step until the induced effect is too small and
we have that beta0Max=beta0(n-1) and beta0Min=alpha0(n). */
%local Stop;
data _ST_MinMaxParms;
	set _ST_InitParms;
	format Ite 3. Stop 1. CurrentBeta0Est 15.10;
	Ite=0;
	CurrentBeta0Est=InitBeta0Est;
	CurrentLogMarEffect=InitLogMarEffect;
	Stop=0;

	Beta0Min=.;
	Beta0Max=.;

	if InitLogMarEffect=TargetLogMarEffect then do;
		Beta0Min=InitBeta0Est;
		Beta0Max=InitBeta0Est;	
		Stop=1;	
	end;

	/* Save value of Stop in a macro variable */
	call symput("Stop",Stop);	
run;

%local k PreviousBeta0Est CurrentBeta0Est CurrentLogMarEffect OptionNotes;

/*Save current value of the Notes option in a macro variable */
%let OptionNotes=%sysfunc(getoption(notes));

option nonotes;
%let k=0;
%put WARNING- Determining minimum and maximum beta0-value:;
%put WARNING- Iteration 0; 
%do %while(&Stop.=0 and &k.<=&MaxIte.);
	%let k=%sysevalf(&k.+1);
	/* Load the parameters from the previous iteration and
	set a new beta0 estimate */
	data _ST_MinMaxParms_ite;
		set _ST_MinMaxParms end=eof;
		if eof;
		ite=&k.;
		/* Save beta0 estimate from last iteration in a 
		macro variable */
		call symput("PreviousBeta0Est",CurrentBeta0Est);
		/* Update the beta0 estimate */
		if InitLogMarEffect<TargetLogMarEffect 
			then CurrentBeta0Est=CurrentBeta0Est+1;
		if InitLogMarEffect>TargetLogMarEffect 
			then CurrentBeta0Est=CurrentBeta0Est-1;
		/* Save the updated beta0 estimate in a macro variable */
		call symput("CurrentBeta0Est",CurrentBeta0Est);
	run;

	data _ST_data;
		set _ST_data;
		TimeToEvent=(-log(u)/(&lambda.*exp(&CurrentBeta0Est.*IndVar+LP_NoInd)))**(1/&eta.);
	run;

	ods listing exclude all;
	proc phreg data=_ST_data;
		class IndVar(ref="0") / param=ref;
		model TimeToEvent=IndVar;
		ods output ParameterEstimates=_ST_InducedMarEffect(keep=Estimate);
	run;
	ods listing select all;

	/* Save the estimate in a macro variable */
	proc sql noprint;
		select Estimate into :CurrentLogMarEffect
			from _ST_InducedMarEffect;
	quit;

	/* Evaluate if a minimum and maximumbeta0-value has been
	found.  */
	data _ST_MinMaxParms_ite;
		set _ST_MinMaxParms_ite;
		CurrentLogMarEffect=&CurrentLogMarEffect.;
		if	(InitLogMarEffect<TargetLogMarEffect 
				and CurrentLogMarEffect>=TargetLogMarEffect) 
		then do;
			Beta0Min=&PreviousBeta0Est.;
			Beta0Max=&CurrentBeta0Est.;
			Stop=1;
		end;
		if	(InitLogMarEffect>TargetLogMarEffect 
				and CurrentLogMarEffect<=TargetLogMarEffect) 
		then do;
			Beta0Min=&CurrentBeta0Est.;
			Beta0Max=&PreviousBeta0Est.;
			Stop=1;
		end;
		/* Save value of Stop in a macro variable */
		call symput("Stop",Stop);
	run;

	/* Add iteration parameters to parameter dataset */
	data _ST_MinMaxParms;
		set _ST_MinMaxParms _ST_MinMaxParms_ite;
	run;
	%put WARNING- Iteration &k.; 
%end; /* End of while-loop */
option &OptionNotes.;


/*******************************************************************************
STEP 3: ESTIMATE BETA0
*******************************************************************************/
/* A bisection approach can now be used to estimate the beta0-value 
that induces the desired marginal effect arbitrarily close 
to the target marginal effect by using the following procedure:
1.	if CurrentLogMarEffect<TargetLogMarEffect then set
		CurrentBeta0Est=CurrentBeta0Est+abs(Beta0Max-Beta0Min)/2
	if CurrentLogMarEffect>TargetLogMarEffect then set
		CurrentBeta0Est=CurrentBetaa0Est-abs(Beta0Max-Beta0Min)/2
2.	Resimulate both potential Time-to-events and estimate the 
	induced marginal effect
3.	If CurrentLogMarEffect>TargetLogMarEffect then set Beta0Max=CurrentBeta0Est.
	If CurrentLogMarEffect<TargetLogMarEffect then set Beta0Min=CurentBeta0Est.
4.	Repeat 1.-3. until Beta0Est induces a marginal effect with an
	absolute difference from the target marginal effect less than &ConvCri.: 
	abs(CurrentLogMarEffect-TargetLogMarEffect)<&ConvCri. */	

/* Create a dataset with initial parameter values */
data _ST_EstParms;
	set _ST_MinMaxParms end=eof;
	format Diff 15.10;
	if eof;
	ite=0;
	Stop=0;
	Diff=.;

	/* If the current induced marginal effect is exactly the desired 
	marginal effect we can set the min and max beta0-values equal to 
	the initial beta0 estimate */
	if CurrentLogMarEffect=TargetLogMarEffect then do;
		Beta0Min=CurrentBeta0Est;
		Beta0Max=CurrentBeta0Est;
		Stop=1;
	end;

	/* Save value of Stop in a macro variable */
	call symput("Stop",Stop);

	drop InitLogMarEffect InitBeta0Est;
run;


%let k=0;
%put WARNING- Estimating beta0 using bisection method:;
%put WARNING- Iteration 0; 

/* Iteratively estimate the beta0-value that induces the desired 
marginal effect  */
options nonotes;
%local Diff;
%do %while (&Stop.=0 and &k.<=&MaxIte.);
	%let k=%sysevalf(&k.+1);
	/* Loading the parameters from the previous iteration and
	calculate a new CondEst estimate from the min and max values
	of the previous iteration */
	data _ST_EstParms_ite;
		set _ST_EstParms end=eof;
		if eof;
		ite=&k.;
		/* Update the Beta0Est estimate */
		if CurrentLogMarEffect<TargetLogMarEffect then do;
			CurrentBeta0Est=CurrentBeta0Est+abs(Beta0Max-Beta0Min)/2;
		end;
		if CurrentLogMarEffect>TargetLogMarEffect then do;
			CurrentBeta0Est=CurrentBeta0Est-abs(Beta0Max-Beta0Min)/2;
		end;		
		/* Save updated beta0 estimate in macro variable */
		call symput("CurrentBeta0Est",CurrentBeta0Est);
	run;

	/* Resimulate TimeToEvent using the new beta0 estimate */
	data _ST_data;
		set _ST_data;
		TimeToEvent=(-log(u)/(&lambda.*exp(&CurrentBeta0Est.*IndVar+LP_NoInd)))**(1/&eta.);
	run;

	/* Estimate the induced marginal effect */
	ods listing exclude all;
	proc phreg data=_ST_data;
		class IndVar(ref="0") / param=ref;
		model TimeToEvent=IndVar;
		ods output ParameterEstimates=_ST_InducedMarEffect(keep=Estimate);
	run;
	ods listing select all;

	/* Save the estimate in a macro variable */
	proc sql noprint;
		select Estimate into :CurrentLogMarEffect
			from _ST_InducedMarEffect;
	quit;

	/* Evaluate if the alpha0 estimate induces a proportion that is
	close enought to the target */
	data _ST_EstParms_ite;
		set _ST_EstParms_ite;
		CurrentLogMarEffect=&CurrentLogMarEffect.;
		if CurrentLogMarEffect<TargetLogMarEffect 
			then Beta0Min=CurrentBeta0Est;
		if CurrentLogMarEffect>TargetLogMarEffect 
			then Beta0Max=CurrentBeta0Est;

		Diff=abs(CurrentLogMarEffect-TargetLogMarEffect);
		if Diff<&ConvCri. then Stop=1;

		/* Save value of Stop in a macro variable */
		call symput("Stop",Stop);
		/* Save value of Diff in a macro variable */
		call symput("Diff",put(Diff,15.10));
	run;

	/* Add iteration parameters to parameter dataset */
	data _ST_EstParms;
		set _ST_EstParms _ST_EstParms_ite;
	run;

	%put WARNING- Iteration: &k.;
	%put WARNING- Diff: &Diff.;

%end; /*End of until-loop*/
options &OptionNotes.;

/* Create dataset with final parameters */
data &OutParms.(rename=(
		CurrentBeta0Est=Beta0Est 
		CurrentLogMarEffect=InducedLogMarEffect
		));
	set _ST_EstParms end=eof;
	if eof;
	InducedMarEffect=exp(CurrentLogMarEffect);
	drop Ite Stop Beta0Min Beta0Max Diff;
run;


%if &DEL=Y %then %do;
	proc datasets nodetails nolist;
		delete _ST_:;
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



%mend SimTimeToEvent;
