/*******************************************************************************
**	DESCRIPTION:
**	Simulate variables using a multivariate normal distribution with 
**	zero mean vector and a covariance matrix with aunit diagonal, eg the 
**	covariance matrix is also the correlation matrix.
**
** PARAMETERS:
**	OutDS:		Name of output dataset with the simulated data	
**	NObs:		Number of observations in each simulated sub dataset	
**	NSim:		Number of simulated sub datasets
**	NVars:		Number of variables to simulate
**	ND:			Number variables that is to be dichotomized 
**	VarCorr:	Speficy the correlation between variables
**	BlockSize:	Number of simulated sub datasets in each block of simulations
**				(see notes below). NSim must be divisible with BlocksSize.
**	Seed:		Specify the seed to be used. Default is 0. Choose a positive 
**				integer as seed to get reproducible simulations.
**	DEL:		Specify if the temporary datasets used in the macro is to be
**				deleted at the end of the macro. Y(Yes)/N(No). Default is Y.
**
**	NOTES:
**	-The output data set contains the variables Dataset and Obs denoting 
**	 the dataset and observation number.
**	-The data is generated in proc IML. When a matrix is created in proc IML
**	 the entire matrix needs to fit in the memory in order for SAS to process it
**	 (unlike datasets). When simulating large amounts of data this can result
**	 in the system running out of memory. To bypass this problem,
**	 the BlockSize parameter can be used to split the simulation of data into 
**	 blocks where each block cointains the specified number of simulated 
**	 subdatasets. These are then combined into a single dataset in the end 
**	 of the macro.
**	-The ND variable can be used to dichotomize some/all of the generated
**	 normal distributed variables. This is done using zero as the cut-off, 
**	 resulting in a prevalence of 50% of the dichotomized variable.
**	-Dichotomzing variables will result in moderate changes in correlations
**	 between the dichomized variables and the rest of the variables.
**
**	EXAMPLE OF USE:
**	%SimVarsMultiNormal(
**		OutDS=test1
**		,NObs=1000
**		,Nsim=1000
**		,NVars=5
**		,ND=3
**		,VarCorr=0.5
**		,BlockSize=1000
**		,seed=4343
**		,DEL=Y
**		);
**
********************************************************************************
**	AUTHOR:			Thomas Bøjer Rasmussen (TBR)
**	CREATION DATE:	2018/01/15
********************************************************************************
*	HISTORY:
*	1.0 MODIFIED:	2018/01/15  BY:TBR
*	- Created.
*******************************************************************************/	
%macro SimVarsMultiNormal(
	OutDS=
	,NObs=
	,Nsim=
	,NVars=
	,ND=0
	,VarCorr=0
	,BlockSize=
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

	%local StartTime EndTime;
	%let StartTime = %sysfunc(datetime(),datetime19.);
	%put WARNING- Start Date: %sysfunc(datepart("&StartTime."dt),date9.);
	%put WARNING- Start Time: %sysfunc(timepart("&StartTime."dt),time.);

	/* Calculating the number of blocks to divide the simulations into. */
	%let NBlocks=%sysevalf(&Nsim./&BlockSize.);
	
	/* Simulation of variables */
	%local i j k;
	proc iml;
		call randseed(&seed.);
		%do i=1 %to &NBlocks.;
			%put WARNING- Simulation Block: &i.;
			/* Specifying the mean vector */
			mu={%do j=1 %to &NVars.; 0 %end;};
			/* Specifying the covariance/corrleation matrix */
			S={
				%do j=1 %to &NVars.;
					%do k=1 %to &NVars.;
						/* Unit diagnoal */
						%if &j.=&k. %then %do;
							1.0
						%end;
						/* Specified correlation elsewhere*/
						%if &j. NE &k. %then %do;
							&VarCorr.
						%end;
						/* Insert commas at the end of each row to indicate
						row changes */
						%if &k.=&NVars. and &j. NE &NVars. %then %do; 
							,
						%end;
					%end;
				%end;
				};

			/* Simulate the block of data */
			DataNormal=randnormal(%sysevalf(&NSim.*&NObs./&NBlocks.),mu,S);

			/* Making an observation vector that denotes the observation
			number in each sub dataset. */
			Obs=repeat(T(do(1,&NObs.,1)),&BlockSize.,1);

			/* Making a dataset vector that is one in the first dataset,
			two in the next data set and so on. */
			Dataset=T(do(%sysevalf(1+(&i.-1)*&BlockSize.),%sysevalf(&i.*&BlockSize.),1))@j(&NObs.,1,1);

			/* Combining the data matrices and the vectors into a single 
			matrix */
			Data=Dataset||Obs||DataNormal;

			/* Saving the combined matrix as a SAS dataset. */
			create _SV_&i. from Data[colname={"Dataset" "Obs" %do j=1 %to &NVars.; "x&j." %end;}];
			append from Data;
			close _SV_&i.;
		%end; /* End of i-loop */
	quit;

	/* A numeric variable with a length of 4 bytes can accurately store an 
	integer up to at least 2,000,000, thus preserving the correct values
	of the Dataset and Obs variable. Using a length of 3 bytes, simulated data 
	from the standard normal distribution seems to have 2-3 correct digits depending on
	the sign of the number. Since the number of correct digits is irrelevant, 
	storing the values in 3 bytes is fine, and will greatly reduce the size of
	the size dataset. Correlations between the correlated variables might be 
	slightly affected? */
	data &OutDS.;
		format 
			Dataset Obs 10. 
			%do i=1 %to &ND.; x&i. 1. %end;
			%do i=%eval(&ND.+1) %to &Nvars.; x&i. 10.2 %end;
			;
		set 
		%do i=1 %to &NBlocks.;
			_SV_&i.
		%end;
		;
		length 
			Dataset Obs 4 
			%do i=1 %to &ND.; x&i. 3 %end;
			%do i=%eval(&ND.+1) %to &Nvars.; x&i. 4 %end;
			;
		%do j=1 %to &ND.;
			format x&j. 1.;
			x&j.=(x&j.>0);
		%end;
	run;

	%if &DEL.=Y %then %do;
		proc datasets nodetails nolist;
			delete _SV_:;
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

%mend SimVarsMultiNormal;
