**************************
*
* P8108 Group Final Project
*
* PGM: Data Cleaning
*
* AUTHORS: Waveley Qiu
*
* DATE: 2022-12-02
*
****************************;

%let fp = C:\Users\14088\Documents\R\RWorkspace\coursework\P8108\p8108_group_project;

proc import file="&fp/data/bladder.csv"
    dbms=csv out=bladder replace;
run;

proc sort data = bladder; by id; run;

* make variable collapsing status = 2 and 3 categories;
data bladder;
	set bladder;
	retained_status = status;
	if status eq 2 or status eq 3 then status = 2;
run;

data bladder_first_id;
	set bladder;
	by id;
	if first.id;
run;

%macro get_bladder_dataset(indf, iddf, status, outdf, keep_non_events = 0);
proc sort data = &indf. out = sorted_end_event_df; where status eq &status.; by id enum; run;

data event_df;
	set sorted_end_event_df;
	by id enum;
	if first.id;
run;

data stacked_bladder;
	set &iddf. event_df;
run;

proc sort data = stacked_bladder noduprec; by id enum; run; 

data &outdf;
	format id treatment number size recur start stop status rtumor rsize enum;
	set stacked_bladder;
	by id enum;
	retain start_time;
	if first.id then start_time = start;
	if last.id;
	drop start_time;
run;
%mend;

* time to first recurrence;
%get_bladder_dataset(bladder, bladder_first_id, 1, ttfr);

* time to death after first recurrence;
/*proc sort data = fin_bladder; by id; run;

data bladder_first_recurrence;
	set ttfr;
	by id;
	start = stop;
	if recur > 0;
run;

data recurrence;
	set bladder;
	if recur > 0;
run;

* for ttdafr: status == 1 indicates that patient did not die, status == 2 indicates that patient died;
%get_bladder_dataset(recurrence, bladder_first_recurrence, 2, ttdafr);
*/

data ttfr;
	set ttfr;
	if treatment eq "placebo" then rx = 0;
	if treatment eq "pyridox" then rx = 1;
	if treatment eq "thiotep" then rx = 2;
	drop treatment;
run;
*Model 1;
proc phreg data=ttfr plots=survival;
    class rx(ref='0');
    model (start,stop)*status(0) = rx size; 
run;

*Model 2;
proc phreg data=ttfr plots=survival;
    class rx(ref='0');
    model (start,stop)*status(0) = rx number; 
run;

*Model 3;
proc phreg data=ttfr plots=survival;
    class rx(ref='0');
    model (start,stop)*status(0) = rx number size; 
run;
/*Use the baseline statement to generate survival plots by group*/

data covar;
    input rx number size;
    datalines;
0 2.025862069 2.0344827586
1 2.025862069 2.0344827586
2 2.025862069 2.0344827586
;
run;

title "Overlaid Survival Curves";

proc phreg data = ttfr plots(overlay)=(survival);
    class rx(ref='0');
    model (start,stop)*status(0) = rx number size;
    baseline covariates=covar out=baseline / rowid=rx;
run;

/*Expanding and interpreting the Cox regression model with interaction terms*/
title "Expanding the Cox regression model with interaction terms";

proc phreg data = ttfr;
    class treatment(ref="placebo");
    model (start,stop)*status(0) = treatment|number number|size;
run;

ods text="Using hazard ratio and graphs to interpret effects, particularly interactions";
proc phreg data = ttfr;
    format treatment treatment.;
    class treatment(ref="placebo");
    model (start,stop)*status(0) = treatment|number number|size;
    hazardratio 'Effect of 1-unit change in size by treatment' size / at(treatment=ALL);
    hazardratio 'Effect of 1-unit change in size across number' size / at(number=1 2 3 4 5 6 7 8);
run;

/*check PH assumptions*/
title "Check PH assumptions using log(-log)";

/*using log(-log)*/
data ttfr;
    set ttfr;
    time = stop - start;
run;

proc lifetest data=ttfr method=km plots=(survival (cl atrisk=0 to 800
by 80) lls) outsurv=survival;
    strata rx;
    time time*status(0);
run;

title "Check PH assumptions by plotting observed vs. fitted";

/*observed vs. fitted*/
data cov;
    length ID $20;
    input treatment id $3-23;
    datalines;
0 placebo
1 pyridoxine
2 thiotepa
;
run;

proc phreg data=ttfr plots(overlay)=survival;
    model time*status(0)=rx/ rl ties = efron;
    baseline covariates=cov out = pred1 survival = _all_ ; 
run;
    
proc lifetest data=ttfr outsurv=km1 plots=s;
    time time*status(0);
    strata rx / adjust = bonferroni;
run;

data trunc_ttfr;
	set ttfr;
	if time > 10;
run;
proc lifetest data=trunc_ttfr outsurv=km1 plots=s;
    time time*status(0);
    strata rx / adjust = bonferroni;
run;

 
proc lifetest data=ttfr outsurv=km1 plots=s;
	where time > 30;
    time time*status(0);
    strata rx;
run;

data all;
    set km1 (in=a) pred1;
    if a and treatment=0 then ID="Observed Placebo";
    if a and treatment=1 then ID="Observed Pyridoxine";
    if a and treatment=2 then ID="Observed Thiotepa";
run;

proc sgplot data=all noborder;
    step x=time y=survival/ group=ID name='s';
    keylegend 's' / linelength=20;
run;

title "Check PH assumptions using Schoenfeld Residuals";
/*Schoenfeld Residuals*/
proc phreg data=ttfr;
    class treatment(ref='0');
    model time*status(0)=treatment size number;
    output out=resid ressch=schtreatment schsize schnumber;
run;

ods text="Plot Schoenfeld residual with variable = number";
proc loess data=resid;
    model schnumber=time / smooth=(0.2 0.4 0.6 0.8);
run;

proc rank data=resid out=resrank ties=mean;
    var time;
    ranks time_rank;
run;

proc corr data=resrank;
    var schnumber;
    with time_rank;
run;

proc reg data=resrank;
    model time_rank=schnumber;
run;

ods text="plot Schoenfeld residual with continuous variable = size";
proc loess data=resid;
    model schsize=time / smooth=(0.2 0.4 0.6 0.8);
run;

proc rank data=resid out=resrank ties=mean;
    var time;
    ranks time_rank;
run;

proc corr data=resrank;
    var schsize;
    with time_rank;
run;

proc reg data=resrank;
    model time_rank=schsize;
run;
