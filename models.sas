/*Create a library for datasets*/

libname proj "/home/u60468981/project";

proc import file="/home/u60468981/project/b1.xlsx"
    dbms=xlsx out=b1 replace;
run;

ods pdf file="/home/u60468981/project/pro.pdf";
ods graphics on;

*Model 1;
proc phreg data=b1 plots=survival;
    class treatment;
    model (start,stop)*status(0) = treatment size; 
run;

*Model 2;
proc phreg data=b1 plots=survival;
    class treatment;
    model (start,stop)*status(0) = treatment number; 
run;

*Model 3;
ods text="A simple Cox regression model 
with treatment as a categorical predictor was fitted through 
proc phreg.";
proc phreg data=b1 plots=survival;
    class treatment;
    model (start,stop)*status(0) = treatment number size; 
run;

/*Use the baseline statement to generate survival plots by group*/
ods text="Expanding and interpreting the Cox regression model 
with interaction terms";

data covar;
    input treatment number size;
    datalines;
0 2.025862069 2.0344827586
1 2.025862069 2.0344827586
2 2.025862069 2.0344827586
;
run;

proc datasets library=work;
   modify covar;
   format treatment BEST.;
quit;

proc phreg data = b1 plots(overlay)=(survival);
    class treatment;
    model (start,stop)*status(0) = treatment number size;
    baseline covariates=covar out=baseline / rowid=treatment;
run;

proc format;
    value treatment 0 = "placebo" 
                    1 = "pyridoxine"
                    2 = "thiotepa";
run;


/*Expanding and interpreting the Cox regression model with interaction terms*/
ods text="Expanding and interpreting the Cox regression model 
with interaction terms";
proc phreg data = b1;
    format treatment treatment.;
    class treatment;
    model (start,stop)*status(0) = treatment|number number|size;
run;

ods text="Using hazard ratio and graphs to interpret effects, particularly interactions";
proc phreg data = b1;
    class treatment;
    model (start,stop)*status(0) = treatment|number number|size;
    hazardratio 'Effect of 1-unit change in size by treatment' size / at(treatment=ALL);
    hazardratio 'Effect of 1-unit change in size across number' size / at(number=1 2 3 4 5 6 7 8);
run;

/*check PH assumptions*/
ods text="check PH assumptions using log(-log)";
/*using log(-log)*/
data b1;
    set b1;
    time = stop - start;
run;

proc lifetest data=b1 method=km plots=(survival (cl atrisk=0 to 800
by 80) lls) outsurv=survival;
    strata treatment;
    time time*status(0);
run;

ods text="Check PH assumptions by plotting observed vs. fitted";
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

proc phreg data=b1 plots(overlay)=survival;
    model time*status(0)=treatment/ rl ties = efron;
    baseline covariates=cov out = pred1 survival = _all_ ; 
run;
    
proc lifetest data=b1 outsurv=km1 plots=s;
    time time*status(0);
    strata treatment;
run;

data all;
    set km1 (in=a) pred1;
    if a and treatment=0 then ID="Observed Placebo";
    if a and treatment=1 then ID="Observed Pyridoxine";
    if a and treatment=2 then ID="Observed Thiotepa";
run;

title 'Observed vs. Fitted';
proc sgplot data=all noborder;
    step x=time y=survival/ group=ID name='s';
    keylegend 's' / linelength=20;
run;

ods text="Check PH assumptions using Schoenfeld Residuals";
/*Schoenfeld Residuals*/
proc phreg data=b1;
    class treatment;
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

ods graphics off;
ods pdf close;


