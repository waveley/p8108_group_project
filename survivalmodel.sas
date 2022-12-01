/*Create a library for datasets*/

libname proj "/home/u60468981/project";

proc import file="/home/u60468981/project/b1.xlsx"
    dbms=xlsx out=b1 replace;
run;

ods pdf file="/home/u60468981/project/models.pdf";
*Model 1;
proc phreg data=b1 plots=survival;
    class treatment(ref='0');
    model (start,stop)*status(0) = treatment size; 
run;

*Model 2;
proc phreg data=b1 plots=survival;
    class treatment(ref='0');
    model (start,stop)*status(0) = treatment number; 
run;

*Model 3;
ods graphics on;
title "Cox regression model with treatment as a categorical predictor";
proc phreg data=b1 plots=survival;
    class treatment(ref='0');
    model (start,stop)*status(0) = treatment number size; 
run;
ods graphics off;

/*Use the baseline statement to generate survival plots by group*/

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

ods graphics on;
title "Overlaid Survival Curves";

proc phreg data = b1 plots(overlay)=(survival);
    class treatment(ref='0');
    model (start,stop)*status(0) = treatment number size;
    baseline covariates=covar out=baseline / rowid=treatment;
run;
ods graphics off;

proc format;
    value treatment 0 = "placebo" 
                    1 = "pyridoxine"
                    2 = "thiotepa";
run;

/*Expanding and interpreting the Cox regression model with interaction terms*/
ods graphics on;
title "Expanding the Cox regression model with interaction terms";

proc phreg data = b1;
    format treatment treatment.;
    class treatment(ref="placebo");
    model (start,stop)*status(0) = treatment|number number|size;
run;

ods text="Using hazard ratio and graphs to interpret effects, particularly interactions";
proc phreg data = b1;
    format treatment treatment.;
    class treatment(ref="placebo");
    model (start,stop)*status(0) = treatment|number number|size;
    hazardratio 'Effect of 1-unit change in size by treatment' size / at(treatment=ALL);
    hazardratio 'Effect of 1-unit change in size across number' size / at(number=1 2 3 4 5 6 7 8);
run;
ods graphics off;

/*check PH assumptions*/
ods graphics on;
title "Check PH assumptions using log(-log)";

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

ods graphics off;

ods graphics on;
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

proc sgplot data=all noborder;
    step x=time y=survival/ group=ID name='s';
    keylegend 's' / linelength=20;
run;
ods graphics off; 

ods graphics on;
title "Check PH assumptions using Schoenfeld Residuals";
/*Schoenfeld Residuals*/
proc phreg data=b1;
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

ods graphics off;

ods graphics on;
title 'Exponential Survival Model';
data b1;
    set b1;
    time = stop - start;
run;

proc lifereg data=b1 noprint;
  model time*status(0) = treatment size number/ distribution=exponential;
  output out=exp cdf=f;
run;

data exp1;
    set exp;
    cox = -log(1-f);
run;

proc lifetest data=exp1 outsurv=surv_exp noprint;
  time cox*status(0);
run;

data surv_exp;
  set surv_exp;
  ls = -log(survival);
run;

goptions reset=all;
axis1 order=(0 to 4 by 1) minor=none label=('Exponential Reg Model Cum Hazard');
axis2 order=(0 to 4 by 1) minor=none label=( a=90 'Kaplan-Meier Cum Hazard');
symbol1 i=l1p  c= blue v=dot h=.4;
symbol2 i = join c = red l = 3;

proc gplot data=surv_exp;
  plot (ls cox)*cox / overlay haxis=axis1 vaxis= axis2;
run;
quit;

/*First we output the estimates of the cumulative distribution function using 
the cdf option in the output statement. 
In the following data step we then calculate the Cox-Snell residual. 
Finally, we use the graphics ability of proc lifetest to plot the graph 
via the plots options in the proc lifetest statement. Furthermore, 
by specifying the Cox-Snell residuals as the time variable 
in the proc lifetest model statement the procedure computes the Kaplan-Meier estimates of 
the cumulative hazard function and graphs it against the Cox-Snell residuals. 
The fitted model is correct if the Cox-Snell residual have an exponential distribution, 
i.e. if the graph is a straight line through the origin and with a slope of 1.*/

ods graphics off;

ods graphics on;
title 'Weibull Survival Model';

proc lifereg data=b1 noprint;
  model time*status(0) = treatment size number / distribution=weibull;
  output out=weibull cdf=f;
run;

data weibull1;
  set weibull;
  cox = -log( 1-f );
run;

proc lifetest data=weibull1 outsurv=surv_wei noprint;
  time cox*status(0);
run;

data surv_wei;
  set surv_wei;
  ls = -log(survival);
run;

goptions reset=all;
axis1 order=(0 to 5 by 1) minor=none label=('Weibull Reg Model Cum Hazard');
axis2 order=(0 to 5 by 1) minor=none label=( a=90 'Kaplan-Meier Cum Hazard');
symbol1 i=l1p  c= blue v=dot h=.4;
symbol2 i = join c = red l = 3;

proc gplot data=surv_wei;
  plot (ls cox)*cox / overlay haxis=axis1 vaxis= axis2;
run;
quit;

ods pdf close;



