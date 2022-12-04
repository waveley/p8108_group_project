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
proc sort data = fin_bladder; by id; run;

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
