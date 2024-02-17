/**********************************************DGP 4 : Linear dependence*******************************/
proc iml; 

/************** AICC ***************************/

N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

n_overfitting_elast = 0;
n_underfitting_elast = 0;
n_fitting_elast = 0;
n_failure_elast = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_step;
        model y=X1-X50 / selection=stepwise (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_elast;
        model y=X1-X50 / selection=ELASTICNET (choose=AICC stop=AICC);
        run;
		proc transpose data=DGP_results_elast out=DGP_results_elast_trans;
		run;
    endsubmit;

/********************************************************************************************************************/

use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set); /* Calculate the cardinality  of 'set' */
cardinal_inter_forward = ncol(inter); /* Calculate the cardinality  of 'inter' */
cardinal_true_forward = nrow(true); /* Calculate the cardinality  of  'true' */

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

use DGP_results_elast_trans; 
read all; 
close DGP_results_elast_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_elast = nrow(set); 
cardinal_inter_elast = ncol(inter); 
cardinal_true_elast = nrow(true); 

*Conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

	if cardinal_set_elast = cardinal_inter_elast & cardinal_inter_elast < cardinal_true_elast then n_underfitting_elast = n_underfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast < cardinal_true_elast & cardinal_inter_elast < cardinal_set_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_fitting_elast = n_fitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_overfitting_elast = n_overfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;

end;

/*calculate the probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

    proba_overfitting_elast = (n_overfitting_elast / N_loop) * 100;
    proba_underfitting_elast = (n_underfitting_elast / N_loop) * 100;
    proba_fitting_elast = (n_fitting_elast / N_loop) * 100;
    proba_failure_elast = (n_failure_elast / N_loop) * 100;

*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;

ELASTICNET = proba_overfitting_elast || proba_underfitting_elast || proba_fitting_elast || proba_failure_elast;
create table_elast from ELASTICNET[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from ELASTICNET;
close table_elast;

submit;

*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso table_elast;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO'
		Col6 = 'ELASTICNET';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=AICC choose=AICC ';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


	/************** BIC ***************************/
	
N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

n_overfitting_elast = 0;
n_underfitting_elast = 0;
n_fitting_elast = 0;
n_failure_elast = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;
 
    submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (stop = BIC choose=BIC);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (stop = BIC choose=BIC);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_step ;
        model y=X1-X50 / selection=stepwise (stop = BIC choose=BIC);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
  	
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (stop=BIC choose=BIC);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (stop = BIC choose=BIC);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;
  
        proc glmselect noprint data=myDataset outdesign = DGP_results_elast;
        model y=X1-X50 / selection=ELASTICNET (stop = BIC choose=BIC);
        run;
		proc transpose data=DGP_results_elast out=DGP_results_elast_trans;
		run;
    endsubmit;
    
/********************************************************************************************************************/

use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set);
cardinal_inter_forward = ncol(inter);
cardinal_true_forward = nrow(true);

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

use DGP_results_elast_trans; 
read all; 
close DGP_results_elast_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_elast = nrow(set); 
cardinal_inter_elast = ncol(inter); 
cardinal_true_elast = nrow(true); 

*Construct conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

	if cardinal_set_elast = cardinal_inter_elast & cardinal_inter_elast < cardinal_true_elast then n_underfitting_elast = n_underfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast < cardinal_true_elast & cardinal_inter_elast < cardinal_set_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_fitting_elast = n_fitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_overfitting_elast = n_overfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;

end;

/*Calculate probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

    proba_overfitting_elast = (n_overfitting_elast / N_loop) * 100;
    proba_underfitting_elast = (n_underfitting_elast / N_loop) * 100;
    proba_fitting_elast = (n_fitting_elast / N_loop) * 100;
    proba_failure_elast = (n_failure_elast / N_loop) * 100;


*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;


ELASTICNET = proba_overfitting_elast || proba_underfitting_elast || proba_fitting_elast || proba_failure_elast;
create table_elast from ELASTICNET[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from ELASTICNET;
close table_elast;

submit;
*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso table_elast;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO'
		Col6 = 'ELASTICNET';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=BIC choose=BIC';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


   /************** Adjusted R-squared ***************************/
  
N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

n_overfitting_elast = 0;
n_underfitting_elast = 0;
n_fitting_elast = 0;
n_failure_elast = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

    submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_step;
        model y=X1-X50 / selection=stepwise (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;

        proc glmselect noprint data=myDataset outdesign = DGP_results_elast;
        model y=X1-X50 / selection=ELASTICNET (stop = ADJRSQ choose=ADJRSQ);
        run;
		proc transpose data=DGP_results_elast out=DGP_results_elast_trans;
		run;
    endsubmit;

/********************************************************************************************************************/
use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set);
cardinal_inter_forward = ncol(inter);
cardinal_true_forward = nrow(true);

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

use DGP_results_elast_trans; 
read all; 
close DGP_results_elast_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_elast = nrow(set); 
cardinal_inter_elast = ncol(inter); 
cardinal_true_elast = nrow(true); 


*Construct conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

	if cardinal_set_elast = cardinal_inter_elast & cardinal_inter_elast < cardinal_true_elast then n_underfitting_elast = n_underfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast < cardinal_true_elast & cardinal_inter_elast < cardinal_set_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_fitting_elast = n_fitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_overfitting_elast = n_overfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;

end;

/*Calculate probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

    proba_overfitting_elast = (n_overfitting_elast / N_loop) * 100;
    proba_underfitting_elast = (n_underfitting_elast / N_loop) * 100;
    proba_fitting_elast = (n_fitting_elast / N_loop) * 100;
    proba_failure_elast = (n_failure_elast / N_loop) * 100;


*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;


ELASTICNET = proba_overfitting_elast || proba_underfitting_elast || proba_fitting_elast || proba_failure_elast;
create table_elast from ELASTICNET[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from ELASTICNET;
close table_elast;

submit;
*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso table_elast;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO'
		Col6 = 'ELASTICNET';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=ADJRSQ choose=ADJRSQ';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;

     /************** CP ***************************/
    
N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

n_overfitting_elast = 0;
n_underfitting_elast = 0;
n_fitting_elast = 0;
n_failure_elast = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

	submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_step;
        model y=X1-X50 / selection=stepwise (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_elast;
        model y=X1-X50 / selection=ELASTICNET (choose=CP stop=CP);
        run;
		proc transpose data=DGP_results_elast out=DGP_results_elast_trans;
		run;
    endsubmit;


use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set); /* Calculate the cardinality  of 'set' */
cardinal_inter_forward = ncol(inter); /* Calculate the cardinality  of 'inter' */
cardinal_true_forward = nrow(true); /* Calculate the cardinality  of  'true' */

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

use DGP_results_elast_trans; 
read all; 
close DGP_results_elast_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_elast = nrow(set); 
cardinal_inter_elast = ncol(inter); 
cardinal_true_elast = nrow(true); 

*Conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

	if cardinal_set_elast = cardinal_inter_elast & cardinal_inter_elast < cardinal_true_elast then n_underfitting_elast = n_underfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast < cardinal_true_elast & cardinal_inter_elast < cardinal_set_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_fitting_elast = n_fitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_overfitting_elast = n_overfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;

end;

/*calculate the probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

    proba_overfitting_elast = (n_overfitting_elast / N_loop) * 100;
    proba_underfitting_elast = (n_underfitting_elast / N_loop) * 100;
    proba_fitting_elast = (n_fitting_elast / N_loop) * 100;
    proba_failure_elast = (n_failure_elast / N_loop) * 100;

*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;

ELASTICNET = proba_overfitting_elast || proba_underfitting_elast || proba_fitting_elast || proba_failure_elast;
create table_elast from ELASTICNET[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from ELASTICNET;
close table_elast;

submit;

*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso table_elast;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO'
		Col6 = 'ELASTICNET';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=CP choose=CP ';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


     /************** Kfold ***************************/
    
N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

n_overfitting_elast = 0;
n_underfitting_elast = 0;
n_fitting_elast = 0;
n_failure_elast = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

    submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  	
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_step;
        model y=X1-X50 / selection=stepwise (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_elast;
        model y=X1-X50 / selection=ELASTICNET (stop = CV choose=CV);
        run;
		proc transpose data=DGP_results_elast out=DGP_results_elast_trans;
		run;
    endsubmit;

use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set);
cardinal_inter_forward = ncol(inter);
cardinal_true_forward = nrow(true);

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

use DGP_results_elast_trans; 
read all; 
close DGP_results_elast_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_elast = nrow(set); 
cardinal_inter_elast = ncol(inter); 
cardinal_true_elast = nrow(true); 


*Construct conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

	if cardinal_set_elast = cardinal_inter_elast & cardinal_inter_elast < cardinal_true_elast then n_underfitting_elast = n_underfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast < cardinal_true_elast & cardinal_inter_elast < cardinal_set_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_fitting_elast = n_fitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast = cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast = cardinal_true_elast then n_overfitting_elast = n_overfitting_elast + 1;
	else tt = 1;
	if cardinal_set_elast > cardinal_true_elast & cardinal_inter_elast < cardinal_true_elast then n_failure_elast = n_failure_elast + 1;
	else tt = 1;

end;

/*Calculate probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

    proba_overfitting_elast = (n_overfitting_elast / N_loop) * 100;
    proba_underfitting_elast = (n_underfitting_elast / N_loop) * 100;
    proba_fitting_elast = (n_fitting_elast / N_loop) * 100;
    proba_failure_elast = (n_failure_elast / N_loop) * 100;


*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;


ELASTICNET = proba_overfitting_elast || proba_underfitting_elast || proba_fitting_elast || proba_failure_elast;
create table_elast from ELASTICNET[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from ELASTICNET;
close table_elast;

submit;
*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso table_elast;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO'
		Col6 = 'ELASTICNET';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=CV choose=CV';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


	/************** Leave one out ***************************/
	/* The STOP=PRESS and CHOOSE=PRESS suboptions cannot be used with SELECTION=LAR, SELECTION=LASSO or SELECTION=ELASTICNET
       unless the LSCOEFFS suboption is also specified*/
      
N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;
n_overfitting_forward = 0;
n_underfitting_forward = 0;
n_fitting_forward = 0;
n_failure_forward = 0;

n_overfitting_backward = 0;
n_underfitting_backward = 0;
n_fitting_backward = 0;
n_failure_backward = 0;

n_overfitting_stepwise = 0;
n_underfitting_stepwise = 0;
n_fitting_stepwise = 0;
n_failure_stepwise = 0;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;


do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

        submit;
        proc glmselect noprint data=myDataset outdesign = DGP_results_forward ;
        model y=X1-X50 / selection=forward (stop=PRESS choose=PRESS);
        run;
		proc transpose data=DGP_results_forward out=DGP_results_for_trans;
		run;
  	
        proc glmselect noprint data=myDataset outdesign = DGP_results_backward;
        model y=X1-X50 / selection=backward (stop=PRESS choose=PRESS);
        run;
		proc transpose data=DGP_results_backward out=DGP_results_back_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_step;
        model y=X1-X50 / selection=stepwise (stop=PRESS choose=PRESS);
        run;
		proc transpose data=DGP_results_step out=DGP_results_step_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (choose=PRESS LSCOEFFS stop=PRESS LSCOEFFS);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
   
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (choose=PRESS LSCOEFFS stop=PRESS LSCOEFFS);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;
    endsubmit;
	
	/*The CHOOSE=PRESS or PRESS LSCOEFFS option is not available for SELECTION=ELASTICNET*/

use DGP_results_for_trans; 
read all; 
close DGP_results_for_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_forward = nrow(set);
cardinal_inter_forward = ncol(inter);
cardinal_true_forward = nrow(true);

use DGP_results_back_trans; 
read all; 
close DGP_results_back_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_backward = nrow(set);
cardinal_inter_backward = ncol(inter);
cardinal_true_backward = nrow(true);

use DGP_results_step_trans; 
read all; 
close DGP_results_step_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_stepwise = nrow(set);
cardinal_inter_stepwise = ncol(inter);
cardinal_true_stepwise = nrow(true);

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);


*Construct conditions;
	if cardinal_set_forward = cardinal_inter_forward & cardinal_inter_forward < cardinal_true_forward then n_underfitting_forward = n_underfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward < cardinal_true_forward & cardinal_inter_forward < cardinal_set_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_fitting_forward = n_fitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward = cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward = cardinal_true_forward then n_overfitting_forward = n_overfitting_forward + 1;
	else tt = 1;
	if cardinal_set_forward > cardinal_true_forward & cardinal_inter_forward < cardinal_true_forward then n_failure_forward = n_failure_forward + 1;
	else tt = 1;

	if cardinal_set_backward = cardinal_inter_backward & cardinal_inter_backward < cardinal_true_backward then n_underfitting_backward = n_underfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward < cardinal_true_backward & cardinal_inter_backward < cardinal_set_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward = cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_fitting_backward = n_fitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward  = cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward = cardinal_true_backward then n_overfitting_backward = n_overfitting_backward + 1;
	else tt = 1;
	if cardinal_set_backward > cardinal_true_backward & cardinal_inter_backward < cardinal_true_backward then n_failure_backward = n_failure_backward + 1;
	else tt = 1;

	if cardinal_set_stepwise = cardinal_inter_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_underfitting_stepwise = n_underfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise < cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_set_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_fitting_stepwise = n_fitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise = cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise = cardinal_true_stepwise then n_overfitting_stepwise = n_overfitting_stepwise + 1;
	else tt = 1;
	if cardinal_set_stepwise > cardinal_true_stepwise & cardinal_inter_stepwise < cardinal_true_stepwise then n_failure_stepwise = n_failure_stepwise + 1;
	else tt = 1;

	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

end;

/*Calculate probabilities*/
    proba_overfitting_forward = (n_overfitting_forward / N_loop) * 100;
    proba_underfitting_forward = (n_underfitting_forward / N_loop) * 100;
    proba_fitting_forward = (n_fitting_forward / N_loop) * 100;
    proba_failure_forward = (n_failure_forward / N_loop) * 100;

    proba_overfitting_backward = (n_overfitting_backward / N_loop) * 100;
    proba_underfitting_backward = (n_underfitting_backward / N_loop) * 100;
    proba_fitting_backward = (n_fitting_backward / N_loop) * 100;
    proba_failure_backward = (n_failure_backward / N_loop) * 100;

	proba_overfitting_stepwise = (n_overfitting_stepwise / N_loop) * 100;
    proba_underfitting_stepwise = (n_underfitting_stepwise / N_loop) * 100;
    proba_fitting_stepwise = (n_fitting_stepwise / N_loop) * 100;
    proba_failure_stepwise = (n_failure_stepwise / N_loop) * 100;

    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

*Create table for each method;
Forward = proba_overfitting_forward || proba_underfitting_forward|| proba_fitting_forward || proba_failure_forward;
create table_for from Forward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Forward;
close table_for;

Backward = proba_overfitting_backward || proba_underfitting_backward || proba_fitting_backward || proba_failure_backward;
create table_back from Backward[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Backward;
close table_back;

Stepwise = proba_overfitting_stepwise || proba_underfitting_stepwise || proba_fitting_stepwise || proba_failure_stepwise;
create table_step from Stepwise[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from Stepwise;
close table_step;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;


submit;
*merge all the tables;
data MergedTable;
  set table_for table_back table_step table_lar table_lasso;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label Col1 = 'Forward'
		Col2 = 'Backward'
		Col3 = 'Stepwise'
		Col4 = 'LAR'
		Col5 = 'LASSO';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=PRESS (LSCOEFFS) choose=PRESS (LSCOEFFS)';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


/************************ Differents choose and stop criterion for machine learning methods *******************************/

/************** Choose=SBC  stop=PRESS LSCOEFFS ***************************/

N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

	submit;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (stop = PRESS LSCOEFFS choose=SBC);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (stop = PRESS LSCOEFFS choose=SBC);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;

    endsubmit;
    
/*************************************************************************/

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

*Construct conditions;
	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	
end;

/*Calculate probabilities*/
  
    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

*Create table for each method;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;

submit;
*merge all the tables;
data MergedTable;
  set table_lar table_lasso;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label	Col1 = 'LAR'
		Col2 = 'LASSO';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=PRESS LSCOEFFS choose=SBC';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;


/************** Choose=PRESS LSCOEFFS  Stop=SBC ***************************/

N_loop = 1000;
N = 500;
nb_var = 50;

*initialization to zero;

n_overfitting_lar = 0;
n_underfitting_lar = 0;
n_fitting_lar = 0;
n_failure_lar = 0;

n_overfitting_lasso = 0;
n_underfitting_lasso = 0;
n_fitting_lasso = 0;
n_failure_lasso = 0;

do i = 1 to N_loop; 
*Mean vector of 5 dependent variables;
mean_5 = j(5,1,0);

*Define a variance-covariance matrix;
cov = toeplitz ({1 0.3 0.7 0.3 0.9});

x_dep = randnormal(N, mean_5, cov);

*Mean vector of 45 independent variables;
mean_45 = j(45,1,0);

x_ind = randnormal(N, mean_45, I(45));

X = x_dep||x_ind;

* Generate random coefficients for first 5 variables et error terms;
beta = {0.4 1.2 0.9 0.2 1.6};
eps = normal(j(nrow(X),1,1))*0.1; 

* Generate response variable with a linear relationship;
y = x[,1:5] * beta` + eps;

* Combine predictors and response into a matrix;
dataMatrix = x || y;
colNames = "X1":"X50" || "y"; /*Create the column names*/

* Ceate a dataset for DPG4;
    create myDataset from datamatrix[colname=colnames];
	append from dataMatrix;
	close myDataset;

	submit;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lar;
        model y=X1-X50 /  selection= LAR (stop = SBC choose=PRESS LSCOEFFS);
        run;
		proc transpose data=DGP_results_lar out=DGP_results_lar_trans;
		run;
    
        proc glmselect noprint data=myDataset outdesign = DGP_results_lasso;
        model y=X1-X50 / selection=LASSO (stop = SBC choose=PRESS LSCOEFFS);
        run;
		proc transpose data=DGP_results_lasso out=DGP_results_lasso_trans;
		run;

    endsubmit;
	/*************************************************************************/

use DGP_results_lar_trans; 
read all; 
close DGP_results_lar_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lar = nrow(set);
cardinal_inter_lar = ncol(inter);
cardinal_true_lar = nrow(true);

use DGP_results_lasso_trans; 
read all; 
close DGP_results_lasso_trans;

set=_LABEL_;
set=set[2:nrow(set)-1];
true=("X1":"X5")`;
inter=Xsect(set,true);
cardinal_set_lasso = nrow(set);
cardinal_inter_lasso = ncol(inter);
cardinal_true_lasso = nrow(true);

*Construct conditions;
	if cardinal_set_lar = cardinal_inter_lar & cardinal_inter_lar < cardinal_true_lar then n_underfitting_lar = n_underfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar < cardinal_true_lar & cardinal_inter_lar < cardinal_set_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_fitting_lar = n_fitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar = cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar = cardinal_true_lar then n_overfitting_lar = n_overfitting_lar + 1;
	else tt = 1;
	if cardinal_set_lar > cardinal_true_lar & cardinal_inter_lar < cardinal_true_lar then n_failure_lar = n_failure_lar + 1;
	else tt = 1;

	if cardinal_set_lasso = cardinal_inter_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_underfitting_lasso = n_underfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso < cardinal_true_lasso & cardinal_inter_lasso < cardinal_set_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_fitting_lasso = n_fitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso = cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso = cardinal_true_lasso then n_overfitting_lasso = n_overfitting_lasso + 1;
	else tt = 1;
	if cardinal_set_lasso > cardinal_true_lasso & cardinal_inter_lasso < cardinal_true_lasso then n_failure_lasso = n_failure_lasso + 1;
	else tt = 1;

end;

/*Calculate probabilities*/
  
    proba_overfitting_lar = (n_overfitting_lar / N_loop) * 100;
    proba_underfitting_lar = (n_underfitting_lar / N_loop) * 100;
    proba_fitting_lar = (n_fitting_lar / N_loop) * 100;
    proba_failure_lar = (n_failure_lar / N_loop) * 100;

    proba_overfitting_lasso = (n_overfitting_lasso / N_loop) * 100;
    proba_underfitting_lasso = (n_underfitting_lasso / N_loop) * 100;
    proba_fitting_lasso = (n_fitting_lasso / N_loop) * 100;
    proba_failure_lasso = (n_failure_lasso / N_loop) * 100;

*Create table for each method;

LAR = proba_overfitting_lar || proba_underfitting_lar || proba_fitting_lar || proba_failure_lar;
create table_lar from LAR[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LAR;
close table_lar;

LASSO = proba_overfitting_lasso || proba_underfitting_lasso || proba_fitting_lasso || proba_failure_lasso ;
create table_lasso from LASSO[colname={'Overfitting','Underfitting','Fitting','Failure'}];
append from LASSO;
close table_lasso;

submit;
*merge all the tables;
data MergedTable;
  set table_lar table_lasso;
run;
proc transpose data=MergedTable out=MergedTable_trans;
		run;
*change the name of all columns;
data MergedTable_trans;
  set MergedTable_trans;
  label	Col1 = 'LAR'
		Col2 = 'LASSO';
run;
proc transpose data=MergedTable_trans out=MergedTable_final;
		run;

*Graphic;
proc template;
define statgraph CombinedGraph;
dynamic __LABEL_ _FITTING _OVERFITTING _UNDERFITTING _FAILURE __LABEL_2;
begingraph / designwidth=1616 designheight=960;
entrytitle halign=center 'stop=SBC choose=PRESS LSCOEFFS';
   layout lattice / rows=2 columns=2 rowgutter=10 columngutter=10;
      layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FITTING / group=__LABEL_2 name='bar_h' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
         discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_OVERFITTING / group=__LABEL_2 name='bar_h2' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_UNDERFITTING / group=__LABEL_2 name='bar_h3' orient=horizontal barwidth=0.85 groupdisplay=Cluster clusterwidth=0.5 grouporder=data BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;

       layout overlay / xaxisopts=( linearopts=( viewmin=1.0 viewmax=100.0)) yaxisopts=(label=('Method') discreteopts=(tickvaluefitpolicy=none));
         barchart category=__LABEL_ response=_FAILURE / group=__LABEL_2 name='bar_h4' orient=horizontal barwidth=1.0 groupdisplay=Cluster clusterwidth=0.85 BARLABEL=TRUE;
		 discretelegend 'bar_h' / opaque=false border=true halign=right valign=top displayclipped=true across=1 order=rowmajor location=inside;
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.MERGEDTABLE_FINAL template=CombinedGraph;
dynamic __LABEL_="'_LABEL_'n" _FITTING="FITTING" _OVERFITTING="OVERFITTING" _UNDERFITTING="UNDERFITTING" _FAILURE="FAILURE" __LABEL_2="'_LABEL_'n";
run;

endsubmit;

