*mediation analysis*;
%macro multiple_mediation(data=,yvar=,avar=,mvar=,cvar=,a0=,a1=,m=,nc=, yreg=,mreg=,
interaction=,casecontrol=false,output=reduced,c=,boot=);

*/data house keeping*/;
data data1;
set &data (keep=&yvar &mvar &avar &cvar);
run;
		%if (&cvar^= & &casecontrol=true) %then %do;
		%LET cvars= &cvar;
		%LET i =1;
		%DO %UNTIL(NOT %LENGTH(%SCAN(&cvars,&i))) ;
proc means noprint data=data1;
where &yvar=0;
var %SCAN(&cvars,&i);
output out=data2&i mean=/autoname;
run;
data data2&i;
set data2&i;
drop _TYPE_ _FREQ_;
run;
proc iml;
use data2&i;
read all into vb;
mean=vb[1,1];
cname1 = {"mean"};
create data2new&i from mean [colname=cname1];
append from mean;
quit;
proc append base=data3 data=data2new&i;
run;
proc sql; 
		%LET i=%EVAL(&i+1);
		%END;
proc iml;
use data3;
read all into vb;
data3=t(vb);
create data2 from data3;
append from data3;
quit;
	%end;

***************************   BOOTSTRAP PROCEDURE   ******************************************************************;

	%if (&boot^= & &boot^=false) %then %do;
	options DMSOUTSIZE=MAX;
	%if &boot=true %then %do;	
%LET n = 1000;
	%end;
	%if &boot^=true %then %do;	
%LET n = &boot;
	%end;

******************* bootstrap samples******************************;
data data1;

do sample = 1 to &n; /* To create b bootstrap replications */
do i = 1 to nobs;
indexbootstrap = round(ranuni(0) * nobs);
set data1
nobs = nobs
point = indexbootstrap;
output;
end;
end;
stop;
run;
  		%do t=1 %to &n;   
   data data1&t;                                                             
   set data1(where=(sample=&t));                              
   run;
		%end;    
	%end;

***************** regression-for bootstrap *************************;
	%if (&boot^= & &boot^=false) %then %do;
		%do t=1 %to &n;

************************************************************************************************************************;
			%if &yreg=logistic  | &yreg=loglinear  |&yreg=poisson | &yreg=negbin %then %do;
************************************************************************************************************************;

				%if &interaction=false & &cvar^= %then %do;
					%if &yreg=logistic %then %do;
proc logistic  data=data1&t descending covout
outest=out0&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  noprint;
model  &yvar=&avar &mvar  &cvar;
run;
					%end;
				%end;
			%end;

************************************************************************************************************************;
			%if  &mreg=linear & &yreg^=linear %then %do;
************************************************************************************************************************;

	    	%LET mvars= &mvar;
        	%LET i =1;
        	%DO %UNTIL(NOT %LENGTH(%SCAN(&mvars,&i))) ;
				%if &cvar^= %then %do;
					%if &casecontrol=true %then %do; 
proc reg data=data1&t covout
outest=out&i&t(drop=_model_ _type_ _name_ _depvar_  %SCAN(&mvars,&i))  noprint;
where &yvar=0;
model  %SCAN(&mvars,&i)=&avar &cvar;
run;
					%end;
				%end;
			%LET i=%EVAL(&i+1);
			%END;
			%end;
		%end;

**********************************************************************************************************************;


***************** regression-for bootstrap 	END *************************;

***************** causal effects for bootstrap  *************************;

 */create objects in which we save the bootstrap samples of causal effects*/;
proc iml;

		%if &mreg=linear & &interaction=false  %then %do;
bootsample=J(&n,3,0);
		%end;

*/compute the causal effects*/;
		%if (&mreg=linear & &interaction=false ) | (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;

			%do t=1 %to &n;
			%LET mvars= &mvar;
        	%LET i =1;
        	%DO %UNTIL(NOT %LENGTH(%SCAN(&mvars,&i))) ;
USE out&i&t;
READ ALL INTO VB;
				%if (&yreg^=linear) %then %do;
beta1_&i=VB[1,3];
				%end;
USE out0&t;
READ ALL INTO VB;
theta1=VB[1,2];
theta2_&i=VB[1,%EVAL(&i+2)];
			%LET i=%EVAL(&i+1);
			%END;

*/cde and nde*/;
				%if (&yreg^=linear & &mreg=linear ) %then %do;
bootsample[&t,1]=exp((theta1)*(&a1-&a0));

bootsample[&t,2]=0;
				%LET mvars= &mvar;
				%LET i =1;
				%DO %UNTIL(NOT %LENGTH(%SCAN(&mvars,&i))) ;
bootsample[&t,2]=((theta2_&i*beta1_&i)*(&a1-&a0))+bootsample[&t,2];
				%LET i=%EVAL(&i+1);
				%END;
bootsample[&t,2]=exp(bootsample[&t,2]);

bootsample[&t,3]=bootsample[&t,1]*bootsample[&t,2];
				%end;
			%end;

x=bootsample;
cname1 = { "boot1" "boot2" "boot3"};
create bootdata from x [colname=cname1];
append from x;
		%end;* noint;

***************** causal effects for bootstrap END  *************************;

***************** causal effects, standard errors and confidence intervals from bootstrap *************************;

*/ no interaction ;
		%if (&mreg=linear & &interaction=false )| (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;

		*effects*;
use bootdata;
read all into bootdata;
effect=J(1,3);
			%do j=1 %to 3;
			
effect[,&j]=sum((bootdata[,&j]))/&n;
		
			%end;
x=(effect);
cname1 = {"effect1" "effect2" "effect3"};
create effect from x [colname=cname1] ;
append from x;
use bootdata;
read all into bootdata;
use effect;
read all into effect;
se=J(3,1);
square=J(&n,3);

*standard errors*;
			%do j=1 %to 3;
				
				%do t=1 %to &n;
square[&t,&j]=((bootdata[&t,&j])-effect[,&j])**2;
				%end;*t loop;

				se[&j,]=(
sqrt(
sum(
(square[,&j])
)
))/sqrt(&n);

			%end;
y=se;
create se from y;
append from y;
quit;
*Percentile confidence intervals*;

			%let alphalev = .05;
			%let a1 = %sysevalf(&alphalev/2*100);
			%let a2 = %sysevalf((1 - &alphalev/2)*100);

			%do j=1 %to 3;
proc univariate data = bootdata alpha = .05;
var boot&j;
output out=pmethod&j mean = effect&j pctlpts=&a1 &a2 pctlpre = p pctlname = _cil&j _ciu&j ;
run;
			%end;

proc iml;
			%do j=1 %to 3;
use pmethod&j;
read all into vb;
cil&j=vb[1,2];
ciu&j=vb[1,3];
			%end;
cil=cil1||cil2||cil3;
ciu=ciu1||ciu2||ciu3;
x= t(cil)||t(ciu) ;
create ci from x;
append from x;
quit;
		%end;*end  noint;
	%end; * boot;

***************************   BOOTSTRAP PROCEDURE -END-  ***************************************************************;

dm "out;clear";

************* REGRESSION TO PRINT **************************************************************************************;
			
************************************************************************************************************************;
	%if &yreg=logistic  | &yreg=loglinear  |&yreg=poisson | &yreg=negbin %then %do;
************************************************************************************************************************;

		%if &interaction=false & &cvar^= %then %do;
			%if &yreg=logistic %then %do;
proc logistic  data=data1 descending covout
  outest=out0(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
  model  &yvar=&avar &mvar  &cvar;
run;
			%end;
		%end;
	%end;

************************************************************************************************************************;
	%if  &mreg=linear & &yreg^=linear %then %do;
************************************************************************************************************************;

	%LET mvars= &mvar;
    %LET i =1;
    %DO %UNTIL(NOT %LENGTH(%SCAN(&mvars,&i))) ;
			%if &casecontrol=true %then %do; 
proc reg data=data1 covout
  outest=out&i(drop=_model_ _type_ _name_ _depvar_  %SCAN(&mvars,&i))  ;
  where &yvar=0;
  model  %SCAN(&mvars,&i)=&avar &cvar;
  proc print;
run;
			%end;
	%LET i=%EVAL(&i+1);
	%END;
	%end;

**********************************************************************************************************************;

***************************   OUTPUT    *************************** ;

*/LINEAR LINEAR no interaction ;
	%if (&mreg=linear & &interaction=false ) | (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;
proc iml;
use effect;
read all into effect;
use se;
read all into se;
		%if (&boot^= & &boot^=false) %then %do;
use ci;
read all into ci;

x= t(effect)||(se)||ci ;
cname1 = { "Estimate"  "s.e." "95% CI lower" "95% CI upper" };
		%END;
create x3 from x [ colname=cname1 ];
append from x;
name='cde=nde' || 'nie'||'total effect' ;
name=t(name);
cname2= {"Effect"};
create x4 from name [colname=cname2];
append from name;
quit;
	%end;

DATA x5; 
MERGE x4 x3; 
drop id;
proc print data=x5;
quit;

***************************   OUTPUT -END-   *************************** ;

data results;
set x5;
run;

proc datasets library=work;
delete data1 data2 data3 out0 out1  x3 x4 x5;
run;

%if &cvar^= %then %do;
 %do i = 1 %to &nc;
 proc datasets library=work;
delete data2&i data2new&i;
run;
 %end;
%end;

%mend;
