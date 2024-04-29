/* Apply a variety of treatment strategies to patients with selected parameter configurations.
   Report the survival times for all strategies.
   Combine optimize_two_drug_responses_trimrates_strategies_index14.c and optimize_three_drug_responses_trimrates_strategies_index4.c.
   Difference from these two files: (1)Operate on multiple drugs, (2)Alter strategy 0, (3)Report the survival times of the best static strategies and the unrealistic settings where all drugs are administered with full dosage.
   Difference from optimize_multidrug_responses_trimrates_strategies_index1.c: fix the bug of converting arguments into drug sensitivities and transition matrices.
   Difference from optimize_multidrug_responses_trimrates_strategies_index1_5.c: correct the error of implementing strategies 2.
   Difference from optimize_multidrug_responses_trimrates_strategies_index1_8.c: (1)Report the parameter configurations, (2)Report the dosage combination sequences over time, (3)Report the population composition sequences over time.
*/


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include "db.h"
# define MAXLINEWIDTH 100000
# define NAMELENGTH 100
# define max(i,j) (i>j ? i : j)
# define min(i,j) (i<j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)
  

extern char * getitemval4(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern double optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);   
extern double dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
extern double dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
extern double dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
extern void num2vec(long index, int ndigits, int *bases, int *numvec);
extern long vec2num(int ndigits, int *bases, int *numvec);
extern long vec2num2(int ndigits, int *bases, int *numvec);
extern double deathtime_multidrug(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates);
extern double optimize_two_drug_responses8(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
extern double global_dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int *truncated);
extern double global_dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);
extern double global_dp_optimize_multidrug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);
extern void recurse_multidrug_response_trimrates2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode);
extern double optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
extern double correct_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
extern double correct_dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);


int nstrategies=3, ndrugs=2, ntypes=4, nvars=17, nsteps=5;
int ntrials=7, totallinewidth=100;
double detected=1e9, mortal=1e13;
int nchar=4, njointchar=16, nnodes=0, initialtime=3;
double epsilon=0.1;
//int nsteps=6, nstates=11, degmode=0;
//int nsteps=6, nstates=21;
int maxnfrontiers=10000;
double timeinterval=45, predictionperiod=45, N0=5e9;
//double mint=0, maxt=7200;   
double mint=0, maxt=1800; 
//int maxncandseqs=500;
int maxncandseqs=50;
double negligible=10, ratiothre=100;

int defaultmaxncandseqs=50, defaultevamode=1;
int highmaxncandseqs=500, highevamode=0;
////int highmaxncandseqs=1000, highevamode=0;
int curevamode=1;
//int highmaxncandseqs=50, highevamode=0;
double bsearchthre=1e-5;

int main(int argc, char *argv[])
{
  // *(argv+1) -- number of drugs.
  // *(argv+2) -- information about parameter interpretations, numbers, configuration values and additional constraints.
  // *(argv+3) -- index of the parameter settings.
  // *(argv+4) -- increment size.
  // *(argv+5) -- output file prefix of the parameter configurations.
  // *(argv+6) -- output file prefix of the stopping times.
  // *(argv+7) -- output file prefix of the dosage combination sequences.
  // *(argv+8) -- output file prefix of the population combination sequences.
    
  int i, j, k, l, m, n, length, nintervals, firstwrite=0, truncated=0;
  int *numvec, *bases, flag, degmode, *binarybases, *resvec, *resvec2, beststatic=-1;
  long index, startindex, endindex, increment, maxnum, maxindex, allowedcnt=0;
  struct variable *vars;
  time_t *t1, *t2;
  double *x0, *T, *g0, *a0, *Sa, *Sg, *t, val=0;
  double *y, *dosage, *dosages, *pops, *stopTs, popthre;
  char ch, *s, *item, *str, *tmpstr, *prefix, *filename;
  char *paramfilename, *stoptimefilename, *dosagefilename, *popfilename;
  FILE *fp, *fout, *fout1, *fout2, *fout3, *fout4;
  double *dstates, lwd, upp;
  int nstates, combind;

  t1=malloc(sizeof(time_t)); t1=time(t1);

  // Load number of drugs.
  ndrugs=atoi(*(argv+1)); 
  ntypes=1;
  for (i=1; i<=ndrugs; i++)
    ntypes*=2;

  // Two drug case: bsearchthre=1e-5.
  // Three drug case: bsearchthre=1e-1.
  //if (ndrugs==2)
  //bsearchthre=1e-5;
  //else
  //bsearchthre=1e-1;
  
  bsearchthre=1e-1;


  // Load variable information.
  if ((fp=fopen(*(argv+2),"r"))==NULL)
    {
      printf("Cannot open %s.\n",*(argv+2));
      return -1;
    } 
  s=malloc(sizeof(char)*MAXLINEWIDTH); fgets(s,MAXLINEWIDTH,fp); 
  nvars=0; vars=malloc(sizeof(struct variable));
  while (fgets(s,MAXLINEWIDTH,fp)!=NULL)
    {
      nvars++; vars=realloc(vars,sizeof(struct variable)*nvars);
      item=getitemval4(s,0,&length,0x09);
      (vars+nvars-1)->name=malloc(sizeof(char)*LONGNAMELENGTH);
      strncpy((vars+nvars-1)->name,item,length); *((vars+nvars-1)->name+length)='\0';
      free(item);
      item=getitemval4(s,1,&length,0x09);
      val=atof2(item,length); free(item);
      (vars+nvars-1)->nvalues=(int)val;
      (vars+nvars-1)->values=malloc(sizeof(double)*(vars+nvars-1)->nvalues);
      for (i=0; i<=((vars+nvars-1)->nvalues-1); i++)
	{
	  item=getitemval4(s,i+2,&length,0x09);
	  val=atof2(item,length); free(item);
	  *((vars+nvars-1)->values+i)=val;
	}
    }
  fclose(fp);
  
  // Load the input index and increment size.
  startindex=atoi(*(argv+3)); increment=atoi(*(argv+4));

  // Initialize parameters.  
  numvec=malloc(sizeof(int)*nvars);
  x0=malloc(sizeof(double)*ntypes);
  g0=malloc(sizeof(double)*ntypes);
  a0=malloc(sizeof(double)*ntypes);
  Sg=malloc(sizeof(double)*ntypes*ndrugs);
  Sa=malloc(sizeof(double)*ntypes*ndrugs);
  T=malloc(sizeof(double)*ntypes*ntypes);  
  //mint=0; maxt=1800; 
  nintervals=ceil((maxt-mint)/timeinterval)+1;
  y=malloc(sizeof(double)*ntypes*nintervals);  
  dosage=malloc(sizeof(double)*ndrugs*nintervals);
  stopTs=malloc(sizeof(double)*nstrategies);
  t=malloc(sizeof(double)*nintervals); *t=mint;
  for (i=1; i<=(nintervals-1); i++)
    *(t+i)=*(t+i-1)+timeinterval;  
  
  dosages=malloc(sizeof(double)*nstrategies*nintervals*ndrugs);
  pops=malloc(sizeof(double)*nstrategies*nintervals*ntypes);
  str=malloc(sizeof(char)*MAXLINEWIDTH);
  tmpstr=malloc(sizeof(char)*MAXLINEWIDTH);
  bases=malloc(sizeof(int)*nvars);
  for (n=0; n<=(nvars-1); n++)
    *(bases+n)=(vars+n)->nvalues;
    
  for (i=0; i<=(nvars-1); i++)
    *(numvec+i)=*(bases+i)-1;
  maxindex=vec2num(nvars,bases,numvec);
  
  
  binarybases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(binarybases+n)=2;
  resvec=malloc(sizeof(int)*ndrugs);
  resvec2=malloc(sizeof(int)*ndrugs);


  maxnum=1;
  for (n=0; n<=(nvars-1); n++)
    maxnum*=(vars+n)->nvalues;
  prefix=malloc(sizeof(char)*MAXLINEWIDTH);

  paramfilename=malloc(sizeof(char)*MAXLINEWIDTH);
  strncpy(prefix,*(argv+5),strlen(*(argv+5))); *(prefix+strlen(*(argv+5)))='\0';  
  sprintf(paramfilename,"%s%ld.txt\0",prefix,startindex);    
  stoptimefilename=malloc(sizeof(char)*MAXLINEWIDTH);
  strncpy(prefix,*(argv+6),strlen(*(argv+6))); *(prefix+strlen(*(argv+6)))='\0';  
  sprintf(stoptimefilename,"%s%ld.txt\0",prefix,startindex);  
  dosagefilename=malloc(sizeof(char)*MAXLINEWIDTH);
  strncpy(prefix,*(argv+7),strlen(*(argv+7))); *(prefix+strlen(*(argv+7)))='\0';  
  sprintf(dosagefilename,"%s%ld.txt\0",prefix,startindex);  
  popfilename=malloc(sizeof(char)*MAXLINEWIDTH);
  strncpy(prefix,*(argv+8),strlen(*(argv+8))); *(prefix+strlen(*(argv+8)))='\0';  
  sprintf(popfilename,"%s%ld.txt\0",prefix,startindex);  

  // Debug
  //startindex=startindex; endindex=increment;

  
  startindex=startindex*increment; endindex=min(maxindex,startindex+increment-1);


  // Enumerate all drug combinations.
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,binarybases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }

  // Debug
  //nstates=2; dstates=malloc(sizeof(double)*nstates*ndrugs);
  //*(dstates+0)=0; *(dstates+1)=1;
  //*(dstates+2)=1; *(dstates+3)=0;


  // Debug
  //startindex=503000; endindex=503000;


  // Debug
  //printf("startindex=%ld, endindex=%ld\n",startindex,endindex);
  
  // Calculate survival times for all strategies on selected parameter configurations.
  for (index=startindex; index<=endindex; index++)
    {
      int allowed=1;
      double *vals;
      long indval;

      /*
      // Debug
      for (i=0; i<=(nvars-1); i++)
	*(numvec+i)=0;
      *(numvec+0)=2; 
      *(numvec+1)=1;
      *(numvec+2)=2; 
      *(numvec+3)=2;
      *(numvec+4)=2;
      *(numvec+5)=0; 
      *(numvec+6)=0;
      *(numvec+7)=0;      
      *(numvec+8)=1;
      *(numvec+9)=1; 
      *(numvec+10)=1;
      *(numvec+11)=0; 
      *(numvec+12)=0; 
      *(numvec+13)=0; 
      *(numvec+14)=2; 
      *(numvec+15)=2; 
      *(numvec+16)=2; 
      indval=vec2num2(nvars,bases,numvec);
      printf("indval=%ld\n",indval);      
      exit(0);
      */

            
      // Convert the index into the parameter configuration.
      num2vec(index,nvars,bases,numvec);
      
      /*
      // Debug
      *(numvec+0)=1; *(numvec+1)=1; *(numvec+2)=1;
      for (i=3; i<=7; i++)
	*(numvec+i)=0;
      for (i=8; i<=13; i++)
	*(numvec+i)=1;
      *(numvec+11)=0;
      //*(numvec+8)=1; *(numvec+9)=2; *(numvec+10)=1;
      //*(numvec+11)=2; *(numvec+12)=1; *(numvec+13)=0;      
      for (i=13; i<=16; i++)
	*(numvec+i)=2;
      index=vec2num(nvars,bases,numvec);
      printf("index=%ld\n",index);
      */

      /*
      // Debug
      printf("numvec= ");
      for (i=0; i<=(nvars-1); i++)
	printf("%d ",*(numvec+i));
      printf("\n");
      exit(0);
      */

      // Convert encoded parameter values into parameter values in the model.
      for (i=0; i<=(ntypes-1); i++)
	{
	  *(x0+i)=0; *(g0+i)=0; *(a0+i)=0;
	  for (j=0; j<=(ntypes-1); j++)
	    *(T+i*ntypes+j)=0;
	  for (j=0; j<=(ndrugs-1); j++)
	    {
	      *(Sa+i*ndrugs+j)=0; *(Sg+i*ndrugs+j)=0; 
	    }
	}

      vals=malloc(sizeof(double)*3);

      for (n=0; n<=(nvars-1); n++)
	{
	  for (i=0; i<=(ndrugs-1); i++)
	    *(resvec+i)=0;
	  k=*(numvec+n); val=*((vars+n)->values+k);
	  
	  // Natural growth rate g0.
	  if (strcmp((vars+n)->name,"g0")==0)
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*(g0+i)=val;
	    }

	  // Initial resistat population x0_Rx.
	  // Obtain the drugs where the cells are resistant.
	  else if (strncmp((vars+n)->name,"x0_R",4)==0)
	    {
	      int starti=4, endi=-1;
	      i=starti;
	      while ((i<strlen((vars+n)->name))&&(endi<0))
		{
		  if (((ch=*((vars+n)->name+i))>='0')&&(ch<='9'))
		    {
		      j=ch-'0'; *(resvec+j-1)=1;
		    }
		  else
		    endi=i;
		  i++;
		}
	      k=vec2num(ndrugs,binarybases,resvec);
	      *(x0+k)=val*N0;
	    }

	  // Drug sensitivity Sa(t1,t2)/g0 or Sa(t1,t2)/Sa(u1,u2).
	  else if (strncmp((vars+n)->name,"Sa(",3)==0)
	    {
	      int typeind1=-1, typeind2=-1, drugind1=-1, drugind2=-1;
	      int starti=3, endi=-1;
	      double refval;
	      i=starti;
	      while ((i<strlen((vars+n)->name))&&(endi<0))
		{
		  if ((ch=*((vars+n)->name+i))==',')
		    endi=i-1;
		  i++;
		}
	      if (strncmp((vars+n)->name+starti,"S",1)==0)
		typeind1=0;
	      else
		{
		  ch=*((vars+n)->name+starti+1);
		  j=ch-'0'; typeind1=j; 
		}
	      starti=endi+2; endi=-1; i=starti;
	      while ((i<strlen((vars+n)->name))&&(endi<0))
		{
		  if ((ch=*((vars+n)->name+i))==')')
		    endi=i-1;
		  i++;
		}
	      ch=*((vars+n)->name+starti+1);
	      j=ch-'0'; drugind1=j-1;
	      starti=endi+3; endi=-1;
	      if (strncmp((vars+n)->name+starti,"Sa",2)==0)
		{
		  starti+=3; endi=-1; i=starti;
		  while ((i<strlen((vars+n)->name))&&(endi<0))
		    {
		      if ((ch=*((vars+n)->name+i))==',')
			endi=i-1;
		      i++;
		    }
		  if (strncmp((vars+n)->name+starti,"S",endi-starti+1)==0)
		    typeind2=0;
		  else
		    {
		      ch=*((vars+n)->name+starti+1);
		      j=ch-'0'; typeind2=j;
		    }
		  starti=endi+2; endi=-1; i=starti;
		  while ((i<strlen((vars+n)->name))&&(endi<0))
		    {
		      if ((ch=*((vars+n)->name+i))==')')
			endi=i-1;
		      i++;
		    }
		  ch=*((vars+n)->name+starti+1);
		  j=ch-'0'; drugind2=j-1;
		}
	      if ((typeind2<0)&&(drugind2<0))
		refval=*(g0+0);
	      else		
		//refval=*(Sa+typeind2*ndrugs+drugind2);
		refval=*(Sa+typeind2*ndrugs+ndrugs-1-drugind2);

	      // If typeind1 is S, then the cell is sensitive to drugind1.
	      // Otherwise the cell is resistant to drugind1.
	      for (i=0; i<=(ntypes-1); i++)
		{
		  num2vec(i,ndrugs,binarybases,resvec);
		  //if ((typeind1==0)&&(*(resvec+drugind1)==0))
		  //*(Sa+i*ndrugs+drugind1)=refval*val;
		  if ((typeind1==0)&&(*(resvec+ndrugs-1-drugind1)==0))
		    *(Sa+i*ndrugs+ndrugs-1-drugind1)=refval*val;
		  //else if ((typeind1>0)&&(*(resvec+drugind1)==1))
		  //*(Sa+i*ndrugs+drugind1)=refval*val;
		  else if ((typeind1>0)&&(*(resvec+ndrugs-1-drugind1)==1))
		    *(Sa+i*ndrugs+ndrugs-1-drugind1)=refval*val;

		  /*
		  // Debug
		  if (i==6)
		    {
		      printf("i=%d, vec[%d]=%d, ind=%d, refval=%f, val=%f\n",i,drugind1,*(resvec+drugind1),i*ndrugs+ndrugs-1-drugind1,refval,val);
		      printf("name=%s\n",(vars+n)->name);
		      printf("resvec: ");
		      for (j=0; j<=(ndrugs-1); j++)
			printf("%d ",*(resvec+j));
		      printf("\n");		      
		    }
		  */
		}

	      // Debug
	      //printf("typeind1=%d, drugind1=%d, typeind2=%d, drugind2=%d, name=%s\n",typeind1,drugind1,typeind2,drugind2,(vars+n)->name);


	      // Write values of Sa ratios to vals for 3-drug cases.
	      if (ndrugs==3)
		{
		  if (strcmp((vars+n)->name,"Sa(S,D1)/g0")==0)
		    *(vals+0)=val;
		  else if (strcmp((vars+n)->name,"Sa(S,D2)/Sa(S,D1)")==0)
		    *(vals+1)=val;
		  else if (strcmp((vars+n)->name,"Sa(S,D3)/Sa(S,D1)")==0)
		    *(vals+2)=val;
		}

	    }

	  // Transition probability T(type1,type2).
	  else if (strncmp((vars+n)->name,"T(",2)==0)
	    {
	      int typeind1=-1, typeind2=-1;
	      int starti=2, endi=-1;	      
	      i=starti;
	      while ((i<strlen((vars+n)->name))&&(endi<0))
		{
		  if ((ch=*((vars+n)->name+i))==',')
		    endi=i-1;
		  i++;
		}
	      if (strncmp((vars+n)->name+starti,"S",endi-starti+1)==0)
		typeind1=0;
	      else
		{
		  ch=*((vars+n)->name+starti+1);
		  j=ch-'0'; typeind1=j;
		}
	      starti=endi+1; endi=-1; i=starti;
	      while ((i<strlen((vars+n)->name))&&(endi<0))
		{
		  if ((ch=*((vars+n)->name+i))==')')
		    endi=i-1;
		  i++;
		}
	      if (strncmp((vars+n)->name+starti,"S",endi-starti+1)==0)
		typeind2=0;
	      else
		{
		  ch=*((vars+n)->name+starti+1);
		  j=ch-'0'; typeind2=j;
		}
	      // For transition from type i to type j, check whether type i is compatible with typeind2, and type j is compatible with typeind1.
	      for (i=0; i<=(ntypes-1); i++)
		{
		  num2vec(i,ndrugs,binarybases,resvec);
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      if (i!=j)
			{
			  num2vec(j,ndrugs,binarybases,resvec2);
			  k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    if (*(resvec+l)!=*(resvec2+l))
			      k++;			  
			  if ((k==1)&&(*(resvec+ndrugs-typeind1)==0)&&(*(resvec2+ndrugs-typeind1)==1))
			      *(T+j*ntypes+i)=val; 
			}
		    }
		}			  	      
	    }
	}

      *(x0+0)=N0;
      for (i=1; i<=(ntypes-1); i++)
	*(x0+0)-=*(x0+i);

      /*
      // Debug      
      printf("g0=%f, x0= ",*(g0+0));
      for (i=0; i<=(ntypes-1); i++)
	printf("%.2e ",*(x0+i));
      printf("\n");
      printf("Sa=\n");
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%f ",*(Sa+i*ndrugs+j));
	  printf("\n");
	}
      printf("T=\n");
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ntypes-1); j++)
	    printf("%.2e ",*(T+i*ntypes+j));
	  printf("\n");
	}     
      //exit(0);
      */

      /*
      // Debug
      {
	double *dosage, *finaly;
	dosage=malloc(sizeof(double)*ndrugs);
	*(dosage+0)=0; *(dosage+1)=1;
	*(x0+0)=0; *(x0+2)=5e9; *(x0+1)=0; *(x0+3)=0;
	finaly=malloc(sizeof(double)*ntypes);
	recurse_multidrug_response_trimrates2_2(x0,T,g0,Sg,a0,Sa,0,45,dosage,finaly,0);
	for (i=0; i<=(ntypes-1); i++)
	  printf("%.4e ",*(finaly+i));
	printf("\n");
	exit(0);
      }
      */



      // Check whether drug sensitities and initial populations satisfy certain constraints.
      // Hard-code those constraints here.
      allowed=1;
      
           
      // Negative initial populations are not allowed.
      i=0; 
      while ((i<=(ntypes-1))&&(allowed==1))
		{
		  if (*(x0+i)<0)
			allowed=0;
		  i++;
		}
	  
      
      // Debug: exclude indices below 3998889.
      //if (index<=3998889)
      //allowed=0;
      
      // Debug: exclude indices below 12167069.
      //if (index<=12167069)
      //allowed=0;


	  //MDM removed for unconstrained run 052819
	  //      // For 3-drug cases, only 7 combinations of drug sensitivities are allowed.
	  //      if ((ndrugs==3)&&(allowed==1))
	  //	{
	  //	  allowed=0;
	  //	  if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-1)<=1e-5)&&(fabs(*(vals+2)-1)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-1)<=1e-5)&&(fabs(*(vals+2)-0.3)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-1)<=1e-5)&&(fabs(*(vals+2)-0.1)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-0.3)<=1e-5)&&(fabs(*(vals+2)-0.3)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-0.3)<=1e-5)&&(fabs(*(vals+2)-0.1)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-5)<=1e-5)&&(fabs(*(vals+1)-0.1)<=1e-5)&&(fabs(*(vals+2)-0.1)<=1e-5))
	  //	    allowed=1;
	  //	  else if ((fabs(*(vals+0)-0.5)<=1e-5)&&(fabs(*(vals+1)-1)<=1e-5)&&(fabs(*(vals+2)-1)<=1e-5))
	  //	    allowed=1;
	  //	}
	  //      
	  //      free(vals);
	  //      
	  //
	  //      // For 3-drug cases, certain combinations of initial population are not allowed.
	  //      if ((ndrugs==3)&&(allowed==1))
	  //	{
	  //	  if ((allowed==1)&&(*(x0+0)<=1e-5))
	  //	    allowed=0;    	  
	  //	  if ((allowed==1)&&(*(x0+1)<=1e-5)&&((*(x0+3)>1)||(*(x0+5)>1)))
	  //	    allowed=0;	 
	  //	  if ((allowed==1)&&(*(x0+2)<=1e-5)&&((*(x0+3)>1)||(*(x0+6)>1)))
	  //	    allowed=0;	   
	  //	  if ((allowed==1)&&(*(x0+4)<=1e-5)&&((*(x0+5)>1)||(*(x0+6)>1)))
	  //	    allowed=0;	  
	  //	  if ((allowed==1)&&((*(x0+1)+*(x0+2)+*(x0+4))<=1e-5))
	  //	    allowed=0;
	  //	}
	  //MDM end removed
	  
	  //MDM removed 052819
      // the total population greatly exceeds N0, then it is not allowed.
	  //	  if (allowed==1)
	  //		{
	  //		  val=0;
	  //		  for (i=0; i<=(ntypes-1); i++)
	  //			val+=*(x0+i);
	  //	  	  val/=(double)N0;
	  //	  	  if (val>=1.05)
	  //	  	    allowed=0;
	  //		}
	  
	        
	  // S population is total population subtracts from the remaining populations.
	  if (allowed==1)
	  	{
	  	  *(x0+0)=N0;
	  	  for (i=1; i<=(ntypes-1); i++)
	  	    *(x0+0)-=*(x0+i);
	  	}    

	      // If the multiple mutants are curable, then the parameters are not allowed.
	  if (allowed==1)
	  	{
	  	  i=0;
	  	  while ((i<=(ndrugs-1))&&(allowed==1))
	  	    {
	  	      if (*(Sa+(ntypes-1)*ndrugs+i)>*(g0+ntypes-1))	       
	  		allowed=0;		  
	  	      i++;
	  	    }
	  	}

	  //MDM 052819 Removed from full run
	  // If some drugs cannot eradicate S cells, then the parameters are not allowed
	  //	  if (allowed==1)
	  //	  	{
	  //	  	  i=0; 
	  //	  	  while ((i<=(ndrugs-1))&&(allowed==1))
	  //	  	    {
	  //	  	      if (*(Sa+0*ndrugs+i)<*(g0+0))
	  //	  		allowed=0;
	  //	  	      i++;
	  //	  	    }
	  //	  	}            

             
      // Debug
      //if (index<=358519445)
      //allowed=0;

      // Debug
      //if (index<=121999381)
      //allowed=0;

      // Debug
      //if (index<=9375220)
      //allowed=0;
      

      // Debug
      //if (index<=11624733)
      //allowed=0;

      // Debug
      //if (index<=12448328)
      //allowed=0;


      // Debug
      //printf("index=%ld, allowed=%d\n",index,allowed);
     
     
      
      // Exclude the trivial cases where all strategies either all cure the patients or yield very similar or identical survival times.
      // To avoid direct computation of survival times of dynamic strategies, we use the results of two static dosages as the lower and upper bounds of survival times.
      // Also, instead of evaluating survival times, we calculate the total population at the max time span (5 years).
      // Lower bound on final population: administer full dosage for all drugs (not allowed).
      // Upper bound on final population: the best legal static strategy.
      // The following cases are excluded.
      // (1)Both upper and lower populations <= negligible values.  All dynamic strategies can cure the patients.
      // (2)Both upper and lower populations are between negligible and mortal values.  Most strategies would allow patients survive beyond the max time span, but cannot eradicate cancer.
      // (3)Both upper and lower populations > mortal values, and upper/lower ratio <= 100.  Patients will die before the max time span, and their survival times are close as the upper and lower bounds are tight.

	  //MDM Removing constraints for full run
	  
	  //      if (allowed==1)
	  //	{
	  //	  double *dosage, *finaly, val1, val2, maxtime=*(t+nintervals-1);
	  //	  dosage=malloc(sizeof(double)*ndrugs);
	  //	  finaly=malloc(sizeof(double)*ntypes);
	  //
	  //	  // Upper bound on final population.  Choose the minimum among the legal strategies.
	  //	  val1=1e500; beststatic=-1;
	  //	  for (i=0; i<=(nstates-1); i++)
	  //	    {
	  //	      for (j=0; j<=(ndrugs-1); j++)
	  //		*(dosage+j)=*(dstates+i*ndrugs+j);
	  //	      recurse_multidrug_response_trimrates2_2(x0,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);  	      
	  //	      //recurse_multidrug_response_trimrates2_2(x0,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1); 
	  //	      val=0;
	  //	      for (j=0; j<=(ntypes-1); j++)
	  //		val+=*(finaly+j);
	  //	      if (val<val1)
	  //		beststatic=i;
	  //	      val1=min(val1,val);
	  //
	  //	      // Debug
	  //	      //printf("dosage ( ");
	  //	      //for (j=0; j<=(ndrugs-1); j++)
	  //	      //printf("%.2f ",*(dosage+j));
	  //	      //printf("), val=%.3e\n",val);
	  //
	  //	    }
	  //
	  //	  // Lower bound on final population.  Administer max dosage to each drug.
	  //	  val2=0;
	  //	  for (i=0; i<=(ndrugs-1); i++)
	  //	    *(dosage+i)=1;
	  //	  recurse_multidrug_response_trimrates2_2(x0,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
	  //	  //recurse_multidrug_response_trimrates2_2(x0,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
	  //	  val=0;
	  //	  for (j=0; j<=(ntypes-1); j++)
	  //	    val+=*(finaly+j);
	  //	  val2=val;
	  //	  
	  //	  // Debug
	  //	  //printf("dosage ( ");
	  //	  //for (j=0; j<=(ndrugs-1); j++)
	  //	  //printf("%.2f ",*(dosage+j));
	  //	  //printf("), val=%.3e\n",val);
	  //
	  //	  free(dosage); free(finaly);
	  //	  
	  //	  // Debug
	  //	  //printf("val1=%.4e, val2=%.4e\n",val1,val2);
	  //
	  // Only consider the parameter configurations where the gap between the lower bound (full dosage for all drugs) and the upper bound (best static therapy) of predicted populations is large enough.
	  //	  if ((val1<=negligible)&&(val2<=negligible))
	  //	    allowed=0;
	  //	  if ((val1>negligible)&&(val1<=mortal)&&(val2>negligible)&&(val2<=mortal))
	  //	    allowed=0;
	  //	  if ((val1>mortal)&&(val2>=mortal)&&(val1<(val2*ratiothre)))
	  //	    allowed=0;
	  //
	  //	  lwd=val2; upp=val1;
	  //
	  //	  // Debug
	  //	  //if (allowed==1)
	  //	  //printf("index=%ld, lwd=%.3e, upp=%.3e, allowed=%d\n",index,lwd,upp,allowed);
	  //	  //exit(0);
	  //
	  //	}            
	  //
	  //      /*
	  //      // Debug
	  //      if (allowed==1)
	  //	{
	  //	  allowedcnt++; allowed=0;	  
	  //	}
	  //      if ((index%1000000)==0)
	  //	printf("index=%ld, allowedcnt=%ld\n",index,allowedcnt);
	  //      */
	  //MDM end of removal      

      /*
      // Debug
      if (allowed==1)
	{
	  sprintf(str,"%ld 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1845.00 0\0",index);
	  printf("%s\n",str);
	  if (firstwrite==0)
	    {
	      firstwrite=1;
	      fout=fopen(filename,"w");	      
	    }
	  fprintf(fout,"%s\n",str);
	  //fprintf(fout,"%ld 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1800.00 1845.00 0\n",index);
	  allowed=0;	  	  
	}
      */

     
      // Debug
      //printf("allowed=%d\n",allowed);

      // If the parameter configuration is allowed, then incur optimization strategies.
      if (allowed==1)
	{
	  //int nsinglestrategies=5;
	  //MDM - reduce number of stratigies tried
	  int nsinglestrategies=3;
	  double lwdT=-1, uppT=-1;

	  // Debug
	  allowedcnt++;

	  // The five previously used heuristics plus a new heuristic of maximizing anticipated death time.
	  nintervals=ceil((maxt-mint)/timeinterval)+1;

//Trying to implement new strategy



	  for (i=0; i<=(nsinglestrategies-1); i++)
	    {
	      int strategyind=i;

/*MDM changed to force only 2.2 and modified 2.2 for first 2 steps
	      if (i==2)
		popthre=1e9;
	      else if (i==3)		
		popthre=1e11; 		
	      else
		popthre=0;
*/
	      if (i==2)
		popthre=1e11;
	      else if (i==1)		
		popthre=1e11; 		
	      else
		popthre=0;

	      if (i>=3)
		strategyind=i-1;
	      
	      *(stopTs+i)=correct_optimize_multidrug_responses2(x0,T,g0,Sg,a0,Sa,t,timeinterval,nintervals,y,strategyind,1,popthre,predictionperiod,dosage);
	      		          

	      // Debug
	      //printf("stopTs[%d]=%f\n",i,*(stopTs+i));

	      for (j=0; j<=(nintervals-1); j++)
		{
		  for (k=0; k<=(ndrugs-1); k++)
		    {
		      if (*(t+j)<=*(stopTs+i))
			*(dosages+i*nintervals*ndrugs+j*ndrugs+k)=*(dosage+k*nintervals+j);
		      else
			*(dosages+i*nintervals*ndrugs+j*ndrugs+k)=-1;
		    }
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      if (*(t+j)<=*(stopTs+i))
			*(pops+i*nintervals*ntypes+j*ntypes+k)=*(y+k*nintervals+j+1);
		      else
			*(pops+i*nintervals*ntypes+j*ntypes+k)=-1;
		    }
		}
	      
	      // Debug
	      //printf("stopT[%d]=%.2f\n",i,*(stopTs+i));

	    }



//MDM commenting out to remove dpm strategies
/*
	  // Dynamic programming with five steps ahead.
	  // The baseline strategies are strategies 1-5.
	  // Do not incur dynamic programming if the single-step strategy already reaches the maximum time limit.
	  for (i=1; i<=(nsinglestrategies-1); i++)	  	  
	    {
	      int strategyind=i;

	      if (i==2)
		popthre=1e9;
	      else if (i==3)		
		popthre=1e11; 		
	      else
		popthre=0;

	      if (i>=3)
		strategyind=i-1;

	      nsteps=5; degmode=0; 
	      
	      for (j=0; j<=(nintervals-1); j++)
		{
		  for (k=0; k<=(ndrugs-1); k++)
		    *(dosages+(nsinglestrategies-1+i)*nintervals*ndrugs+j*ndrugs+k)=*(dosages+i*nintervals*ndrugs+j*ndrugs+k);
		  for (k=0; k<=(ntypes-1); k++)
		    *(pops+(nsinglestrategies-1+i)*nintervals*ntypes+j*ntypes+k)=*(pops+i*nintervals*ntypes+j*ntypes+k);
		}

	      if (*(stopTs+i)>=(maxt-0.1))
		*(stopTs+nsinglestrategies-1+i)=*(stopTs+i);
	      else	      
		{		
		  *(stopTs+nsinglestrategies-1+i)=correct_dp_optimize_multidrug_responses2_2(x0,T,g0,Sg,a0,Sa,t,timeinterval,nintervals,y,nsteps,dosage,degmode,strategyind,popthre);
		  
		  for (j=0; j<=(nintervals-1); j++)
		    {
		      for (k=0; k<=(ndrugs-1); k++)
			{
			  if (*(t+j)<=*(stopTs+nsinglestrategies-1+i))
			    *(dosages+(nsinglestrategies-1+i)*nintervals*ndrugs+j*ndrugs+k)=*(dosage+k*nintervals+j);
			  else
			    *(dosages+(nsinglestrategies-1+i)*nintervals*ndrugs+j*ndrugs+k)=-1;
			}
		      for (k=0; k<=(ntypes-1); k++)
			{
			  if (*(t+j)<=*(stopTs+nsinglestrategies-1+i))
			    *(pops+(nsinglestrategies-1+i)*nintervals*ntypes+j*ntypes+k)=*(y+k*nintervals+j+1);
			  else
			    *(pops+(nsinglestrategies-1+i)*nintervals*ntypes+j*ntypes+k)=-1;
			}
		    }
		}	      

	      // Debug
	      //printf("stopT[%d]=%.2f\n",i,*(stopTs+i));

	    }
*/

/* MDM removed	    	 
	  if (ndrugs==2)
	    nsteps=5;	
	  else
	    nsteps=3;
*/


	  // Debug
	  //nsteps=5;

	  /*
	  // If all but strategy 0 report maxtime, then set the global dp stoptime to maxtime.  Do not need to spend time on calculation.
	  flag=1; i=1;
	  while ((i<=(nstrategies-2))&&(flag==1))
	    {
	      if (*(stopTs+i)<(*(t+nintervals-1)+timeinterval-0.1))
		flag=0;
	      i++;
	    }
	  */

/* MDM removed	    	 
	  for (j=0; j<=(nintervals-1); j++)
	    {
	      for (k=0; k<=(ndrugs-1); k++)
		*(dosages+(nstrategies-1)*nintervals*ndrugs+j*ndrugs+k)=-1;
	      for (k=0; k<=(ntypes-1); k++)
		*(pops+(nstrategies-1)*nintervals*ntypes+j*ntypes+k)=-1;
	    }


	  // If the majority of strategies report maxtime, then set the global dp stoptime to maxtime.  Do not need to spend time on calculation.
	  k=0; flag=0;
	  for (i=0; i<=(nstrategies-2); i++)
	    if (*(stopTs+i)>=(*(t+nintervals-1)+timeinterval-0.1))
	      k++;
	  if ((k*2)>=(nstrategies-1))
	    flag=1;
	  

	  // Debug
	  flag=0;


	  // If all strategies are curable, then set the global dp stoptime is curable.  Do not need to spend time on calculation.
	  k=0; flag=0;
	  for (i=0; i<=(nstrategies-2); i++)
	    if (*(stopTs+i)>=(*(t+nintervals-1)+timeinterval-0.1))
	      k++;
	  if (k==(nstrategies-1))
	    {
	      flag=1; *(stopTs+nstrategies-1)=*(t+nintervals-1)+timeinterval;
	    }

	  
	  // Debug
	  // If the max survival time of all strategies >= 1800, then set the global dp stoptime to the max survival time.
	  val=1800;
	  for (i=0; i<=(nstrategies-2); i++)
	    if (*(stopTs+i)<val)
	      val=*(stopTs+i);
	  if (val>=1800)
	    {
	      flag=1; *(stopTs+nstrategies-1)=val;
	    }
MDM removed */ 	 

	  /*
	  // If all strategies report maxtime, then set the global dp stoptime maxtime.  Do not need to spend time on calculation.  
	  if (flag==0)
	    {
	      k=0;
	      for (i=0; i<=(nstrategies-2); i++)
		if (fabs(*(stopTs+i)-*(t+nintervals-1))<0.1)
		  k++;
	      if (k==(nstrategies-1))
		{
		  flag=1; *(stopTs+nstrategies-1)=*(t+nintervals-1);
		}
	    }
	  */

	  // Debug
	  //printf("index=%ld\n",index);

	  //if (flag==1)
	  //*(stopTs+nstrategies-1)=*(t+nintervals-1)+timeinterval;


	  /*
	  // Debug: skip if all other strategies report survival time 1800.
	  allowed=0; i=0;
	  while ((i<=(nstrategies-2))&&(allowed==0))
	    {
	      if (fabs(*(stopTs+i)-*(t+nintervals-1))>0.1)
		allowed=1;
	      i++;
	    }
	  */
	 

	  //if (allowed==1)
	  //{


	  // Debug
	  //flag=0;
	  //printf("flag=%d\n",flag);
	  
	  //else
/* MDM removed	    	 
	  if (flag==0)
	    {	  
	      
	      // First set maxncandseqs to 50 and curevamode to 1.  These parameters restrict the max parallel candidates for recursion to 50 and apply direct matrix exponential to calculate the bounds of survival times.  Turnouts are faster.
	      maxncandseqs=defaultmaxncandseqs; curevamode=defaultevamode;

	      // Debug
	      //maxncandseqs=highmaxncandseqs; 
	      //curevamode=highevamode;

	      // Debug
	      //curevamode=highevamode; maxncandseqs=highmaxncandseqs; maxncandseqs=300;

	      // Debug
	      //*(stopTs+nstrategies-1)=0;
	      

	      *(stopTs+nstrategies-1)=global_dp_optimize_multidrug_responses2(x0,T,g0,Sg,a0,Sa,t,timeinterval,nintervals,y,nsteps,dosage,&truncated);

	      // Debug
	      //printf("t1=%.2f\n",*(stopTs+nstrategies-1));

	      // Debug
	      //*(stopTs+nstrategies-1)=0;	      
	      
	      // If the survival time from the global optimization is smaller than any heuristic, then change maxncandseqs to highmaxncandseqs and evamode to highevamode, so that the results are more accurate.
	      // But if ndrugs>2, then do not increase maxncandseqs as it will substantially slow down the speed.
	      i=0;
	      while ((i<=(nstrategies-2))&&(*(stopTs+nstrategies-1)>=*(stopTs+i)))
		i++;
	      if (i<(nstrategies-1))
		{
		  if (ndrugs==2)
		    maxncandseqs=highmaxncandseqs; 
		  curevamode=highevamode;

		  // Debug
		  //maxncandseqs=1000; nsteps=5;
		  //maxncandseqs=500;

		  // Debug
		  //maxncandseqs=highmaxncandseqs;

		  // Debug
		  //nsteps=7; //maxncandseqs=1000;

		  *(stopTs+nstrategies-1)=global_dp_optimize_multidrug_responses2(x0,T,g0,Sg,a0,Sa,t,timeinterval,nintervals,y,nsteps,dosage,&truncated);

		  // Debug
		  //printf("nsteps=%d, t2=%.2f\n",nsteps,*(stopTs+nstrategies-1));

		}
	      
	      for (j=0; j<=(nintervals-1); j++)
		{
		  for (k=0; k<=(ndrugs-1); k++)
		    {
		      if (*(t+j)<=*(stopTs+nstrategies-1))
			*(dosages+(nstrategies-1)*nintervals*ndrugs+j*ndrugs+k)=*(dosage+k*nintervals+j);
		      else
			*(dosages+(nstrategies-1)*nintervals*ndrugs+j*ndrugs+k)=-1;
		    }
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      if (*(t+j)<=*(stopTs+nstrategies-1))
			*(pops+(nstrategies-1)*nintervals*ntypes+j*ntypes+k)=*(y+k*nintervals+j+1);
		      else
			*(pops+(nstrategies-1)*nintervals*ntypes+j*ntypes+k)=-1;
		    }
		}

	    }


	  // Debug
	  //exit(0);

	  // Debug
	  //*(stopTs+nstrategies-1)=0;

	  // Debug
	  //printf("stopT[%d]=%.2f\n",nstrategies-1,*(stopTs+nstrategies-1));
	  
MDM removed */	  
 	  	      
	  // Report simulation outcomes.

	  if (firstwrite==0)
	    {
	      firstwrite=1;
	      fout1=fopen(paramfilename,"w");
	      fout2=fopen(stoptimefilename,"w");
	      fout3=fopen(dosagefilename,"w");
	      fout4=fopen(popfilename,"w");
	    }

	  // Report parameter combination values: x0, g0, Sa, T.
	  sprintf(str,"%ld \0",index);
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      sprintf(tmpstr,"%.4e \0",*(x0+i));
	      strncat(str,tmpstr,strlen(tmpstr));
	    }
	  sprintf(tmpstr,"%.4e \0",*(g0+0));
	  strncat(str,tmpstr,strlen(tmpstr));
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ndrugs-1); j++)
		{
		  sprintf(tmpstr,"%.4e \0",*(Sa+i*ndrugs+j));
		  strncat(str,tmpstr,strlen(tmpstr));
		}
	    }
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ntypes-1); j++)
		{
		  sprintf(tmpstr,"%.4e \0",*(T+i*ntypes+j));
		  strncat(str,tmpstr,strlen(tmpstr));
		}
	    }
	  fprintf(fout1,"%s\n",str);


	  // Report stoptimes: each line reports the parameter configuration index and stoptimes for all strategies.
	  sprintf(str,"%ld \0",index);
	  for (i=0; i<=(nstrategies-1); i++)
	    {
	      if (i<(nstrategies-1))
		sprintf(tmpstr,"%.2f \0",*(stopTs+i));
	      else
		sprintf(tmpstr,"%.2f\0",*(stopTs+i));	      
	      strncat(str,tmpstr,strlen(tmpstr));
	    }
	  fprintf(fout2,"%s\n",str);


	  // Report dosages: each line reports the dosage combination sequence of each strategy for each parameter configuration.
	  // parameter index -- strategy index -- t=0 (dosage1,dosage2,...) -- t=1 (dosage1,dosage2,...) -- ....
	  for (i=0; i<=(nstrategies-1); i++)
	    {
	      sprintf(str,"%ld %d \0",index,i);
	      //for (j=0; j<=(nintervals-1); j++)
	      for (j=0; j<=(nintervals-2); j++)
		{
		  for (k=0; k<=(ndrugs-1); k++)
		    {
		      val=*(dosages+i*nintervals*ndrugs+j*ndrugs+k);
		      if ((j==(nintervals-1))&&(k==(ndrugs-1)))
			sprintf(tmpstr,"%.2f\0",val);
		      else
			sprintf(tmpstr,"%.2f \0",val);
		      strncat(str,tmpstr,strlen(tmpstr));
		    }
		}
	      
	      fprintf(fout3,"%s\n",str);
	    }


	  // Report pops: each line reports the population composition sequence of each strategy for each parameter configuration.
	  // parameter index -- strategy index -- t=0 (pop1,pop2,...) -- t=1 (pop1,pop2,...) -- ....
	  for (i=0; i<=(nstrategies-1); i++)
	    {
	      sprintf(str,"%ld %d \0",index,i);
	      //for (j=0; j<=(nintervals-1); j++)
	      for (j=0; j<=(nintervals-2); j++)
		{
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      val=*(pops+i*nintervals*ntypes+j*ntypes+k);
		      if ((j==(nintervals-1))&&(k==(ntypes-1)))
			sprintf(tmpstr,"%.4e\0",val);
		      else
			sprintf(tmpstr,"%.4e \0",val);
		      strncat(str,tmpstr,strlen(tmpstr));
		    }
		}	      
	      fprintf(fout4,"%s\n",str);
	    }
		
	}

      // Debug
      //if ((index%1000)==0)
      //printf("index=%ld, allowedcnt=%ld\n",index,allowedcnt);

    }

  
  if (firstwrite==1)
    {
      fclose(fout1);
      fclose(fout2);
      fclose(fout3);
      fclose(fout4);
    }

  t2=malloc(sizeof(time_t)); t2=time(t2);

  val=difftime(t2,t1);


  // Debug
  //printf("dt=%f\n",val);

  return 0;
}

