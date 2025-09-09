// Utility functions of optimization strategies for multiple drugs.

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include "db.h"
# define max(i,j) (i>j ? i : j)
# define min(i,j) (i<j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)

/*
struct pair2 {
  int i;
  double s;
};
*/

extern int ndrugs, ntypes, initialtime;
extern double detected, mortal;
extern int maxnfrontiers;
extern double bsearchthre;

extern void padm(int ndim, double *A, double t, int p, double *E);
extern char * getitemval3(char *s, int ind, int *length, char sepch);
extern char * getitemval4(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
extern void diagonal(int ndim, double *vec, double *mat);
extern int pair2_cmp(struct pair2 *p1, struct pair2 *p2);
extern void drugmat_expm(int ntypes, double *A, double *expA);

double optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
void recurse_multidrug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
void recurse_multidrug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
void recurse_multidrug_response_trimrates2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode);
void recurse_multidrug_response_trimrates3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode);
double dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
double dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
double dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
struct seqnode * recurse_sequence_multidrug_responses(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);
struct seqnode * recurse_sequence_multidrug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);
double deathtime_multidrug(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates);
void num2vec(long index, int ndigits, int *bases, int *numvec);
long vec2num(int ndigits, int *bases, int *numvec);
long vec2num2(int ndigits, int *bases, int *numvec);
double optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
void recurse_multidrug_response_trimrates2_3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timespan, double *dosage, double *responses, int evamode);
double correct_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
double correct_dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre);
int power(int base, int exponent);


// Implement optimization heuristics for multiple drugs.
// Adapt optimize_two_drug_responses8.
// Difference from optimize_two_drug_responses8: exercise fewer possible drug combinations.
// For k drugs, exhaust combinations x1+...+xk=1.
// For all the xk!=0, equalize their dosages.
// There are totally 2^{k}-1 drug combinations.
double optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage)
{
  int i, j, k, l, m, n, cnt, *bases, combind, *resvec;
  double *d, *mininds, *mininds2, *mininds3, minval, minval2, minval3, *mask, maxval;
  double *F, *template, *constF, *y0, *vec, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, val2, pop, *drugrecords;

  d=malloc(sizeof(double)*ndrugs); 
  mininds=malloc(sizeof(double)*ndrugs);
  mininds2=malloc(sizeof(double)*ndrugs);
  mininds3=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
  mask=malloc(sizeof(double)*ntypes*ntypes);
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (i=0; i<=(ndrugs-1); i++)
    *(bases+i)=2;  

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
 
  drugrecords=malloc(sizeof(double)*(ndrugs+1+ntypes)*ntypes);
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(y0+i)=*(x0+i);
	      *(y+i*nintervals+n)=*(y0+i);
	    }
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);


      // Strategy 0 does not need prediction.      
      if (mode==0)
	{	  
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;	  

	  
	  // Initially find the drug that kills the highest number of cancer cells.
	  // Count the total populations killable by each drug.
	  // Exclude the sensitive clones.
	  if (n==0)
	    {
	      double killed[ndrugs];

	      for (i=0; i<=(ndrugs-1); i++)
		killed[i]=0;
	      for (i=1; i<=(ntypes-1); i++)
		{
		  num2vec(i,ndrugs,bases,resvec);
		  for (j=0; j<=(ndrugs-1); j++)
		    if (*(resvec+j)==0)
		      killed[j]+=*(y0+i);
		}	  
	      k=-1; maxval=-1;
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  if (killed[i]>maxval)
		    {
		      maxval=killed[i]; k=i;
		    }
		}
	      for (i=0; i<=(ndrugs-1); i++)
		*(mininds+i)=0;
	      *(mininds+k)=1;
	    }	      
	  
	  // Switch the drug when the total population reaches the double the nadir or the population reemerges.
	  else
	    {
	      int switching=0;
	      double nadir=1e50;

	      // Find the nadir of the current regimen.
	      // A nadir is the local minimum of the population under the current regimen.
	      // Start from the current time step, trace back until the change of regimen.
	      i=n-1; k=1;
	      while ((i>=0)&&(k==1))
		{
		  int p, q;

		  // Check whether the dosages at time step i equals the dosages at time step n-1.
		  p=0; q=1;
		  while ((p<=(ndrugs-1))&&(q==1))
		    {
		      if (fabs(*(dosage+p*nintervals+i)-*(dosage+p*nintervals+n-1))>1e-12)
			q=0;
		      p++;
		    }
		  // If so then update nadir.
		  if (q==1)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  if (i==0)
			    val+=*(x0+j);
			  else
			    val+=*(y+j*nintervals+i);
			}
		      nadir=min(nadir,val);
		    }
		  else
		    k=0;
		  i--;
		}

     
	      // If the current population >= 2 x nadir population, then switch.
	      if (pop>=(nadir*2))
		switching=1;

	      // If the current population reemerges from the level below the detection threshold, then switch.
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  if ((n-1)==0)
		    val+=*(x0+i);
		  else
		    val+=*(y+i*nintervals+n-1);
		}

	      if ((val<detected)&&(pop>=detected))
		switching=1;
	      

	      // Find the most effective drug other than the current choice.
	      if (switching==1)
		{
		  double killed[ndrugs];

		  // Count the total populations killable by each drug.
		  // Exclude the sensitive clones.
		  for (i=0; i<=(ndrugs-1); i++)
		    killed[i]=0;
		  for (i=1; i<=(ntypes-1); i++)
		    {
		      num2vec(i,ndrugs,bases,resvec);
		      for (j=0; j<=(ndrugs-1); j++)
			if (*(resvec+j)==0)
			  killed[j]+=*(y0+i);
		    }
		  j=-1; maxval=0;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if (*(dosage+i*nintervals+n-1)>maxval)
			{
			  maxval=*(dosage+i*nintervals+n-1); j=i;
			}
		    }
		  k=-1; maxval=-1;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if ((i!=j)&&(killed[i]>maxval))
			{
			  maxval=killed[i]; k=i;
			}
		    }
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=0;
		  *(mininds+k)=1;
		}
	      else
		{
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=*(dosage+i*nintervals+n-1);
		}	      
	    }
	  
	  for (i=0; i<=(ntypes-1); i++)
	    *(vec2+i)=*(y0+i);
	}

      // For other strategies, calculate the predictions.
      else
	{
	  // Mask matrix.
	  // Mask the terms involved in the initial populations below 1.
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ntypes-1); j++)
		{
		  if (*(y0+j)<1)
		    *(mask+i*ntypes+j)=0;
		  else
		    *(mask+i*ntypes+j)=1;
		}
	    }   	  
	  
	  // Change the drug combinations to optimize certain criterion.
	  // Here the predicted responses are within the designated time interval (predictionperiod) instead of the sampling time interval.

	  // mininds2 and minval2: the drug dosage that minimizes the total population.
	  // mininds3 and minval3: the drug dosage that minimizes the R12 population.

	  // Consider the drug combinations (d1,d2,...,dk), where di's are either 0 or 1/m, where m is the number of nonzero di's.
	  // Exhaust all combinations of drugs except (0,0,...0).
	  
	  for (i=0; i<=(ndrugs-1); i++)
	    {
	      *(mininds+i)=0; *(mininds2+i)=0; *(mininds3+i)=0;
	    }

	  minval=1e100; minval2=1e100; minval3=1e100; cnt=0; 
	  for (combind=1; combind<=(ntypes-1); combind++)
	    {	     
	      num2vec(combind,ndrugs,bases,resvec);
	      k=0;
	      for (i=0; i<=(ndrugs-1); i++)	       
		k+=*(resvec+i);
	      val=(double)1/(double)k;
	      
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  *(vec+i)=*(resvec+i)*val;
		}
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(y0+i);
	  	  
		  
	      // compmode=0: use F to estimate the population.
	      if (compmode==0)
		{		  
		  mat_multiply(ntypes,ndrugs,1,Sg,vec,vec2);
		  diagonal(ntypes,vec2,template);	      
		  for (i=0; i<=(ntypes-1); i++)
		    for (j=0; j<=(ntypes-1); j++)
		      *(template2+i*ntypes+j)=0;
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      *(template2+i*ntypes+i)=1;
		      for (j=0; j<=(ntypes-1); j++)	
			*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
		    }		  
		  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);		  
		  mat_multiply(ntypes,ndrugs,1,Sa,vec,vec2);		  
		  diagonal(ntypes,vec2,template);		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)-=*(template+i*ntypes+j); 		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*predictionperiod;
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);		  
		}

	      // compmode=1: incur the recursive function to estimate the population.
	      else if (compmode==1)
		recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,predictionperiod,vec,vec2);
	     
	      /*
	      // Debug
	      printf("( ");
	      for (i=0; i<=(ndrugs-1); i++)
		printf("%.2f ",*(vec+i));
	      printf(") ( ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.2e ",*(vec2+i));
	      printf(")\n");
	      */
	      
	      // vec2 is the estimated population vector at the end of the time interval.
	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      
	      // Record the following quantities for each drug combination.
	      // (1)dosages of each drug, (2)total population, (3)population of each subclone.
	      cnt++; 
	      for (i=0; i<=(ndrugs-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+i)=*(vec+i);
	      *(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs)=g;
	      for (i=0; i<=(ntypes-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs+1+i)=*(vec2+i);
	    	
	      /*
	      // Debug
	      if ((mode==2)&&(popthre>1e10)&&(n==2))
		{
		  printf("cand %d, dosage (%.2f %.2f %.2f), pop (",cnt-1,*(vec+0),*(vec+1),*(vec+2));
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.2e ",*(vec2+i));
		  printf(")\n");
		}
	      */
	    }
		  
	  // Select the optimal dosage combination.
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;
	  	  
	  
	  // Strategy 1: Select the dosage combination to minimize total population.
	  if (mode==1)
	    {
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
		    {
		      for (j=0; j<=(ndrugs-1); j++)
			*(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
		      minval=val;		      
		    }
		}
	    }
	  
	  // Strategy 2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
	  else if (mode==2)
	    {
	      minval=1e20;

	      /*
	      // Debug
	      if ((popthre>1e10)&&(n>=5)&&(n<=9))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      printf("%d dosage=( ",i);
		      for (j=0; j<=(ndrugs-1); j++)
			printf("%.2f ",*(drugrecords+i*(ndrugs+1+ntypes)+j));
		      printf("), y=( ");
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.2e ",*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j));
		      printf(") %.2e\n",*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs));
		    }
		}	       
	      */	  

	      // Total population is above the threshold: minimize total population.
	      if (pop>=popthre)	      
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}
	      
	      // Total population is below the threshold: minimize the risk of developing incurable cells.
	      // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
	      // prob((1,...1) \rightarrow (1,...,1))=1.
	      // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
	      else
		{		  
		  for (i=0; i<=(cnt-1); i++)
		    {
		      double risk=0;

		      for (j=0; j<=(ntypes-1); j++)
			{
			  val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);			  
			  num2vec(j,ndrugs,bases,resvec);			  
			  for (k=0; k<=(ndrugs-1); k++)			    
			    if (*(resvec+k)==0)
			      val*=*(T+(2^k)*ntypes+0);
			  //val*=*(T+(k+1)*ntypes+0);
			  risk+=val;
			}		     

		      if (risk<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=risk;
			}
		    }
		}	
	    }

	  // Strategy 3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
	  else if (mode==3)
	    {
	      minval=1e20;

	      // If incurable cells > 1, then minimize total population.
	      if (*(y0+ntypes-1)>(1+1e-10))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}

	      // Otherwise minimize the predicted total multiple resistant populations.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      val=0;
		      for (j=1; j<=(ntypes-1); j++)
			{
			  num2vec(j,ndrugs,bases,resvec); k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    k+=*(resvec+l);
			  if (k>1)
			    val+=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);
			}			  		      
		      if (val<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;	
			}
		    }
		}
	    }
	}
	  	     	                  
      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);

      /*
      // Debug
      if ((mode==2)&&(popthre>1e10)&&(n==2))
	{
	  printf("y0=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.2e ",*(y0+i));
	  printf(")\n");
	  printf("T=\n");
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.2e ",*(T+i*ntypes+j));
	      printf("\n");
	    }
	  printf("g0=%f\n",*(g0+0));
	  printf("Sa=\n");
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%f ",*(Sa+i*ndrugs+j));
	      printf("\n");
	    }
	  printf("timeinterval=%f\n",timeinterval);
	  printf("mininds=( ");
	  for (i=0; i<=(ndrugs-1); i++)
	    printf("%.2f ",*(mininds+i));
	  printf(")\n");
	  printf("vec2=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.2e ",*(vec2+i));
	  printf(")\n");
	  exit(0);
	}
      */
      
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=*(mininds+i);              
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
              
     
      /*
      // Debug
      if ((mode==0)&&(n==0))
      //if ((mode==2)&&(popthre>1e10)&&(n==0))
      //if ((mode==1)&&(n==0))
      //if ((mode==3)&&(n==0))
	{
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
	}
      */

      /*
      // Debug
      //if (mode==0)
      //if ((mode==2)&&(popthre>1e10))
      if (mode==1)
      //if (mode==3)
	{	  
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(y+i*nintervals+n+1);
	  printf("t(%d)=%.2f, dosage=( ",n,*(t+n));
	  for (i=0; i<=(ndrugs-1); i++)
	    printf("%.2f ",*(dosage+i*nintervals+n));
	  printf("), y=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*(y+i*nintervals+n+1));
	  printf(") %.4e\n",val);	  
	  //printf("y0=( "); val=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //{
	  //  printf("%.4e ",*(y0+i)); val+=*(y0+i);
	  //}
	  //printf(") %.4e\n",val);
	}      
      */
           
 
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
           
      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask); free(resvec); free(bases);
  
  for (i=0; i<=(ndrugs-1); i++)
    *(dosage+i*nintervals+n)=*(dosage+i*nintervals+n-1);
  

  // If stop<0, then decide whether the patient recovers or not.
  // Check if each final subpopulation < 1.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {      
      k=-1; i=0;
      while ((i<=(ntypes-1))&&(k<0))
	{
	  if (*(y+i*nintervals+nintervals-2)>1)
	    k=i;
	  i++;
	}
      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }
  

  
  // If stopT<0 then the patient recovers.
  //if (stopT<0)
  //stopT=*(t+nintervals-1);   
  
 
  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}

  
// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
// Difference from recurse_two_drug_response_trimrates: multi-drug version.
void recurse_multidrug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;
  
  F=malloc(sizeof(double)*ntypes*ntypes);
  F2=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);
  vec3=malloc(sizeof(double)*ntypes);
  yn=malloc(sizeof(double)*ntypes); 
  mask=malloc(sizeof(double)*ntypes*ntypes);

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
  
  // Set initial conditions.
  y0=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *(y0+i)=*(x0+i);    

  // Mask matrix.
  // Mask the terms involved in the initial populations below 1.
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(y0+j)<1)
	    *(mask+i*ntypes+j)=0;
	  else
	    *(mask+i*ntypes+j)=1;
	}
    }   

  // Get the masked differential equation rates.
  mat_multiply(ntypes,ndrugs,1,Sg,dosage,vec2);
  diagonal(ntypes,vec2,template);  
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(template2+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(template2+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      
      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);
  
  mat_multiply(ntypes,ndrugs,1,Sa,dosage,vec2);
  diagonal(ntypes,vec2,template);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)		  
      *(F+i*ntypes+j)-=*(template+i*ntypes+j);
   
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));


  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template); 

  mat_multiply(ntypes,ntypes,1,template,y0,yn);
      
  alters=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(alters+i)=0;
      if ((*(y0+i)<1)&&(*(yn+i)>=1))
	*(alters+i)=1;
      else if ((*(y0+i)>=1)&&(*(yn+i)<1))
	*(alters+i)=-1;
    }

  k=0;
  for (i=0; i<=(ntypes-1); i++)
    if (*(alters+i)!=0)
      k++;
  

  // If there is no boundary crossing between startt and endt, then use the same differential equation throughout the time interval.
  if (k==0)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(responses+i)=*(yn+i);

      // Release memory.
      free(F); free(F2); free(vec2); free(vec3); free(constF); free(tempF); free(mask);
      free(y0); free(yn); free(template); free(template2); free(template3); free(alters);

    }
  // Otherwise find the time point where the first boundary crossing occurs.
  else
    {
      double t1=startt, t2=endt, curt=(t1+t2)/2;      
      for (i=0; i<=(ntypes-1); i++)
	{
	  *(vec2+i)=*(y0+i);
	}
      while ((t2-t1)>1e-5)
	{
	  curt=(t1+t2)/2; 

	  for (i=0; i<=(ntypes-1); i++)
	    for (j=0; j<=(ntypes-1); j++)
	      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(curt-t1);
	  padm(ntypes,tempF,1,7,template);

	  //padm(ntypes,F2,curt-t1,7,template);	      
	  mat_multiply(ntypes,ntypes,1,template,vec2,vec3);	  
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(alters+i)=0;

	      
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>0.01)&&((*(vec3+i)-1)>0.01))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>0.01)&&((1-*(vec3+i))>0.01))
		*(alters+i)=1;
		
	      //
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      //if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
	      //*(alters+i)=1;
	      //else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
	      //*(alters+i)=1;
	      //

	      //
	      //if ((*(vec2+i)<1)&&(*(vec3+i)>1))
	      //*(alters+i)=1;
	      //else if ((*(vec2+i)>1)&&(*(vec3+i)<1))
	      //*(alters+i)=-1;
	      //
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;
	
	  if (k==0)
	    {
	      t1=curt;
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(vec3+i);
	    }
	  else
	    {
	      t2=curt;
	    }	  
	}
   
      
      // Release memory.
      free(F); free(F2); free(vec2); free(constF); free(tempF); free(mask);
      free(yn); free(template); free(template2); free(template3); free(alters);
      
      if (max(t1,t2)>=endt)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(responses+i)=*(vec3+i);
	}

      // Modify the differential equation terms and run the equations again.      
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(vec3+i);
	  recurse_multidrug_response_trimrates(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
	}      

      // Release memory.
      free(y0); free(vec3);
    }
	    

  return;
}



// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
// Difference from recurse_two_drug_response_trimrates: reduce the tolerance error range for determining boundary crossing.  Can correct some errors from previous simulations.
// Difference from recurse_two_drug_response_trimrates2: multi-drug version.
void recurse_multidrug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;
  
  F=malloc(sizeof(double)*ntypes*ntypes);
  F2=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);
  vec3=malloc(sizeof(double)*ntypes);
  yn=malloc(sizeof(double)*ntypes); 
  mask=malloc(sizeof(double)*ntypes*ntypes);

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
  
  // Set initial conditions.
  y0=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *(y0+i)=*(x0+i);    

  // Mask matrix.
  // Mask the terms involved in the initial populations below 1.
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(y0+j)<1)
	    *(mask+i*ntypes+j)=0;
	  else
	    *(mask+i*ntypes+j)=1;
	}
    }   

  // Get the masked differential equation rates.
  mat_multiply(ntypes,ndrugs,1,Sg,dosage,vec2);
  diagonal(ntypes,vec2,template);  
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(template2+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(template2+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      
      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);
  
  mat_multiply(ntypes,ndrugs,1,Sa,dosage,vec2);
  diagonal(ntypes,vec2,template);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)		  
      *(F+i*ntypes+j)-=*(template+i*ntypes+j);
  
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));

  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template);

  //padm(ntypes,F2,endt-startt,7,template); 
 
  mat_multiply(ntypes,ntypes,1,template,y0,yn); 

  /*
  // Debug
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.2e ",*(y0+i));
  printf("\n");
  printf("template:\n");
  for (i=0; i<=(ntypes-1); i++)	
    {	
      for (j=0; j<=(ntypes-1); j++)
	printf("%.2e ",*(template+i*ntypes+j));
      printf("\n");
    }
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.2e ",*(yn+i));
  printf("\n");
  */

  alters=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(alters+i)=0;
      if ((*(y0+i)<1)&&(*(yn+i)>=1))
	*(alters+i)=1;
      else if ((*(y0+i)>=1)&&(*(yn+i)<1))
	*(alters+i)=-1;
    }

  k=0;
  for (i=0; i<=(ntypes-1); i++)
    if (*(alters+i)!=0)
      k++;

  /*
  // Debug
  printf("alters= ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%d ",*(alters+i));
  printf("\n");
  printf("k=%d\n",k);
  */

  // If there is no boundary crossing between startt and endt, then use the same differential equation throughout the time interval.
  if (k==0)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(responses+i)=*(yn+i);

      /*
      // Debug
      printf("responses: ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.2e ",*(responses+i));
      printf("\n");
      */

      // Release memory.
      free(F); free(F2); free(vec2); free(vec3); free(constF); free(tempF); free(mask);
      free(y0); free(yn); free(template); free(template2); free(template3); free(alters);

    }
  // Otherwise find the time point where the first boundary crossing occurs.
  else
    {
      double t1=startt, t2=endt, curt=(t1+t2)/2;      
      for (i=0; i<=(ntypes-1); i++)
	{
	  *(vec2+i)=*(y0+i);
	}
      while ((t2-t1)>1e-5)
	{
	  curt=(t1+t2)/2; 

	  for (i=0; i<=(ntypes-1); i++)
	    for (j=0; j<=(ntypes-1); j++)
	      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(curt-t1);
	  padm(ntypes,tempF,1,7,template);

	  //padm(ntypes,F2,curt-t1,7,template);	      
	  mat_multiply(ntypes,ntypes,1,template,vec2,vec3);	  
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(alters+i)=0;

	      /*
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>0.01)&&((*(vec3+i)-1)>0.01))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>0.01)&&((1-*(vec3+i))>0.01))
		*(alters+i)=1;
		*/

	      
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
		*(alters+i)=1;
				

              //
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      //if (((1-*(vec2+i))>1e-6)&&((*(vec3+i)-1)>1e-6))
              //*(alters+i)=1;
	      //else if (((*(vec2+i)-1)>1e-6)&&((1-*(vec3+i))>1e-6))
              //*(alters+i)=1;
	      //

	      //
	      //if ((*(vec2+i)<1)&&(*(vec3+i)>1))
              //*(alters+i)=1;
	      //else if ((*(vec2+i)>1)&&(*(vec3+i)<1))
              //*(alters+i)=-1;
	      //
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;
	  if (k==0)
	    {
	      t1=curt;
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(vec3+i);
	    }
	  else
	    {
	      t2=curt;
	    }	  
	}

      
      // Release memory.
      free(F); free(F2); free(vec2); free(constF); free(tempF); free(mask);
      free(yn); free(template); free(template2); free(template3); free(alters);

      /*
      // Debug
      printf("startt=%.2f curtt=%.2f endt=%.2f, t1=%.2f, t2=%.2f\n",startt,curt,endt,t1,t2);
      printf("vec3= ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.2e ",*(vec3+i));
      printf("\n");      
      */

      if (max(t1,t2)>=endt)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(responses+i)=*(vec3+i);
	}

      // Modify the differential equation terms and run the equations again.      
      else
	{
	  // Debug
	  //printf("here\n");

	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(vec3+i);
	  recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
	}      

      // Release memory.
      free(y0); free(vec3);
    }
	    

  return;
}


// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
// Difference from recurse_two_drug_response_trimrates: reduce the tolerance error range for determining boundary crossing.  Can correct some errors from previous simulations.
// Difference from recurse_two_drug_response_trimrates2: multi-drug version.
// Difference from recurse_multidrug_response_trimrates2: do not apply U(x-1) in the equation if evamode=1.
// evamode=0: same as recurse_multidrug_response_trimrates2.
// evamode=1: direct calculation of matrix exponential.
void recurse_multidrug_response_trimrates2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;

  F=malloc(sizeof(double)*ntypes*ntypes);
  F2=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);
  vec3=malloc(sizeof(double)*ntypes);
  yn=malloc(sizeof(double)*ntypes); 
  mask=malloc(sizeof(double)*ntypes*ntypes);

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
  
  // Set initial conditions.
  y0=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *(y0+i)=*(x0+i);    

  // Mask matrix.
  // Mask the terms involved in the initial populations below 1.
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(y0+j)<1)
	    *(mask+i*ntypes+j)=0;
	  else
	    *(mask+i*ntypes+j)=1;
	  
	  if (evamode==1)
	    *(mask+i*ntypes+j)=1;
	}
    }   

  // Get the masked differential equation rates.
  mat_multiply(ntypes,ndrugs,1,Sg,dosage,vec2);
  diagonal(ntypes,vec2,template);  
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(template2+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(template2+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      
      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);
  
  mat_multiply(ntypes,ndrugs,1,Sa,dosage,vec2);
  diagonal(ntypes,vec2,template);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)		  
      *(F+i*ntypes+j)-=*(template+i*ntypes+j);
  
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));

  /*
  // Debug
  printf("startt=%.4f, endt=%.4f, F2=\n",startt,endt);
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.4e ",*(F2+i*ntypes+j));
      printf("\n");
    }
  */

  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template);

  //padm(ntypes,F2,endt-startt,7,template); 
 
  mat_multiply(ntypes,ntypes,1,template,y0,yn); 

  /*
  // Debug
  printf("startt=%.3e, endt=%.3e\n",startt,endt);
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.4e ",*(y0+i));
  printf("\n");
  
  printf("template:\n");
  for (i=0; i<=(ntypes-1); i++)	
    {	
      for (j=0; j<=(ntypes-1); j++)
	printf("%.2e ",*(template+i*ntypes+j));
      printf("\n");
    }
  
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.4e ",*(yn+i));
  printf("\n");
  */

  alters=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(alters+i)=0;
      if ((*(y0+i)<1)&&(*(yn+i)>=1))
	*(alters+i)=1;
      else if ((*(y0+i)>=1)&&(*(yn+i)<1))
	*(alters+i)=-1;
    }

  k=0;
  for (i=0; i<=(ntypes-1); i++)
    if (*(alters+i)!=0)
      k++;


  if (evamode==1)
    k=0;

  /*
  // Debug
  printf("alters= ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%d ",*(alters+i));
  printf("\n");
  printf("k=%d\n",k);
  */

  // If there is no boundary crossing between startt and endt, then use the same differential equation throughout the time interval.
  if (k==0)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(responses+i)=*(yn+i);

      /*
      // Debug
      printf("responses: ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.2e ",*(responses+i));
      printf("\n");
      */

      // Release memory.
      free(F); free(F2); free(vec2); free(vec3); free(constF); free(tempF); free(mask);
      free(y0); free(yn); free(template); free(template2); free(template3); free(alters);

    }
  // Otherwise find the time point where the first boundary crossing occurs.
  else
    {
      double t1=startt, t2=endt, curt=(t1+t2)/2;      
      for (i=0; i<=(ntypes-1); i++)
	{
	  *(vec2+i)=*(y0+i);
	}
      //while ((t2-t1)>1e-5)      
      //while ((t2-t1)>1e-1)
      while ((t2-t1)>bsearchthre)
	{
	  // Debug
	  //printf("t1=%.3e, t2=%.3e\n",t1,t2);

	  curt=(t1+t2)/2; 	 
	  for (i=0; i<=(ntypes-1); i++)
	    for (j=0; j<=(ntypes-1); j++)
	      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(curt-t1);
	  padm(ntypes,tempF,1,7,template);

	  //padm(ntypes,F2,curt-t1,7,template);	      
	  mat_multiply(ntypes,ntypes,1,template,vec2,vec3);	  

	  /*
	  // Debug
	  printf("vec2: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(vec2+i));
	  printf(", vec3: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(vec3+i));
	  printf("\n");
	  */


	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(alters+i)=0;

	      /*
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>0.01)&&((*(vec3+i)-1)>0.01))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>0.01)&&((1-*(vec3+i))>0.01))
		*(alters+i)=1;
		*/

	      
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
		*(alters+i)=1;
				

              //
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      //if (((1-*(vec2+i))>1e-6)&&((*(vec3+i)-1)>1e-6))
              //*(alters+i)=1;
	      //else if (((*(vec2+i)-1)>1e-6)&&((1-*(vec3+i))>1e-6))
              //*(alters+i)=1;
	      //

	      //
	      //if ((*(vec2+i)<1)&&(*(vec3+i)>1))
              //*(alters+i)=1;
	      //else if ((*(vec2+i)>1)&&(*(vec3+i)<1))
              //*(alters+i)=-1;
	      //
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;
	  if (k==0)
	    {
	      t1=curt;
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(vec3+i);
	    }
	  else
	    {
	      t2=curt;
	    }	  
	}

      
      // Release memory.
      free(F); free(F2); free(vec2); free(constF); free(tempF); free(mask);
      free(yn); free(template); free(template2); free(template3); free(alters);
      
      if (max(t1,t2)>=endt)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(responses+i)=*(vec3+i);
	}

      // Modify the differential equation terms and run the equations again.      
      else
	{	  
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(vec3+i);
	  //recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
	  recurse_multidrug_response_trimrates2_2(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses,evamode);
	}      

      // Release memory.
      free(y0); free(vec3);
    }
	    

  return;
}


// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
// Difference from dp_optimize_two_drug_responses4: apply branch and bound to trim unnecessary nodes in recursion.
// Difference from dp_optimize_two_drug_responses5: (1)multi-drug version, (2)make the format of y and dosage compatible with optimize_multidrug_responses1, (2)can choose one of multiple strategies as baseline for dp optimization.
// mode=0: monotherapy switching, benchmark.
// mode=1: minimize total population.
// mode=2: select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
// mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
// mode=4: maximize the anticipated death time for the best monotonous treatment.
double dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, *resvec, combind, *bases;
  long curind, *seq;
  double *dstates, val, stopT, *frontiers, pop, incurablepop;
  struct seqnode *seqnodes;  
    
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }
  

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }

  /*
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+n*ndrugs+i)=-1;
    }    

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);
    */

  seq=malloc(sizeof(long)*nsteps); n=0;
  nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

      // Debug
      //printf("n=%d\n",n);

      // (n%nsteps)=0: incur dynamic programming.
      if ((n%nsteps)==0)
	{
	  if (n>0)
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  free((seqnodes+i)->children);
		  free((seqnodes+i)->dosage);
		  free((seqnodes+i)->y);
		}
	      free(seqnodes);	     	      
	    }
	  nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
	  (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
	  (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
	  (seqnodes+0)->children=malloc(sizeof(long));
	  (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
	  (seqnodes+0)->y=malloc(sizeof(double)*ntypes);
	  if (n==0)
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+i);
	      
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+0);
	    }
	  else
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);

	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+n-1);
	    }

	  for (i=0; i<=nsteps; i++)
	    *(nfrontiers+i)=0;

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  	  	 
	  // Find the max depth of the seqnodes.
	  maxdepth=-1;
	  for (i=0; i<=(nseqnodes-1); i++)
	    maxdepth=max(maxdepth,(seqnodes+i)->depth);
	  	  

	  // Find the candidate terminal seqnodes.
	  ncandidates=0; candidates=malloc(sizeof(long));
	  
	  // Find the terminal nodes with zero population.
	  // If so then pick up the one that minimizes the total population.
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      if ((seqnodes+i)->nchildren==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  if (val<1)
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	    }
	  
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);
	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  
	  
	  // If the tumor cannot be eradicated within nsteps, then find the candidate sequences.
	  // A candidate sequence has the max depth, terminates at a leaf, and is not inferior to any other sequences.
	  // A sequence s1 is inferior to another sequence s2 if the total population and R12 population of s1 are both greater than those of s2.
	  if (ncandidates==0)
	    {
	      struct pair2 *pairs1, *pairs2;
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates, allmortal;
	      double p1, p2, p, *deathtimes;

	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0;
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  if ((c1==1)&&(c2==1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }		  
		}	     

	      // Sort candidate sequences by total population and by incurable population.
	      pairs1=malloc(sizeof(struct pair2)*ncandidates);
	      pairs2=malloc(sizeof(struct pair2)*ncandidates);
	      lookup1=malloc(sizeof(int)*ncandidates);
	      lookup2=malloc(sizeof(int)*ncandidates);
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  (pairs1+i)->i=i;
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  (pairs1+i)->s=val;
		  (pairs2+i)->i=i;
		  val=*((seqnodes+*(candidates+i))->y+ntypes-1);
		  (pairs2+i)->s=val;
		}

	      qsort(pairs1,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      qsort(pairs2,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  j=(pairs1+i)->i; *(lookup1+j)=i;
		  j=(pairs2+i)->i; *(lookup2+j)=i;
		}

	      
	      // Find the candidates inferior to others.
	      // Direct comparison of all candidate pairs take n^2 time, n=ncandidates.
	      // To reduce computing time use probabilistic arguments.
	      // Suppose a candidate's total population > n1 other candidates, and its incurable population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and incurable population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in incurable population) =
	      // 1-Pr(there is no overlap ...).
	      // Pr(there is no overlap ...) = 
	      // Pr(place n2 candidates in (n-n1) among n slots) or
	      // Pr(place n1 candidates in (n-n2) among n slots).
	      // If (n1+n2)>=n, then Pr(there is an overlap ...)=1.
	      // If (n1+n2)<n, then
	      // P1=((n-n1)/n)^n2, P2=((n-n2)/n)^n1. 
	      // Pr(the candidate is inferior) = 1 - max(P1,P2).
	      inferior=malloc(sizeof(int)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{		  
		  *(inferior+i)=0; p=0;
		  j=*(lookup1+i); n1=max(j-1,0);
		  j=*(lookup2+i); n2=max(j-1,0);
		  if ((n1+n2)>=ncandidates)
		    p=1;
		  else
		    {
		      val=(double)(ncandidates-n1)/(double)ncandidates;
		      p1=n2*log(val); p1=exp(p1);
		      val=(double)(ncandidates-n2)/(double)ncandidates;
		      p2=n1*log(val); p2=exp(p2);
		      p=1-max(p1,p2);
		    }

		  if (p>=(1-1e-4))
		    *(inferior+i)=1;
		  else
		    {
		      // Compare the superior list in two values.
		      j=0; k=-1;
		      while ((j<=n1)&&(k<0))
			{
			  l=(pairs1+j)->i; l=*(lookup2+l);
			  if (((pairs1+*(lookup1+i))->s>(pairs1+j)->s)&&((pairs2+*(lookup2+i))->s>(pairs2+l)->s))
			    //if (l<=n2)
			    k=l;
			  j++;
			}
		      if (k>=0)
			*(inferior+i)=1;
		    }
		}

	      // Check whether all candidates lead to mortality.
	      allmortal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allmortal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<mortal)
		    allmortal=0;
		  i++;
		}
	      
	      // Debug
	      //printf("ncandidates=%d, before\n",ncandidates);

	      // Choose a candidate according to the selected strategy.

	      // mode=1: Minimize the total population. 
	      if (mode==1)
		{		  
		  minind=-1; minval=1e100;	      
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			val+=*((seqnodes+*(candidates+i))->y+j);
		      if (val<minval)
			{
			  minind=*(candidates+i); minval=val;
			}
		    }	      
		}	      

	      // mode=2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==2)
		{		
		  minind=-1; minval=1e20;		  		  
		  
		  // Total population is above the threshold: minimize total population.
		  // Total population is above the threshold or all drug combinations lead to mortality: minimize total population.
		  //if (pop>=popthre)	      		  
		  if ((pop>=popthre)||(allmortal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }		  	     		
	      
		  // Total population is below the threshold: minimize the risk of developing incurable cells.
		  // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
		  // prob((1,...1) \rightarrow (1,...,1))=1.
		  // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
		  // Discard the combinations that lead to mortality.
		  else
		    {		  		      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  val*=*(T+(2^k)*ntypes+0);
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);			      			     
			    }

			  //if (risk<minval)
			  if ((risk<minval)&&(val2<mortal))
			    {			      
			      minind=*(candidates+i); minval=val;			      
			    }			  			  
			}		    
		    }
		}	
	    
	      // mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==3)
		{
		  minval=1e20;

		  // If incurable cells > 1, then minimize total population.
		  //if (incurablepop>=(1+1e-10))
		  if ((incurablepop>=(1+1e-10))||(allmortal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }	

		  // Otherwise minimize the predicted total multiple resistant populations.
		  // Discard the combinations that lead to mortality.
		  else
		    {
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val2;

			  val=0; val2=0;
			  for (j=1; j<=(ntypes-1); j++)
			    {
			      num2vec(j,ndrugs,bases,resvec); k=0;
			      for (l=0; l<=(ndrugs-1); l++)
				k+=*(resvec+l);
			      if (k>1)
				val+=*((seqnodes+*(candidates+i))->y+j);
			      val2+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  //if (val<minval)
			  if ((val<minval)&&(val2<mortal))
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }
		}

	      // mode=4: Maximize the anticipated death time based on static therapy.
	      else if (mode==4)
		{
	      	      	     	      
		  // Evaluate the expected death time for each non-inferior candidate.
		  // Apply the best monotherapy.
		  deathtimes=malloc(sizeof(double)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      *(deathtimes+i)=-1;
		      if (*(inferior+i)==0)
			{
			  val=deathtime_multidrug((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
			  
			  *(deathtimes+i)=val;		      		      
			}
		    }	
		  
		  // Debug
		  //printf("ncandidates=%d, after\n",ncandidates);
		  
		  
		  // Choose the candidate with the max death time.
		  minind=-1; minval=-1;
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if (*(deathtimes+i)>minval)
			{
			  minind=*(candidates+i); minval=*(deathtimes+i);
			}
		    }
		  
		  // If there are degenerate death times, then choose the one with minimum incurable population.
		  ndegenerates=0; degenerates=malloc(sizeof(int)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if ((minval-*(deathtimes+i))<=1e-5)
			{
			  ndegenerates++; *(degenerates+ndegenerates-1)=*(candidates+i);
			}
		    }
		  if (ndegenerates>1)
		    {
		      minind=-1; minval=1e100;
		      for (i=0; i<=(ndegenerates-1); i++)
			{		      		      
			  if (degmode==0)
			    {
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				val+=*((seqnodes+*(degenerates+i))->y+j);
			    }
			  else
			    val=*((seqnodes+*(degenerates+i))->y+ntypes-1);
			  

			  if (val<minval)
			    {
			      minind=*(degenerates+i); minval=val;
			    }
			}
		    }
		  free(degenerates); free(deathtimes);
		}

	 	      	      	    	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); 
	      
	    }
	     
	  free(candidates);
	  
	  for (i=0; i<=(nsteps-1); i++)
	    *(seq+i)=-1;

	  if (minind>=0)
	    {
	      k=(seqnodes+minind)->depth-1;
	      *(seq+k)=minind; i=k-1;
	      //*(seq+nsteps-1)=minind; i=nsteps-2;
	      while (i>=0)
		{
		  *(seq+i)=(seqnodes+*(seq+i+1))->pa;
		  i--;
		}	  
	    }	  	 
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  /*
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
	    */

	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j*nintervals+n)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+i*nintervals+n)=*((seqnodes+l)->y+i);



	  maxval=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //maxval+=*(y+n*ntypes+i);
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+i*nintervals+n);
	
	  pop=maxval; incurablepop=*(y+(ntypes-1)*nintervals+n);

	  // Debug
	  //printf("maxval=%.2e\n",maxval);

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);	  

	  /*
	  // Debug
	  if ((mode==2)&&(popthre>1e10))
	    {
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		val+=*(y+i*nintervals+n);
	      printf("dp t(%d)=%.2f, dosage=( ",n,*(t+n));
	      for (i=0; i<=(ndrugs-1); i++)
		printf("%.2f ",*(dosage+i*nintervals+n));
	      printf("), y=( ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.4e ",*(y+i*nintervals+n));
	      printf(") %.4e\n",val);
	    }
	  */
	}

      // Debug
      //printf("n=%d, stopT=%.2f\n",n,stopT);

      n++;
    }

  if (stopT<0)
    stopT=*(t+nintervals-1);

  // Release memory.
  free(seq); free(dstates);
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); free(nfrontiers); free(frontiers);
  free(resvec); free(bases);
  
  return stopT;
}



// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
// Difference from recurse_sequence_two_drug_responses3: multi-drug version.
struct seqnode * recurse_sequence_multidrug_responses(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt)
{
  int i, j, k, l, m, n, inferior, depth;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2, val1, val2;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
    

  // Return when the node depth>=nsteps.
  if ((seqnodes+curind)->depth>=nsteps)
    return seqnodes;

  // Return when the total population exceeds the mortal threshold.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val+=*((seqnodes+curind)->y+i);
  
 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;


  // Return when the total and incurable populations of the current node are inferior to the frontiers at the current level.
  depth=(seqnodes+curind)->depth; inferior=0; n=0; 
  val1=0;
  for (i=0; i<=(ntypes-1); i++)
    val1+=*((seqnodes+curind)->y+i);
  val2=*((seqnodes+curind)->y+ntypes-1); 
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // Check whether both total and incurable populations are larger than a frontier.           
      if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))&&(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	inferior=1;
      n++;
    }


  if (inferior==1)
    return seqnodes;

  // Augment/update the frontier points if the current node is not inferior to any of them.
  else
    {
      int nleft;
      double *left;      
      nleft=0; left=malloc(sizeof(double)*maxnfrontiers*2);
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  if ((val1<=*(frontiers+depth*maxnfrontiers*2+n*2+0))||(val2<=*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	    {
	      nleft++; 
	      *(left+(nleft-1)*2+0)=*(frontiers+depth*maxnfrontiers*2+n*2+0);
	      *(left+(nleft-1)*2+1)=*(frontiers+depth*maxnfrontiers*2+n*2+1);
	    }
	}
      for (n=0; n<=(nleft-1); n++)
	{
	  *(frontiers+depth*maxnfrontiers*2+n*2+0)=*(left+n*2+0);
	  *(frontiers+depth*maxnfrontiers*2+n*2+1)=*(left+n*2+1);
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+0)=val1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+1)=val2;
	}
      free(left);

    }

  
  // Spin off each dosage combination.
  seqnodes=realloc(seqnodes,sizeof(struct seqnode)*(nseqnodes+nstates));
  vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);


  (seqnodes+curind)->nchildren=nstates;
  (seqnodes+curind)->children=malloc(sizeof(long)*nstates);
  for (n=0; n<=(nstates-1); n++)
    {
      *((seqnodes+curind)->children+n)=nseqnodes+n;
      (seqnodes+nseqnodes+n)->index=nseqnodes+n;
      (seqnodes+nseqnodes+n)->depth=(seqnodes+curind)->depth+1;
      (seqnodes+nseqnodes+n)->pa=curind;
      (seqnodes+nseqnodes+n)->nchildren=0;
      (seqnodes+nseqnodes+n)->children=malloc(sizeof(long));
      (seqnodes+nseqnodes+n)->dosage=malloc(sizeof(double)*ndrugs);
      (seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ndrugs-1); i++)
	*(vec+i)=*(dstates+n*ndrugs+i);      
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);
     
      recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
           
      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=recurse_sequence_multidrug_responses(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind2);
    }

  return seqnodes;
}



// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
// Difference from recurse_sequence_two_drug_responses3: multi-drug version.
// Difference from recurse_sequence_multidrug_responses: remove the frontiers which are no longer frontiers.
struct seqnode * recurse_sequence_multidrug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt)
{
  int i, j, k, l, m, n, inferior, depth;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2, val1, val2;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
    

  // Return when the node depth>=nsteps.
  if ((seqnodes+curind)->depth>=nsteps)
    return seqnodes;

  // Return when the total population exceeds the mortal threshold.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val+=*((seqnodes+curind)->y+i);
  
 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;


  // Return when the total and incurable populations of the current node are inferior to the frontiers at the current level.
  depth=(seqnodes+curind)->depth; inferior=0; n=0; 
  val1=0;
  for (i=0; i<=(ntypes-1); i++)
    val1+=*((seqnodes+curind)->y+i);
  val2=*((seqnodes+curind)->y+ntypes-1); 
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // Check whether both total and incurable populations are larger than a frontier.           
      if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))&&(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	inferior=1;
      n++;
    }

  /*
  // Debug
  if ((seqnodes+curind)->depth==1)
    {
      //printf("nfrontiers=%d, frontiers= ",*(nfrontiers+depth));
      //for (n=0; n<*(nfrontiers+depth); n++)
      //printf("(%.3e %.3e) ",*(frontiers+depth*maxnfrontiers*2+n*2+0),*(frontiers+depth*maxnfrontiers*2+n*2+1));
      //printf("\n");
      printf("curind=%d, dosage ( ",curind);
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*((seqnodes+curind)->dosage+i));
      printf("), pop ( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	{
	  printf("%.3e ",*((seqnodes+curind)->y+i));
	  val+=*((seqnodes+curind)->y+i);
	}
      printf(" %.3e ), inferior=%d\n",val,inferior);
    }
  */


  if (inferior==1)
    return seqnodes;

  
  // Augment/update the frontier points if the current node is not inferior to any of them.
  // Also remove the frontier points inferior to the current node.
  else
    {
      int nleft;
      double *left;      
      nleft=0; left=malloc(sizeof(double)*maxnfrontiers*2);
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))||(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	    {
	      nleft++; 
	      *(left+(nleft-1)*2+0)=*(frontiers+depth*maxnfrontiers*2+n*2+0);
	      *(left+(nleft-1)*2+1)=*(frontiers+depth*maxnfrontiers*2+n*2+1);
	    }
	}
      for (n=0; n<=(nleft-1); n++)
	{
	  *(frontiers+depth*maxnfrontiers*2+n*2+0)=*(left+n*2+0);
	  *(frontiers+depth*maxnfrontiers*2+n*2+1)=*(left+n*2+1);
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+0)=val1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+1)=val2;
	}
      free(left);
    }
  

  /*
  else
    {
      int nleft;
      double *left;      
      nleft=0; left=malloc(sizeof(double)*maxnfrontiers*2);
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  if ((val1<=*(frontiers+depth*maxnfrontiers*2+n*2+0))||(val2<=*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	    {
	      nleft++; 
	      *(left+(nleft-1)*2+0)=*(frontiers+depth*maxnfrontiers*2+n*2+0);
	      *(left+(nleft-1)*2+1)=*(frontiers+depth*maxnfrontiers*2+n*2+1);
	    }
	}
      for (n=0; n<=(nleft-1); n++)
	{
	  *(frontiers+depth*maxnfrontiers*2+n*2+0)=*(left+n*2+0);
	  *(frontiers+depth*maxnfrontiers*2+n*2+1)=*(left+n*2+1);
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+0)=val1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+1)=val2;
	}
      free(left);
    }
  */

  
  // Debug
  //printf("depth=%d, nfrontiers=%d\n",depth,*(nfrontiers+depth));
  //for (i=0; i<=(*(nfrontiers+depth)-1); i++)
  //printf("%.3e %.3e\n",*(frontiers+depth*maxnfrontiers*2+i*2+0),*(frontiers+depth*maxnfrontiers*2+i*2+1));
  

  
  // Spin off each dosage combination.
  seqnodes=realloc(seqnodes,sizeof(struct seqnode)*(nseqnodes+nstates));
  vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);

  // Debug
  //printf("curind=%d, nseqnodes=%d, nstates=%d, nn=%d\n",curind,nseqnodes,nstates,nseqnodes+nstates);


  (seqnodes+curind)->nchildren=nstates;
  (seqnodes+curind)->children=malloc(sizeof(long)*nstates);
  for (n=0; n<=(nstates-1); n++)
    {
      *((seqnodes+curind)->children+n)=nseqnodes+n;
      (seqnodes+nseqnodes+n)->index=nseqnodes+n;
      (seqnodes+nseqnodes+n)->depth=(seqnodes+curind)->depth+1;
      (seqnodes+nseqnodes+n)->pa=curind;
      (seqnodes+nseqnodes+n)->nchildren=0;
      (seqnodes+nseqnodes+n)->children=malloc(sizeof(long));
      (seqnodes+nseqnodes+n)->dosage=malloc(sizeof(double)*ndrugs);
      (seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ndrugs-1); i++)
	*(vec+i)=*(dstates+n*ndrugs+i);      
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);
     
      recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
           
      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=recurse_sequence_multidrug_responses2(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind2);
    }

  return seqnodes;
}



// Estimate the upper bound of death time given the vector of initial population.
// Impose max dosage to all drugs.
// Use simple matrix product to estimate the population at time t.
// Apply binary search to find the time step prior to death.
// Difference from deathtime: apply nstates drug combinations and choose the max survival time among them.
// Difference from deathtime2: apply recurse_two_drug_response_trimrates to estimate population size.
// Using matrix exponential will amplify the effect of R12 fractional population.
// Difference from deathtime3: multi-drug version.
double deathtime_multidrug(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates)
{
  int i, j, k, l, m, n, cnt, lwd, upp, mid, ind, tind, maxind=-1;
  double tspan, *vec;
  double stopT=-1, val, val2, *dstates, *y;
  int *resvec, combind, *bases;
  
  
  // If each subpopulation size < 1, then death time is max time interval.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val=max(val,*(x0+i));
  if (val<1)
    {
      stopT=(nintervals-1)*timeinterval;
      return stopT;
    }

  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }  
  y=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs);

  // For each drug combination, find the survival time.
  for (tind=0; tind<=(nstates-1); tind++)
    {
      n=0;   
           
      // Fix the drug combination to the designated one throughout the whole episode.
      for (i=0; i<=(ndrugs-1); i++)
	*(vec+i)=*(dstates+tind*ndrugs+i);            
      
      // Apply binary search to find the time to death.
      // Use matrix exponential to estimate population.
      // Impose max dosage to all drugs.
      lwd=0; upp=nintervals-1; ind=-1;
      while ((lwd<upp)&&(ind<0))
	{
	  mid=(int)((lwd+upp)/2); tspan=timeinterval*mid;	  
	  recurse_multidrug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,tspan,vec,y);
	  	  
	  // Check if the total population >= mortal.
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(y+i);

	  /*
	  // Debug
	  printf("dosage ( ");
	  for (i=0; i<=(ndrugs-1); i++)
	    printf("%.2f ",*(vec+i));
	  printf("), time=%.2f, totalpop=%.2e, pop=( ",tspan,val);
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.2e ",*(y+i));
	  printf(")\n");
	  */

	  if (val>=mortal)
	    upp=mid;
	  else
	    lwd=mid;
	  
	  if (lwd==upp)
	    ind=mid;
	  else if (lwd==(upp-1))
	    ind=mid;
	}
    val=timeinterval*ind;
    
    /*
    // Debug
    printf("dosage=( ");
    for (i=0; i<=(ndrugs-1); i++)
      printf("%.2f ",*(vec+i));
    printf("), deathtime=%f\n",val);
    */


    // Debug
    if (val>stopT)
      maxind=tind;

    stopT=max(val,stopT);
    
   }

  free(vec); free(dstates); free(y);

  return stopT;
}


// Convert a number into a vector.
void num2vec(long index, int ndigits, int *bases, int *numvec)
{
  int i=ndigits-1, j;
  long residual=index;

  while (i>=0)
    {
      j=residual%(*(bases+i));
      *(numvec+i)=j;
      residual=(residual-j)/(*(bases+i));
      i--;
    }
  return;
}

		     
// Convert a vector into a number.
// Read the digits from end to beginning.
long vec2num(int ndigits, int *bases, int *numvec)
{
  int i=ndigits-1;
  long index=0;
  
  
  while (i>=0)
    {
      index=*(bases+i)*index+*(numvec+i);      
      i--;
    }
  

  /*
  i=0; index=0;
  while (i<=(ndigits-1))
    {
      index=index+*(numvec+i);
      if (i<(ndigits-1))
	index=*(bases+i+1)*index;
      i++;
    }
  */      


  return index;
}


// Convert a vector into a number.
// Read the digits from beginning to end.
long vec2num2(int ndigits, int *bases, int *numvec)
{
  int i=ndigits-1;
  long index=0;
  
  /*
  while (i>=0)
    {
      index=*(bases+i)*index+*(numvec+i);      
      i--;
    }
  */

  
  i=0; index=0;
  while (i<=(ndigits-1))
    {
      index=index+*(numvec+i);
      if (i<(ndigits-1))
	index=*(bases+i+1)*index;
      i++;
    }        

  return index;
}



// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
// Difference from dp_optimize_two_drug_responses4: apply branch and bound to trim unnecessary nodes in recursion.
// Difference from dp_optimize_two_drug_responses5: (1)multi-drug version, (2)make the format of y and dosage compatible with optimize_multidrug_responses1, (2)can choose one of multiple strategies as baseline for dp optimization.
// mode=0: monotherapy switching, benchmark.
// mode=1: minimize total population.
// mode=2: select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
// mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
// mode=4: maximize the anticipated death time for the best monotonous treatment.
// Difference from dp_optimize_multidrug_responses1: modify the criteria to improve strategy 2.  Ensure the predicted total population in multiple steps does not exceed the threshold.  If so then minimize the total population.
double dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, *resvec, combind, *bases;
  long curind, *seq;
  double *dstates, val, stopT, *frontiers, pop, incurablepop;
  struct seqnode *seqnodes;  
  
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }
  

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }

  /*
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+n*ndrugs+i)=-1;
    }    

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);
    */

  seq=malloc(sizeof(long)*nsteps); n=0;
  nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

      // Debug
      //printf("n=%d\n",n);

      // (n%nsteps)=0: incur dynamic programming.
      if ((n%nsteps)==0)
	{
	  if (n>0)
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  free((seqnodes+i)->children);
		  free((seqnodes+i)->dosage);
		  free((seqnodes+i)->y);
		}
	      free(seqnodes);	     	      
	    }
	  nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
	  (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
	  (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
	  (seqnodes+0)->children=malloc(sizeof(long));
	  (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
	  (seqnodes+0)->y=malloc(sizeof(double)*ntypes);
	  if (n==0)
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+i);
	      
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+0);
	    }
	  else
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);

	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+n-1);
	    }

	  for (i=0; i<=nsteps; i++)
	    *(nfrontiers+i)=0;

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  	  	 
	  // Find the max depth of the seqnodes.
	  maxdepth=-1;
	  for (i=0; i<=(nseqnodes-1); i++)
	    maxdepth=max(maxdepth,(seqnodes+i)->depth);
	  	  
	  /*
	  // Debug
	  //if ((mode==2)&&(popthre>1e10)&&(n==5))
	  if ((mode==2)&&(popthre>1e10)&&(n==0))
	    {	      
	      int ii;
	      double val;

	      printf("level 1 frontiers: ");
	      for (i=0; i<=(*(nfrontiers+1)-1); i++)
		printf("(%.2e %.2e) ",*(frontiers+1*maxnfrontiers*2+i*2+0),*(frontiers+1*maxnfrontiers*2+i*2+1));
	      printf("\n");

	      printf("nseqnodes=%d\n",nseqnodes);
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  printf("ind %d, pa %d, depth %d, dosage ( ",i,(seqnodes+i)->pa,(seqnodes+i)->depth);
		  for (j=0; j<=(ndrugs-1); j++)
		    printf("%.2f ",*((seqnodes+i)->dosage+j));
		  val=0;
		  printf("), y ( ");
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      printf("%.2e ",*((seqnodes+i)->y+j)); val+=*((seqnodes+i)->y+j);
		    }
		  printf(") %.2e, children ",val);
		  for (j=0; j<=((seqnodes+i)->nchildren-1); j++)
		    printf("%d ",*((seqnodes+i)->children+j));
		  printf("\n");
		}


	      printf("trace seqnodes\n");
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  if ((seqnodes+i)->depth>=maxdepth)
		    {
		      ii=i;
		      printf("node %d ",ii);
		      while ((seqnodes+ii)->pa>=0)
			{
			  printf("( ");
			  for (j=0; j<=(ndrugs-1); j++)
			    printf("%.2f ",*((seqnodes+ii)->dosage+j));
			  printf(") ");
			  ii=(seqnodes+ii)->pa;
			}
		      printf("\n");
		    }
		}
	      exit(0);
	    }
	  */



	  // Find the candidate terminal seqnodes.
	  ncandidates=0; candidates=malloc(sizeof(long));
	  
	  // Find the terminal nodes with zero population.
	  // If so then pick up the one that minimizes the total population.
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      if ((seqnodes+i)->nchildren==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  if (val<1)
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	    }	  
	  
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);
	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  
	  
	  // If the tumor cannot be eradicated within nsteps, then find the candidate sequences.
	  // A candidate sequence has the max depth, terminates at a leaf, and is not inferior to any other sequences.
	  // A sequence s1 is inferior to another sequence s2 if the total population and R12 population of s1 are both greater than those of s2.
	  if (ncandidates==0)
	    {
	      struct pair2 *pairs1, *pairs2;
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates, allhightotal, allmortal;
	      double p1, p2, p, *deathtimes;

	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0;
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  if ((c1==1)&&(c2==1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }		  
		}	     

	      // Sort candidate sequences by total population and by incurable population.
	      pairs1=malloc(sizeof(struct pair2)*ncandidates);
	      pairs2=malloc(sizeof(struct pair2)*ncandidates);
	      lookup1=malloc(sizeof(int)*ncandidates);
	      lookup2=malloc(sizeof(int)*ncandidates);
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  (pairs1+i)->i=i;
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  (pairs1+i)->s=val;
		  (pairs2+i)->i=i;
		  val=*((seqnodes+*(candidates+i))->y+ntypes-1);
		  (pairs2+i)->s=val;
		}

	      qsort(pairs1,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      qsort(pairs2,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  j=(pairs1+i)->i; *(lookup1+j)=i;
		  j=(pairs2+i)->i; *(lookup2+j)=i;
		}


	      /*
	      // Debug
	      //if ((mode==2)&(popthre>1e10)&&(n==1))
	      //if ((mode==2)&(popthre>1e10)&&(n==2))
	      //if ((mode==2)&(popthre>1e10)&&(n==5))
	      if ((mode==2)&(popthre>1e10)&&(n==0))
		{
		  printf("ncandidates=%d, nseqnodes=%d\n",ncandidates,nseqnodes);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      int ii;

		      printf("%d %d dosage=( ",i,*(candidates+i));
		      for (j=0; j<=(ndrugs-1); j++)
			printf("%.2f ",*((seqnodes+*(candidates+i))->dosage+j));
		      printf("), ");		      
		      printf("y=( ");
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  printf("%.2e ",*((seqnodes+*(candidates+i))->y+j));
			  val+=*((seqnodes+*(candidates+i))->y+j);
			}
		      printf(") %.2e\n",val);

		      printf("history ");
		      ii=*(candidates+i);
		      while ((seqnodes+ii)->pa>=0)
			{
			  printf("( ");
			  for (j=0; j<=(ndrugs-1); j++)
			    printf("%.2f ",*((seqnodes+ii)->dosage+j));
			  printf(") ");
			  ii=(seqnodes+ii)->pa;
			}
		      printf("\n");

		    }
		}		      
	      */


	      // Find the candidates inferior to others.
	      // Direct comparison of all candidate pairs take n^2 time, n=ncandidates.
	      // To reduce computing time use probabilistic arguments.
	      // Suppose a candidate's total population > n1 other candidates, and its incurable population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and incurable population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in incurable population) =
	      // 1-Pr(there is no overlap ...).
	      // Pr(there is no overlap ...) = 
	      // Pr(place n2 candidates in (n-n1) among n slots) or
	      // Pr(place n1 candidates in (n-n2) among n slots).
	      // If (n1+n2)>=n, then Pr(there is an overlap ...)=1.
	      // If (n1+n2)<n, then
	      // P1=((n-n1)/n)^n2, P2=((n-n2)/n)^n1. 
	      // Pr(the candidate is inferior) = 1 - max(P1,P2).
	      inferior=malloc(sizeof(int)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{		  
		  *(inferior+i)=0; p=0;
		  j=*(lookup1+i); n1=max(j-1,0);
		  j=*(lookup2+i); n2=max(j-1,0);
		  if ((n1+n2)>=ncandidates)
		    p=1;
		  else
		    {
		      val=(double)(ncandidates-n1)/(double)ncandidates;
		      p1=n2*log(val); p1=exp(p1);
		      val=(double)(ncandidates-n2)/(double)ncandidates;
		      p2=n1*log(val); p2=exp(p2);
		      p=1-max(p1,p2);
		    }

		  if (p>=(1-1e-4))
		    *(inferior+i)=1;
		  else
		    {
		      // Compare the superior list in two values.
		      j=0; k=-1;
		      while ((j<=n1)&&(k<0))
			{
			  l=(pairs1+j)->i; l=*(lookup2+l);
			  if (((pairs1+*(lookup1+i))->s>(pairs1+j)->s)&&((pairs2+*(lookup2+i))->s>(pairs2+l)->s))
			    //if (l<=n2)
			    k=l;
			  j++;
			}
		      if (k>=0)
			*(inferior+i)=1;
		    }
		}

	      // Check whether all candidates lead to mortality.
	      allmortal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allmortal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<mortal)
		    allmortal=0;
		  i++;
		}

	      // Check whether all candidates lead to high total population.
	      allhightotal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allhightotal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<popthre)
		    allhightotal=0;
		  i++;
		}
	      
	      // Debug
	      //printf("ncandidates=%d, before\n",ncandidates);

	      // Choose a candidate according to the selected strategy.

	      // mode=1: Minimize the total population. 
	      if (mode==1)
		{		  
		  minind=-1; minval=1e100;	      
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			val+=*((seqnodes+*(candidates+i))->y+j);
		      if (val<minval)
			{
			  minind=*(candidates+i); minval=val;
			}
		    }	      
		}	      

	      // mode=2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==2)
		{		
		  minind=-1; minval=1e20;		  		  
		  
		  // Total population is above the threshold: minimize total population.
		  // Total population is above the threshold or all drug combinations lead to mortality: minimize total population.
		  //if (pop>=popthre)	      		  
		  if ((pop>=popthre)||(allhightotal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }		  	     		
	      
		  // Total population is below the threshold: minimize the risk of developing incurable cells.
		  // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
		  // prob((1,...1) \rightarrow (1,...,1))=1.
		  // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
		  // Discard the combinations that lead to mortality.
		  // Discard the combinations that have high total population.
		  else
		    {		  		      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  val*=*(T+(2^k)*ntypes+0);
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);				      
			    }

			  //if (risk<minval)
			  //if ((risk<minval)&&(val2<mortal))
			  if ((risk<minval)&&(val2<popthre))
			    {			      
			      minind=*(candidates+i); minval=risk;

			      // Debug
			      //if ((popthre>1e10)&&(n==5))
			      //printf("minind=%d, minval=%.2e, val2=%.2e, mortal=%.2e\n",minind,minval,val2,mortal);

			    }			  			  
			}		    
		    }
		}	
	    
	      // mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==3)
		{
		  minval=1e20;

		  // If incurable cells > 1, then minimize total population.
		  //if (incurablepop>=(1+1e-10))
		  if ((incurablepop>=(1+1e-10))||(allmortal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }	

		  // Otherwise minimize the predicted total multiple resistant populations.
		  // Discard the combinations that lead to mortality.
		  else
		    {
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val2;

			  val=0; val2=0;
			  for (j=1; j<=(ntypes-1); j++)
			    {
			      num2vec(j,ndrugs,bases,resvec); k=0;
			      for (l=0; l<=(ndrugs-1); l++)
				k+=*(resvec+l);
			      if (k>1)
				val+=*((seqnodes+*(candidates+i))->y+j);
			      val2+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  //if (val<minval)
			  if ((val<minval)&&(val2<mortal))
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }
		}

	      // mode=4: Maximize the anticipated death time based on static therapy.
	      else if (mode==4)
		{
	      	      	     	      
		  // Evaluate the expected death time for each non-inferior candidate.
		  // Apply the best monotherapy.
		  deathtimes=malloc(sizeof(double)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      *(deathtimes+i)=-1;
		      if (*(inferior+i)==0)
			{
			  val=deathtime_multidrug((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
			  
			  *(deathtimes+i)=val;		      		      
			}
		    }	
		  
		  // Debug
		  //printf("ncandidates=%d, after\n",ncandidates);
		  
		  
		  // Choose the candidate with the max death time.
		  minind=-1; minval=-1;
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if (*(deathtimes+i)>minval)
			{
			  minind=*(candidates+i); minval=*(deathtimes+i);
			}
		    }
		  
		  // If there are degenerate death times, then choose the one with minimum incurable population.
		  ndegenerates=0; degenerates=malloc(sizeof(int)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if ((minval-*(deathtimes+i))<=1e-5)
			{
			  ndegenerates++; *(degenerates+ndegenerates-1)=*(candidates+i);
			}
		    }
		  if (ndegenerates>1)
		    {
		      minind=-1; minval=1e100;
		      for (i=0; i<=(ndegenerates-1); i++)
			{		      		      
			  if (degmode==0)
			    {
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				val+=*((seqnodes+*(degenerates+i))->y+j);
			    }
			  else
			    val=*((seqnodes+*(degenerates+i))->y+ntypes-1);
			  

			  if (val<minval)
			    {
			      minind=*(degenerates+i); minval=val;
			    }
			}
		    }
		  free(degenerates); free(deathtimes);
		}

	 	      	      	    	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); 
	      
	    }
	     
	  free(candidates);
	  
	  for (i=0; i<=(nsteps-1); i++)
	    *(seq+i)=-1;

	  if (minind>=0)
	    {
	      k=(seqnodes+minind)->depth-1;
	      *(seq+k)=minind; i=k-1;
	      //*(seq+nsteps-1)=minind; i=nsteps-2;
	      while (i>=0)
		{
		  *(seq+i)=(seqnodes+*(seq+i+1))->pa;
		  i--;
		}	  
	    }	  	 
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);

      // Debug
      //printf("n=%d, k=%d, l=%d\n",n,k,l);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  /*
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
	    */

	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j*nintervals+n)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+i*nintervals+n)=*((seqnodes+l)->y+i);



	  maxval=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //maxval+=*(y+n*ntypes+i);
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+i*nintervals+n);
	
	  pop=maxval; incurablepop=*(y+(ntypes-1)*nintervals+n);

	  // Debug
	  //printf("maxval=%.2e\n",maxval);

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);	  

	  
	  /*
	  // Debug
	  //if ((mode==2)&&(popthre>1e10))
	  //if (mode==1)
	  if (mode==3)
	    {
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		val+=*(y+i*nintervals+n);
	      printf("dp t(%d)=%.2f, dosage=( ",n,*(t+n));
	      for (i=0; i<=(ndrugs-1); i++)
		printf("%.2f ",*(dosage+i*nintervals+n));
	      printf("), y=( ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.4e ",*(y+i*nintervals+n));
	      printf(") %.4e\n",val);
	    }
	  */


	}

      // Debug
      //printf("n=%d, stopT=%.2f\n",n,stopT);

      n++;
    }

  if (stopT<0)
    stopT=*(t+nintervals-1);

  // Release memory.
  free(seq); free(dstates);
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); free(nfrontiers); free(frontiers);
  free(resvec); free(bases);
  
  return stopT;
}



// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
// Difference from dp_optimize_two_drug_responses4: apply branch and bound to trim unnecessary nodes in recursion.
// Difference from dp_optimize_two_drug_responses5: (1)multi-drug version, (2)make the format of y and dosage compatible with optimize_multidrug_responses1, (2)can choose one of multiple strategies as baseline for dp optimization.
// mode=0: monotherapy switching, benchmark.
// mode=1: minimize total population.
// mode=2: select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
// mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
// mode=4: maximize the anticipated death time for the best monotonous treatment.
// Difference from dp_optimize_multidrug_responses1: modify the criteria to improve strategy 2.  Ensure the predicted total population in multiple steps does not exceed the threshold.  If so then minimize the total population.
// Difference from dp_optimize_multidrug_responses2: modify recurse_sequence_multidrug_responses to create valid boundaries.
double dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, *resvec, combind, *bases;
  long curind, *seq;
  double *dstates, val, stopT, *frontiers, pop, incurablepop;
  struct seqnode *seqnodes;  
 
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }
  

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }

  /*
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+n*ndrugs+i)=-1;
    }    

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);
    */

  seq=malloc(sizeof(long)*nsteps); n=0;
  nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

      // Debug
      //printf("n=%d\n",n);

      // (n%nsteps)=0: incur dynamic programming.
      if ((n%nsteps)==0)
	{
	  if (n>0)
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  free((seqnodes+i)->children);
		  free((seqnodes+i)->dosage);
		  free((seqnodes+i)->y);
		}
	      free(seqnodes);	     	      
	    }
	  nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
	  (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
	  (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
	  (seqnodes+0)->children=malloc(sizeof(long));
	  (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
	  (seqnodes+0)->y=malloc(sizeof(double)*ntypes);
	  if (n==0)
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+i);
	      
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+0);
	    }
	  else
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);

	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+n-1);
	    }

	  for (i=0; i<=nsteps; i++)
	    *(nfrontiers+i)=0;

	  curind=0;	  	  	      	      	 
	  //seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  seqnodes=recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  	  	 
	  // Find the max depth of the seqnodes.
	  maxdepth=-1;
	  for (i=0; i<=(nseqnodes-1); i++)
	    maxdepth=max(maxdepth,(seqnodes+i)->depth);	  	  	  


	  // Find the candidate terminal seqnodes.
	  ncandidates=0; candidates=malloc(sizeof(long));
	  
	  // Find the terminal nodes with zero population.
	  // If so then pick up the one that minimizes the total population.
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      if ((seqnodes+i)->nchildren==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  if (val<1)
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	    }	  
	  
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);
	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  
	  
	  // If the tumor cannot be eradicated within nsteps, then find the candidate sequences.
	  // A candidate sequence has the max depth, terminates at a leaf, and is not inferior to any other sequences.
	  // A sequence s1 is inferior to another sequence s2 if the total population and R12 population of s1 are both greater than those of s2.
	  if (ncandidates==0)
	    {
	      struct pair2 *pairs1, *pairs2;
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates, allhightotal, allmortal;
	      double p1, p2, p, *deathtimes;

	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0;
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  if ((c1==1)&&(c2==1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }		  
		}	     

	      // Sort candidate sequences by total population and by incurable population.
	      pairs1=malloc(sizeof(struct pair2)*ncandidates);
	      pairs2=malloc(sizeof(struct pair2)*ncandidates);
	      lookup1=malloc(sizeof(int)*ncandidates);
	      lookup2=malloc(sizeof(int)*ncandidates);
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  (pairs1+i)->i=i;
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  (pairs1+i)->s=val;
		  (pairs2+i)->i=i;
		  val=*((seqnodes+*(candidates+i))->y+ntypes-1);
		  (pairs2+i)->s=val;
		}

	      qsort(pairs1,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      qsort(pairs2,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  j=(pairs1+i)->i; *(lookup1+j)=i;
		  j=(pairs2+i)->i; *(lookup2+j)=i;
		}


	      /*
	      // Debug
	      //if ((mode==2)&(popthre>1e10)&&(n==1))
	      //if ((mode==2)&(popthre>1e10)&&(n==2))
	      //if ((mode==2)&(popthre>1e10)&&(n==5))
	      if ((mode==2)&(popthre>1e10)&&(n==0))
		{
		  printf("ncandidates=%d, nseqnodes=%d\n",ncandidates,nseqnodes);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      int ii;

		      printf("%d %d dosage=( ",i,*(candidates+i));
		      for (j=0; j<=(ndrugs-1); j++)
			printf("%.2f ",*((seqnodes+*(candidates+i))->dosage+j));
		      printf("), ");		      
		      printf("y=( ");
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  printf("%.2e ",*((seqnodes+*(candidates+i))->y+j));
			  val+=*((seqnodes+*(candidates+i))->y+j);
			}
		      printf(") %.2e\n",val);

		      printf("history ");
		      ii=*(candidates+i);
		      while ((seqnodes+ii)->pa>=0)
			{
			  printf("( ");
			  for (j=0; j<=(ndrugs-1); j++)
			    printf("%.2f ",*((seqnodes+ii)->dosage+j));
			  printf(") ");
			  ii=(seqnodes+ii)->pa;
			}
		      printf("\n");

		    }
		}		      
	      */


	      // Find the candidates inferior to others.
	      // Direct comparison of all candidate pairs take n^2 time, n=ncandidates.
	      // To reduce computing time use probabilistic arguments.
	      // Suppose a candidate's total population > n1 other candidates, and its incurable population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and incurable population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in incurable population) =
	      // 1-Pr(there is no overlap ...).
	      // Pr(there is no overlap ...) = 
	      // Pr(place n2 candidates in (n-n1) among n slots) or
	      // Pr(place n1 candidates in (n-n2) among n slots).
	      // If (n1+n2)>=n, then Pr(there is an overlap ...)=1.
	      // If (n1+n2)<n, then
	      // P1=((n-n1)/n)^n2, P2=((n-n2)/n)^n1. 
	      // Pr(the candidate is inferior) = 1 - max(P1,P2).
	      inferior=malloc(sizeof(int)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{		  
		  *(inferior+i)=0; p=0;
		  j=*(lookup1+i); n1=max(j-1,0);
		  j=*(lookup2+i); n2=max(j-1,0);
		  if ((n1+n2)>=ncandidates)
		    p=1;
		  else
		    {
		      val=(double)(ncandidates-n1)/(double)ncandidates;
		      p1=n2*log(val); p1=exp(p1);
		      val=(double)(ncandidates-n2)/(double)ncandidates;
		      p2=n1*log(val); p2=exp(p2);
		      p=1-max(p1,p2);
		    }

		  if (p>=(1-1e-4))
		    *(inferior+i)=1;
		  else
		    {
		      // Compare the superior list in two values.
		      j=0; k=-1;
		      while ((j<=n1)&&(k<0))
			{
			  l=(pairs1+j)->i; l=*(lookup2+l);
			  if (((pairs1+*(lookup1+i))->s>(pairs1+j)->s)&&((pairs2+*(lookup2+i))->s>(pairs2+l)->s))
			    //if (l<=n2)
			    k=l;
			  j++;
			}
		      if (k>=0)
			*(inferior+i)=1;
		    }
		}

	      // Check whether all candidates lead to mortality.
	      allmortal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allmortal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<mortal)
		    allmortal=0;
		  i++;
		}

	      // Check whether all candidates lead to high total population.
	      allhightotal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allhightotal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<popthre)
		    allhightotal=0;
		  i++;
		}
	      
	      // Debug
	      //printf("ncandidates=%d, before\n",ncandidates);

	      // Choose a candidate according to the selected strategy.

	      // mode=1: Minimize the total population. 
	      if (mode==1)
		{		  
		  minind=-1; minval=1e100;	      
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			val+=*((seqnodes+*(candidates+i))->y+j);
		      if (val<minval)
			{
			  minind=*(candidates+i); minval=val;
			}
		    }	      
		}	      

	      // mode=2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==2)
		{		
		  minind=-1; minval=1e20;		  		 	  
		  
		  // Total population is above the threshold: minimize total population.
		  // Total population is above the threshold or all drug combinations lead to mortality: minimize total population.
		  //if (pop>=popthre)	      		  
		  if ((pop>=popthre)||(allhightotal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }		  	     		
	      
		  // Total population is below the threshold: minimize the risk of developing incurable cells.
		  // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
		  // prob((1,...1) \rightarrow (1,...,1))=1.
		  // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
		  // Discard the combinations that lead to mortality.
		  // Discard the combinations that have high total population.
		  else
		    {		  		      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  val*=*(T+(2^k)*ntypes+0);				  
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);	
			    }

			  //if (risk<minval)
			  //if ((risk<minval)&&(val2<mortal))
			  if ((risk<minval)&&(val2<popthre))
			    {			      
			      minind=*(candidates+i); minval=risk;
			    }		

			  // Debug
			  //printf("i=%d, cand=%d, risk=%.4e, val2=%.4e\n",i,*(candidates+i),risk,val2);
	  			  
			}		    
		    }

		  /*
		  // Debug
		  if ((popthre>1e10)&&(n==0))
		    {
		      printf("n=%d\n",n);
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  val*=*(T+(2^(ndrugs-1-k))*ntypes+0);
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);	

			      // Debug
			      if ((i==318)||(i==402))
				{
				  printf("( ");
				  for (k=0; k<=(ndrugs-1); k++)
				    printf("%d ",*(resvec+k));
				  printf("%.3e ",val);
				  
				  if (j==3)
				    {
				      printf("((%.3e ",*((seqnodes+*(candidates+i))->y+j));
				      for (k=0; k<=(ndrugs-1); k++)			    
					if (*(resvec+k)==0)
					  printf("%.3e ",*(T+(2^(ndrugs-1-k))*ntypes+0));
				      printf("))");
				    }

				  //for (k=0; k<=(ndrugs-1); k++)
				  //printf("k=%d, resvec=%d, transition=%.3e ",k,*(resvec+k),*(T+(2^k)*ntypes+0));
				  printf(")");
				}
			    }

			  //
			  //printf("candidate %d %d, ( ",i,*(candidates+i));
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  val=*((seqnodes+*(candidates+i))->y+j);
			  //  printf("%.4e ",val);
			  //}			  		  
			  //printf(")\n");
			  //

			  printf("candidate %d %d, risk=%.4e, ",i,*(candidates+i),risk);
			  k=*(candidates+i);

			  if ((k==620)||(k==767))
			    {
			      printf("( ");
			      for (j=0; j<=(ntypes-1); j++)
				{
				  val=*((seqnodes+*(candidates+i))->y+j);
				  printf("%.4e ",val);
				}
			      printf(")");
			    }		

			  while (k>=0)
			    {
			      printf("( ");
			      for (j=0; j<=(ndrugs-1); j++)
				{
				  val=*((seqnodes+k)->dosage+j);
				  printf("%.2f ",val);
				}
			      printf(")");
			      k=(seqnodes+k)->pa;
			    }
			  printf("\n");
			}
		      printf("minind=%d\n",minind);
		    }
		  */
		  

		}	
	    
	      // mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==3)
		{
		  minval=1e20;

		  // If incurable cells > 1, then minimize total population.
		  //if (incurablepop>=(1+1e-10))
		  if ((incurablepop>=(1+1e-10))||(allmortal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }	

		  // Otherwise minimize the predicted total multiple resistant populations.
		  // Discard the combinations that lead to mortality.
		  else
		    {
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val2;

			  val=0; val2=0;
			  for (j=1; j<=(ntypes-1); j++)
			    {
			      num2vec(j,ndrugs,bases,resvec); k=0;
			      for (l=0; l<=(ndrugs-1); l++)
				k+=*(resvec+l);
			      if (k>1)
				val+=*((seqnodes+*(candidates+i))->y+j);
			      val2+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  //if (val<minval)
			  if ((val<minval)&&(val2<mortal))
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }
		}

	      // mode=4: Maximize the anticipated death time based on static therapy.
	      else if (mode==4)
		{
	      	      	     	      
		  // Evaluate the expected death time for each non-inferior candidate.
		  // Apply the best monotherapy.
		  deathtimes=malloc(sizeof(double)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      *(deathtimes+i)=-1;
		      if (*(inferior+i)==0)
			{
			  val=deathtime_multidrug((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
			  
			  *(deathtimes+i)=val;		      		      
			}
		    }	
		  
		  // Debug
		  //printf("ncandidates=%d, after\n",ncandidates);
		  
		  
		  // Choose the candidate with the max death time.
		  minind=-1; minval=-1;
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if (*(deathtimes+i)>minval)
			{
			  minind=*(candidates+i); minval=*(deathtimes+i);
			}
		    }
		  
		  // If there are degenerate death times, then choose the one with minimum incurable population.
		  ndegenerates=0; degenerates=malloc(sizeof(int)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if ((minval-*(deathtimes+i))<=1e-5)
			{
			  ndegenerates++; *(degenerates+ndegenerates-1)=*(candidates+i);
			}
		    }
		  if (ndegenerates>1)
		    {
		      minind=-1; minval=1e100;
		      for (i=0; i<=(ndegenerates-1); i++)
			{		      		      
			  if (degmode==0)
			    {
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				val+=*((seqnodes+*(degenerates+i))->y+j);
			    }
			  else
			    val=*((seqnodes+*(degenerates+i))->y+ntypes-1);
			  

			  if (val<minval)
			    {
			      minind=*(degenerates+i); minval=val;
			    }
			}
		    }
		  free(degenerates); free(deathtimes);
		}

	 	      	      	    	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); 
	      
	    }
	     
	  free(candidates);
	  
	  for (i=0; i<=(nsteps-1); i++)
	    *(seq+i)=-1;

	  if (minind>=0)
	    {
	      k=(seqnodes+minind)->depth-1;
	      *(seq+k)=minind; i=k-1;
	      //*(seq+nsteps-1)=minind; i=nsteps-2;
	      while (i>=0)
		{
		  *(seq+i)=(seqnodes+*(seq+i+1))->pa;
		  i--;
		}	  
	    }	  	 
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);

      // Debug
      //printf("n=%d, k=%d, l=%d\n",n,k,l);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  /*
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
	    */

	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j*nintervals+n)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+i*nintervals+n)=*((seqnodes+l)->y+i);



	  maxval=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //maxval+=*(y+n*ntypes+i);
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+i*nintervals+n);
	
	  pop=maxval; incurablepop=*(y+(ntypes-1)*nintervals+n);

	  // Debug
	  //printf("maxval=%.2e\n",maxval);

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  // Report the survival time to the max time span + timeinterval.
	  // This distinguishes from the cases where the patients survival for the max time span but are not cured.
	  //if (maxval<=0)	
	  if (maxval<1)	    
	    //stopT=*(t+nintervals-1);	  
	    stopT=*(t+nintervals-1)+timeinterval;

	  
	  /*
	  // Debug	  
	  if ((mode==2)&&(popthre>1e10))
	  //if (mode==1)
	  //if (mode==3)
	    {
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		val+=*(y+i*nintervals+n);
	      printf("dp t(%d)=%.2f, dosage=( ",n,*(t+n));
	      for (i=0; i<=(ndrugs-1); i++)
		printf("%.2f ",*(dosage+i*nintervals+n));
	      printf("), y=( ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.4e ",*(y+i*nintervals+n));
	      printf(") %.4e\n",val);
	    }
	  */
	  

	}

      // Debug
      //printf("n=%d, stopT=%.2f\n",n,stopT);

      n++;
    }

  /*
  // If stop<0, then decide whether the patient recovers or not.
  // Check if each final subpopulation < 1.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      k=-1; i=0;
      while ((i<=(ntypes-1))&&(k<0))
	{
	  //if (*(y+i*nintervals+nintervals-2)>1)
	  if (*(y+i*nintervals+nintervals-2)>10)
	    k=i;
	  i++;
	}
      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }
  */

  //if (stopT<0)
  //stopT=*(t+nintervals-1);
    
  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *localdosage, *inity, *finaly, maxtime=nintervals*timeinterval;
      localdosage=malloc(sizeof(double)*ndrugs);
      finaly=malloc(sizeof(double)*ntypes);
      inity=malloc(sizeof(double)*ntypes);
      
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(localdosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,localdosage,finaly,0);
	  flag=1; m=0;
	  while ((m<=(ntypes-1))&&(flag==1))
	    {
	      if (*(finaly+m)>1)
		flag=0;
	      m++;
	    }
	  l++;
	}
      if (flag==1)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
      free(localdosage); free(inity); free(finaly);
    }


  // Release memory.
  free(seq); free(dstates);
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); free(nfrontiers); free(frontiers);
  free(resvec); free(bases);
  
  return stopT;
}



// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
// Difference from recurse_sequence_two_drug_responses3: multi-drug version.
// Difference from recurse_sequence_multidrug_responses: remove the frontiers which are no longer frontiers.
struct seqnode * deb_recurse_sequence_multidrug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt)
{
  int i, j, k, l, m, n, inferior, depth;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2, val1, val2;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
    
  // Debug
  //printf("step 1\n");


  // Return when the node depth>=nsteps.
  if ((seqnodes+curind)->depth>=nsteps)
    return seqnodes;

  // Return when the total population exceeds the mortal threshold.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val+=*((seqnodes+curind)->y+i);
  
 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;

  // Debug
  //printf("step 2\n");



  // Return when the total and incurable populations of the current node are inferior to the frontiers at the current level.
  depth=(seqnodes+curind)->depth; inferior=0; n=0; 
  val1=0;
  for (i=0; i<=(ntypes-1); i++)
    val1+=*((seqnodes+curind)->y+i);
  val2=*((seqnodes+curind)->y+ntypes-1); 
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // Check whether both total and incurable populations are larger than a frontier.           
      if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))&&(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	inferior=1;
      n++;
    }

  // Debug
  //printf("step 3\n");


  if (inferior==1)
    return seqnodes;
  
  // Augment/update the frontier points if the current node is not inferior to any of them.
  // Also remove the frontier points inferior to the current node.
  else
    {
      int nleft;
      double *left;      
      nleft=0; left=malloc(sizeof(double)*maxnfrontiers*2);
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))||(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	    {
	      nleft++; 
	      *(left+(nleft-1)*2+0)=*(frontiers+depth*maxnfrontiers*2+n*2+0);
	      *(left+(nleft-1)*2+1)=*(frontiers+depth*maxnfrontiers*2+n*2+1);
	    }
	}
      for (n=0; n<=(nleft-1); n++)
	{
	  *(frontiers+depth*maxnfrontiers*2+n*2+0)=*(left+n*2+0);
	  *(frontiers+depth*maxnfrontiers*2+n*2+1)=*(left+n*2+1);
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+0)=val1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+1)=val2;
	}
      free(left);
    }
  

  /*
  else
    {
      int nleft;
      double *left;      
      nleft=0; left=malloc(sizeof(double)*maxnfrontiers*2);
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  if ((val1<=*(frontiers+depth*maxnfrontiers*2+n*2+0))||(val2<=*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	    {
	      nleft++; 
	      *(left+(nleft-1)*2+0)=*(frontiers+depth*maxnfrontiers*2+n*2+0);
	      *(left+(nleft-1)*2+1)=*(frontiers+depth*maxnfrontiers*2+n*2+1);
	    }
	}
      for (n=0; n<=(nleft-1); n++)
	{
	  *(frontiers+depth*maxnfrontiers*2+n*2+0)=*(left+n*2+0);
	  *(frontiers+depth*maxnfrontiers*2+n*2+1)=*(left+n*2+1);
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+0)=val1;
	  *(frontiers+depth*maxnfrontiers*2+(*(nfrontiers+depth)-1)*2+1)=val2;
	}
      free(left);
    }
  */

  
  // Debug
  //printf("depth=%d, nfrontiers=%d\n",depth,*(nfrontiers+depth));
  //for (i=0; i<=(*(nfrontiers+depth)-1); i++)
  //printf("%.3e %.3e\n",*(frontiers+depth*maxnfrontiers*2+i*2+0),*(frontiers+depth*maxnfrontiers*2+i*2+1));
  

  // Debug
  //printf("step 4\n");


  
  // Spin off each dosage combination.
  seqnodes=realloc(seqnodes,sizeof(struct seqnode)*(nseqnodes+nstates));
  vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);


  (seqnodes+curind)->nchildren=nstates;
  (seqnodes+curind)->children=malloc(sizeof(long)*nstates);
  for (n=0; n<=(nstates-1); n++)
    {
      *((seqnodes+curind)->children+n)=nseqnodes+n;
      (seqnodes+nseqnodes+n)->index=nseqnodes+n;
      (seqnodes+nseqnodes+n)->depth=(seqnodes+curind)->depth+1;
      (seqnodes+nseqnodes+n)->pa=curind;
      (seqnodes+nseqnodes+n)->nchildren=0;
      (seqnodes+nseqnodes+n)->children=malloc(sizeof(long));
      (seqnodes+nseqnodes+n)->dosage=malloc(sizeof(double)*ndrugs);
      (seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ndrugs-1); i++)
	*(vec+i)=*(dstates+n*ndrugs+i);      
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);
     
      recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
           
      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=deb_recurse_sequence_multidrug_responses2(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind2);
    }

  // Debug
  //printf("step 5\n");


  return seqnodes;
}



// Implement optimization heuristics for multiple drugs.
// Adapt optimize_two_drug_responses8.
// Difference from optimize_two_drug_responses8: exercise fewer possible drug combinations.
// For k drugs, exhaust combinations x1+...+xk=1.
// For all the xk!=0, equalize their dosages.
// There are totally 2^{k}-1 drug combinations.
// Difference from optimize_multidrug_responses1: Alter strategy 0.  Rank strategies by their rates of reducing the dominant subclone.
double optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage)
{
  int i, j, k, l, m, n, cnt, *bases, combind, *resvec, nstates;
  double *d, *mininds, *mininds2, *mininds3, minval, minval2, minval3, *mask, maxval;
  double *F, *template, *constF, *y0, *vec, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, val2, pop, *drugrecords, *dstates;

  d=malloc(sizeof(double)*ndrugs); 
  mininds=malloc(sizeof(double)*ndrugs);
  mininds2=malloc(sizeof(double)*ndrugs);
  mininds3=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
  mask=malloc(sizeof(double)*ntypes*ntypes);
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (i=0; i<=(ndrugs-1); i++)
    *(bases+i)=2;  
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
 
  drugrecords=malloc(sizeof(double)*(ndrugs+1+ntypes)*ntypes);
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(y0+i)=*(x0+i);
	      *(y+i*nintervals+n)=*(y0+i);
	    }
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);


      // Strategy 0 does not need prediction.      
      if (mode==0)
	{	  
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;	  

	  // Initially find the drug that yields the highest rate of reducing the dominant population.
	  // The growth rate of type i is g0(1+\sum_{j!=i}T(j->i))-\sum_{k}Sa(i,k)dk.
	  if (n==0)
	    {
	      int dominant=-1;
	      double rates[ndrugs];	      
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  if (*(y0+i)>val)
		    {
		      dominant=i; val=*(y0+i);
		    }
		}
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  rates[i]=*(g0+dominant);
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      if (j!=dominant)
			{
			  val=*(T+dominant*ntypes+j)*(*(g0+dominant));
			  rates[i]+=val;
			}
		    }
		  val=*(Sa+dominant*ndrugs+i);
		  rates[i]-=val;
		}
	      k=-1; val=1e500;
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  if (rates[i]<val)
		    {
		      k=i; val=rates[i];
		    }
		}
	       for (i=0; i<=(ndrugs-1); i++)
		*(mininds+i)=0;
	      *(mininds+k)=1;
	    }

	  // Switch the drug when the total population reaches the double the nadir or the population reemerges.
	  else
	    {
	      int switching=0;
	      double nadir=1e50;
	      
	      // Find the nadir of the current regimen.
	      // A nadir is the local minimum of the population under the current regimen.
	      // Start from the current time step, trace back until the change of regimen.
	      i=n-1; k=1;
	      while ((i>=0)&&(k==1))
		{
		  int p, q;

		  // Check whether the dosages at time step i equals the dosages at time step n-1.
		  p=0; q=1;
		  while ((p<=(ndrugs-1))&&(q==1))
		    {
		      if (fabs(*(dosage+p*nintervals+i)-*(dosage+p*nintervals+n-1))>1e-12)
			q=0;
		      p++;
		    }
		  // If so then update nadir.
		  if (q==1)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  if (i==0)
			    val+=*(x0+j);
			  else
			    val+=*(y+j*nintervals+i);
			}
		      nadir=min(nadir,val);
		    }
		  else
		    k=0;
		  i--;
		}
     
	      // If the current population >= 2 x nadir population, then switch.
	      if (pop>=(nadir*2))
		switching=1;

	      // If the current population reemerges from the level below the detection threshold, then switch.
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  if ((n-1)==0)
		    val+=*(x0+i);
		  else
		    val+=*(y+i*nintervals+n-1);
		}

	      if ((val<detected)&&(pop>=detected))
		switching=1;
	      
	      // Find the most effective drug other than the current choice.
	      if (switching==1)
		{
		  int dominant=-1;
		  double rates[ndrugs];	      
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      if (*(y0+i)>val)
			{
			  dominant=i; val=*(y0+i);
			}
		    }
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      rates[i]=*(g0+dominant);
		      for (j=0; j<=(ntypes-1); j++)
			{
			  if (j!=dominant)
			    {
			      val=*(T+dominant*ntypes+j)*(*(g0+dominant));
			      rates[i]+=val;
			    }
			}
		      val=*(Sa+dominant*ndrugs+i);
		      rates[i]-=val;
		    }
		  j=-1; maxval=0;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if (*(dosage+i*nintervals+n-1)>maxval)
			{
			  maxval=*(dosage+i*nintervals+n-1); j=i;
			}
		    }
		  k=-1; val=1e500;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if ((i!=j)&&(rates[i]<val))
			{
			  k=i; val=rates[i];
			}
		    }
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=0;
		  *(mininds+k)=1;
		}
	      else
		{
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=*(dosage+i*nintervals+n-1);
		}
	    }	  
	  
	  for (i=0; i<=(ntypes-1); i++)
	    *(vec2+i)=*(y0+i);
	}

      // For other strategies, calculate the predictions.
      else
	{
	  // Mask matrix.
	  // Mask the terms involved in the initial populations below 1.
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ntypes-1); j++)
		{
		  if (*(y0+j)<1)
		    *(mask+i*ntypes+j)=0;
		  else
		    *(mask+i*ntypes+j)=1;
		}
	    }   	  
	  
	  // Change the drug combinations to optimize certain criterion.
	  // Here the predicted responses are within the designated time interval (predictionperiod) instead of the sampling time interval.

	  // mininds2 and minval2: the drug dosage that minimizes the total population.
	  // mininds3 and minval3: the drug dosage that minimizes the R12 population.

	  // Consider the drug combinations (d1,d2,...,dk), where di's are either 0 or 1/m, where m is the number of nonzero di's.
	  // Exhaust all combinations of drugs except (0,0,...0).
	  
	  for (i=0; i<=(ndrugs-1); i++)
	    {
	      *(mininds+i)=0; *(mininds2+i)=0; *(mininds3+i)=0;
	    }

	  minval=1e100; minval2=1e100; minval3=1e100; cnt=0; 
	  for (combind=1; combind<=(ntypes-1); combind++)
	    {	     
	      num2vec(combind,ndrugs,bases,resvec);
	      k=0;
	      for (i=0; i<=(ndrugs-1); i++)	       
		k+=*(resvec+i);
	      val=(double)1/(double)k;
	      
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  *(vec+i)=*(resvec+i)*val;
		}
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(y0+i);
	  	  
		  
	      // compmode=0: use F to estimate the population.
	      if (compmode==0)
		{		  
		  mat_multiply(ntypes,ndrugs,1,Sg,vec,vec2);
		  diagonal(ntypes,vec2,template);	      
		  for (i=0; i<=(ntypes-1); i++)
		    for (j=0; j<=(ntypes-1); j++)
		      *(template2+i*ntypes+j)=0;
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      *(template2+i*ntypes+i)=1;
		      for (j=0; j<=(ntypes-1); j++)	
			*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
		    }		  
		  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);		  
		  mat_multiply(ntypes,ndrugs,1,Sa,vec,vec2);		  
		  diagonal(ntypes,vec2,template);		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)-=*(template+i*ntypes+j); 		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*predictionperiod;
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);		  
		}

	      // compmode=1: incur the recursive function to estimate the population.
	      else if (compmode==1)
		recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,predictionperiod,vec,vec2);
	     	  
	      // vec2 is the estimated population vector at the end of the time interval.
	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      
	      // Record the following quantities for each drug combination.
	      // (1)dosages of each drug, (2)total population, (3)population of each subclone.
	      cnt++; 
	      for (i=0; i<=(ndrugs-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+i)=*(vec+i);
	      *(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs)=g;
	      for (i=0; i<=(ntypes-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs+1+i)=*(vec2+i);	    		     
	    }
		  
	  // Select the optimal dosage combination.
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;
	  	  
	  
	  // Strategy 1: Select the dosage combination to minimize total population.
	  if (mode==1)
	    {
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
		    {
		      for (j=0; j<=(ndrugs-1); j++)
			*(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
		      minval=val;		      
		    }
		}
	    }
	  
	  // Strategy 2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
	  else if (mode==2)
	    {
	      minval=1e20;	    

	      // Total population is above the threshold: minimize total population.
	      if (pop>=popthre)	      
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}
	      
	      // Total population is below the threshold: minimize the risk of developing incurable cells.
	      // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
	      // prob((1,...1) \rightarrow (1,...,1))=1.
	      // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
	      else
		{		  
		  for (i=0; i<=(cnt-1); i++)
		    {
		      double risk=0;

		      for (j=0; j<=(ntypes-1); j++)
			{
			  val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);			  
			  num2vec(j,ndrugs,bases,resvec);			  
			  for (k=0; k<=(ndrugs-1); k++)			    
			    if (*(resvec+k)==0)
			      val*=*(T+(2^k)*ntypes+0);
			  //val*=*(T+(k+1)*ntypes+0);
			  risk+=val;
			}		     

		      if (risk<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=risk;
			}
		    }
		}	
	    }

	  // Strategy 3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
	  else if (mode==3)
	    {
	      minval=1e20;

	      // If incurable cells > 1, then minimize total population.
	      if (*(y0+ntypes-1)>(1+1e-10))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}

	      // Otherwise minimize the predicted total multiple resistant populations.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      val=0;
		      for (j=1; j<=(ntypes-1); j++)
			{
			  num2vec(j,ndrugs,bases,resvec); k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    k+=*(resvec+l);
			  if (k>1)
			    val+=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);
			}			  		      
		      if (val<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;	
			}
		    }
		}
	    }
	}
	  	     	                  
      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
     
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=*(mininds+i);              
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
              
     
      /*
      // Debug
      //if ((mode==0)&&(n==0))
      if ((mode==2)&&(popthre>1e10)&&(n==0))
      //if ((mode==1)&&(n==0))
      //if ((mode==3)&&(n==0))
	{
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
	}
      */

      /*
      // Debug
      //if (mode==0)
      if ((mode==2)&&(popthre>1e10))
      //if (mode==1)
      //if (mode==3)
      //if ((mode==0)||((mode==2)&&(popthre>1e10)))
	{	  
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(y+i*nintervals+n+1);
	  printf("t(%d)=%.2f, dosage=( ",n,*(t+n));
	  for (i=0; i<=(ndrugs-1); i++)
	    printf("%.2f ",*(dosage+i*nintervals+n));
	  printf("), y=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*(y+i*nintervals+n+1));
	  printf(") %.4e\n",val);	  
	  //printf("y0=( "); val=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //{
	  //  printf("%.4e ",*(y0+i)); val+=*(y0+i);
	  //}
	  //printf(") %.4e\n",val);
	}      
      */
           
 
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
           
      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask); free(resvec); free(bases);
  
  for (i=0; i<=(ndrugs-1); i++)
    *(dosage+i*nintervals+n)=*(dosage+i*nintervals+n-1);
  
  /*
  // If stop<0, then decide whether the patient recovers or not.
  // Check if each final subpopulation < 1.
  // If the patient recovers, then report survival time to max time span + timeinterval.  
  if (stopT<0)
    {
      k=-1; i=0;
      while ((i<=(ntypes-1))&&(k<0))
	{
	  //if (*(y+i*nintervals+nintervals-1)>1)
	  if (*(y+i*nintervals+nintervals-1)>10)
	    k=i;
	  i++;
	}

      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }
  */

  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *localdosage, *inity, *finaly, maxtime=nintervals*timeinterval;
      localdosage=malloc(sizeof(double)*ndrugs);
      finaly=malloc(sizeof(double)*ntypes);
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(localdosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,localdosage,finaly,0);
	  flag=1; m=0;
	  while ((m<=(ntypes-1))&&(flag==1))
	    {
	      if (*(finaly+m)>1)
		flag=0;
	      m++;
	    }
	  l++;
	}
      if (flag==1)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);

      // Debug
      //printf("finaly: ( ");
      //for (i=0; i<=(ntypes-1); i++)
      //printf("%.3e ",*(finaly+i));
      //printf(")\n");

      free(localdosage); free(inity); free(finaly);
    }

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF); free(dstates);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}


// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
// Difference from recurse_two_drug_response_trimrates: reduce the tolerance error range for determining boundary crossing.  Can correct some errors from previous simulations.
// Difference from recurse_two_drug_response_trimrates2: multi-drug version.
// Difference from recurse_multidrug_response_trimrates2: do not apply U(x-1) in the equation if evamode=1.
// evamode=0: same as recurse_multidrug_response_trimrates2.
// evamode=1: direct calculation of matrix exponential.
// Difference from recurse_multidrug_response_trimrates2_2: speed up matrix exponentiation for growth matrix with a specific format.
void recurse_multidrug_response_trimrates3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;

  F=malloc(sizeof(double)*ntypes*ntypes);
  F2=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);
  vec3=malloc(sizeof(double)*ntypes);
  yn=malloc(sizeof(double)*ntypes); 
  mask=malloc(sizeof(double)*ntypes*ntypes);

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
  
  // Set initial conditions.
  y0=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *(y0+i)=*(x0+i);    

  // Mask matrix.
  // Mask the terms involved in the initial populations below 1.
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(y0+j)<1)
	    *(mask+i*ntypes+j)=0;
	  else
	    *(mask+i*ntypes+j)=1;
	  
	  if (evamode==1)
	    *(mask+i*ntypes+j)=1;
	}
    }   

  // Get the masked differential equation rates.
  mat_multiply(ntypes,ndrugs,1,Sg,dosage,vec2);
  diagonal(ntypes,vec2,template);  
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(template2+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(template2+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      
      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);
  
  mat_multiply(ntypes,ndrugs,1,Sa,dosage,vec2);
  diagonal(ntypes,vec2,template);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)		  
      *(F+i*ntypes+j)-=*(template+i*ntypes+j);
  
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));

  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);

  drugmat_expm(ntypes,tempF,template);

  mat_multiply(ntypes,ntypes,1,template,y0,yn); 

  alters=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(alters+i)=0;
      if ((*(y0+i)<1)&&(*(yn+i)>=1))
	*(alters+i)=1;
      else if ((*(y0+i)>=1)&&(*(yn+i)<1))
	*(alters+i)=-1;
    }

  k=0;
  for (i=0; i<=(ntypes-1); i++)
    if (*(alters+i)!=0)
      k++;

  if (evamode==1)
    k=0;
  
  // If there is no boundary crossing between startt and endt, then use the same differential equation throughout the time interval.
  if (k==0)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(responses+i)=*(yn+i);

      // Release memory.
      free(F); free(F2); free(vec2); free(vec3); free(constF); free(tempF); free(mask);
      free(y0); free(yn); free(template); free(template2); free(template3); free(alters);

    }
  // Otherwise find the time point where the first boundary crossing occurs.
  else
    {
      double t1=startt, t2=endt, curt=(t1+t2)/2;      
      for (i=0; i<=(ntypes-1); i++)
	{
	  *(vec2+i)=*(y0+i);
	}
      //while ((t2-t1)>1e-5)      
      //while ((t2-t1)>1e-1)
      while ((t2-t1)>bsearchthre)
	{
	  curt=(t1+t2)/2; 	 
	  for (i=0; i<=(ntypes-1); i++)
	    for (j=0; j<=(ntypes-1); j++)
	      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(curt-t1);

	  drugmat_expm(ntypes,tempF,template);	 
	  mat_multiply(ntypes,ntypes,1,template,vec2,vec3);	  
	 
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(alters+i)=0;

	      /*
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>0.01)&&((*(vec3+i)-1)>0.01))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>0.01)&&((1-*(vec3+i))>0.01))
		*(alters+i)=1;
		*/
	      
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
		*(alters+i)=1;				
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;
	  if (k==0)
	    {
	      t1=curt;
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(vec3+i);
	    }
	  else
	    {
	      t2=curt;
	    }	  
	}

      
      // Release memory.
      free(F); free(F2); free(vec2); free(constF); free(tempF); free(mask);
      free(yn); free(template); free(template2); free(template3); free(alters);
      
      if (max(t1,t2)>=endt)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(responses+i)=*(vec3+i);
	}

      // Modify the differential equation terms and run the equations again.      
      else
	{	  
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(vec3+i);
	  //recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
	  recurse_multidrug_response_trimrates3(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses,evamode);
	}      

      // Release memory.
      free(y0); free(vec3);
    }
	    
  return;
}

  
/*
// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
// Difference from recurse_two_drug_response_trimrates: reduce the tolerance error range for determining boundary crossing.  Can correct some errors from previous simulations.
// Difference from recurse_two_drug_response_trimrates2: multi-drug version.
// Difference from recurse_multidrug_response_trimrates2: do not apply U(x-1) in the equation if evamode=1.
// evamode=0: same as recurse_multidrug_response_trimrates2.
// evamode=1: direct calculation of matrix exponential.
// Difference from recurse_multidrug_response_trimrates2_2: do not incur recursion, cut duplicated matrix exponential calculations.
void recurse_multidrug_response_trimrates2_3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timespan, double *dosage, double *responses, int evamode)
{
  int i, j, k, l, m, n, *alters, nalters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val, startt, endt;

  startt=0; endt=timespan;
  F=malloc(sizeof(double)*ntypes*ntypes);
  F2=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);
  vec3=malloc(sizeof(double)*ntypes);
  yn=malloc(sizeof(double)*ntypes); 
  mask=malloc(sizeof(double)*ntypes*ntypes);

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
  
  // Set initial conditions.
  y0=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *(y0+i)=*(x0+i);    

  // Mask matrix.
  // Mask the terms involved in the initial populations below 1.
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(y0+j)<1)
	    *(mask+i*ntypes+j)=0;
	  else
	    *(mask+i*ntypes+j)=1;	  
	  if (evamode==1)
	    *(mask+i*ntypes+j)=1;
	}
    }   

  // Get the masked differential equation rates.
  mat_multiply(ntypes,ndrugs,1,Sg,dosage,vec2);
  diagonal(ntypes,vec2,template);  
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(template2+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(template2+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)      
      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);
  
  mat_multiply(ntypes,ndrugs,1,Sa,dosage,vec2);
  diagonal(ntypes,vec2,template);
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)		  
      *(F+i*ntypes+j)-=*(template+i*ntypes+j);  
  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template);
  mat_multiply(ntypes,ntypes,1,template,y0,yn); 
 
  alters=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(alters+i)=0;
      if ((*(y0+i)<1)&&(*(yn+i)>=1))
	*(alters+i)=1;
      else if ((*(y0+i)>=1)&&(*(yn+i)<1))
	*(alters+i)=-1;
    }

  nalters=0;
  for (i=0; i<=(ntypes-1); i++)
    if (*(alters+i)!=0)
      nalters++;

  if (evamode==1)
    nalters=0;

  
  // If there is no boundary crossing between startt and endt, then use the same differential equation throughout the time interval.
  if (nalters==0)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(responses+i)=*(yn+i);
    }

  // Otherwise proceed matrix exponential until either nalters=0 or it reachs the terminal time.
  else
    {
      while ((startt<endt)&&(nalters>0))
	{
	  double t1=startt, t2=endt, curt=(t1+t2)/2;
  
	  // Find the time point where the first boundary crossing occurs.
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(vec2+i)=*(y0+i);
	    }
	  while ((t2-t1)>1e-5)
	    {
	      curt=(t1+t2)/2; 
	      for (i=0; i<=(ntypes-1); i++)
		for (j=0; j<=(ntypes-1); j++)
		  *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(curt-t1);
	      padm(ntypes,tempF,1,7,template);
	      mat_multiply(ntypes,ntypes,1,template,vec2,vec3);	  
	      for (i=0; i<=(ntypes-1); i++)
		{
		  *(alters+i)=0;	      
		  // Avoid the deadlock at boundary 1.
		  if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
		    *(alters+i)=1;
		  else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
		    *(alters+i)=1;				             
		}
	      k=0;
	      for (i=0; i<=(ntypes-1); i++)
		if (*(alters+i)!=0)
		  k++;
	      if (k==0)
		{
		  t1=curt;
		  for (i=0; i<=(ntypes-1); i++)
		    *(vec2+i)=*(vec3+i);
		}
	      else
		{
		  t2=curt;
		}	  
	    }
	  
	  // If the boundary crossing occurs at the right end, then stop.
	  if (max(t1,t2)>=endt)
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*(responses+i)=*(vec3+i);
	      nalters=0;
	    }

	  // Otherwise modify the differential equation terms and run the equations again.  
	  else
	    {	      
	      // Update the initial time.
	      startt=curt;
	      // Update the initial condition.
	      for (i=0; i<=(ntypes-1); i++)
		*(y0+i)=*(vec3+i);
	      // Update the mask matrix.
	      // Mask the terms involved in the initial populations below 1.
	      for (i=0; i<=(ntypes-1); i++)
		{
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      if (*(y0+j)<1)
			*(mask+i*ntypes+j)=0;
		      else
			*(mask+i*ntypes+j)=1;	  
		      if (evamode==1)
			*(mask+i*ntypes+j)=1;
		    }
		}   
	      for (i=0; i<=(ntypes-1); i++)		
		for (j=0; j<=(ntypes-1); j++)
		  *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));
	      for (i=0; i<=(ntypes-1); i++)
		for (j=0; j<=(ntypes-1); j++)
		  *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
	      padm(ntypes,tempF,1,7,template);
	      mat_multiply(ntypes,ntypes,1,template,y0,yn); 
	    
	      // Count the number of boundary crosses.
	      for (i=0; i<=(ntypes-1); i++)
		{
		  *(alters+i)=0;
		  if ((*(y0+i)<1)&&(*(yn+i)>=1))
		    *(alters+i)=1;
		  else if ((*(y0+i)>=1)&&(*(yn+i)<1))
		    *(alters+i)=-1;
		}	      
	      nalters=0;
	      for (i=0; i<=(ntypes-1); i++)
		if (*(alters+i)!=0)
		  nalters++;
	      if (nalters==0)
		{		  		  
		  for (i=0; i<=(ntypes-1); i++)
		    *(responses+i)=*(yn+i);		  
		}
	    }
	}
    }    
    
  // Release memory.
  free(F); free(F2); free(vec2); free(constF); free(tempF); free(mask);
  free(yn); free(template); free(template2); free(template3); free(alters);
  free(y0); free(vec3);	    

  return;
}
*/



// Implement optimization heuristics for multiple drugs.
// Adapt optimize_two_drug_responses8.
// Difference from optimize_two_drug_responses8: exercise fewer possible drug combinations.
// For k drugs, exhaust combinations x1+...+xk=1.
// For all the xk!=0, equalize their dosages.
// There are totally 2^{k}-1 drug combinations.
// Difference from optimize_multidrug_responses1: Alter strategy 0.  Rank strategies by their rates of reducing the dominant subclone.
// Difference from optimize_multidrug_responses2: correct the bug of evaluating risk for strategies 2.1 and 2.2.
// The order of drugs in the cell types is (d3,d2,d1).  But the order of drugs in the transition matrix is (d1,d2,d3).  Thus need to invert the order when evaluating risk.
double correct_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage)
{
  int i, j, k, l, m, n, cnt, *bases, combind, *resvec, nstates;
  double *d, *mininds, *mininds2, *mininds3, minval, minval2, minval3, *mask, maxval;
  double *F, *template, *constF, *y0, *vec, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, val2, pop, *drugrecords, *dstates;

  d=malloc(sizeof(double)*ndrugs); 
  mininds=malloc(sizeof(double)*ndrugs);
  mininds2=malloc(sizeof(double)*ndrugs);
  mininds3=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
  mask=malloc(sizeof(double)*ntypes*ntypes);
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (i=0; i<=(ndrugs-1); i++)
    *(bases+i)=2;  
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }

  // Constant matrix.
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)=0;
  for (i=0; i<=(ntypes-1); i++)
    {
      *(constF+i*ntypes+i)=1;
      for (j=0; j<=(ntypes-1); j++)	
	*(constF+i*ntypes+j)+=*(T+i*ntypes+j);
    }
  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(constF+i*ntypes+j)*=*(g0+j);
  for (i=0; i<=(ntypes-1); i++)
    *(constF+i*ntypes+i)-=*(a0+i);
 
  drugrecords=malloc(sizeof(double)*(ndrugs+1+ntypes)*ntypes);
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      *(y0+i)=*(x0+i);
	      *(y+i*nintervals+n)=*(y0+i);
	    }
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);


      // Strategy 0 does not need prediction.      
      if (mode==0)
	{	  
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;	  

	  // Initially find the drug that yields the highest rate of reducing the dominant population.
	  // The growth rate of type i is g0(1+\sum_{j!=i}T(j->i))-\sum_{k}Sa(i,k)dk.
	  if (n==0)
	    {
	      int dominant=-1;
	      double rates[ndrugs];	      
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  if (*(y0+i)>val)
		    {
		      dominant=i; val=*(y0+i);
		    }
		}
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  rates[i]=*(g0+dominant);
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      if (j!=dominant)
			{
			  val=*(T+dominant*ntypes+j)*(*(g0+dominant));
			  rates[i]+=val;
			}
		    }
		  val=*(Sa+dominant*ndrugs+i);
		  rates[i]-=val;
		}
	      k=-1; val=1e500;
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  if (rates[i]<val)
		    {
		      k=i; val=rates[i];
		    }
		}
	       for (i=0; i<=(ndrugs-1); i++)
		*(mininds+i)=0;
	      *(mininds+k)=1;
	    }

	  // Switch the drug when the total population reaches the double the nadir or the population reemerges.
	  else
	    {
	      int switching=0;
	      double nadir=1e50;
	      
	      // Find the nadir of the current regimen.
	      // A nadir is the local minimum of the population under the current regimen.
	      // Start from the current time step, trace back until the change of regimen.
	      i=n-1; k=1;
	      while ((i>=0)&&(k==1))
		{
		  int p, q;

		  // Check whether the dosages at time step i equals the dosages at time step n-1.
		  p=0; q=1;
		  while ((p<=(ndrugs-1))&&(q==1))
		    {
		      if (fabs(*(dosage+p*nintervals+i)-*(dosage+p*nintervals+n-1))>1e-12)
			q=0;
		      p++;
		    }
		  // If so then update nadir.
		  if (q==1)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  if (i==0)
			    val+=*(x0+j);
			  else
			    val+=*(y+j*nintervals+i);
			}
		      nadir=min(nadir,val);
		    }
		  else
		    k=0;
		  i--;
		}
     
	      // If the current population >= 2 x nadir population, then switch.
	      if (pop>=(nadir*2))
		switching=1;

	      // If the current population reemerges from the level below the detection threshold, then switch.
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  if ((n-1)==0)
		    val+=*(x0+i);
		  else
		    val+=*(y+i*nintervals+n-1);
		}

	      if ((val<detected)&&(pop>=detected))
		switching=1;
	      
	      // Find the most effective drug other than the current choice.
	      if (switching==1)
		{
		  int dominant=-1;
		  double rates[ndrugs];	      
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      if (*(y0+i)>val)
			{
			  dominant=i; val=*(y0+i);
			}
		    }
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      rates[i]=*(g0+dominant);
		      for (j=0; j<=(ntypes-1); j++)
			{
			  if (j!=dominant)
			    {
			      val=*(T+dominant*ntypes+j)*(*(g0+dominant));
			      rates[i]+=val;
			    }
			}
		      val=*(Sa+dominant*ndrugs+i);
		      rates[i]-=val;
		    }
		  j=-1; maxval=0;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if (*(dosage+i*nintervals+n-1)>maxval)
			{
			  maxval=*(dosage+i*nintervals+n-1); j=i;
			}
		    }
		  k=-1; val=1e500;
		  for (i=0; i<=(ndrugs-1); i++)
		    {
		      if ((i!=j)&&(rates[i]<val))
			{
			  k=i; val=rates[i];
			}
		    }
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=0;
		  *(mininds+k)=1;
		}
	      else
		{
		  for (i=0; i<=(ndrugs-1); i++)
		    *(mininds+i)=*(dosage+i*nintervals+n-1);
		}
	    }	  
	  
	  for (i=0; i<=(ntypes-1); i++)
	    *(vec2+i)=*(y0+i);
	}

      // For other strategies, calculate the predictions.
      else
	{
	  // Mask matrix.
	  // Mask the terms involved in the initial populations below 1.
	  for (i=0; i<=(ntypes-1); i++)
	    {
	      for (j=0; j<=(ntypes-1); j++)
		{
		  if (*(y0+j)<1)
		    *(mask+i*ntypes+j)=0;
		  else
		    *(mask+i*ntypes+j)=1;
		}
	    }   	  
	  
	  // Change the drug combinations to optimize certain criterion.
	  // Here the predicted responses are within the designated time interval (predictionperiod) instead of the sampling time interval.

	  // mininds2 and minval2: the drug dosage that minimizes the total population.
	  // mininds3 and minval3: the drug dosage that minimizes the R12 population.

	  // Consider the drug combinations (d1,d2,...,dk), where di's are either 0 or 1/m, where m is the number of nonzero di's.
	  // Exhaust all combinations of drugs except (0,0,...0).
	  
	  for (i=0; i<=(ndrugs-1); i++)
	    {
	      *(mininds+i)=0; *(mininds2+i)=0; *(mininds3+i)=0;
	    }

	  minval=1e100; minval2=1e100; minval3=1e100; cnt=0; 
	  for (combind=1; combind<=(ntypes-1); combind++)
	    {	     
	      num2vec(combind,ndrugs,bases,resvec);
	      k=0;
	      for (i=0; i<=(ndrugs-1); i++)	       
		k+=*(resvec+i);
	      val=(double)1/(double)k;
	      
	      for (i=0; i<=(ndrugs-1); i++)
		{
		  *(vec+i)=*(resvec+i)*val;
		}
	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(y0+i);
	  	  
		  
	      // compmode=0: use F to estimate the population.
	      if (compmode==0)
		{		  
		  mat_multiply(ntypes,ndrugs,1,Sg,vec,vec2);
		  diagonal(ntypes,vec2,template);	      
		  for (i=0; i<=(ntypes-1); i++)
		    for (j=0; j<=(ntypes-1); j++)
		      *(template2+i*ntypes+j)=0;
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      *(template2+i*ntypes+i)=1;
		      for (j=0; j<=(ntypes-1); j++)	
			*(template2+i*ntypes+j)+=*(T+i*ntypes+j);
		    }		  
		  mat_multiply(ntypes,ntypes,ntypes,template2,template,template3);
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)=*(constF+i*ntypes+j)-*(template3+i*ntypes+j);		  
		  mat_multiply(ntypes,ndrugs,1,Sa,vec,vec2);		  
		  diagonal(ntypes,vec2,template);		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(F+i*ntypes+j)-=*(template+i*ntypes+j); 		  
		  for (i=0; i<=(ntypes-1); i++)		
		    for (j=0; j<=(ntypes-1); j++)		  
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*predictionperiod;
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);		  
		}

	      // compmode=1: incur the recursive function to estimate the population.
	      else if (compmode==1)
		recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,predictionperiod,vec,vec2);
	     	  
	      // vec2 is the estimated population vector at the end of the time interval.
	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      
	      // Record the following quantities for each drug combination.
	      // (1)dosages of each drug, (2)total population, (3)population of each subclone.
	      cnt++; 
	      for (i=0; i<=(ndrugs-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+i)=*(vec+i);
	      *(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs)=g;
	      for (i=0; i<=(ntypes-1); i++)
		*(drugrecords+(cnt-1)*(ndrugs+1+ntypes)+ndrugs+1+i)=*(vec2+i);	    		     
	    }
		  
	  // Select the optimal dosage combination.
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;
	  	  
	  //MDM modifityin for mode 1 to be DPM 2.2 and mode 2 to be DPM 2.2 for first 2 steps, then switch to mode 0
	  if (mode==1)
	    {
	      minval=1e20;

	      // If incurable cells > 1, then minimize total population.
	      if (*(y0+ntypes-1)>(1+1e-10))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}

	      // Otherwise minimize the predicted total multiple resistant populations.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      val=0;
		      for (j=1; j<=(ntypes-1); j++)
			{
			  num2vec(j,ndrugs,bases,resvec); k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    k+=*(resvec+l);
			  if (k>1)
			    val+=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);
			}			  		      
		      if (val<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;	
			}
		    }
		}
	    }
	  //MDM DPM 2.2 for first 2 steps then strategy 0
	  else if (mode==2)
	    {
	      minval=1e20;

	      // If incurable cells > 1, then minimize total population.
	      if (*(y0+ntypes-1)>(1+1e-10))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}

	      // Otherwise minimize the predicted total multiple resistant populations.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      val=0;
		      for (j=1; j<=(ntypes-1); j++)
			{
			  num2vec(j,ndrugs,bases,resvec); k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    k+=*(resvec+l);
			  if (k>1)
			    val+=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);
			}			  		      
		      if (val<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;	
			}
		    }
		}
//MDM switch to mode 0 after second timestep
		  if (n>1)
		mode=0;
	    }
	  // Strategy 3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
	  else if (mode==3)
	    {
	      minval=1e20;

	      // If incurable cells > 1, then minimize total population.
	      if (*(y0+ntypes-1)>(1+1e-10))
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs))<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;		      
			}
		    }
		}

	      // Otherwise minimize the predicted total multiple resistant populations.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      val=0;
		      for (j=1; j<=(ntypes-1); j++)
			{
			  num2vec(j,ndrugs,bases,resvec); k=0;
			  for (l=0; l<=(ndrugs-1); l++)
			    k+=*(resvec+l);
			  if (k>1)
			    val+=*(drugrecords+i*(ndrugs+1+ntypes)+ndrugs+1+j);
			}			  		      
		      if (val<minval)
			{
			  for (j=0; j<=(ndrugs-1); j++)
			    *(mininds+j)=*(drugrecords+i*(ndrugs+1+ntypes)+j);
			  minval=val;	
			}
		    }
		}
	    }
	}
	  	     	                  
      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_multidrug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
     
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=*(mininds+i);              
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
              
     
      /*
      // Debug
      //if ((mode==0)&&(n==0))
      if ((mode==2)&&(popthre>1e10)&&(n==0))
      //if ((mode==1)&&(n==0))
      //if ((mode==3)&&(n==0))
	{
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
	}
      */

      /*
      // Debug
      //if (mode==0)
      if ((mode==2)&&(popthre>1e10))
      //if (mode==1)
      //if (mode==3)
      //if ((mode==0)||((mode==2)&&(popthre>1e10)))
	{	  
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(y+i*nintervals+n+1);
	  printf("t(%d)=%.2f, dosage=( ",n,*(t+n));
	  for (i=0; i<=(ndrugs-1); i++)
	    printf("%.2f ",*(dosage+i*nintervals+n));
	  printf("), y=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*(y+i*nintervals+n+1));
	  printf(") %.4e\n",val);	  
	  //printf("y0=( "); val=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //{
	  //  printf("%.4e ",*(y0+i)); val+=*(y0+i);
	  //}
	  //printf(") %.4e\n",val);
	}      
      */
           
 
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
           
      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask); free(resvec); free(bases);
  
  for (i=0; i<=(ndrugs-1); i++)
    *(dosage+i*nintervals+n)=*(dosage+i*nintervals+n-1);
  
  /*
  // If stop<0, then decide whether the patient recovers or not.
  // Check if each final subpopulation < 1.
  // If the patient recovers, then report survival time to max time span + timeinterval.  
  if (stopT<0)
    {
      k=-1; i=0;
      while ((i<=(ntypes-1))&&(k<0))
	{
	  //if (*(y+i*nintervals+nintervals-1)>1)
	  if (*(y+i*nintervals+nintervals-1)>10)
	    k=i;
	  i++;
	}

      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }
  */

  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *localdosage, *inity, *finaly, maxtime=nintervals*timeinterval;
      localdosage=malloc(sizeof(double)*ndrugs);
      finaly=malloc(sizeof(double)*ntypes);
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(localdosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,localdosage,finaly,0);
	  flag=1; m=0;
	  while ((m<=(ntypes-1))&&(flag==1))
	    {
	      if (*(finaly+m)>1)
		flag=0;
	      m++;
	    }
	  l++;
	}
      if (flag==1)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);

      // Debug
      //printf("finaly: ( ");
      //for (i=0; i<=(ntypes-1); i++)
      //printf("%.3e ",*(finaly+i));
      //printf(")\n");

      free(localdosage); free(inity); free(finaly);
    }

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF); free(dstates);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}


// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
// Difference from dp_optimize_two_drug_responses4: apply branch and bound to trim unnecessary nodes in recursion.
// Difference from dp_optimize_two_drug_responses5: (1)multi-drug version, (2)make the format of y and dosage compatible with optimize_multidrug_responses1, (2)can choose one of multiple strategies as baseline for dp optimization.
// mode=0: monotherapy switching, benchmark.
// mode=1: minimize total population.
// mode=2: select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.
// mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.
// mode=4: maximize the anticipated death time for the best monotonous treatment.
// Difference from dp_optimize_multidrug_responses1: modify the criteria to improve strategy 2.  Ensure the predicted total population in multiple steps does not exceed the threshold.  If so then minimize the total population.
// Difference from dp_optimize_multidrug_responses2: modify recurse_sequence_multidrug_responses to create valid boundaries.
// Difference from dp_optimize_multidrug_responses2_2: correct the bug of evaluating risk for strategies 2.1 and 2.2.
// The order of drugs in the cell types is (d3,d2,d1).  But the order of drugs in the transition matrix is (d1,d2,d3).  Thus need to invert the order when evaluating risk.
double correct_dp_optimize_multidrug_responses2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int degmode, int mode, double popthre)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, *resvec, combind, *bases;
  long curind, *seq;
  double *dstates, val, stopT, *frontiers, pop, incurablepop;
  struct seqnode *seqnodes;  
 
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;
  nstates=ntypes-1; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (combind=1; combind<=(ntypes-1); combind++)
    {	     
      num2vec(combind,ndrugs,bases,resvec);
      k=0;
      for (i=0; i<=(ndrugs-1); i++)
	if (*(resvec+i)==1)
	  k++;
      val=(double)1/(double)k;
      for (i=0; i<=(ndrugs-1); i++)
	*(dstates+(combind-1)*ndrugs+i)=*(resvec+i)*val;
    }
  

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }

  /*
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+n*ndrugs+i)=-1;
    }    

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);
    */

  seq=malloc(sizeof(long)*nsteps); n=0;
  nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

      // Debug
      //printf("n=%d\n",n);

      // (n%nsteps)=0: incur dynamic programming.
      if ((n%nsteps)==0)
	{
	  if (n>0)
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  free((seqnodes+i)->children);
		  free((seqnodes+i)->dosage);
		  free((seqnodes+i)->y);
		}
	      free(seqnodes);	     	      
	    }
	  nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
	  (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
	  (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
	  (seqnodes+0)->children=malloc(sizeof(long));
	  (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
	  (seqnodes+0)->y=malloc(sizeof(double)*ntypes);
	  if (n==0)
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+i);
	      
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+0);
	    }
	  else
	    {
	      //for (i=0; i<=(ntypes-1); i++)
	      //*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);

	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i*nintervals+n-1);
	    }

	  for (i=0; i<=nsteps; i++)
	    *(nfrontiers+i)=0;

	  curind=0;	  	  	      	      	 
	  //seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  seqnodes=recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  	  	 
	  // Find the max depth of the seqnodes.
	  maxdepth=-1;
	  for (i=0; i<=(nseqnodes-1); i++)
	    maxdepth=max(maxdepth,(seqnodes+i)->depth);	  	  	  
	  
	  // Debug
	  //if (n==5)
	  //exit(0);


	  /*
	  // Debug
	  if ((mode==2)&&(popthre>1e10)&&(n==5))
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
		    {
		      int seq[5];
		      seq[0]=i;		      
		      for (j=1; j<=4; j++)
			seq[j]=(seqnodes+seq[j-1])->pa;
		      printf("node %d, dosage ",i);
		      for (l=4; l>=0; l--)
			{
			  printf("( ");
			  for (j=0; j<=(ndrugs-1); j++)
			    printf("%.2f ",*((seqnodes+seq[l])->dosage+j));
			  printf(") ");
			}
		      printf("\n");
		      printf("pop ( ");
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.3e ",*((seqnodes+i)->y+j));
		      printf(")\n");
		    }
		}
	    }
	  */
	      


	  // Find the candidate terminal seqnodes.
	  ncandidates=0; candidates=malloc(sizeof(long));
	  
	  // Find the terminal nodes with zero population.
	  // If so then pick up the one that minimizes the total population.
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      if ((seqnodes+i)->nchildren==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  if (val<1)
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	    }	  
	  
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);
	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  
	  
	  // If the tumor cannot be eradicated within nsteps, then find the candidate sequences.
	  // A candidate sequence has the max depth, terminates at a leaf, and is not inferior to any other sequences.
	  // A sequence s1 is inferior to another sequence s2 if the total population and R12 population of s1 are both greater than those of s2.
	  if (ncandidates==0)
	    {
	      struct pair2 *pairs1, *pairs2;
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates, allhightotal, allmortal;
	      double p1, p2, p, *deathtimes;

	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0;
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  if ((c1==1)&&(c2==1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }		  
		}	     

	      // Sort candidate sequences by total population and by incurable population.
	      pairs1=malloc(sizeof(struct pair2)*ncandidates);
	      pairs2=malloc(sizeof(struct pair2)*ncandidates);
	      lookup1=malloc(sizeof(int)*ncandidates);
	      lookup2=malloc(sizeof(int)*ncandidates);
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  (pairs1+i)->i=i;
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  (pairs1+i)->s=val;
		  (pairs2+i)->i=i;
		  val=*((seqnodes+*(candidates+i))->y+ntypes-1);
		  (pairs2+i)->s=val;
		}

	      qsort(pairs1,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      qsort(pairs2,ncandidates,sizeof(struct pair2),(* pair2_cmp));
	      
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  j=(pairs1+i)->i; *(lookup1+j)=i;
		  j=(pairs2+i)->i; *(lookup2+j)=i;
		}


	      /*
	      // Debug
	      //if ((mode==2)&(popthre>1e10)&&(n==1))
	      //if ((mode==2)&(popthre>1e10)&&(n==2))
	      //if ((mode==2)&(popthre>1e10)&&(n==5))
	      if ((mode==2)&(popthre>1e10)&&(n==0))
		{
		  printf("ncandidates=%d, nseqnodes=%d\n",ncandidates,nseqnodes);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      int ii;

		      printf("%d %d dosage=( ",i,*(candidates+i));
		      for (j=0; j<=(ndrugs-1); j++)
			printf("%.2f ",*((seqnodes+*(candidates+i))->dosage+j));
		      printf("), ");		      
		      printf("y=( ");
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			{
			  printf("%.2e ",*((seqnodes+*(candidates+i))->y+j));
			  val+=*((seqnodes+*(candidates+i))->y+j);
			}
		      printf(") %.2e\n",val);

		      printf("history ");
		      ii=*(candidates+i);
		      while ((seqnodes+ii)->pa>=0)
			{
			  printf("( ");
			  for (j=0; j<=(ndrugs-1); j++)
			    printf("%.2f ",*((seqnodes+ii)->dosage+j));
			  printf(") ");
			  ii=(seqnodes+ii)->pa;
			}
		      printf("\n");

		    }
		}		      
	      */


	      // Find the candidates inferior to others.
	      // Direct comparison of all candidate pairs take n^2 time, n=ncandidates.
	      // To reduce computing time use probabilistic arguments.
	      // Suppose a candidate's total population > n1 other candidates, and its incurable population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and incurable population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in incurable population) =
	      // 1-Pr(there is no overlap ...).
	      // Pr(there is no overlap ...) = 
	      // Pr(place n2 candidates in (n-n1) among n slots) or
	      // Pr(place n1 candidates in (n-n2) among n slots).
	      // If (n1+n2)>=n, then Pr(there is an overlap ...)=1.
	      // If (n1+n2)<n, then
	      // P1=((n-n1)/n)^n2, P2=((n-n2)/n)^n1. 
	      // Pr(the candidate is inferior) = 1 - max(P1,P2).
	      inferior=malloc(sizeof(int)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{		  
		  *(inferior+i)=0; p=0;
		  j=*(lookup1+i); n1=max(j-1,0);
		  j=*(lookup2+i); n2=max(j-1,0);
		  if ((n1+n2)>=ncandidates)
		    p=1;
		  else
		    {
		      val=(double)(ncandidates-n1)/(double)ncandidates;
		      p1=n2*log(val); p1=exp(p1);
		      val=(double)(ncandidates-n2)/(double)ncandidates;
		      p2=n1*log(val); p2=exp(p2);
		      p=1-max(p1,p2);
		    }

		  if (p>=(1-1e-4))
		    *(inferior+i)=1;
		  else
		    {
		      // Compare the superior list in two values.
		      j=0; k=-1;
		      while ((j<=n1)&&(k<0))
			{
			  l=(pairs1+j)->i; l=*(lookup2+l);
			  if (((pairs1+*(lookup1+i))->s>(pairs1+j)->s)&&((pairs2+*(lookup2+i))->s>(pairs2+l)->s))
			    //if (l<=n2)
			    k=l;
			  j++;
			}
		      if (k>=0)
			*(inferior+i)=1;
		    }
		}

	      // Check whether all candidates lead to mortality.
	      allmortal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allmortal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<mortal)
		    allmortal=0;
		  i++;
		}

	      // Check whether all candidates lead to high total population.
	      allhightotal=1; i=0;
	      while ((i<=(ncandidates-1))&&(allhightotal==1))
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*((seqnodes+*(candidates+i))->y+j);
		  if (val<popthre)
		    allhightotal=0;
		  i++;
		}
	      
	      // Debug
	      //printf("ncandidates=%d, before\n",ncandidates);

	      // Choose a candidate according to the selected strategy.

	      // mode=1: Minimize the total population. 
	      if (mode==1)
		{		  
		  minind=-1; minval=1e100;	      
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      val=0;
		      for (j=0; j<=(ntypes-1); j++)
			val+=*((seqnodes+*(candidates+i))->y+j);
		      if (val<minval)
			{
			  minind=*(candidates+i); minval=val;
			}
		    }	      
		}	      

	      // mode=2: Select the dosage combination to minimize the risk of incurable cells developing unless the total population exceeds popthre.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==2)
		{		
		  minind=-1; minval=1e20;		  		 	  
		  
		  // Total population is above the threshold: minimize total population.
		  // Total population is above the threshold or all drug combinations lead to mortality: minimize total population.
		  //if (pop>=popthre)	      		  
		  if ((pop>=popthre)||(allhightotal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }		  	     		
	      
		  // Total population is below the threshold: minimize the risk of developing incurable cells.
		  // Risk: \sum_{type(i1,i2,...,ik)} population(type(i1,i2,...,ik)) x prob(type(i1,i2,...,ik) \rightarrow type(1,1,...1)).
		  // prob((1,...1) \rightarrow (1,...,1))=1.
		  // prob((i1,...,ik) \rightarrow (1,...,1))=\prod_{ij=0}prob(S \rightarrow R_{j}).
		  // Discard the combinations that lead to mortality.
		  // Discard the combinations that have high total population.
		  else
		    {		  		      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  //val*=*(T+(2^(ndrugs-1-k))*ntypes+0);				  
				  val*=*(T+power(2,(ndrugs-1-k))*ntypes+0);
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);	
			    }

			  /*
			  // Debug
			  //if ((mode==2)&&(popthre>1e10)&&(n==0))
			  if ((mode==2)&&(popthre>1e10)&&(n==5))
			    {
			      int seq[5];
			      seq[0]=*(candidates+i);
			      for (j=1; j<=4; j++)
				seq[j]=(seqnodes+seq[j-1])->pa;
			      
			      printf("dp, candidate %d, dosage ",*(candidates+i));
			      for (l=4; l>=0; l--)
				{
				  printf("( ");
				  for (j=0; j<=(ndrugs-1); j++)
				    printf("%.2f ",*((seqnodes+seq[l])->dosage+j));
				  printf(") ");
				}
			      printf("risk=%.3e\n",risk);

			      //printf("dp, dosage ( ");
			      //for (j=0; j<=(ndrugs-1); j++)
			      //printf("%.2f ",*((seqnodes+*(candidates+i))->dosages+j));
			      //printf("), risk=%.3e\n",risk);
			      for (j=0; j<=(ntypes-1); j++)
				{
				  val=*((seqnodes+*(candidates+i))->y+j);
				  num2vec(j,ndrugs,bases,resvec);
				  printf("j=%d, vec=(%d %d %d), pop=%.3e, ",j,*(resvec+0),*(resvec+1),*(resvec+2),val);
				  for (k=0; k<=(ndrugs-1); k++)
				    {			    
				      if (*(resvec+k)==0)
					{
					  //val*=*(T+(2^(ndrugs-1-k))*ntypes+0);
					  val*=*(T+power(2,(ndrugs-1-k))*ntypes+0);
					  //printf("(%d %d %.3e) ",k,(2^(ndrugs-1-k)),*(T+(2^(ndrugs-1-k))*ntypes+0));
					  printf("(%d %d %.3e) ",k,power(2,(ndrugs-1-k)),*(T+power(2,(ndrugs-1-k))*ntypes+0));
					}
				    }
				  printf("contribution %.3e",val);
				  printf("\n");
				  //risk2+=val;
				}
			      //printf("risk2=%.3e\n",risk2);
			    }
			  */
			  

			  //if (risk<minval)
			  //if ((risk<minval)&&(val2<mortal))
			  if ((risk<minval)&&(val2<popthre))
			    {			      
			      minind=*(candidates+i); minval=risk;
			    }		

			  // Debug
			  //printf("i=%d, cand=%d, risk=%.4e, val2=%.4e\n",i,*(candidates+i),risk,val2);
	  			  
			}

		      // Debug
		      //if ((mode==2)&&(popthre>1e10)&&(n==5))
		      //printf("minind=%d\n",minind);

		    }

		  /*
		  // Debug
		  if ((popthre>1e10)&&(n==0))
		    {
		      printf("n=%d\n",n);
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double risk=0, val2=0; 
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val=*((seqnodes+*(candidates+i))->y+j);			      
			      num2vec(j,ndrugs,bases,resvec);			  
			      for (k=0; k<=(ndrugs-1); k++)			    
				if (*(resvec+k)==0)
				  val*=*(T+(2^(ndrugs-1-k))*ntypes+0);
			      //val*=*(T+(k+1)*ntypes+0);
			      risk+=val;
			      val2+=*((seqnodes+*(candidates+i))->y+j);	

			      // Debug
			      if ((i==318)||(i==402))
				{
				  printf("( ");
				  for (k=0; k<=(ndrugs-1); k++)
				    printf("%d ",*(resvec+k));
				  printf("%.3e ",val);
				  
				  if (j==3)
				    {
				      printf("((%.3e ",*((seqnodes+*(candidates+i))->y+j));
				      for (k=0; k<=(ndrugs-1); k++)			    
					if (*(resvec+k)==0)
					  printf("%.3e ",*(T+(2^(ndrugs-1-k))*ntypes+0));
				      printf("))");
				    }

				  //for (k=0; k<=(ndrugs-1); k++)
				  //printf("k=%d, resvec=%d, transition=%.3e ",k,*(resvec+k),*(T+(2^k)*ntypes+0));
				  printf(")");
				}
			    }

			  //
			  //printf("candidate %d %d, ( ",i,*(candidates+i));
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  val=*((seqnodes+*(candidates+i))->y+j);
			  //  printf("%.4e ",val);
			  //}			  		  
			  //printf(")\n");
			  //

			  printf("candidate %d %d, risk=%.4e, ",i,*(candidates+i),risk);
			  k=*(candidates+i);

			  if ((k==620)||(k==767))
			    {
			      printf("( ");
			      for (j=0; j<=(ntypes-1); j++)
				{
				  val=*((seqnodes+*(candidates+i))->y+j);
				  printf("%.4e ",val);
				}
			      printf(")");
			    }		

			  while (k>=0)
			    {
			      printf("( ");
			      for (j=0; j<=(ndrugs-1); j++)
				{
				  val=*((seqnodes+k)->dosage+j);
				  printf("%.2f ",val);
				}
			      printf(")");
			      k=(seqnodes+k)->pa;
			    }
			  printf("\n");
			}
		      printf("minind=%d\n",minind);
		    }
		  */
		  

		}	
	    
	      // mode=3: Minimize total population unless predicted multiple resistant population>=1.  In that case minimize total predicted multiple resistant population.  If in the current population incurable cells > 1, then minimize total population.  But when some (but not all) candidate combinations lead to mortality, then discard these combinations.
	      else if (mode==3)
		{
		  minval=1e20;

		  // If incurable cells > 1, then minimize total population.
		  //if (incurablepop>=(1+1e-10))
		  if ((incurablepop>=(1+1e-10))||(allmortal==1))
		    {
		      minind=-1; minval=1e100;	      
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    val+=*((seqnodes+*(candidates+i))->y+j);
			  if (val<minval)
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }	

		  // Otherwise minimize the predicted total multiple resistant populations.
		  // Discard the combinations that lead to mortality.
		  else
		    {
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val2;

			  val=0; val2=0;
			  for (j=1; j<=(ntypes-1); j++)
			    {
			      num2vec(j,ndrugs,bases,resvec); k=0;
			      for (l=0; l<=(ndrugs-1); l++)
				k+=*(resvec+l);
			      if (k>1)
				val+=*((seqnodes+*(candidates+i))->y+j);
			      val2+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  //if (val<minval)
			  if ((val<minval)&&(val2<mortal))
			    {
			      minind=*(candidates+i); minval=val;
			    }
			}
		    }
		}

	      // mode=4: Maximize the anticipated death time based on static therapy.
	      else if (mode==4)
		{
	      	      	     	      
		  // Evaluate the expected death time for each non-inferior candidate.
		  // Apply the best monotherapy.
		  deathtimes=malloc(sizeof(double)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      *(deathtimes+i)=-1;
		      if (*(inferior+i)==0)
			{
			  val=deathtime_multidrug((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
			  
			  *(deathtimes+i)=val;		      		      
			}
		    }	
		  
		  // Debug
		  //printf("ncandidates=%d, after\n",ncandidates);
		  
		  
		  // Choose the candidate with the max death time.
		  minind=-1; minval=-1;
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if (*(deathtimes+i)>minval)
			{
			  minind=*(candidates+i); minval=*(deathtimes+i);
			}
		    }
		  
		  // If there are degenerate death times, then choose the one with minimum incurable population.
		  ndegenerates=0; degenerates=malloc(sizeof(int)*ncandidates);
		  for (i=0; i<=(ncandidates-1); i++)
		    {
		      if ((minval-*(deathtimes+i))<=1e-5)
			{
			  ndegenerates++; *(degenerates+ndegenerates-1)=*(candidates+i);
			}
		    }
		  if (ndegenerates>1)
		    {
		      minind=-1; minval=1e100;
		      for (i=0; i<=(ndegenerates-1); i++)
			{		      		      
			  if (degmode==0)
			    {
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				val+=*((seqnodes+*(degenerates+i))->y+j);
			    }
			  else
			    val=*((seqnodes+*(degenerates+i))->y+ntypes-1);
			  

			  if (val<minval)
			    {
			      minind=*(degenerates+i); minval=val;
			    }
			}
		    }
		  free(degenerates); free(deathtimes);
		}

	 	      	      	    	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); 
	      
	    }
	     
	  free(candidates);
	  
	  for (i=0; i<=(nsteps-1); i++)
	    *(seq+i)=-1;

	  if (minind>=0)
	    {
	      k=(seqnodes+minind)->depth-1;
	      *(seq+k)=minind; i=k-1;
	      //*(seq+nsteps-1)=minind; i=nsteps-2;
	      while (i>=0)
		{
		  *(seq+i)=(seqnodes+*(seq+i+1))->pa;
		  i--;
		}	  
	    }

	  /*
	  // Debug
	  if ((mode==2)&&(popthre>1e10)&&(n==5))
	    {
	      printf("seq\n");
	      for (i=0; i<=4; i++)
		{
		  printf("%d ( ",*(seq+i));
		  for (j=0; j<=(ntypes-1); j++)
		    printf("%.3e ",*((seqnodes+*(seq+i))->y+j));
		  printf(") ");
		  printf("( ");
		  for (j=0; j<=(ndrugs-1); j++)
		    printf("%.2f ",*((seqnodes+*(seq+i))->dosage+j));
		  printf(")\n");
		}
	      exit(0);
	    }
	  */	 
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);

      // Debug
      //printf("n=%d, k=%d, l=%d\n",n,k,l);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  /*
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
	    */

	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j*nintervals+n)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+i*nintervals+n)=*((seqnodes+l)->y+i);



	  maxval=0;
	  //for (i=0; i<=(ntypes-1); i++)
	  //maxval+=*(y+n*ntypes+i);
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+i*nintervals+n);
	
	  pop=maxval; incurablepop=*(y+(ntypes-1)*nintervals+n);

	  // Debug
	  //printf("maxval=%.2e\n",maxval);

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  // Report the survival time to the max time span + timeinterval.
	  // This distinguishes from the cases where the patients survival for the max time span but are not cured.
	  //if (maxval<=0)	
	  if (maxval<1)	    
	    //stopT=*(t+nintervals-1);	  
	    stopT=*(t+nintervals-1)+timeinterval;

	  
	  /*
	  // Debug	  
	  if ((mode==2)&&(popthre>1e10))
	  //if (mode==1)
	  //if (mode==3)
	    {
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		val+=*(y+i*nintervals+n);
	      printf("dp t(%d)=%.2f, dosage=( ",n,*(t+n));
	      for (i=0; i<=(ndrugs-1); i++)
		printf("%.2f ",*(dosage+i*nintervals+n));
	      printf("), y=( ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.4e ",*(y+i*nintervals+n));
	      printf(") %.4e\n",val);
	    }
	  */
	  

	}

      // Debug
      //printf("n=%d, stopT=%.2f\n",n,stopT);

      n++;
    }

  /*
  // If stop<0, then decide whether the patient recovers or not.
  // Check if each final subpopulation < 1.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      k=-1; i=0;
      while ((i<=(ntypes-1))&&(k<0))
	{
	  //if (*(y+i*nintervals+nintervals-2)>1)
	  if (*(y+i*nintervals+nintervals-2)>10)
	    k=i;
	  i++;
	}
      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }
  */

  //if (stopT<0)
  //stopT=*(t+nintervals-1);
    
  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *localdosage, *inity, *finaly, maxtime=nintervals*timeinterval;
      localdosage=malloc(sizeof(double)*ndrugs);
      finaly=malloc(sizeof(double)*ntypes);
      inity=malloc(sizeof(double)*ntypes);
      
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(localdosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,localdosage,finaly,0);
	  flag=1; m=0;
	  while ((m<=(ntypes-1))&&(flag==1))
	    {
	      if (*(finaly+m)>1)
		flag=0;
	      m++;
	    }
	  l++;
	}
      if (flag==1)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
      free(localdosage); free(inity); free(finaly);
    }


  // Release memory.
  free(seq); free(dstates);
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); free(nfrontiers); free(frontiers);
  free(resvec); free(bases);
  
  return stopT;
}


// Compute the power of a number.
int power(int base, int exponent)
{
  int i, val=1;
  for (i=1; i<=exponent; i++)
    val*=base;
  return val;
}
