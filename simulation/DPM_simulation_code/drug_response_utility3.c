/* Utility functions involved in drug responses and optimization.
   Apply dynamic programming to find the optimal sequence(s) of dosages over several steps.
   Exhaust all possible sequences up to a given number of lookahead steps.
   The cost function is the survival time, which is expressed as the sum of the step function values for the total population.
*/

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

extern void padm(int ndim, double *A, double t, int p, double *E);
extern char * getitemval3(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
extern void diagonal(int ndim, double *vec, double *mat);
extern void recurse_two_drug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
extern void recurse_two_drug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);




extern int ndrugs, ntypes, initialtime;
extern double detected, mortal;
extern int maxnfrontiers;


double dp_optimize_two_drug_responses(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage);
struct seqnode * recurse_sequence_two_drug_responses(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, long *curind_pnt);
double deathtime(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals);
double dp_optimize_two_drug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode);
double dp_optimize_two_drug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode);
int pair2_cmp(struct pair2 *p1, struct pair2 *p2);
double deathtime2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates);
double deb_deathtime2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates);
double deathtime3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates);
double dp_optimize_two_drug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode);
double dp_optimize_two_drug_responses5(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode);
struct seqnode * recurse_sequence_two_drug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, long *curind_pnt);
struct seqnode * recurse_sequence_two_drug_responses3(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);



// Find an optimal treatment sequence up to nsteps using dynamic programming.
double dp_optimize_two_drug_responses(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage)
{
  int i, j, k, l, m, n, nseqnodes=0;
  long curind, *seq;
  double *dstates, val, stopT;
  struct seqnode *seqnodes;

  /*
  // Debug
  *(dosage+0)=1; *(dosage+1)=0;
  recurse_two_drug_response_trimrates(x0,T,g0,Sg,a0,Sa,0,timeinterval,dosage,y);
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y+i));
  printf("\n");
  exit(0);
  */
  
  
  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  /*
  // Debug
  for (i=0; i<=(nstates-1); i++)
    printf("%.3f ",*(dstates+i));
  printf("\n");
  */ 

  // Initialize the records.
  stopT=-1;
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      *(dosage+n)=-1;
    }

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);

  seq=malloc(sizeof(long)*nsteps); n=0;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

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
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i);
	    }
	  else
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);
	    }

	  curind=0;	  	  	      	      

	  /*
	  // Debug
	  {
	    double *vec, *vec2;
	    for (i=0; i<=(ntypes-1); i++)
	      *(y+i)=*(x0+i);
	    vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);
	    *(vec+0)=1; *(vec+1)=0;
	    recurse_two_drug_response_trimrates(y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	    for (i=0; i<=(ntypes-1); i++)
	      printf("%.3e ",*(vec2+i));
	    printf("\n");
	    exit(0);
	  }
	  */

	  seqnodes=recurse_sequence_two_drug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind);
	  	  

	  /*
	  // Debug
	  printf("nseqnodes=%d\n",nseqnodes);
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      printf("seqnode[%d], index=%d, pa=%d, depth=%d, dosage=( ",i,(seqnodes+i)->index,(seqnodes+i)->pa,(seqnodes+i)->depth);
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%.3f ",*((seqnodes+i)->dosage+j));
	      printf("), pop=( ");
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.3e ",*((seqnodes+i)->y+j));
	      printf("), totalpop= ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((seqnodes+i)->y+j);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */

	  // Find the max depth of the seqnodes.
	  maxdepth=-1;
	  for (i=0; i<=(nseqnodes-1); i++)
	    maxdepth=max(maxdepth,(seqnodes+i)->depth);
	  
	  // Debug
	  //printf("maxdepth=%d\n",maxdepth);
	  

	  // Find the candidate terminal seqnodes.
	  ncandidates=0; candidates=malloc(sizeof(long));
	  

	  // Find the terminal nodes with zero population.
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      if ((seqnodes+i)->nchildren==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  //if (val<=1e-10)
		  if (val<1)
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	    }
	  

	  /*
	  // If no sequences can eradicate the total population, then choose the ones with the longest surviving time.
	  if (ncandidates==0)
	    {
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0, c3=0, c4=0;
		  
		  // Find the terminal nodes where R12 does not pose a threat.
		  // This means either R12<0.1 or (R12<1 and current R12 <= 1.01 * previous R12).	      
		  // Minimize total population.
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  if (val<mortal)
		    c4=1;
		  if  (*((seqnodes+i)->y+3)<0.1)
		    c3=1;
		  else if ((*((seqnodes+i)->y+3)>=0.1)&&(*((seqnodes+i)->y+3)<=(*((seqnodes+0)->y+3)*1.01)))
		    c3=1;
		  if ((c1==1)&&(c2==1)&&(c3==1)&&(c4==1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		      opt=0;
		    }
		}
	      
	      // If no candidates satisfy these criteria, then discard this requirement.  Minimize R12.
	      if (ncandidates==0)
		{
		  for (i=0; i<=(nseqnodes-1); i++)
		    {
		      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			{
			  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
			  *(candidates+ncandidates-1)=i;
			}
		    }
		  opt=1;
		}
	    }

	  if (opt<0)
	    opt=0;


	  // Debug
	  printf("opt=%d\n",opt);

	  // Choose the optimal sequence that minimizes the terminal population size.
	  // opt=0: minimize total population.
	  // opt=1: minimize R12 population.
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      if (opt==0)
		{
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+*(candidates+i))->y+j);
		}	      
	      else
		val=*((seqnodes+*(candidates+i))->y+3);

	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  */


	  /*
	  // If no sequences can eradicate the total population, then choose the ones with the longest surviving time.
	  if (ncandidates==0)
	    {
	      // Rule out the terminal nodes with total population >= mortal or R12 size > 1 or R12 size >= 0.5.
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  int c1=0, c2=0, c3=0, c4=0;

		  val=0;
		  for (j=0; j<=(ntypes-1); j++)		
		    val+=*((seqnodes+i)->y+j);
		  
		  if ((seqnodes+i)->nchildren==0)
		    c1=1;
		  if ((seqnodes+i)->depth==maxdepth)
		    c2=1;
		  //if (*((seqnodes+i)->y+3)<0.2)
		  //c3=1;
		  //else if ((*((seqnodes+i)->y+3)>=0.2)&&(*((seqnodes+i)->y+3)<=(*((seqnodes+0)->y+3)*1.2)))
		  //c3=1;
		  c3=1;
		  if (val<mortal)
		    c4=1;

		  
		  if ((c1==1)&&(c2==1)&&(c3==1)&&(c4==1))
		  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth)&&(*((seqnodes+i)->y+3)<1))		  
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		  

		  //
		  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth)&&(*((seqnodes+i)->y+3)<0.5)&&(val<mortal))
		    //{
		      //ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      //*(candidates+ncandidates-1)=i;
		    //}
		  //
		}
	      // If no candidates satisfy these criteria, then discard this requirement.
	      if (ncandidates==0)
		{
		  for (i=0; i<=(nseqnodes-1); i++)
		    {
		      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			{
			  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
			  *(candidates+ncandidates-1)=i;
			}
		    }
		}
	    }

	  // Choose the optimal sequence that minimizes the terminal population size.
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);
	      
	      val=*((seqnodes+*(candidates+i))->y+3);

	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  */

	  
	  
	  // If no sequences can eradicate the total population, then choose the ones with the longest surviving time.
	  if (ncandidates==0)
	    {
	      // Rule out the terminal nodes with R12 size > 1.	  
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth)&&(*((seqnodes+i)->y+3)<=1e-10))
		  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth)&&(*((seqnodes+i)->y+3)<1))
		    {
		      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
		      *(candidates+ncandidates-1)=i;
		    }
		}
	      
	      // If all candidates have R12 > 0, then discard this requirement.
	      if (ncandidates==0)
		{
		  for (i=0; i<=(nseqnodes-1); i++)
		    {
		      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			{
			  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
			  *(candidates+ncandidates-1)=i;
			}
		    }
		}
	    }

	  // Debug
	  //printf("ncandidates=%d\n",ncandidates);
	  

	  // Choose the optimal sequence that minimizes the terminal population size.
	  minind=-1; minval=1e100;
	  for (i=0; i<=(ncandidates-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		val+=*((seqnodes+*(candidates+i))->y+j);

	      // Debug
	      //val=*((seqnodes+*(candidates+i))->y+3);


	      if (val<minval)
		{
		  minind=*(candidates+i); minval=val;
		}
	    }
	  


	  // Debug
	  //printf("minind=%d\n",minind);
	  //exit(0);

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
	  for (i=0; i<=(nsteps-1); i++)
	    {
	      j=*(seq+i);
	      printf("seqnode[%d], index=%ld, pa=%d, depth=%d, dosage=( ",j,(seqnodes+j)->index,(seqnodes+j)->pa,(seqnodes+j)->depth);
	      for (k=0; k<=(ndrugs-1); k++)
		printf("%.3f ",*((seqnodes+j)->dosage+k));
	      printf("), pop=( ");
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((seqnodes+j)->y+k));
	      printf("), totalpop=");
	      val=0;
	      for (k=0; k<=(ntypes-1); k++)
		val+=*((seqnodes+j)->y+k);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  /*
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf(")\n");
	  */

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
      
	  /*
	  // Debug
	  printf("pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf(")\n");
	  */

	  maxval=0;
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+n*ntypes+i);

	  // Debug
	  //printf("totalpop=%.3e\n",maxval);
	  

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);

	  /*
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf("), pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf("), totalpop=%.3e\n",maxval);
	  //exit(0);
	  */
	}

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
  free(seqnodes);  
  
  return stopT;
}


// Recursively exhaust drug responses over all sequences up to depth=nsteps.
struct seqnode * recurse_sequence_two_drug_responses(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, long *curind_pnt)
{
  int i, j, k, l, m, n;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
  
  // Return when the node depth>=nsteps.
  if ((seqnodes+curind)->depth>=nsteps)
    return seqnodes;

  // Return when the total population exceeds the mortal threshold.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val+=*((seqnodes+curind)->y+i);

  // Debug
  //printf("curind=%ld, depth=%d, pop=%.4e\n",curind,(seqnodes+curind)->depth,val);
  //for (i=0; i<=(ntypes-1); i++)
  //printf("%.4e ",*((seqnodes+curind)->y+i));
  //printf("\n");

 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;

  
  // Spin off each dosage combination.
  seqnodes=realloc(seqnodes,sizeof(struct seqnode)*(nseqnodes+nstates));
  vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);

  /*
  // Debug
  *(vec+0)=1; *(vec+1)=0;
  recurse_two_drug_response_trimrates((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(vec2+i));
  printf("\n");
  exit(0);
  */


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
      *(vec+0)=*(dstates+n); *(vec+1)=1-*(vec+0);
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);

      // Debug
      //printf("before\n");

      recurse_two_drug_response_trimrates((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);

      // Debug
      //for (i=0; i<=(ntypes-1); i++)
      //*(vec2+i)=1e5;


      /*
      // Debug
      printf("dosage: ");
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.3f ",*(vec+i));
      printf("\n");
      printf("pop: ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.3e ",*(vec2+i));
      printf("\n");
      exit(0);
      */


      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=recurse_sequence_two_drug_responses(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind2);
    }

  return seqnodes;
}



// Estimate the upper bound of death time given the vector of initial population.
// Impose max dosage to all drugs.
// Use simple matrix product to estimate the population at time t.
// Apply binary search to find the time step prior to death.
double deathtime(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals)
{
  int i, j, k, l, m, n, cnt, lwd, upp, mid, ind;
  double tspan, *F, *vec, *template, *constF, *vec2, *template2, *template3, *tempF, *vec3;
  double stopT, val, val2;

  // If each subpopulation size < 1, then death time is max time interval.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val=max(val,*(x0+i));
  if (val<1)
    {
      stopT=(nintervals-1)*timeinterval;
      return stopT;
    }

  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);

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

  n=0; 
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  for (i=0; i<=(ndrugs-1); i++)
    *(vec+i)=1.0;

  // Apply binary search to find the time to death.
  // Use matrix exponential to estimate population.
  // Impose max dosage to all drugs.
  lwd=0; upp=nintervals-1; ind=-1;
  while ((lwd<upp)&&(ind<0))
    {
      mid=(int)((lwd+upp)/2); tspan=timeinterval*mid;
    
      // Debug
      //mid=10; tspan=timeinterval*mid;
      
  
      // Estimate the population at time mid*timeinterval.
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
	  *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*tspan;
		   		  		  		  		      
      padm(ntypes,tempF,1,7,template);	
      mat_multiply(ntypes,ntypes,1,template,x0,vec2);
      

      // Check if the total population >= mortal.
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(vec2+i);

      
      // Debug
      //printf("lwd=%d, mid=%d, upp=%d, vec2: ",lwd,mid,upp);
      //for (i=0; i<=(ntypes-1); i++)
      //printf("%.4e ",*(vec2+i));
      //printf(", total %.4e\n",val);

      
      if (val>=mortal)
	upp=mid;
      else
	lwd=mid;
      
      if (lwd==upp)
	ind=lwd;
      else if (lwd==(upp-1))
	ind=lwd;
    }

  stopT=timeinterval*ind;

  free(F); free(tempF); free(template); free(template2); 
  free(template3); free(constF); free(vec2);

  return stopT;
}


  
// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
double dp_optimize_two_drug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode)
{
  int i, j, k, l, m, n, nseqnodes=0;
  long curind, *seq;
  double *dstates, val, stopT;
  struct seqnode *seqnodes;
    
  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  // Initialize the records.
  stopT=-1;
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      *(dosage+n)=-1;
    }

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);


  /*
  // Debug
  *(y+0)=2.441e+03; *(y+1)=5.321e+04; *(y+2)=1.137e+04; *(y+3)=1.218e+00;
  *(y+0)=7.529e+11; *(y+1)=4.074e+04; *(y+2)=1.615e+07; *(y+3)=9.263e-02;
  //val=deathtime(y,T,g0,Sg,a0,Sa,timeinterval,nintervals);
  val=deathtime(y,T,g0,Sg,a0,Sa,1,nintervals*45);
  printf("t1=%f\n",val);
  *(y+0)=2.576e+00; *(y+1)=6.779e+05; *(y+2)=9.336e-01; *(y+3)=8.292e-01;
  *(y+0)=7.995e+11; *(y+1)=3.999e+04; *(y+2)=1.739e+07; *(y+3)=9.142e-02;
  //val=deathtime(y,T,g0,Sg,a0,Sa,timeinterval,nintervals);
  val=deathtime(y,T,g0,Sg,a0,Sa,1,nintervals*45);
  printf("t2=%f\n",val);
  exit(0);
  */


  seq=malloc(sizeof(long)*nsteps); n=0;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

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
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i);
	    }
	  else
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);
	    }

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_two_drug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind);
	  	  

	  /*
	  // Debug
	  printf("nseqnodes=%d\n",nseqnodes);
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      printf("seqnode[%d], index=%d, pa=%d, depth=%d, dosage=( ",i,(seqnodes+i)->index,(seqnodes+i)->pa,(seqnodes+i)->depth);
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%.3f ",*((seqnodes+i)->dosage+j));
	      printf("), pop=( ");
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.3e ",*((seqnodes+i)->y+j));
	      printf("), totalpop= ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((seqnodes+i)->y+j);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */

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
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates;
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

	      // Debug
	      //printf("ncandidates=%d\n",ncandidates); 

	      // Debug
	      //minind=*(candidates+0);

	      
	      // Sort candidate sequences by total population and by R12 population.
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


	      // Debug
	      //printf("ind=%d, pop=%.4e\n",(pairs1+0)->i,(pairs1+0)->s);

	      
	      // Debug
	      //printf("l1=%d, pair1=(%d,%.4e), l2=%d, pair2=(%d,%.4e)\n",*(lookup1+0),(pairs1+*(lookup1+0))->i,(pairs1+*(lookup1+0))->s,*(lookup2+0),(pairs2+*(lookup2+0))->i,(pairs2+*(lookup2+0))->s);
		  
	      // Find the candidates inferior to others.
	      // Direct comparison of all candidate pairs take n^2 time, n=ncandidates.
	      // To reduce computing time use probabilistic arguments.
	      // Suppose a candidate's total population > n1 other candidates, and its R12 population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and R12 population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in R12 population) =
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

		  // Debug
		  //val=0;
		  //for (j=0; j<=(ntypes-1); j++)		
		  //val+=*((seqnodes+*(candidates+i))->y+j);
		  //printf("n1=%d, n2=%d, p=%f, v1=%.4e, v2=%.4e\n",n1,n2,p,val,*((seqnodes+*(candidates+i))->y+3)); exit(0);

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
	      

	      // Debug
	      //printf("inferior[%d]=%d, inferior[%d]=%d\n",161051,*(inferior+161051),102486,*(inferior+102486));
		
	     	      
	      // Evaluate the expected death time for each non-inferior candidate.
	      deathtimes=malloc(sizeof(double)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  *(deathtimes+i)=-1;
		  if (*(inferior+i)==0)
		    {
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,timeinterval,nintervals);

		      val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval);

		      //val=deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);


		      *(deathtimes+i)=val;		      		      
		    }
		}
	
	      // Debug
	      //printf("deathtime[%d]=%f, deathtime[%d]=%f\n",161051,*(deathtimes+161051),102486,*(deathtimes+102486));
	
	      //printf("%.4e %.4e %.4e %.4e\n",*((seqnodes+*(candidates+161051))->y+0),*((seqnodes+*(candidates+161051))->y+1),*((seqnodes+*(candidates+161051))->y+2),*((seqnodes+*(candidates+161051))->y+3));
	      
	      //printf("%.4e %.4e %.4e %.4e\n",*((seqnodes+*(candidates+102486))->y+0),*((seqnodes+*(candidates+102486))->y+1),*((seqnodes+*(candidates+102486))->y+2),*((seqnodes+*(candidates+102486))->y+3));
	      


	      // Choose the candidate with the max death time.
	      minind=-1; minval=-1;
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(deathtimes+i)>minval)
		    {
		      minind=*(candidates+i); minval=*(deathtimes+i);
		    }
		}

	      // If there are degenerate death times, then choose the one with minimum R12 population.	      
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
	      free(degenerates);

	      
	      // Debug
	      //for (i=0; i<=(ncandidates-1); i++)
	      //if (*(candidates+i)==minind)
	      //printf("minind=%d\n",i);


	      /*
	      // Debug
	      printf("ndegenerates=%d\n",ndegenerates);
	      printf("( ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		{
		  val+=*((seqnodes+minind)->y+j);
		  printf("%.4e ",*((seqnodes+minind)->y+j));
		}
	      printf("), total %.4e\n",val);
	      */
	      
	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); free(deathtimes);
	      
	    }
	     

	  // Debug
	  //printf("minind=%d\n",minind);
	  //exit(0);

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
	  for (i=0; i<=(nsteps-1); i++)
	    {
	      j=*(seq+i);
	      printf("seqnode[%d], index=%ld, pa=%d, depth=%d, dosage=( ",j,(seqnodes+j)->index,(seqnodes+j)->pa,(seqnodes+j)->depth);
	      for (k=0; k<=(ndrugs-1); k++)
		printf("%.3f ",*((seqnodes+j)->dosage+k));
	      printf("), pop=( ");
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((seqnodes+j)->y+k));
	      printf("), totalpop=");
	      val=0;
	      for (k=0; k<=(ntypes-1); k++)
		val+=*((seqnodes+j)->y+k);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
      
	  maxval=0;
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+n*ntypes+i);
	

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);

	  
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf("), pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf("), totalpop=%.3e\n",maxval);
	  //exit(0);
	  
	}

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
  free(seqnodes);  
  
  return stopT;
}


// Compare two pair2 structures.
int pair2_cmp(struct pair2 *p1, struct pair2 *p2)
{
  int retval=0;

  if (p1->s<p2->s)
    retval=-1;
  else if (p1->s>p2->s)
    retval=1;
  return retval;
}


// Estimate the upper bound of death time given the vector of initial population.
// Impose max dosage to all drugs.
// Use simple matrix product to estimate the population at time t.
// Apply binary search to find the time step prior to death.
// Difference from deathtime: apply nstates drug combinations and choose the max survival time among them.
double deathtime2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates)
{
  int i, j, k, l, m, n, cnt, lwd, upp, mid, ind, tind;
  double tspan, *F, *vec, *template, *constF, *vec2, *template2, *template3, *tempF, *vec3;
  double stopT, val, val2, *dstates;
  
  // If each subpopulation size < 1, then death time is max time interval.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val=max(val,*(x0+i));
  if (val<1)
    {
      stopT=(nintervals-1)*timeinterval;
      return stopT;
    }

  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);

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

  vec=malloc(sizeof(double)*ndrugs); stopT=-1;

  
  // For each drug combination, find the survival time.
  for (tind=0; tind<=(nstates-1); tind++)
    {
      n=0;   
      
      //for (i=0; i<=(ndrugs-1); i++)
      //*(vec+i)=1.0;

      // Vary the drug combination.
      *(vec+0)=*(dstates+tind); *(vec+1)=1-*(vec+0);

      // Apply binary search to find the time to death.
      // Use matrix exponential to estimate population.
      // Impose max dosage to all drugs.
      lwd=0; upp=nintervals-1; ind=-1;
      while ((lwd<upp)&&(ind<0))
	{
	  mid=(int)((lwd+upp)/2); tspan=timeinterval*mid;
	  
	  // Debug
	  //mid=10; tspan=timeinterval*mid;
	  
	  
	  // Estimate the population at time mid*timeinterval.
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
	      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*tspan;
	  
	  padm(ntypes,tempF,1,7,template);	
	  mat_multiply(ntypes,ntypes,1,template,x0,vec2);
	  
	  
	  // Check if the total population >= mortal.
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(vec2+i);
	  
      
	  // Debug
	  //printf("lwd=%d, mid=%d, upp=%d, vec2: ",lwd,mid,upp);
	  //for (i=0; i<=(ntypes-1); i++)
	  //printf("%.4e ",*(vec2+i));
	  //printf(", total %.4e\n",val);
      
	  if (val>=mortal)
	    upp=mid;
	  else
	    lwd=mid;
	  
	  if (lwd==upp)
	    ind=lwd;
	  else if (lwd==(upp-1))
	    ind=lwd;
	}

      //stopT=timeinterval*ind;

      val=timeinterval*ind;
      stopT=max(val,stopT);
    }


  free(F); free(tempF); free(template); free(template2); 
  free(template3); free(constF); free(vec2); free(dstates);

  return stopT;
}



// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
double dp_optimize_two_drug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode)
{
  int i, j, k, l, m, n, nseqnodes=0;
  long curind, *seq;
  double *dstates, val, stopT;
  struct seqnode *seqnodes;
    
  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  // Initialize the records.
  stopT=-1;
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      *(dosage+n)=-1;
    }

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);

  seq=malloc(sizeof(long)*nsteps); n=0;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

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
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i);
	    }
	  else
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);
	    }

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_two_drug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind);
	  	  

	  /*
	  // Debug
	  printf("nseqnodes=%d\n",nseqnodes);
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      printf("seqnode[%d], index=%d, pa=%d, depth=%d, dosage=( ",i,(seqnodes+i)->index,(seqnodes+i)->pa,(seqnodes+i)->depth);
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%.3f ",*((seqnodes+i)->dosage+j));
	      printf("), pop=( ");
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.3e ",*((seqnodes+i)->y+j));
	      printf("), totalpop= ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((seqnodes+i)->y+j);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */

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
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates;
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

	      // Debug
	      //printf("ncandidates=%d\n",ncandidates); 

	      /*
	      // Debug
	      for (i=0; i<=(nseqnodes-1); i++)
		{
		  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
		    {
		      if ((fabs(*((seqnodes+i)->dosage+0)-0)<0.01)&&(fabs(*((seqnodes+i)->dosage+0)-1)<0.01))
			{
			  k=(seqnodes+i)->pa;
			  if ((fabs(*((seqnodes+k)->dosage+0)-0)<0.01)&&(fabs(*((seqnodes+k)->dosage+0)-1)<0.01))
			    {
			      k=(seqnodes+i)->pa;
			      if ((fabs(*((seqnodes+k)->dosage+0)-0)<0.01)&&(fabs(*((seqnodes+k)->dosage+0)-1)<0.01))
				{
				  k=(seqnodes+i)->pa;
				  if ((fabs(*((seqnodes+k)->dosage+0)-0)<0.01)&&(fabs(*((seqnodes+k)->dosage+0)-1)<0.01))
				    {
				      k=(seqnodes+i)->pa;
				      if ((fabs(*((seqnodes+k)->dosage+0)-0.75)<0.01)&&(fabs(*((seqnodes+k)->dosage+0.25)-1)<0.01))
	      */	


	     
	      
	      // Sort candidate sequences by total population and by R12 population.
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
	      // Suppose a candidate's total population > n1 other candidates, and its R12 population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and R12 population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in R12 population) =
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
	      	
	      /*
	      // Debug
	      k=0;
	      for (i=0; i<=(ncandidates-1); i++)
		if (*(inferior+i)==0)
		  k++;
	      printf("ncandidates=%d\n",k);
	      */
	     	      
	      // Evaluate the expected death time for each non-inferior candidate.
	      deathtimes=malloc(sizeof(double)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  *(deathtimes+i)=-1;
		  if (*(inferior+i)==0)
		    {
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,timeinterval,nintervals);
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval);
		      
		      val=deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		      
		      //val=deathtime3((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);


		      *(deathtimes+i)=val;		      		      
		    }
		}	

	      /*
	      // Debug
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(inferior+i)==0)
		    {
		      printf("ind %d, pop ",i);
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.4e ",*((seqnodes+*(candidates+i))->y+j));
		      printf(", deathtime=%.0f\n",*(deathtimes+i));
		    }
		}
	      //exit(0);
	      */

	      

	      // Choose the candidate with the max death time.
	      minind=-1; minval=-1;
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(deathtimes+i)>minval)
		    {
		      minind=*(candidates+i); minval=*(deathtimes+i);
		    }
		}

	      // If there are degenerate death times, then choose the one with minimum R12 population.	      
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
	      free(degenerates);
	 
	      
	      // Debug
	      //if (n==0)
	      //minind=*(candidates+1376869);

	      /*
	      // Debug
	      if (n==0)
		{
		  i=1376869;
		  i=316206;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		  i=1376886;
		  i=797819;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		}
	      */

	      /*
	      // Debug
	      if (n==0)
		{
		  double *vec, *vec2, *y0;
		  y0=malloc(sizeof(double)*ntypes);
		  vec=malloc(sizeof(double)*ndrugs);
		  vec2=malloc(sizeof(double)*ntypes);
		  *(vec+0)=0.75; *(vec+1)=0.2;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  printf("\n\n");

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2); 
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		  *(vec+0)=0.75; *(vec+1)=0.25;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  exit(0);

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		}
	      */

	      /*
	      // Debug
	      printf("ndegenerates=%d\n",ndegenerates);
	      printf("( ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		{
		  val+=*((seqnodes+minind)->y+j);
		  printf("%.4e ",*((seqnodes+minind)->y+j));
		}
	      printf("), total %.4e\n",val);
	      */
	      
	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); free(deathtimes);
	      
	    }
	     

	  // Debug
	  //printf("minind=%d\n",minind);
	  //exit(0);

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
	  for (i=0; i<=(nsteps-1); i++)
	    {
	      j=*(seq+i);
	      printf("seqnode[%d], index=%ld, pa=%d, depth=%d, dosage=( ",j,(seqnodes+j)->index,(seqnodes+j)->pa,(seqnodes+j)->depth);
	      for (k=0; k<=(ndrugs-1); k++)
		printf("%.3f ",*((seqnodes+j)->dosage+k));
	      printf("), pop=( ");
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((seqnodes+j)->y+k));
	      printf("), totalpop=");
	      val=0;
	      for (k=0; k<=(ntypes-1); k++)
		val+=*((seqnodes+j)->y+k);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
      
	  maxval=0;
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+n*ntypes+i);
	

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);

	  /*
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf("), pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf("), totalpop=%.3e\n",maxval);
	  //exit(0);
	  */
	}

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
  free(seqnodes);  
  
  return stopT;
}



// Estimate the upper bound of death time given the vector of initial population.
// Impose max dosage to all drugs.
// Use simple matrix product to estimate the population at time t.
// Apply binary search to find the time step prior to death.
// Difference from deathtime: apply nstates drug combinations and choose the max survival time among them.
// Difference from deathtime2: debug version.
double deb_deathtime2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates)
{
  int i, j, k, l, m, n, cnt, lwd, upp, mid, ind, tind;
  double tspan, *F, *vec, *template, *constF, *vec2, *template2, *template3, *tempF, *vec3;
  double stopT, val, val2, *dstates;
  
  // If each subpopulation size < 1, then death time is max time interval.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val=max(val,*(x0+i));
  if (val<1)
    {
      stopT=(nintervals-1)*timeinterval;
      return stopT;
    }

  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);

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

  vec=malloc(sizeof(double)*ndrugs); stopT=-1;

  
  // For each drug combination, find the survival time.
  for (tind=0; tind<=(nstates-1); tind++)
    {
      n=0;   
      
      //for (i=0; i<=(ndrugs-1); i++)
      //*(vec+i)=1.0;

      // Vary the drug combination.
      *(vec+0)=*(dstates+tind); *(vec+1)=1-*(vec+0);

      // Apply binary search to find the time to death.
      // Use matrix exponential to estimate population.
      // Impose max dosage to all drugs.
      lwd=0; upp=nintervals-1; ind=-1;
      while ((lwd<upp)&&(ind<0))
	{
	  mid=(int)((lwd+upp)/2); tspan=timeinterval*mid;
	  
	  // Debug
	  //mid=10; tspan=timeinterval*mid;
	  
	  
	  // Estimate the population at time mid*timeinterval.
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
	      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*tspan;
	  
	  padm(ntypes,tempF,1,7,template);	
	  mat_multiply(ntypes,ntypes,1,template,x0,vec2);
	  
	  
	  // Check if the total population >= mortal.
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(vec2+i);
	  
      
	  // Debug
	  //printf("lwd=%d, mid=%d, upp=%d, vec2: ",lwd,mid,upp);
	  //for (i=0; i<=(ntypes-1); i++)
	  //printf("%.4e ",*(vec2+i));
	  //printf(", total %.4e\n",val);
      
	  if (val>=mortal)
	    upp=mid;
	  else
	    lwd=mid;
	  
	  if (lwd==upp)
	    ind=lwd;
	  else if (lwd==(upp-1))
	    ind=lwd;
	}

      //stopT=timeinterval*ind;

      val=timeinterval*ind;

      // Debug
      //printf("%.2f %.2f %.0f\n",*(vec+0),*(vec+1),val);

      stopT=max(val,stopT);
    }


  free(F); free(tempF); free(template); free(template2); free(vec);
  free(template3); free(constF); free(vec2); free(dstates);

  return stopT;
}



// Estimate the upper bound of death time given the vector of initial population.
// Impose max dosage to all drugs.
// Use simple matrix product to estimate the population at time t.
// Apply binary search to find the time step prior to death.
// Difference from deathtime: apply nstates drug combinations and choose the max survival time among them.
// Difference from deathtime2: apply recurse_two_drug_response_trimrates to estimate population size.
// Using matrix exponential will amplify the effect of R12 fractional population.
double deathtime3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int nintervals, int nstates)
{
  int i, j, k, l, m, n, cnt, lwd, upp, mid, ind, tind, maxind=-1;
  double tspan, *vec;
  double stopT=-1, val, val2, *dstates, *y;
  
  // If each subpopulation size < 1, then death time is max time interval.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val=max(val,*(x0+i));
  if (val<1)
    {
      stopT=(nintervals-1)*timeinterval;
      return stopT;
    }

  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }
  y=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs);

  // For each drug combination, find the survival time.
  for (tind=0; tind<=(nstates-1); tind++)
    {
      n=0;   
      
      //for (i=0; i<=(ndrugs-1); i++)
      //*(vec+i)=1.0;

      // Vary the drug combination.
      *(vec+0)=*(dstates+tind); *(vec+1)=1-*(vec+0);

      // Debug
      //printf("dosage (%.2f %.2f)\n",*(vec+0),*(vec+1));

      // Apply binary search to find the time to death.
      // Use matrix exponential to estimate population.
      // Impose max dosage to all drugs.
      lwd=0; upp=nintervals-1; ind=-1;
      while ((lwd<upp)&&(ind<0))
	{
	  mid=(int)((lwd+upp)/2); tspan=timeinterval*mid;	  
	  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,tspan,vec,y);
	  
	  /*
	  // Debug
	  if (*(vec+1)>=0.99)
	    {
	      printf("%d %d %d tspan=%.0f, y: ",lwd,mid,upp,tspan);
	      val=0;
	      for (i=0; i<=(ntypes-1); i++)
		{
		  val+=*(y+i);
		  printf("%.4e ",*(y+i));
		}
	      printf("%.4e ",val);
	      printf("\n");
	    }
	  */	  


	  // Check if the total population >= mortal.
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(y+i);
	  
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

    // Debug
    //printf("%.2f %.2f %.2f\n",*(vec+0),*(vec+1),val);


    // Debug
    if (val>stopT)
      maxind=tind;

    stopT=max(val,stopT);
    
   }

  // Debug
  //printf("max dosage=(%.2f,%.2f), stopT=%.2f\n",*(dstates+maxind),1-*(dstates+maxind),stopT);

  free(vec); free(dstates); free(y);

  return stopT;
}



// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
double dp_optimize_two_drug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode)
{
  int i, j, k, l, m, n, nseqnodes=0;
  long curind, *seq;
  double *dstates, val, stopT;
  struct seqnode *seqnodes;

  /*
  // Debug
  {
    double *dosage, *responses, val;

    dosage=malloc(sizeof(double)*ndrugs);
    responses=malloc(sizeof(double)*ntypes);
    *(dosage+0)=0.85; *(dosage+1)=0.05;
    recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,45,dosage,responses);
    val=0;
    printf("dosage=( ");
    for (i=0; i<=(ndrugs-1); i++)
      printf("%.2f ",*(dosage+i));
    printf("), pop=( ");
    for (i=0; i<=(ntypes-1); i++)
      {
	val+=*(responses+i);
	printf("%.3e ",*(responses+i));
      }
    printf("), totalpop=%.3e\n",val);
    *(dosage+0)=0.4; *(dosage+1)=0.6;
    recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,45,dosage,responses);
    val=0;
    printf("dosage=( ");
    for (i=0; i<=(ndrugs-1); i++)
      printf("%.2f ",*(dosage+i));
    printf("), pop=( ");
    for (i=0; i<=(ntypes-1); i++)
      {
	val+=*(responses+i);
	printf("%.3e ",*(responses+i));
      }
    printf("), totalpop=%.3e\n",val);
    exit(0);
  }
  */

    
  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  // Initialize the records.
  stopT=-1;
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      *(dosage+n)=-1;
    }

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);

  seq=malloc(sizeof(long)*nsteps); n=0;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

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
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i);
	    }
	  else
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);
	    }

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_two_drug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind);
	  	  

	  /*
	  // Debug
	  printf("nseqnodes=%d\n",nseqnodes);
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      printf("seqnode[%d], index=%d, pa=%d, depth=%d, dosage=( ",i,(seqnodes+i)->index,(seqnodes+i)->pa,(seqnodes+i)->depth);
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%.3f ",*((seqnodes+i)->dosage+j));
	      printf("), pop=( ");
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.3e ",*((seqnodes+i)->y+j));
	      printf("), totalpop= ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((seqnodes+i)->y+j);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */

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
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates;
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

	      // Debug
	      //printf("ncandidates=%d\n",ncandidates); 
	     	     
	      
	      // Sort candidate sequences by total population and by R12 population.
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
	      // Suppose a candidate's total population > n1 other candidates, and its R12 population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and R12 population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in R12 population) =
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

	      /*
	      // Debug
	      k=0;
	      for (i=0; i<=(ncandidates-1); i++)
		if (*(inferior+i)==0)
		  k++;
	      printf("ncandidates=%d\n",k);
	      */
	     	      
	      // Evaluate the expected death time for each non-inferior candidate.
	      deathtimes=malloc(sizeof(double)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  *(deathtimes+i)=-1;
		  if (*(inferior+i)==0)
		    {
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,timeinterval,nintervals);
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval);
		      
		      //val=deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		      
		      val=deathtime3((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);


		      *(deathtimes+i)=val;		      		      
		    }
		}	

	      /*
	      // Debug
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(inferior+i)==0)
		    {
		      printf("ind %d, pop ",i);
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.4e ",*((seqnodes+*(candidates+i))->y+j));
		      printf(", deathtime=%.0f\n",*(deathtimes+i));
		    }
		}
	      //exit(0);
	      */
	      

	      // Choose the candidate with the max death time.
	      minind=-1; minval=-1;
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(deathtimes+i)>minval)
		    {
		      minind=*(candidates+i); minval=*(deathtimes+i);
		    }
		}

	      // If there are degenerate death times, then choose the one with minimum R12 population.	      
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
	      free(degenerates);
	 
	      
	      // Debug
	      //if (n==0)
	      //minind=*(candidates+1376869);

	      /*
	      // Debug
	      if (n==0)
		{
		  i=1376869;
		  i=316206;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		  i=1376886;
		  i=797819;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		}
	      */

	      /*
	      // Debug
	      if (n==0)
		{
		  double *vec, *vec2, *y0;
		  y0=malloc(sizeof(double)*ntypes);
		  vec=malloc(sizeof(double)*ndrugs);
		  vec2=malloc(sizeof(double)*ntypes);
		  *(vec+0)=0.75; *(vec+1)=0.2;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  printf("\n\n");

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2); 
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		  *(vec+0)=0.75; *(vec+1)=0.25;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  exit(0);

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		}
	      */

	      /*
	      // Debug
	      printf("ndegenerates=%d\n",ndegenerates);
	      printf("( ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		{
		  val+=*((seqnodes+minind)->y+j);
		  printf("%.4e ",*((seqnodes+minind)->y+j));
		}
	      printf("), total %.4e\n",val);
	      */
	      
	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); free(deathtimes);
	      
	    }
	     

	  // Debug
	  //printf("minind=%d\n",minind);
	  //exit(0);

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
	  for (i=0; i<=(nsteps-1); i++)
	    {
	      j=*(seq+i);
	      printf("seqnode[%d], index=%ld, pa=%d, depth=%d, dosage=( ",j,(seqnodes+j)->index,(seqnodes+j)->pa,(seqnodes+j)->depth);
	      for (k=0; k<=(ndrugs-1); k++)
		printf("%.3f ",*((seqnodes+j)->dosage+k));
	      printf("), pop=( ");
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((seqnodes+j)->y+k));
	      printf("), totalpop=");
	      val=0;
	      for (k=0; k<=(ntypes-1); k++)
		val+=*((seqnodes+j)->y+k);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
      
	  maxval=0;
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+n*ntypes+i);
	

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);

	  /*
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf("), pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf("), totalpop=%.3e\n",maxval);
	  //exit(0);
	  */
	}

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
  free(seqnodes);  
  
  return stopT;
}


// Find an optimal treatment sequence up to nsteps using dynamic programming.
// Difference from dp_optimize_two_drug_responses: maximize the estimated death time.
// degmode=0: choose a degenerate dosage to minimize total population.
// degmode=1: choose a degenerate dosage to minimize R12 population.
// Difference from dp_optimize_two_drug_responses2: choose a better estimation of death times.
// Difference from dp_optimize_two_drug_responses3: (1)reduce errors in recurse_two_drug_response_trimrates, (2)choose an even better estimation of death times.
// Difference from dp_optimize_two_drug_responses4: apply branch and bound to trim unnecessary nodes in recursion.
double dp_optimize_two_drug_responses5(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, int nstates, double *dosage, int degmode)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers;
  long curind, *seq;
  double *dstates, val, stopT, *frontiers;
  struct seqnode *seqnodes;

  /*
  // Debug
  {
    double *dosage, *responses, val;

    dosage=malloc(sizeof(double)*ndrugs);
    responses=malloc(sizeof(double)*ntypes);
    *(dosage+0)=0.85; *(dosage+1)=0.05;
    recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,45,dosage,responses);
    val=0;
    printf("dosage=( ");
    for (i=0; i<=(ndrugs-1); i++)
      printf("%.2f ",*(dosage+i));
    printf("), pop=( ");
    for (i=0; i<=(ntypes-1); i++)
      {
	val+=*(responses+i);
	printf("%.3e ",*(responses+i));
      }
    printf("), totalpop=%.3e\n",val);
    *(dosage+0)=0.4; *(dosage+1)=0.6;
    recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,45,dosage,responses);
    val=0;
    printf("dosage=( ");
    for (i=0; i<=(ndrugs-1); i++)
      printf("%.2f ",*(dosage+i));
    printf("), pop=( ");
    for (i=0; i<=(ntypes-1); i++)
      {
	val+=*(responses+i);
	printf("%.3e ",*(responses+i));
      }
    printf("), totalpop=%.3e\n",val);
    exit(0);
  }
  */

    
  dstates=malloc(sizeof(double)*nstates); val=(double)1/(double)(nstates-1);
  *(dstates+0)=0; 
  for (i=1; i<=(nstates-1); i++)
    {
      *(dstates+i)=*(dstates+i-1)+val;      
    }

  // Initialize the records.
  stopT=-1;
  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+n*ntypes+i)=0;
      *(dosage+n)=-1;
    }

  for (i=0; i<=(ntypes-1); i++)
    *(y+i)=*(x0+i);

  seq=malloc(sizeof(long)*nsteps); n=0;
  nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  while ((n<(nintervals-1))&&(stopT<0))
    {
      int maxdepth, maxind, minind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1;

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
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+i);
	    }
	  else
	    {
	      for (i=0; i<=(ntypes-1); i++)
		*((seqnodes+0)->y+i)=*(y+(n-1)*ntypes+i);
	    }

	  for (i=0; i<=nsteps; i++)
	    *(nfrontiers+i)=0;

	  curind=0;	  	  	      	      	 
	  seqnodes=recurse_sequence_two_drug_responses3(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
	  	  

	  /*
	  // Debug
	  printf("nseqnodes=%d\n",nseqnodes);
	  for (i=0; i<=(nseqnodes-1); i++)
	    {
	      printf("seqnode[%d], index=%d, pa=%d, depth=%d, dosage=( ",i,(seqnodes+i)->index,(seqnodes+i)->pa,(seqnodes+i)->depth);
	      for (j=0; j<=(ndrugs-1); j++)
		printf("%.3f ",*((seqnodes+i)->dosage+j));
	      printf("), pop=( ");
	      for (j=0; j<=(ntypes-1); j++)
		printf("%.3e ",*((seqnodes+i)->y+j));
	      printf("), totalpop= ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((seqnodes+i)->y+j);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */

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
	      int *inferior, n1, n2, *lookup1, *lookup2, ndegenerates, *degenerates;
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

	      // Debug
	      //printf("ncandidates=%d\n",ncandidates); 
	     	     
	      
	      // Sort candidate sequences by total population and by R12 population.
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
	      // Suppose a candidate's total population > n1 other candidates, and its R12 population > n2 other candidates.
	      // Pr(the candidate is inferior) = 
	      // Pr(there is at least one other candidate with both total and R12 population < the target candidate) =
	      // Pr(there is an overlap of n1 superior candidates in total population and n2 superior candidates in R12 population) =
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

	      /*
	      // Debug
	      k=0;
	      for (i=0; i<=(ncandidates-1); i++)
		if (*(inferior+i)==0)
		  k++;
	      printf("ncandidates=%d\n",k);
	      */
	     	      
	      // Evaluate the expected death time for each non-inferior candidate.
	      deathtimes=malloc(sizeof(double)*ncandidates);
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  *(deathtimes+i)=-1;
		  if (*(inferior+i)==0)
		    {
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,timeinterval,nintervals);
		      //val=deathtime((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval);
		      
		      //val=deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		      
		      val=deathtime3((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);


		      *(deathtimes+i)=val;		      		      
		    }
		}	

	      /*
	      // Debug
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(inferior+i)==0)
		    {
		      printf("ind %d, pop ",i);
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.4e ",*((seqnodes+*(candidates+i))->y+j));
		      printf(", deathtime=%.0f\n",*(deathtimes+i));
		    }
		}
	      //exit(0);
	      */
	      

	      // Choose the candidate with the max death time.
	      minind=-1; minval=-1;
	      for (i=0; i<=(ncandidates-1); i++)
		{
		  if (*(deathtimes+i)>minval)
		    {
		      minind=*(candidates+i); minval=*(deathtimes+i);
		    }
		}

	      // If there are degenerate death times, then choose the one with minimum R12 population.	      
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
	      free(degenerates);
	 
	      
	      // Debug
	      //if (n==0)
	      //minind=*(candidates+1376869);

	      /*
	      // Debug
	      if (n==0)
		{
		  i=1376869;
		  i=316206;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		  i=1376886;
		  i=797819;
		  val=deb_deathtime2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,1,nintervals*timeinterval,nstates);
		}
	      */

	      /*
	      // Debug
	      if (n==0)
		{
		  double *vec, *vec2, *y0;
		  y0=malloc(sizeof(double)*ntypes);
		  vec=malloc(sizeof(double)*ndrugs);
		  vec2=malloc(sizeof(double)*ntypes);
		  *(vec+0)=0.75; *(vec+1)=0.2;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  printf("\n\n");

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2); 
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		  *(vec+0)=0.75; *(vec+1)=0.25;
		  recurse_two_drug_response_trimrates2(x0,T,g0,Sg,a0,Sa,0,timeinterval,vec,y0);
		  printf("t0: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(y0+i));
		  printf("\n");

		  exit(0);

		  //*(vec+0)=0; *(vec+1)=1;
		  //recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		  //printf("t1: ");
		  //for (i=0; i<=(ntypes-1); i++)
		  //printf("%.3e ",*(vec2+i));
		  //printf("\n");
		}
	      */

	      /*
	      // Debug
	      printf("ndegenerates=%d\n",ndegenerates);
	      printf("( ");
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)		
		{
		  val+=*((seqnodes+minind)->y+j);
		  printf("%.4e ",*((seqnodes+minind)->y+j));
		}
	      printf("), total %.4e\n",val);
	      */
	      
	      
	      // Release memory.
	      free(inferior);
	      free(pairs1); free(pairs2); 
	      free(lookup1); free(lookup2); free(deathtimes);
	      
	    }
	     

	  // Debug
	  //printf("minind=%d\n",minind);
	  //exit(0);

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
	  for (i=0; i<=(nsteps-1); i++)
	    {
	      j=*(seq+i);
	      printf("seqnode[%d], index=%ld, pa=%d, depth=%d, dosage=( ",j,(seqnodes+j)->index,(seqnodes+j)->pa,(seqnodes+j)->depth);
	      for (k=0; k<=(ndrugs-1); k++)
		printf("%.3f ",*((seqnodes+j)->dosage+k));
	      printf("), pop=( ");
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((seqnodes+j)->y+k));
	      printf("), totalpop=");
	      val=0;
	      for (k=0; k<=(ntypes-1); k++)
		val+=*((seqnodes+j)->y+k);
	      printf("%.3e\n",val);
	    }
	  //exit(0);
	  */
	}


      // Transcribe the dosages and responses determined by dynamic programming.
      k=(n%nsteps); l=*(seq+k);
      
      if (l<0)
	stopT=*(t+n);

      else
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+n*ndrugs+j)=*((seqnodes+l)->dosage+j);

	  for (i=0; i<=(ntypes-1); i++)
	    *(y+n*ntypes+i)=*((seqnodes+l)->y+i);
      
	  maxval=0;
	  for (i=0; i<=(ntypes-1); i++)
	    maxval+=*(y+n*ntypes+i);
	

	  // Stop when the total population >= mortal.
	  if (maxval>=mortal)
	    stopT=*(t+n);

	  // Stop when the total population < 1.
	  //if (maxval<=0)	
	  if (maxval<1)
	    stopT=*(t+nintervals-1);

	  /*
	  // Debug
	  printf("n=%d, t=%.3f, dosage=( ",n,*(t+n));
	  for (j=0; j<=(ndrugs-1); j++)
	    printf("%.3f ",*(dosage+n*ndrugs+j));
	  printf("), pop=( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(y+n*ntypes+i));
	  printf("), totalpop=%.3e\n",maxval);
	  //exit(0);
	  */

	}

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
  
  return stopT;
}



// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
struct seqnode * recurse_sequence_two_drug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, long *curind_pnt)
{
  int i, j, k, l, m, n;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
  
  // Return when the node depth>=nsteps.
  if ((seqnodes+curind)->depth>=nsteps)
    return seqnodes;

  // Return when the total population exceeds the mortal threshold.
  val=0;
  for (i=0; i<=(ntypes-1); i++)
    val+=*((seqnodes+curind)->y+i);

  // Debug
  //printf("curind=%ld, depth=%d, pop=%.4e\n",curind,(seqnodes+curind)->depth,val);
  //for (i=0; i<=(ntypes-1); i++)
  //printf("%.4e ",*((seqnodes+curind)->y+i));
  //printf("\n");

 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;
  
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
      *(vec+0)=*(dstates+n); *(vec+1)=1-*(vec+0);
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);

     
      recurse_two_drug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
      
      /*
      // Debug
      printf("dosage: ");
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.3f ",*(vec+i));
      printf("\n");
      printf("pop: ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.3e ",*(vec2+i));
      printf("\n");
      exit(0);
      */


      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=recurse_sequence_two_drug_responses2(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,&curind2);
    }

  return seqnodes;
}



// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
struct seqnode * recurse_sequence_two_drug_responses3(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt)
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

  // Debug
  //printf("curind=%ld, depth=%d, pop=%.4e\n",curind,(seqnodes+curind)->depth,val);
  //for (i=0; i<=(ntypes-1); i++)
  //printf("%.4e ",*((seqnodes+curind)->y+i));
  //printf("\n");

 
  if (val>=mortal)
    return seqnodes;

  // Return when the total population is below 0.
  //if (val<=1e-10)
  if (val<1)
    return seqnodes;


  // Return when the total and R12 populations of the current node are inferior to the frontiers at the current level.
  depth=(seqnodes+curind)->depth; inferior=0; n=0; 
  val1=0;
  for (i=0; i<=(ntypes-1); i++)
    val1+=*((seqnodes+curind)->y+i);
  val2=*((seqnodes+curind)->y+ntypes-1); 
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // Check whether both total and R12 populations are larger than a frontier.           
      if ((val1>*(frontiers+depth*maxnfrontiers*2+n*2+0))&&(val2>*(frontiers+depth*maxnfrontiers*2+n*2+1)))
	inferior=1;
      n++;
    }

  // Debug
  //printf("curind=%d, inferior=%d\n",curind,inferior);

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
      *(vec+0)=*(dstates+n); *(vec+1)=1-*(vec+0);
      for (i=0; i<=(ndrugs-1); i++)
	*((seqnodes+nseqnodes+n)->dosage+i)=*(vec+i);

     
      recurse_two_drug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
      
      /*
      // Debug
      printf("dosage: ");
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.3f ",*(vec+i));
      printf("\n");
      printf("pop: ");
      for (i=0; i<=(ntypes-1); i++)
	printf("%.3e ",*(vec2+i));
      printf("\n");
      exit(0);
      */


      for (i=0; i<=(ntypes-1); i++)
	*((seqnodes+nseqnodes+n)->y+i)=*(vec2+i);
    }
  free(vec); free(vec2);

  *nseqnodes_pnt=nseqnodes+nstates;

  for (n=0; n<=(nstates-1); n++)
    {
      curind2=nseqnodes+n;
      seqnodes=recurse_sequence_two_drug_responses3(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind2);
    }

  return seqnodes;
}



