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


extern int ndrugs, ntypes, initialtime;
extern double detected, mortal;
extern int maxnfrontiers;

extern void padm(int ndim, double *A, double t, int p, double *E);
extern char * getitemval3(char *s, int ind, int *length, char sepch);
extern char * getitemval4(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
extern void diagonal(int ndim, double *vec, double *mat);
extern int pair2_cmp(struct pair2 *p1, struct pair2 *p2);
extern int pair_cmp_ascend(struct pair *p1, struct pair *p2);
extern void recurse_multidrug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
extern void recurse_multidrug_response_trimrates2_2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode);
extern void recurse_multidrug_response_trimrates3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses, int evamode);

extern struct seqnode * recurse_sequence_multidrug_responses(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);

extern struct seqnode * recurse_sequence_multidrug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);

extern struct seqnode * deb_recurse_sequence_multidrug_responses2(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers, long *curind_pnt);


//int maxncandseqs=500;
//int maxncandseqs=50;

extern int maxncandseqs;
extern int curevamode;

int frontier2flag=1;
double valthre=0.01, eqvalthre=1e-10;
//double vallowthre=0.01, valuppthre=0.05;
//double vallowthre=0.01, valuppthre=1.0;
double vallowthre=0.01, valuppthre=100.0;


double global_dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int *truncated);

struct seqnode * recurse_sequence_multidrug_responses3(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers1, double *frontiers2, double *frontiers3, double maxtime, long *curind_pnt);

struct seqnode * recurse_sequence_multidrug_responses4(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers1, double *frontiers2, double *frontiers3, double maxtime, long *curind_pnt);


double global_dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);

struct seqnode * deb_recurse_sequence_multidrug_responses3(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers1, double *frontiers2, double *frontiers3, double maxtime, long *curind_pnt);

int vectorbelow(int ndim, double *vec1, double *vec2);
int vectorbeloweq(int ndim, double *vec1, double *vec2);
int vectoreq(int ndim, double *vec1, double *vec2);
int deb_vectorbelow(int ndim, double *vec1, double *vec2);
int deb_vectorbeloweq(int ndim, double *vec1, double *vec2);

double global_dp_optimize_multidrug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);

double global_dp_optimize_multidrug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);

double global_dp_optimize_multidrug_monotherapy_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated);


// Find a global optimal treatment sequence using dynamic programming.
// Traverse the decision tree of all treatment sequences and pick the best one.
// Apply branch and bound to trim inferior branches.
// Branch and bound cannot retract previously established branches.
// Stop every nsteps and keep only superior sequences.  Continue traversing.
// The ultimate criterion is the survival time.

double global_dp_optimize_multidrug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosage, int *truncated)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, combind, *resvec, *bases;
  long curind, ncandseqs, nnewcandseqs, minind;
  double *dstates, val, stopT, *frontiers, pop, incurablepop, finaltreatmenttime=0, minval;
  struct seqnode *seqnodes; 
  struct sequence *candseqs, *newcandseqs;
  time_t *tp;

  *truncated=0;
  tp=malloc(sizeof(time_t)); srand((int)(time(tp)));
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


  n=0; nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*2);
  ncandseqs=1; candseqs=malloc(sizeof(struct sequence));
  (candseqs+ncandseqs-1)->depth=0;
  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs);
  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    *((candseqs+ncandseqs-1)->ys+i)=*(y+i*nintervals+0);
  nnewcandseqs=0; newcandseqs=malloc(sizeof(struct sequence));

  // Consider all treatment sequences throughout nintervals.
  while ((n<(nintervals-1))&&(stopT<0))
    {      
      int maxdepth, maxind, minind, candind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1, extinction=0, *valid;

      // (n%nsteps)=0: incur dynamic programming.
      // For each candidate sequence, run the recursive program to generate all descendant nodes of depth nsteps.
      // Discard the nodes by branch-and-bound.
      if ((n%nsteps)==0)
	{	  
	  nnewcandseqs=0;
	  valid=malloc(sizeof(int));

	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);	  	  
	  
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      printf("candseq %d, ( ",i);
	      for (j=0; j<=(ntypes-1); j++)
		{
		  printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j));
		  val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j);
		}
	      printf("), %.3e\n",val);	      
	    }
	  */


	  // For each candidate sequence, run recurse_sequence_multidrug_responses to exhaust all valid following sequences up to nsteps.
	  // If the candidate sequence already leads to mortality, then do not proceed.
	  for (candind=0; candind<=(ncandseqs-1); candind++)
	    {
	      // If the prior treatment sequences lead to zero population, then don't consider new ones.
	      if (extinction==0)
		{
		  // Do not proceed recursion if the current candidate already leads to mortality.
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);

		  // Debug
		  //printf("candind=%d, val=%.3e\n",candind,val);


		  if (val>=mortal)
		    {		      
		      nnewcandseqs++;
		      valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
		      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);		      
		      (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth;
		      (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
		      (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
		      for (j=0; j<=((candseqs+candind)->depth-1); j++)
			{
			  for (k=0; k<=(ndrugs-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			}
		      for (j=0; j<=(candseqs+candind)->depth; j++)
			{
			  for (k=0; k<=(ntypes-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			}		      
		    }
		  
		  else
		    {
		      /*
		      // Clear up and initialize memory for seqnodes.
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
		      */

		      // Initialize memory for seqnodes.
		      nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
		      (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
		      (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
		      (seqnodes+0)->children=malloc(sizeof(long));
		      (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
		      (seqnodes+0)->y=malloc(sizeof(double)*ntypes);
		      for (i=0; i<=(ntypes-1); i++)
			*((seqnodes+0)->y+i)=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);		    
		      for (i=0; i<=nsteps; i++)
			*(nfrontiers+i)=0;		  
		      curind=0;	  	  	
		      
		      /*
		      // Debug
		      if ((n==5)&&(candind==0))
			{
			  seqnodes=deb_recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
			  //for (i=0; i<=(10-1); i++)
			  //printf("(%d %d %d)\n",i,(seqnodes+i)->nchildren,(seqnodes+i)->depth); 
			  maxdepth=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    if ((seqnodes+i)->nchildren==0)
			      maxdepth=max(maxdepth,(seqnodes+i)->depth);
			  printf("maxdepth=%d\n",maxdepth);
			  exit(0);
			}
		      */

		    
		      //seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
		      seqnodes=recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);
		      	      
		      /*
		      // Debug
		      printf("n=%d, candind=%d, nseqnodes=%d\n",n,candind,nseqnodes);
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  printf("i=%d, pa=%d, nchildren=%d, depth=%d, ( ",i,(seqnodes+i)->pa,(seqnodes+i)->nchildren,(seqnodes+i)->depth);
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      printf("%.3e ",*((seqnodes+i)->y+j));
			      val+=*((seqnodes+i)->y+j);
			    }
			  printf(") %.3e\n",val);
			}
		      //exit(0);
		      */
		      
		      maxdepth=0;
		      for (i=0; i<=(nseqnodes-1); i++)
			if ((seqnodes+i)->nchildren==0)
			  maxdepth=max(maxdepth,(seqnodes+i)->depth);


		      // Find the candidate terminal seqnoces that yield almost zero total population.
		      // Each subpopulation < 1, and this trend lasts for maxdepth-1 steps.
		      ncandidates=0; candidates=malloc(sizeof(long));
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			    {	
			      int flag=1, curind=i;

			      // If maxdepth=1, then check the population composition at the current level.
			      // If maxdepth>1, then check the population composition up to (maxdepth-1) levels.
			      if (maxdepth==1)
				m=maxdepth;
			      else
				m=maxdepth-1;				  
			      l=0;
			      //while ((l<(maxdepth-1))&&(flag==1)&&(curind>=0))
			      while ((l<m)&&(flag==1)&&(curind>=0))
				{
				  j=0;
				  while ((j<=(ntypes-1))&&(flag==1))
				    {				      
				      if (*((seqnodes+curind)->y+j)>=1)
				      //if (*((seqnodes+curind)->y+j)>1)
					flag=0;
				      j++;
				    }
				  //if (*((seqnodes+curind)->y+ntypes-1)>=0.1)
				  //flag=0;
				  curind=(seqnodes+curind)->pa; l++;
				}

			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}

			      /*
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)		
				val+=*((seqnodes+i)->y+j);
			      if (val<1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;
				}
			      */
			    }
			}
		      

		      // Debug
		      //printf("ncandidates=%d\n",ncandidates);
		      

		      // If the candidate sequences lead to zero population, then pick up any instance sequence and stop iteration.  Nullify the sequences derived from other candidates.
		      if (ncandidates>0)
			{
			  extinction=1;
			  for (i=0; i<=(nnewcandseqs-1); i++)
			    *(valid+i)=0;			  
			}
		      
		      // If no seqnodes yield zero total population, then keep the ones not inferior to others.
		      else
			{
			  int nnewfrontiers, nleft, *visited;
			  double *newfrontiers, *template;
			  
			  // Establish frontiers of terminal nodes.			  
			  nnewfrontiers=0; newfrontiers=malloc(sizeof(double)*2);			  
			  template=malloc(sizeof(double)*2);
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1=0, val2=0;
				  for (j=0; j<=(ntypes-1); j++)
				    val1+=*((seqnodes+i)->y+j);
				  val2=*((seqnodes+i)->y+ntypes-1);				  
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      if ((val1>=*(newfrontiers+j*2+0))&&(val2>=*(newfrontiers+j*2+1)))
					k=j;
				      j++;
				    }
				  if (k<0)
				    {				      
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*2);
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  if ((*(newfrontiers+j*2+0)<=val1)||(*(newfrontiers+j*2+1)<=val2))
					    {
					      nleft++; *(template+(nleft-1)*2+0)=*(newfrontiers+j*2+0);
					      *(template+(nleft-1)*2+1)=*(newfrontiers+j*2+1);
					    }
					}
				      nleft++;
				      *(template+(nleft-1)*2+0)=val1;
				      *(template+(nleft-1)*2+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers=realloc(newfrontiers,sizeof(double)*nnewfrontiers*2);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  *(newfrontiers+j*2+0)=*(template+j*2+0);
					  *(newfrontiers+j*2+1)=*(template+j*2+1);
					}
				    }
				}
			    }				      
			  
			  free(template);
	  			  			  
			  /*
			  // Establish frontiers of terminal nodes.
			  nnewfrontiers=0; newfrontiers=malloc(sizeof(double)*2);
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
				{
				  double val1=0, val2=0;
				  for (j=0; j<=(ntypes-1); j++)
				    val1+=*((seqnodes+i)->y+j);
				  val2=*((seqnodes+i)->y+ntypes-1);
				  				  
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      if ((val1>=*(newfrontiers+j*2+0))&&(val2>=*(newfrontiers+j*2+1)))
					k=j;
				      j++;
				    }
				  if (k<0)
				    {
				      nnewfrontiers++;
				      newfrontiers=realloc(newfrontiers,sizeof(double)*nnewfrontiers*2);
				      *(newfrontiers+(nnewfrontiers-1)*2+0)=val1;
				      *(newfrontiers+(nnewfrontiers-1)*2+1)=val2;
				    }

				  // Debug
				  if (k<0)
				    printf("*");
				  printf("%d %.3e %.3e\n",i,val1,val2);

				}
			    }
			  */			

			  
			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1=0, val2=0;
				  for (j=0; j<=(ntypes-1); j++)
				    val1+=*((seqnodes+i)->y+j);
				  val2=*((seqnodes+i)->y+ntypes-1);
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      if ((val1>*(newfrontiers+j*2+0))&&(val2>*(newfrontiers+j*2+1)))
					k=j;
				      j++;
				    }
				 				  
				  
				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  if ((val1==*(newfrontiers+j*2+0))&&(val2==*(newfrontiers+j*2+1)))
					    k=j;
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  *(visited+k)=1; k=-1;
					}
				    }
				  					    

				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;
				    }
				}
			    }
			  
			  free(visited);

			  /*
			  // Find the candidates on the new frontiers.
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1=0, val2=0;
				  for (j=0; j<=(ntypes-1); j++)
				    val1+=*((seqnodes+i)->y+j);
				  val2=*((seqnodes+i)->y+ntypes-1);
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      if ((val1>*(newfrontiers+j*2+0))&&(val2>*(newfrontiers+j*2+1)))
					k=j;
				      j++;
				    }

				  //
				  //while ((j<=*(nfrontiers+nsteps))&&(k<0))
				  //{
				  //  if ((val1>*(frontiers+nsteps*maxnfrontiers*2+j*2+0))&&(val2>*(frontiers+nsteps*maxnfrontiers*2+j*2+1)))
				  //  k=j;
				  //  j++;
				  //}
				  //

				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;
				    }
				}
			    }
			  */

			  free(newfrontiers);
			}			    
		    	      		    
		      /*
		      // Debug
		      printf("n=%d, candind=%d, ncandidates=%d, nnewcandseqs=%d\n",n,candind,ncandidates,nnewcandseqs);
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  printf("%d ",*(candidates+i));
			}
		      printf("\n");
		      exit(0);
		      //m=nnewcandseqs;
		      */

		      // Append the candidate sequences to the existing ones.
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  nnewcandseqs++;
			  valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
			  newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);
			  //(newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+nsteps;
			  (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+maxdepth;
			  (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
			  (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
			  for (j=0; j<=((candseqs+candind)->depth-1); j++)
			    {
			      for (k=0; k<=(ndrugs-1); k++)
				*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			    }
			  for (j=0; j<=(candseqs+candind)->depth; j++)
			    {
			      for (k=0; k<=(ntypes-1); k++)
				*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			    }
			  curind=*(candidates+i); l=(newcandseqs+nnewcandseqs-1)->depth;
			  while (curind>0)
			    {
			      for (j=0; j<=(ndrugs-1); j++)			    
				*((newcandseqs+nnewcandseqs-1)->dosages+(l-1)*ndrugs+j)=*((seqnodes+curind)->dosage+j);
			      for (j=0; j<=(ntypes-1); j++)
				*((newcandseqs+nnewcandseqs-1)->ys+l*ntypes+j)=*((seqnodes+curind)->y+j);
			      curind=(seqnodes+curind)->pa; l--;
			    }			 
			}

		      // Debug
		      //printf("after appending\n");

		      /*
		      // Debug
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  printf("candidate %d, ( ",i); val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      printf("%.3e ",*((newcandseqs+m+i)->ys+((newcandseqs+m+i)->depth-1)*ntypes+j));
			      val+=*((newcandseqs+m+i)->ys+((newcandseqs+m+i)->depth-1)*ntypes+j);
			    }
			  printf("), %.3e\n",val);
			}
		      */      
		
		      free(candidates);

		      // Clear up memory for seqnodes.
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  free((seqnodes+i)->children);
			  free((seqnodes+i)->dosage);
			  free((seqnodes+i)->y);
			}
		      free(seqnodes);	     

		    }		 
		}
	    }

	  /*
	  // Debug
	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  printf("n=%d, nnewcandseqs=%d, maxdepth=%d, extinction=%d\n",n,nnewcandseqs,maxdepth,extinction);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      //if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  printf("newcand %d, valid=%d, depth=%d, ( ",i,*(valid+i),(newcandseqs+i)->depth);
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      printf("%.3e ",*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j));
		      val+=*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j);
		    }
		  printf(") %.3e\n",val);
		}
	    }
	  */


	  // Distinct treatment sequences may lead to identical or inferior terminal population composition.  Establish frontiers again on all treatment sequences and filter out inferior/identical ones.
	  if (nnewcandseqs>0)
	    {
	      int nnewfrontiers, nleft, *visited;
	      double *newfrontiers, *template;
	
	      maxdepth=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		maxdepth=max(maxdepth,(newcandseqs+i)->depth);

	      // Establish frontiers of terminal nodes.			  
	      nnewfrontiers=0; newfrontiers=malloc(sizeof(double)*2);			  
	      template=malloc(sizeof(double)*2);
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1=0, val2=0;
		      for (j=0; j<=(ntypes-1); j++)
			val1+=*((newcandseqs+i)->ys+maxdepth*ntypes+j);
		      val2=*((newcandseqs+i)->ys+maxdepth*ntypes+ntypes-1);
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  if ((val1>=*(newfrontiers+j*2+0))&&(val2>=*(newfrontiers+j*2+1)))
			    k=j;
			  j++;
			}
		      if (k<0)
			{				      
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*2);
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      if ((*(newfrontiers+j*2+0)<=val1)||(*(newfrontiers+j*2+1)<=val2))
				{
				  nleft++; *(template+(nleft-1)*2+0)=*(newfrontiers+j*2+0);
				  *(template+(nleft-1)*2+1)=*(newfrontiers+j*2+1);
				}
			    }
			  nleft++;
			  *(template+(nleft-1)*2+0)=val1;
			  *(template+(nleft-1)*2+1)=val2;
			  nnewfrontiers=nleft;
			  newfrontiers=realloc(newfrontiers,sizeof(double)*nnewfrontiers*2);
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      *(newfrontiers+j*2+0)=*(template+j*2+0);
			      *(newfrontiers+j*2+1)=*(template+j*2+1);
			    }
			}
		    }
		}				      
	      
	      free(template);

	      /*
	      // Debug
	      printf("nnewfrontiers=%d\n",nnewfrontiers);
	      
	      for (i=0; i<=(nnewfrontiers-1); i++)
		printf("%d %f %f\n",i,*(newfrontiers+i*2+0),*(newfrontiers+i*2+1));
	      */

		  
	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1=0, val2=0;
		      for (j=0; j<=(ntypes-1); j++)
			val1+=*((newcandseqs+i)->ys+maxdepth*ntypes+j);
		      val2=*((newcandseqs+i)->ys+maxdepth*ntypes+ntypes-1);		      
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  if ((val1>*(newfrontiers+j*2+0))&&(val2>*(newfrontiers+j*2+1)))
			  //if ((val1>=*(newfrontiers+j*2+0))&&(val2>=*(newfrontiers+j*2+1))&&(~((val1==*(newfrontiers+j*2+0))&&(val2==*(newfrontiers+j*2+1)))))
			    k=j;
			  else if ((val1>*(newfrontiers+j*2+0))&&(val2==*(newfrontiers+j*2+1)))
			    k=j;
			  else if ((val1==*(newfrontiers+j*2+0))&&(val2>*(newfrontiers+j*2+1)))
			    k=j;
			  j++;
			}
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      if ((val1==*(newfrontiers+j*2+0))&&(val2==*(newfrontiers+j*2+1)))
				k=j;
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;
			    }
			}

		      if (k>=0)
			*(valid+i)=0;
		    }
		}
	      free(visited); free(newfrontiers);
	    }	      
	  

	  	  
	  /*
	  // Discard the candidates that yield identical or inferior terminal population composition.
	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  for (j=(i+1); j<=(nnewcandseqs-1); j++)
		    {
		      if ((*(valid+j)==1)&&((newcandseqs+j)->depth==maxdepth))
			{
			  int flag=1; 			 

			  k=0;
			  while ((k<=(ntypes-1))&&(flag==1))
			    {
			      if (fabs(*((newcandseqs+i)->ys+(maxdepth)*ntypes+k)-*((newcandseqs+j)->ys+(maxdepth)*ntypes+k))>0.01)
				flag=0;
			      k++;
			    }
			  if (flag==1)
			    *(valid+j)=0;
			}
		    }
		}
	    }
	  	  
	  */
	
	  /*
	  // Debug
	  printf("nnewcandseqs=%d\n",nnewcandseqs);	  	  	  
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      printf("newcandseq %d, depth %d\n",i,(newcandseqs+i)->depth);
	      printf("dosages: ");
	      for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((newcandseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(newcandseqs+i)->depth; j++)
		{
		  printf("(%d ",j);
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((newcandseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  //exit(0);
	  //if (n==10)
	  //exit(0);
	  */


	  // Update the candidate sequences.
	  // Remove the new candidate sequeces that are inferior to others.
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      free((candseqs+i)->dosages); free((candseqs+i)->ys);
	    }	  
	  ncandseqs=0;

	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  ncandseqs++; candseqs=realloc(candseqs,sizeof(struct sequence)*ncandseqs);
		  (candseqs+ncandseqs-1)->depth=(newcandseqs+i)->depth;
		  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(candseqs+ncandseqs-1)->depth);
		  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*((candseqs+ncandseqs-1)->depth+1));
		  for (j=0; j<=((candseqs+ncandseqs-1)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((candseqs+ncandseqs-1)->dosages+j*ndrugs+k)=*((newcandseqs+i)->dosages+j*ndrugs+k);
		  for (j=0; j<=(candseqs+ncandseqs-1)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((candseqs+ncandseqs-1)->ys+j*ntypes+k)=*((newcandseqs+i)->ys+j*ntypes+k);		  
		}
	    }
	  
	 

	  // Clear up memory for new candidates.
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      free((newcandseqs+i)->dosages);
	      free((newcandseqs+i)->ys);
	    }
	  free(valid); 

	  // Debug
	  //printf("ncandseqs=%d\n",ncandseqs);
	  
	  /*
	  // Debug
	  if (n==4)
	    {
	      printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  val=0;
		  printf("candseq %d, depth %d, y=( ",i,(candseqs+i)->depth);
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k));
		      val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k);
		    }
		  printf("), %.3e\n",val);
		}
	    }
	  */

	  
	  // If ncandseqs exceeds the maximum number, then only keep maxncandseqs of them and report truncation.
	  if (ncandseqs>=maxncandseqs)
	    {
	      struct pair *pairs;	      
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  k=rand(); (pairs+i)->i=i; (pairs+i)->s=(double)k;
		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs;
	      free(pairs); *truncated=1;	    
	    }
	  

	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      printf("candseq %d, depth %d\n",i,(candseqs+i)->depth);
	      printf("dosages: ");
	      for (j=0; j<=((candseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((candseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(candseqs+i)->depth; j++)
		{
		  printf("( ");
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  */

	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      printf("candseq %d, depth %d, y=( ",i,(candseqs+i)->depth);
	      for (k=0; k<=(ntypes-1); k++)
		printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k));
	      printf(")\n");
	    }
	  */

	  // Debug
	  //printf("n=%d, stopT=%f, extinction=%d\n",n,stopT,extinction);


	  // Stop when some sequences yield zero population.
	  // Report the survival time to be max time span + timeinterval.
	  // This number distinguishes the cases where the patients are not cured but can survive for the max time span.
	  if (extinction==1)
	    {
	      //stopT=*(t+nintervals-1);
	      stopT=*(t+nintervals-1)+timeinterval;
	      finaltreatmenttime=*(t+maxdepth);
	    }

	  // Stop when the best candidate treatment sequences lead to mortality.
	  minval=1e5000;
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((candseqs+i)->ys+ntypes*(candseqs+i)->depth+j);
	      minval=min(minval,val);
	    }
	
	  //if (val>=mortal)
	  if (minval>=mortal)
	    {
	      if (maxdepth>=1)
		stopT=*(t+maxdepth-1);
	      else
		stopT=0;
	      finaltreatmenttime=stopT;
	    }	  
	  
	}

      n++;
    }
  
  if (stopT<0)
    stopT=*(t+nintervals-1);

  
  // Debug
  //printf("stopT=%f\n",stopT);
  


  // Extract the best sequence from candidates.
  // In principle, all the final candidates are equivalent.
  // Pick up the one with minimum total population.
  n=0; minind=-1; minval=1e5000;
  for (n=0; n<=(ncandseqs-1); n++)
    {
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*((candseqs+n)->ys+ntypes*(candseqs+n)->depth+i);
      if (val<minval)
	{
	  minind=n; minval=val;
	}
    }
  
  // Write the selected treatment sequence to dosage and population to y.
  for (n=0; n<=((candseqs+minind)->depth-1); n++)    
    for (i=0; i<=(ndrugs-1); i++)	
      *(dosage+i*nintervals+n)=*((candseqs+minind)->dosages+n*ndrugs+i);
  for (n=0; n<=(candseqs+minind)->depth; n++)
    for (i=0; i<=(ntypes-1); i++)
      *(y+i*nintervals+n)=*((candseqs+minind)->ys+n*ntypes+i);
  if ((minval>=mortal)||(minval<1))
    {
      for (n=(candseqs+minind)->depth; n<=(nintervals-1); n++)
	for (i=0; i<=(ndrugs-1); i++)
	  *(dosage+i*nintervals+n)=0;
      for (n=((candseqs+minind)->depth+1); n<=(nintervals-1); n++)
	for (i=0; i<=(ntypes-1); i++)
	  *(y+i*nintervals+n)=*((candseqs+minind)->ys+(candseqs+minind)->depth*ntypes+i);
    }

  /*
  // Debug  
  for (n=0; n<=((candseqs+minind)->depth-1); n++)
    {
      printf("gb t(%d)=%.2f, dosage=( ",n,*(t+n));
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*(dosage+i*nintervals+n));
      printf("), y=( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)	
	{
	  printf("%.4e ",*(y+i*nintervals+n));
	  val+=*(y+i*nintervals+n);
	}
      printf(") %.4e\n",val);
    }
  */

  
  // Release memory.
  free(dstates);
  /*
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); 
  */
  free(nfrontiers); free(frontiers);
  free(resvec); free(bases);
  for (n=0; n<=(ncandseqs-1); n++)
    {
      free((candseqs+n)->dosages);
      free((candseqs+n)->ys);
    }
  free(candseqs); free(tp);   
  return stopT;
}
		
		      
// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
// Difference from recurse_sequence_two_drug_responses3: multi-drug version.
// Difference from recurse_sequence_multidrug_responses: remove the frontiers which are no longer frontiers.
// Difference from recurse_sequence_multidrug_responses2: change the criteria for frontiers.
// Old criteria are total population and multiple-resistant population.
// New criteria are population composition and time to mortality.
// Time to mortality cannot be directly assessed for adaptive treatment sequences.  Thus consider the hypothetically worst and best cases.
// Worst case: stick to one drug throughout the whole cycle.
// Best case: administer max dosage for all drugs throughout the whole cycle. 
// Instead of evaluating time to death, calculate total population at maxt under both cases.
// Since population dynamics here are monotic functions, the one with smaller final population should have a longer survival time.
// If the worst case of parameter A has smaller final population than the best case of parameter B, then A is superior to B.
// frontiers1: population composition.
// frontiers2: worst case final total population.
// frontiers3: best case final total population.
// Also return val1 and val2 of each seqnode, so that there is no need to re-evaluate them.
// Store them in ->y fields after the values of all subpopulations.
struct seqnode * recurse_sequence_multidrug_responses3(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers1, double *frontiers2, double *frontiers3, double maxtime, long *curind_pnt)
{
  int i, j, k, l, m, n, inferior, depth, flag;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2, val1, val2;
  double *dosage, *finaly;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
        
  
  // Evaluate population composition at maxtime.
  // Worst case: the highest outcomes of single drugs: val1.
  // Best case: max dosage for each drug: val2.
  // Here the worst case is actually a lower bound for the final population achievable from the current population composition.
  // Thus choose the best time-invariant therapy (monotherapy or combination) and evaluate its final population.
  dosage=malloc(sizeof(double)*ndrugs); val1=0;
  finaly=malloc(sizeof(double)*ntypes);

  
  if (frontier2flag==1)
    {
      val1=1e5000;
      for (n=0; n<=(nstates-1); n++)
	{
	  for (i=0; i<=(ndrugs-1); i++)
	    *(dosage+i)=*(dstates+n*ndrugs+i);
	  //recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
	  //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
	  //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
	  recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(finaly+i);
	  val1=min(val1,val);      
	}
    }
  

  /*
  val1=0;
  for (n=0; n<=(ndrugs-1); n++)
    {
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i)=0;
      *(dosage+n)=1;
      recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(finaly+i);
      val1=max(val1,val);      
    }
  */

  if (frontier2flag==1)
    {      
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i)=1;
      //recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
      //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);      
      //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
      recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(finaly+i);
      val2=val; val2=min(val2,val1);                 
    }

  free(dosage); free(finaly);

  // Debug
  //if ((curind>=85)&&(curind<=87))
  //if ((seqnodes+curind)->pa==2)
  //printf("curind=%d, val1=%.4e, val2=%.4e\n",curind,val1,val2);


  // Debug
  //if ((curind==87)||(curind==94))
  //printf("curind=%d, val1=%.4e, val2=%.4e\n",curind,val1,val2);


  if (frontier2flag==1)
    {
      *((seqnodes+curind)->y+ntypes)=val1;
      *((seqnodes+curind)->y+ntypes+1)=val2;
    }
      

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

  /*
  // Debug
  if ((curind==87)||(curind==94))
    {
      depth=(seqnodes+curind)->depth;
      printf("curind=%d, depth=%d, nfrontiers=%d, ",curind,depth,*(nfrontiers+depth));
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  printf("( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i));
	  printf("%.4e %.4e ) ",*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n));
	}
      printf("\n");
    }
  */
 
  // Check whether the population composition of the current node is inferior to frontiers1.
  // A is inferior to B if all the components of A >= the corresponding components of B, and A B are not identical.
  // If A is not inferior in population composition, then check whether the final total population in the best case of A is larger than that of the worst case of B.  If so then A is inferior to B.
  depth=(seqnodes+curind)->depth; inferior=0; n=0;
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // The current node is inferior if the frontier is below the current node.
      inferior=vectorbelow(ntypes,frontiers1+depth*maxnfrontiers*ntypes+n*ntypes,(seqnodes+curind)->y);

      /*
      // Debug
      //if (curind==94)
      if (curind==34)
	{
	  printf("n=%d, inferior1=%d, ",n,inferior);
	  printf("cury ( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*((seqnodes+curind)->y+i));
	  printf("), fronty ( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.4e ",*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i));
	  printf(")\n");
	  //if (n==2)
	  //k=deb_vectorbelow(ntypes,frontiers1+depth*maxnfrontiers*ntypes+n*ntypes,(seqnodes+curind)->y);
	}
      */


      /*
      // Debug
      if (curind==86)
	{
	  printf("n=%d, inferior1=%d, ft= ",n,inferior);
	  for (j=0; j<=(ntypes-1); j++)
	    printf("%.4e ",*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+j));
	  printf("\n");
	}
      */

      /*
      // Check whether each subpopulation of the current node >= the frontier.
      inferior=1; i=0;
      while ((i<=(ntypes-1))&&(inferior==1))
	{
	  if (*((seqnodes+curind)->y+i)<*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i))
	    inferior=0;
	  i++;
	}  

      // Rule out the cases when the two population compositions are identical.
      if (inferior==1)
	{
	  i=0; inferior=0; 
	  while (i<=(ntypes-1)&&(inferior==0))
	    {
	      if (fabs(*((seqnodes+curind)->y+i)-*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i))>0.01)
		inferior=1;
	      i++;
	    }
	}
      */

      // Consider the cases when the final population of the best case in the current node is worse than the frontier.
      if (frontier2flag==1)
	{
	  if ((inferior==0)&&(val2>*(frontiers2+depth*maxnfrontiers+n)))
	    inferior=1;
	}     

      // Debug
      //if (curind==94)
      //if (curind==34)
      //printf("n=%d, inferior2=%d, val1=%f, val2=%f, ft2=%f, ft3=%f\n",n,inferior,val1,val2,*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n));



      // Debug
      //if (curind==86)
      //printf("n=%d, inferior2=%d, val1=%.4e, val2=%.4e, ft=(%.4e %.4e)\n",n,inferior,val1,val2,*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n));


      n++;
    }  

  // Debug
  //if (curind==86)
  //if (curind==34)
  //printf("inferior=%d\n",inferior);

  

  // Return when the current node is inferior to the frontier.
  if (inferior==1)
    return seqnodes;
  
  // Augment/update the frontier points if the current node is not inferior to any of them.
  // Also remove the frontier points inferior to the current node.
  else
    {
      int nleft;
      double *left1, *left2, *left3;      
      nleft=0; left1=malloc(sizeof(double)*maxnfrontiers*ntypes);
      left2=malloc(sizeof(double)*maxnfrontiers);
      left3=malloc(sizeof(double)*maxnfrontiers);

      /*
      // Debug
      printf("before curind=%d, depth=%d, nfrontiers=%d\n",curind,depth,*(nfrontiers+depth));
      for (n=0; n<=(*(nfrontiers+depth)-1); n++)
	{
	  printf("( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i));
	  printf(") %.3e %.3e ",*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n));
	}
      printf("\n");
      */

      
      // Remove the frontiers that are inferior or equal to the current node.
      for (n=0; n<*(nfrontiers+depth); n++)
	{	  
	  flag=vectorbeloweq(ntypes,(seqnodes+curind)->y,frontiers1+depth*maxnfrontiers*ntypes+n*ntypes);
	  if (frontier2flag==1)
	    {
	      //if (val2<=*(frontiers2+depth*maxnfrontiers+n))
	      if ((flag==0)&&(*(frontiers3+depth*maxnfrontiers+n)>val1))
		flag=1;
	    }

	  if (flag==0)
	    {
	      nleft++; 
	      for (i=0; i<=(ntypes-1); i++)		
		*(left1+(nleft-1)*ntypes+i)=*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i);
	      *(left2+(nleft-1))=val1; *(left3+(nleft-1))=val2;	    	      	      
	    }
	}
      

      /*
      for (n=0; n<*(nfrontiers+depth); n++)
	{
	  i=0; flag=1;
	  while ((i<=(ntypes-1))&&(flag==1))
	    {
	      if (*((seqnodes+curind)->y+i)>=*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i))
		flag=0;
	      i++;
	    }

	  
	  i=0;
	  if ((i<=(ntypes-1))&&(flag==0))
	    {
	      if (fabs(*((seqnodes+curind)->y+i)-*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i))>0.01)
		flag=1;
	      i++;
	    }

	  //if (val2<=*(frontiers2+depth*maxnfrontiers+n))
	  //flag=0;

	  if (flag==0)
	    {
	      nleft++; 
	      for (i=0; i<=(ntypes-1); i++)		
		*(left1+(nleft-1)*ntypes+i)=*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i);
	      *(left2+(nleft-1))=val1; *(left3+(nleft-1))=val2;	    	      	      
	    }
	}
      */

      /*
      // Debug
      printf("middle curind=%d, depth=%d, nleft=%d\n",curind,depth,nleft);
      for (n=0; n<=(nleft-1); n++)
	{
	  printf("( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(left1+n*ntypes+i));
	  printf(") %.3e %.3e ",*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n));	  
	}
      printf("\n");
      */

      for (n=0; n<=(nleft-1); n++)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i)=*(left1+n*ntypes+i);
	  *(frontiers2+depth*maxnfrontiers+n)=*(left2+n);
	  *(frontiers3+depth*maxnfrontiers+n)=*(left3+n);	  
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  n=*(nfrontiers+depth)-1;
	  for (i=0; i<=(ntypes-1); i++)
	    *(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i)=*((seqnodes+curind)->y+i);
	  *(frontiers2+depth*maxnfrontiers+n)=val1;
	  *(frontiers3+depth*maxnfrontiers+n)=val2;	  
	}
      free(left1); free(left2); free(left3);
    
      /*
      // Debug
      printf("after curind=%d, depth=%d, nfrontiers=%d\n",curind,depth,*(nfrontiers+depth));
      for (n=0; n<=(*(nfrontiers+depth)-1); n++)
	{
	  printf("( ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i));
	  printf(") %.3e %.3e ",*(frontiers2+depth*maxnfrontiers+n),*(frontiers3+depth*maxnfrontiers+n)); 	  
	}
      printf("\n");
      */

    }
  
  
  // Spin off each dosage combination.
  seqnodes=realloc(seqnodes,sizeof(struct seqnode)*(nseqnodes+nstates));
  vec=malloc(sizeof(double)*ndrugs); vec2=malloc(sizeof(double)*ntypes);

  // Debug
  //printf("curind=%d, depth=%d, nchildren=%d\n",curind,(seqnodes+curind)->depth,nstates);

  
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
      //(seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*ntypes);
      
      // Also store val1 and val2 under ->y.
      (seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*(ntypes+2));
      
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
      seqnodes=recurse_sequence_multidrug_responses3(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind2);
    }

  return seqnodes;
}




// Find a global optimal treatment sequence using dynamic programming.
// Traverse the decision tree of all treatment sequences and pick the best one.
// Apply branch and bound to trim inferior branches.
// Branch and bound cannot retract previously established branches.
// Stop every nsteps and keep only superior sequences.  Continue traversing.
// The ultimate criterion is the survival time.
// Difference from global_dp_optimize_multidrug_responses1: apply more rigorous criteria of population composition and worst-case final population to set frontiers.
// One criterion is the predicted final population.  Since evaluating population dynamics is time-consuming, make sure that each population composition only evaluates val1 and val2 once.
// Set candseqs to have the following forms.
// At each time step, dosage specifies the treatment dosage through the interval, population specifies the initial population composition.
// Thus the number of population compositions = the number of dosages + 1.

double global_dp_optimize_multidrug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, combind, *resvec, *bases;
  long curind, ncandseqs, nnewcandseqs, minind;
  double *dstates, val, stopT, *frontiers1, *frontiers2, *frontiers3, pop, incurablepop, finaltreatmenttime=0, minval, maxtime;
  double *gvals1, *gvals2, *cvals1, *cvals2, *dosage, *finaly;  
  struct seqnode *seqnodes; 
  struct sequence *candseqs, *newcandseqs;
  time_t *tp, *t1, *t2;  
        

  *truncated=0; maxtime=nintervals*timeinterval;
  tp=malloc(sizeof(time_t)); srand((int)(time(tp)));
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

  /*
  // Debug
  nstates=2; dstates=malloc(sizeof(double)*nstates*ndrugs);
  *(dstates+0)=0; *(dstates+1)=1;
  *(dstates+2)=1; *(dstates+3)=0;
  */

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosages+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }


  n=0; nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers1=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*ntypes);
  frontiers2=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  frontiers3=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  ncandseqs=1; candseqs=malloc(sizeof(struct sequence));
  (candseqs+ncandseqs-1)->depth=0;
  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs);
  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*2);
  for (i=0; i<=(ntypes-1); i++)
    *((candseqs+ncandseqs-1)->ys+i)=*(y+i*nintervals+0);
  nnewcandseqs=0; newcandseqs=malloc(sizeof(struct sequence));
  gvals1=malloc(sizeof(double)); gvals2=malloc(sizeof(double));
  cvals1=malloc(sizeof(double)); cvals2=malloc(sizeof(double));
  dosage=malloc(sizeof(double)*ndrugs);
  finaly=malloc(sizeof(double)*ntypes);


  // Consider all treatment sequences throughout nintervals.
  while ((n<(nintervals-1))&&(stopT<0))  
    {      
      int maxdepth, maxind, minind, candind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1, extinction=0, *valid;

      // (n%nsteps)=0: incur dynamic programming.
      // For each candidate sequence, run the recursive program to generate all descendant nodes of depth nsteps.
      // Discard the nodes by branch-and-bound.
      if ((n%nsteps)==0)
	{		  
  
	  nnewcandseqs=0;
	  valid=malloc(sizeof(int));
	  
	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);	  
	  
	  for (i=0; i<=(ncandseqs-1); i++)
	    {	      
	      minval=1e500;
	      for (l=0; l<=(nstates-1); l++)
		{
		  for (j=0; j<=(ndrugs-1); j++)
		    *(dosage+j)=*(dstates+l*ndrugs+j);
		  recurse_multidrug_response_trimrates2_2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*(finaly+j);
		  minval=min(minval,val);
		}
	      
	      val=0;
	      printf("candseq %d, ( ",i);
	      for (j=0; j<=(ntypes-1); j++)
		{
		  printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j));
		  val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j);
		}
	      printf("), %.3e, %.3e\n",val,minval);	      	      	      
	    }
	  */


	  // For each candidate sequence, run recurse_sequence_multidrug_responses to exhaust all valid following sequences up to nsteps.
	  // If the candidate sequence already leads to mortality, then do not proceed.
	  for (candind=0; candind<=(ncandseqs-1); candind++)
	    {
	      // If the prior treatment sequences lead to zero population, then don't consider new ones.
	      if (extinction==0)
		{
		  // Do not proceed recursion if the current candidate already leads to mortality.
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		  
		  if (val>=mortal)
		    {		      
		      nnewcandseqs++;
		      valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
		      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);		      
		      (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth;
		      (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
		      (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
		      for (j=0; j<=((candseqs+candind)->depth-1); j++)
			{
			  for (k=0; k<=(ndrugs-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			}
		      for (j=0; j<=(candseqs+candind)->depth; j++)
			{
			  for (k=0; k<=(ntypes-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			}		      
		    }
				  
		  else
		    {
		      double *lvals1, *lvals2;

		      /*
		      // Clear up and initialize memory for seqnodes.
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
		      */
		      
		      // Initialize memory for seqnodes.
		      nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
		      (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
		      (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
		      (seqnodes+0)->children=malloc(sizeof(long));
		      (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
		      // Also store val1 and val2 under ->y.
		      //(seqnodes+0)->y=malloc(sizeof(double)*ntypes);
		      (seqnodes+0)->y=malloc(sizeof(double)*(ntypes+2));
		      for (i=0; i<=(ntypes-1); i++)
			*((seqnodes+0)->y+i)=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i); 

		      for (i=0; i<=nsteps; i++)
			*(nfrontiers+i)=0;		  
		      curind=0;	  	  	
		      		    
		      //seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);

		      //seqnodes=recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);

		      // Debug
		      //t1=malloc(sizeof(time_t));  
		      //time(t1);

		      seqnodes=recurse_sequence_multidrug_responses3(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind);		      

		      // Debug
		      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
		      //printf("candind=%d, td0=%.4f\n",candind,val);
		      //free(t1); free(t2);


		      // Debug
		      //printf("candind=%d, nseqnodes=%d\n",candind,nseqnodes);

		      /*
		      // Debug		      
		      //if (n==5)
		      if (n==0)
		      //if (n==5)
		      //if (n==10)
		      //if ((n==5)&&(candind==3))
		      //if ((n==10)||(n==15))	
		      //if ((n==5)&&((candind==0)||(candind==5)))
			{
			  printf("n=%d, candind=%d, nseqnodes=%d\n",n,candind,nseqnodes); 		      		      
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      printf("i=%d, pa=%d, nchildren=%d, depth=%d, ( ",i,(seqnodes+i)->pa,(seqnodes+i)->nchildren,(seqnodes+i)->depth);
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				{
				  printf("%.4e ",*((seqnodes+i)->y+j));
				  val+=*((seqnodes+i)->y+j);
				}
			      printf(") %.4e, ",val);
			      printf("dosage: ( ");
			      for (j=0; j<=(ndrugs-1); j++)
				printf("%.2f ",*((seqnodes+i)->dosage+j));
			      printf(")\n");
			    }
			  //exit(0);
			}
		      */      
      

		      /*
		      // Debug
		      //if ((n==0)&&(candind==0))
		      //if (n==5)
		      //if ((n==5)&&((candind==0)||(candind==5)))
		      //if (n==0)
			{
			  printf("n=%d, candind=%d\n",n,candind);

			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if ((seqnodes+i)->nchildren==0)
				{
				  int curind;

				  printf("i=%d, depth=%d, pa=%d, dosage sequence is ",i,(seqnodes+i)->depth,(seqnodes+i)->pa);
				  curind=i;
				  while (curind>=0)
				    {
				      printf("( ");
				      for (j=0; j<=(ndrugs-1); j++)
					printf("%.2f ",*((seqnodes+curind)->dosage+j));
				      printf(") ");
				      curind=(seqnodes+curind)->pa;
				    }
				  val=0;
				  printf("pop is ( ");
				  for (j=0; j<=(ntypes-1); j++)
				    {
				      val+=*((seqnodes+i)->y+j);
				      printf("%.3e ",*((seqnodes+i)->y+j));
				    }
				  printf(") %.3e, ",val);
				  printf("bd is (%.3e, %.3e)\n",*((seqnodes+i)->y+ntypes),*((seqnodes+i)->y+ntypes+1));
				}
			    }
			  //exit(0);
			}
		      */



		      lvals1=malloc(sizeof(double)*nseqnodes);
		      lvals2=malloc(sizeof(double)*nseqnodes);
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  *(lvals1+i)=0; *(lvals2+i)=0;
			}


		      maxdepth=0;
		      for (i=0; i<=(nseqnodes-1); i++)
			if ((seqnodes+i)->nchildren==0)
			  maxdepth=max(maxdepth,(seqnodes+i)->depth);

		    
		      // Find the candidate terminal seqnoces that yield almost zero total population.
		      // Each subpopulation < 1, and this trend lasts for maxdepth-1 steps.		      
		      ncandidates=0; candidates=malloc(sizeof(long));
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  // When the population dynamics lasts for maxdepth steps.
			  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			    {	
			      int flag=1, curind=i;

			      // If maxdepth=1, then check the population composition at the current level.
			      // If maxdepth>1, then check the population composition up to (maxdepth-1) levels.
			      if (maxdepth==1)
				m=maxdepth;
			      else
				m=maxdepth-1;				  
			      l=0;
			      //while ((l<(maxdepth-1))&&(flag==1)&&(curind>=0))
			      while ((l<m)&&(flag==1)&&(curind>=0))
				{
				  j=0;
				  while ((j<=(ntypes-1))&&(flag==1))
				    {				      
				      if (*((seqnodes+curind)->y+j)>=1)
				      //if (*((seqnodes+curind)->y+j)>1)
					flag=0;
				      j++;
				    }
				  //if (*((seqnodes+curind)->y+ntypes-1)>=0.1)
				  //flag=0;
				  curind=(seqnodes+curind)->pa; l++;
				}

			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}

			      /*
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)		
				val+=*((seqnodes+i)->y+j);
			      if (val<1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;
				}
			      */
			    }

			  
			  // When the population dynamics ends early.  Add a candidate if each subpopulation < 1.
			  else if ((seqnodes+i)->nchildren==0)
			    {
			      int flag=1; 
			      j=0;
			      while ((j<=(ntypes-1))&&(flag==1))
				{				      
				  if (*((seqnodes+i)->y+j)>=1)
				    flag=0;
				  j++;
				}
			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }
			    			  
			}

		      
		      // Debug
		      //printf("n=%d, candind=%d, ncandidates=%d\n",n,candind,ncandidates);


		      
		      // If the no candidates present negligible population, then check whether static therapies can cure the patients.  If so, then label them as curable.
		      if (ncandidates<=0)
			{

			  // Debug
			  //t1=malloc(sizeof(time_t));  
			  //time(t1);

			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  int flag=0;	
				  l=0; 
				  while ((l<=(nstates-1))&&(flag==0))
				    {
				      for (j=0; j<=(ndrugs-1); j++)
					*(dosage+j)=*(dstates+l*ndrugs+j);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
				      
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
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;				 

				      // Debug
				      //printf("finaly: ( ");
				      //for (l=0; l<=(ntypes-1); l++)
				      //printf("%.3e ",*(finaly+l));
				      //printf(")\n");

				    }
				}
			    }

			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td1=%.4f\n",val);
			  //free(t1); free(t2);
			}				      				  
		      

		      
		      // Debug
		      //ncandidates=0;


		      // If the candidate sequences lead to zero population, then pick up any instance sequence and stop iteration.  Nullify the sequences derived from other candidates.
		      if (ncandidates>0)
			{
			  extinction=1;
			  for (i=0; i<=(nnewcandseqs-1); i++)
			    *(valid+i)=0;
			  maxdepth=nsteps;
			}
		      
		      // If no seqnodes yield zero total population, then keep the ones not inferior to others.
		      else
			{
			  int nnewfrontiers, nleft, *visited;
			  double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;			  
			  
			  // Debug
			  //t1=malloc(sizeof(time_t)); time(t1);

			  // Establish frontiers of terminal nodes.			  
			  nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
			  newfrontiers2=malloc(sizeof(double));
			  newfrontiers3=malloc(sizeof(double));			  
			  template=malloc(sizeof(double)*ntypes);			  
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome for static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  // val1 and val2 are already evaluated and stored in *((seqnodes+i)->y+ntypes) and *((seqnodes+i)->y+ntypes+1).
				  
				  if (frontier2flag==1)
				    {
				      /*
				      val1=1e5000;
				      for (j=0; j<=(nstates-1); j++)
					{
					  for (l=0; l<=(ndrugs-1); l++)
					    *(dosage+l)=*(dstates+j*ndrugs+l);				    
					  //recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
					  recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
					  val=0;
					  for (l=0; l<=(ntypes-1); l++)
					    val+=*(finaly+l);
					  val1=min(val1,val);
					}
				      */

				      val1=*((seqnodes+i)->y+ntypes);

				      *(lvals1+i)=val1;				      				      
				    }
				  
				  /*
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  */

				  if (frontier2flag==1)
				    {
				      /*
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=1;
				      //recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val2=val; val2=min(val2,val1);
				      */

				      val2=*((seqnodes+i)->y+ntypes+1);

				      *(lvals2+i)=val2;
				    }

				  //free(dosage); free(finaly);
				
				  // Debug
				  //printf("n=%d, candind=%d, i=%d, val1=%.3e, val2=%.3e\n",n,candind,i,val1,val2);
				  

				  // Debug
				  //printf("i=%d, nnewfrontiers=%d\n",i,nnewfrontiers);
				  
				  // Check population composition and predicted final population.
				  // If the current seqnode lies below every new frontier, then it is a new frontier.
				  // If the best of the current seqnode is larger than the worst of any new frontier, then it is not a new frontier.
				  // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
				      
				      /*
				      // Debug
				      //if ((n==10)&&(candind==14)&&(i==79))
				      //if ((n==5)&&(candind==3)&&(i==100)&&(j==2))
				      if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
					{
					  printf("node %d val (%f %f %f %f), ",i,*((seqnodes+i)->y+0),*((seqnodes+i)->y+1),*((seqnodes+i)->y+2),*((seqnodes+i)->y+3));
					  printf("ft %d val (%f %f %f %f), ",j,*(newfrontiers1+j*ntypes+0),*(newfrontiers1+j*ntypes+1),*(newfrontiers1+j*ntypes+2),*(newfrontiers1+j*ntypes+3));
					  m=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  printf("beloweq=%d\n",m);
					}
				      */


				      if (l==1)
					k=j;
				      
				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;

					  /*
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    {
					      //if (*(newfrontiers2+j)>*(newfrontiers3+j))
					      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
						k=j;
					    }
					  */
					}					  
				      j++;
				    }
				
				  // Debug
				  //if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
				  //printf("node %d, k=%d, nnewfrontiers=%d\n",i,k,nnewfrontiers);

				  /*
				  // Check population composition and predicted final population.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=0; k=j;
				      while ((l<=(ntypes-1))&&(k>=0))
					{
					  if (*((seqnodes+i)->y+l)<=*(newfrontiers1+j*ntypes+l))
					    k=-1;
					  l++;
					}
				      if (k<0)
					{
					  l=0; m=-1;
					  while ((l<=(ntypes-1))&&(m<0))
					    {
					      if (fabs(*((seqnodes+i)->y+l)-*(newfrontiers1+j*ntypes+l))>0.01)
						m=l;
					      l++;
					    }
					  if (m<0)
					    k=j;
					}
				      //if ((k<0)&&(val2>*(newfrontiers2+j)))
				      //k=j;					  
				      j++;
				    }
				  */
				
				  
				  // Debug
				  //if ((n==10)&&(candind==14)&&(i==79))
				  //printf("k=%d\n",k);


				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      // An unqualified frontier lies above the current node.
				      // The best of an unqualified frontier is larger than the worst of the current node.
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;

					  l=vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);
					  if (l==1)
					    qualified=0;
					 
					  // Debug
					  //if ((n==10)&&(candind==14)&&(i==79)&&(j==0))
					  //l=deb_vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);

					  /*
					  // Debug
					  if ((n==0)&&(candind==0)&&(i==31))
					    {
					      printf("j=%d, ft=( ",j);
					      for (l=0; l<=(ntypes-1); l++)
						printf("%f ",*(newfrontiers1+j*ntypes+l));
					      printf("), qualified=%d\n",qualified);
					    }
					  */

					  if (frontier2flag==1)
					    {
					      //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						qualified=0;
					      
					      /*
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						{
						  // If the upper and lower bounds are identical, then don't discard the old frontiers.
						  //if (*(newfrontiers2+j)>*(newfrontiers3+j))
						  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(*(newfrontiers3+j)>=10))
						    qualified=0;
						}
					      */
					    }
					  					 
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}

				      /*
				      // Debug
				      if ((n==5)&&(candind==0))
					{
					  printf("n=%d, candind=%d, i=%d, nnewfrontiers=%d, ft1=( ",n,candind,i,nnewfrontiers);
					  for (l=0; l<=(ntypes-1); l++)
					    printf("%.3e ",*((seqnodes+i)->y+l));
					  printf(")\n");
					}
				      */
				      
				      /*
				      // Debug
				      if ((n==10)&&(candind==14)&&(i==79))
					{
					  printf("nnewfrontiers=%d: ",nnewfrontiers);
					  for (j=0; j<=(nnewfrontiers-1); j++)
					    {
					      printf("( ");
					      for (l=0; l<=(ntypes-1); l++)
						printf("%f ",*(newfrontiers1+j*ntypes+l));
					      printf("%.4e %.4e)\n",*(newfrontiers2+j),*(newfrontiers3+j));
					    }
					  //exit(0);
					}
				      */
				    }


				  /*
				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;
					  l=0;
					  while ((l<=(ntypes-1))&&(qualified==1))
					    {
					      if (*(newfrontiers1+j*ntypes+l)>*((seqnodes+i)->y+l))
						qualified=0;
					      l++;
					    }
					  //if ((qualified==1)&&(val2>*(newfrontiers2+j)))
					  //qualified=0;
					  					  
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}
				      
				      //
				      // Debug
				      //printf("nnewfrontiers=%d: ",nnewfrontiers);
				      //for (j=0; j<=(nnewfrontiers-1); j++)
				      //{
				      //  printf("( ");
				      //  for (l=0; l<=(ntypes-1); l++)
				      //    printf("%.3e ",*(newfrontiers1+j*ntypes+l));
				      //  printf("%.3e %.3e)\n",*(newfrontiers2+j),*(newfrontiers3+j));
				      //}
				      //exit(0);
				      //
				    }
				  */
				}
			    }				      
			
			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td2=%.4f\n",val);
			  //free(t1); free(t2);


			  free(template);
	  		
			  /*
			  // Debug
			  if (n==0)
			    {
			      for (i=0; i<=(nseqnodes-1); i++)
				{
				  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				    printf("i=%d, val1=%.3e, val2=%.3e\n",i,*(lvals1+i),*(lvals2+i));
				}
			    }
			  */	  			  
			 
			  /*
			  // Debug
			  printf("nnewfrontiers=%d\n",nnewfrontiers);			  
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    printf("%.3e %.3e\n",*(newfrontiers+i*2+0),*(newfrontiers+i*2+1));
			  exit(0);
			  */

			  // Debug
			  //printf("n=%d, candind=%d, nnewfrontiers=%d\n",n,candind,nnewfrontiers);

			  /*
			  // Debug
			  if ((n==0)&&(candind==0))
			    {
			      printf("nnewfrontiers=%d\n",nnewfrontiers);
			      
			      for (i=0; i<=(nnewfrontiers-1); i++)
				{
				  printf("(%d ",i);
				  for (j=0; j<=(ntypes-1); j++)
				    printf("%f ",*(newfrontiers1+i*ntypes+j));
				  printf(") ");
				  printf("%f %f ",*(newfrontiers2+i),*(newfrontiers3+i));
				  printf("\n");
				}
			    }			  			  
			  */
			
			  // Debug
			  //t1=malloc(sizeof(time_t));  
			  //time(t1);

			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome of static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  if (frontier2flag==1)
				    {
				      /*
				      val1=1e5000;
				      for (j=0; j<=(nstates-1); j++)
					{
					  for (l=0; l<=(ndrugs-1); l++)
					    *(dosage+l)=*(dstates+j*ndrugs+l);				      
					  recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
					  val=0;
					  for (l=0; l<=(ntypes-1); l++)
					    val+=*(finaly+l);
					  val1=min(val1,val);  					  
					}
				      */
				      val1=*(lvals1+i);
				    }

				  /*
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  */

				  if (frontier2flag==1)
				    {
				      /*
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val2=val; val2=min(val2,val1); 
				      */
				      val2=*(lvals2+i);
				    }

				  //free(dosage); free(finaly);
				
				  // Check whether the current node lies on the new frontier.
				  // If the best of the current node is larger than the worst of the new frontier, then it is not a nwe frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y); 
				      if (l==1)
					k=j;

				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;

					  /*
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    {
					      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
						k=j;
					    }
					  */
					}	
				      j++;
				    }				  

				  // Debug
				  //if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
				  //printf("i=%d, k=%d\n",i,k);

				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=vectoreq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  if (l==1)
					    k=j;					  
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i; 				      
				    }
				}
			    }

			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td3=%.4f\n",val);
			  //free(t1); free(t2);


			  /*
			  // Debug
			  //if ((n==10)&&(candind==14))
			  //if ((n==5)&&(candind==3))
			  //if ((n==0)&&(candind==0))
			  if ((n==5)&&(candind==0))
			    {
			      printf("n=%d, candind=%d, ncandidates=%d\n",n,candind,ncandidates);
			      for (i=0; i<=(ncandidates-1); i++)
				printf("(%d %d) ",i,*(candidates+i));
			      printf("\n");
			    }
			  */


			  /*
			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  dosage=malloc(sizeof(double)*ndrugs); 
				  finaly=malloc(sizeof(double)*ntypes);
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  for (l=0; l<=(ndrugs-1); l++)
				    *(dosage+l)=1;
				  recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				  val=0;
				  for (l=0; l<=(ntypes-1); l++)
				    val+=*(finaly+l);
				  val2=val; 
				  free(dosage); free(finaly);
				  
				  // Check population composition and predicted final population.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=0; k=j;
				      while ((l<=(ntypes-1))&&(k>=0))
					{
					  if (*((seqnodes+i)->y+l)<=*(newfrontiers1+j*ntypes+l))
					    k=-1;
					  l++;
					}				      
				      //if ((k<0)&&(val2>*(newfrontiers2+j)))
				      //k=j;					  
				      j++;
				    }

				  
				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=0; k=j;
					  while ((l<=(ntypes-1))&&(k>=0))
					    {
					      if (fabs(*((seqnodes+i)->y+l)-*(newfrontiers1+j*ntypes+l))>0.01)
						k=-1;
					      l++;
					    }
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  // Debug
					  printf("i=%d, visited[%d]=1\n",i,k);

					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;
				      				      
				    }
				}
			    }
			  */
			  			  
			  free(visited); 			  			  
			  free(newfrontiers1); free(newfrontiers2); free(newfrontiers3); 

			  // Debug
			  //printf("n=%d, candind=%d, ncandidates=%d, nnewfrontiers=%d\n",n,candind,ncandidates,nnewfrontiers);

			  
			}			    

		      /*
		      // Debug
		      if ((n==0)&&(candind==0))
			{
			  printf("ncandidates=%d, candidates: ",ncandidates);
			  for (i=0; i<=(ncandidates-1); i++)
			    printf("%d ",*(candidates+i));
			  printf("\n");
			}
		      */
      		    
		      /*
		      // Debug
		      printf("n=%d, candind=%d, ncandidates=%d, nnewcandseqs=%d\n",n,candind,ncandidates,nnewcandseqs);
		      
		      //for (i=0; i<=(ncandidates-1); i++)
		      //{
		      //  printf("%d ",*(candidates+i));
		      //}
		      //printf("\n");

		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val1, val2, tt;
			  double *dosage, *finaly;

			  // Evaluate worst and best case final population.
			  // Here the worst case is the best outcome of static therapies.
			  dosage=malloc(sizeof(double)*ndrugs); 
			  finaly=malloc(sizeof(double)*ntypes);

			  val1=1e5000;
			  for (j=0; j<=(nstates-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l);			    
			      recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  
			  //
			  //val1=0;
			  //for (j=0; j<=(ndrugs-1); j++)
			  //{
			  //  for (l=0; l<=(ndrugs-1); l++)
			  //	*(dosage+l)=0;
			  //  *(dosage+j)=1;
			  //  recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  //  val=0;
			  //    for (l=0; l<=(ntypes-1); l++)
			  //	val+=*(finaly+l);
			  //  val1=max(val1,val);
			  //}
			  //

			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val;
			  //free(dosage); free(finaly);


			  printf("node %d, ( ",*(candidates+i)); val=0;			  
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val+=*((seqnodes+*(candidates+i))->y+j);
			      printf("%.3e ",*((seqnodes+*(candidates+i))->y+j));
			    }
			  printf(") %.3e  %.3e %.3e\n",val,val1,val2);

			  //
			  //tt=maxtime;
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0;
			  //*(dosage+0)=1;
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage1: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0;
			  //*(dosage+1)=1;
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage2: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0.5;			  
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage12: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //
			}
			  
		      //exit(0);
		      //m=nnewcandseqs;
		      */

		      /*
		      // Debug
		      if ((n==5)&&(candind==3))
			{
			  printf("ncandidates=%d: ",ncandidates);
			  for (i=0; i<=(ncandidates-1); i++)
			    printf("(%d %d)",*(candidates+i),(seqnodes+*(candidates+i))->depth);
			  printf("\n");
			}
		      */

		      // Debug
		      //t1=malloc(sizeof(time_t)); time(t1);
		      
		      
		      // Append the candidate sequences to the existing ones.
		      // Also record the val1 and val2 of each population composition.
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  nnewcandseqs++;

			  // Debug
			  //if ((n==5)&&(candind==0)&&((i==0)||(i==148)))
			  //printf("nnewcandseqs-1=%d\n",nnewcandseqs-1);
			  
			  // Write the candidates at the prior level.
			  valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
			  newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);
			  //(newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+nsteps;
			  (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+maxdepth;
			  (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
			  (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
			  for (j=0; j<=((candseqs+candind)->depth-1); j++)
			    {
			      for (k=0; k<=(ndrugs-1); k++)
				*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			    }
			  for (j=0; j<=(candseqs+candind)->depth; j++)
			    {
			      for (k=0; k<=(ntypes-1); k++)
				*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			    }

			  // Append the current candidates to the existing ones.
			  curind=*(candidates+i); l=(newcandseqs+nnewcandseqs-1)->depth;
			  while (curind>0)
			    {
			      for (j=0; j<=(ndrugs-1); j++)			    
				*((newcandseqs+nnewcandseqs-1)->dosages+(l-1)*ndrugs+j)=*((seqnodes+curind)->dosage+j);
			      for (j=0; j<=(ntypes-1); j++)
				*((newcandseqs+nnewcandseqs-1)->ys+l*ntypes+j)=*((seqnodes+curind)->y+j);

			      /*
			      // Debug
			      if ((n==5)&&(candind==0)&&(i==0))
				{
				  printf("curind=%d, dosage=( ",curind);
				  for (j=0; j<=(ndrugs-1); j++)
				    printf("%.2f ",*((seqnodes+curind)->dosage+j));
				  printf(")\n");
				}
			      */

			      /*
			      // Debug
			      if (n==0)
				{
				  printf("candind %d, candidate %d, curind=%d, y ",candind,i,curind);
				  for (j=0; j<=(ntypes-1); j++)
				    printf("%.4e ",*((seqnodes+curind)->y+j));
				  printf("\n");
				}
			      */			      			      

			      curind=(seqnodes+curind)->pa; l--;
			    }

			  gvals1=realloc(gvals1,sizeof(double)*nnewcandseqs);			  
			  *(gvals1+nnewcandseqs-1)=*(lvals1+*(candidates+i));
			  gvals2=realloc(gvals2,sizeof(double)*nnewcandseqs);
			  *(gvals2+nnewcandseqs-1)=*(lvals2+*(candidates+i));

			  /*
			  // Debug
			  //if ((n==5)&&((candind==0)||(candind==12)))
			  if (n==5)
			    {
			      printf("n=%d, candind=%d, i=%d, curind=%d, index=%d, depth=%d\n",n,candind,i,*(candidates+i),nnewcandseqs-1,(newcandseqs+nnewcandseqs-1)->depth);
			      //printf("curdosage: ");
			      //for (k=0; k<=(ndrugs-1); k++)
			      //printf("%.2f ",*((seqnodes+*(candidates+i))->dosage+k));
			      //printf("\n");
			      printf("dosage sequence ");
			      for (j=0; j<=((newcandseqs+nnewcandseqs-1)->depth-1); j++)
				{
				  printf("( ");
				  for (k=0; k<=(ndrugs-1); k++)
				    printf("%.2f ",*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k));
				  printf(") ");
				}
			      printf("\n");
			      printf("pop sequence ");
			      for (j=0; j<=((newcandseqs+nnewcandseqs-1)->depth); j++)
				{
				  printf("( ");
				  for (k=0; k<=(ntypes-1); k++)
				    printf("%.3e ",*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k));
				  printf(") ");
				}
			      printf("\n\n");
			    }
			  */
			  

			  /*
			  // Debug
			  printf("n=%d, candind=%d, i=%d, depth=%d, ( ",n,candind,i,(seqnodes+*(candidates+i))->depth);
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      printf("%.3e ",*((seqnodes+*(candidates+i))->y+j));
			      val+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  printf("), %.3e\n",val);
			  */
			  

			}

		      // Debug
		      //printf("after appending\n");
		      

		      // Debug
		      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
		      //printf("td4=%.4f\n",val);
		      //free(t1); free(t2);


		      free(candidates);
		      free(lvals1); free(lvals2);

		      // Clear up memory for seqnodes.
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  free((seqnodes+i)->children);
			  free((seqnodes+i)->dosage);
			  free((seqnodes+i)->y);
			}
		      free(seqnodes);	     	      
		       
		    }
		  
		  /*
		  // Debug		  
		  printf("n=%d, candind=%d, depth=%d, extinction=%d, nnewcandseqs=%d\n",n,candind,(candseqs+candind)->depth,extinction,nnewcandseqs);
		  val=0; printf("( ");
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      printf("%.3e ",*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i));
		      val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		    }
		  printf(") %.3e\n",val);
		  */		
		}
	    }

	  /*
	  // Debug	  
	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  printf("n=%d, nnewcandseqs=%d, maxdepth=%d, extinction=%d\n",n,nnewcandseqs,maxdepth,extinction);	  	  
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      //if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  printf("newcand %d, valid=%d, depth=%d, ( ",i,*(valid+i),(newcandseqs+i)->depth);
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      printf("%.3e ",*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j));
		      val+=*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j);
		    }
		  printf(") %.3e\n",val);
		}
	    }
	  */

	  // Debug
	  //printf("n=%d, nnewcandseqs=%d\n",n,nnewcandseqs);
	
		
	  // Distinct treatment sequences may lead to identical or inferior terminal population composition.  Establish frontiers again on all treatment sequences and filter out inferior/identical ones.
	  if (nnewcandseqs>0)
	    {
	      int nnewfrontiers, nleft, *visited;
	      double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;
	
	      maxdepth=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		maxdepth=max(maxdepth,(newcandseqs+i)->depth);

	      // Debug
	      //t1=malloc(sizeof(time_t)); time(t1);

	      // Establish frontiers of terminal nodes.			  
	      nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
	      newfrontiers2=malloc(sizeof(double));
	      newfrontiers3=malloc(sizeof(double));
	      template=malloc(sizeof(double)*(ntypes+2));
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{		  		  
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      //double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      // val1 and val2 are already calculated and stored in gvals1 and gvals2.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);
		      
		      if (frontier2flag==1)
			{
			  /*
			  val1=1e5000;
			  for (j=0; j<=(ndrugs-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l);
			      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  */
			  val1=*(gvals1+i);
			}

		      /*
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      */

		      if (frontier2flag==1)
			{
			  /*
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val; val2=min(val2,val1);
			  */
			  val2=*(gvals2+i);
			}

		      //free(dosage); free(finaly);
		    

		      // Check population composition and predicted final population.
		      // If the candidate sequence lies below every new frontier, then it is a new frontier.
		      // If the best of the candidate sequence is larger than the worst of a frontier, then it is not a frontier.
		      // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  //l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			  l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))			      
				k=j;

			      /*
			      if ((k<0)&&(val2>*(newfrontiers2+j)))	
				{
				  //if (*(newfrontiers2+j)>*(newfrontiers3+j))
				  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
				    k=j;
				}
			      */
			    }					  
			  j++;
			}

		      /*
		      // Debug
		      if ((n==5)&&(i==31))
			{
			  j=0;
			  printf("flag1=%d, flag2=%d, flag3=%d\n",vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes),(val2>*(newfrontiers2+j)),*(newfrontiers2+j)>*(newfrontiers3+j));
			}
		      */

		      /*
		      // Debug
		      if ((n==5)&&(i==10))
			{
			  double *dosage, *finaly;

			  dosage=malloc(sizeof(double)*ndrugs);
			  finaly=malloc(sizeof(double)*ntypes);
			  printf("n=%d, i=%d, k=%d\n",n,i,k);
			  printf("ys: ( ");
			  for (l=0; l<=(ntypes-1); l++)
			    //printf("%.3e ",*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l));
			    printf("%.3e ",*((newcandseqs+i)->ys+maxdepth*ntypes+l));
			  printf("), val1=%.3e, val2=%.3e\n",val1,val2);
			  printf("ft1: ( ");
			  for (l=0; l<=(ntypes-1); l++)
			    printf("%.3e ",*(newfrontiers1+k*ntypes+l));
			  printf("), ft2=%.3e, ft3=%.3e\n",*(newfrontiers2+k),*(newfrontiers3+k));
			  //printf("flag1=%d, flag2=%d\n",vectorbeloweq(ntypes,newfrontiers1+k*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes),val2>*(newfrontiers2+k));
			  printf("flag1=%d, flag2=%d\n",vectorbeloweq(ntypes,newfrontiers1+k*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes),val2>*(newfrontiers2+k));
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1.0;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_1=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+0)=1;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_2=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+1)=1;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_3=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.5;			  
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_4=%.3e\n",val);

			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1.0;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_1=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+0)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_2=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+1)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_3=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.5;
			  *(dosage+0)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_4=%.3e\n",val);
			  
			}
		      */
			
		      
		      /*
		      // Check population composition and predicted final population.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=0; k=j;
			  while ((l<=(ntypes-1))&&(k>=0))
			    {
			      if (*((newcandseqs+i)->ys+maxdepth*ntypes+l)<*(newfrontiers1+j*ntypes+l))			      
				k=-1;
			      l++;
			    }				      
			  //if ((k<0)&&(val2>*(newfrontiers2+j)))
			  //k=j;					  
			  j++;
			}
		      */

		      // The current candidate sequence is a new frontier.
		      if (k<0)
			{
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("before nnewfrontiers=%d\n", nnewfrontiers);


			  // Remove unqualified old frontiers.				      
			  // An unqualified frontier lies above the current candidate sequence.
			  // The best of an unqualified frontier is larger than the worst of the current node.
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      
			      //l=vectorbelow(ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes,newfrontiers1+j*ntypes);
			      l=vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes);	
			      if (l==1)
				qualified=0;
			      
			      if (frontier2flag==1)
				{
				  //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))				  
				    qualified=0;
				  /*
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))
				    {
				      // If the upper and lower bounds are identical, then don't discard the old frontiers.
				      //if (*(newfrontiers2+j)>*(newfrontiers3+j))				      
				      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(*(newfrontiers3+j)>=10))
					qualified=0;
				    }
				  */
				}


			      // Debug
			      //if ((n==5)&&(i==29)&&(j==(nnewfrontiers-1)))
			      //printf("qualified=%d, flag1=%d, flag2=%d, flag3=%d, diff=%f\n",qualified,vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes),*(newfrontiers3+j)>val1,*(newfrontiers2+j)>*(newfrontiers3+j),*(newfrontiers2+j)-*(newfrontiers3+j));

			      
			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}					      
			    }
		       
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("middle, nleft=%d\n", nleft);
			  

			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l);
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l);
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;

			  // Debug
			  //if (n==5)
			  //{
			  //  printf("n=%d, i=%d, y=( ",n,i);
			  //  for (l=0; l<=(ntypes-1); l++)
			  //	printf("%.3e ",*((newcandseqs+i)->ys+maxdepth*ntypes+l));
			  //  printf("), val1=%.3e, val2=%.3e\n",val1,val2);
			  //}

			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }

			  /*
			  // Debug
			  if (n==5)
			    {
			      //printf("n=%d, i=%d, newfrontier[%d]=%d, maxdepth=%d: ( ",n,i,nnewfrontiers-1,i,maxdepth);
			      printf("n=%d, i=%d, nnewfrontiers=%d, maxdepth=%d\n",n,i,nnewfrontiers,maxdepth);
			      for (j=0; j<=(nnewfrontiers-1); j++)
				{
				  printf("ft[%d]: ",j);
				  for (l=0; l<=(ntypes-1); l++)
				    printf("%f ",*(newfrontiers1+j*ntypes+l));
				  printf("%f %f\n",*(newfrontiers2+j),*(newfrontiers3+j));
				}
			      //for (l=0; l<=(ntypes-1); l++)
			      //printf("%.3e ",*(template+(nnewfrontiers-1)*(ntypes+2)+l));
			      //printf(")\n");
			    }
			  */
			  
			}


		      /*
		      if (k<0)
			{	
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      l=0;
			      while ((l<=(ntypes-1))&&(qualified==1))
				{
				  if (*(newfrontiers1+j*ntypes+l)>*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l))
				    qualified=0;
				  l++;
				}
			      //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
			      //qualified=0;

			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}	
			    }
			  
			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l); 
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l); 
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers); 			 
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }					 
			}
		      */
		    }
		}				      
	    
	      free(template);

	      // Debug
	      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
	      //printf("td5=%.4f\n",val);
	      //free(t1); free(t2);

	      /*
	      // Debug
	      printf("nnewfrontiers=%d\n",nnewfrontiers);	      	      
	      for (i=0; i<=(nnewfrontiers-1); i++)
		printf("%d %f %f\n",i,*(newfrontiers+i*2+0),*(newfrontiers+i*2+1));
	      */

	      
	      /*
	      // Debug
	      if (n==5)
		{
		  printf("nnewfrontiers=%d\n",nnewfrontiers);	      
		  for (i=0; i<=(nnewfrontiers-1); i++)
		    {
		      printf("%d ",i);
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.3e ",*(newfrontiers1+i*ntypes+j));
		      printf("%.3e %.3e\n",*(newfrontiers2+i),*(newfrontiers3+i));
		    }
		  printf("nnewcandseqs=%d, maxdepth=%d\n",nnewcandseqs,maxdepth);
		  for (i=0; i<=(nnewcandseqs-1); i++)
		    if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		      printf("%d ",i);
		  printf("\n");

		}
	      */
	      	      
	      // Debug
	      //t1=malloc(sizeof(time_t)); time(t1);


	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))		
		    {
		      double val1, val2;
		      //double *dosage, *finaly;

		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);

		      if (frontier2flag==1)
			{
			  /*
			  val1=1e5000;
			  for (j=0; j<=(ndrugs-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l); 			  
			      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  */
			  val1=*(gvals1+i);
			}

		      /*
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      */

		      if (frontier2flag==1)
			{
			  /*
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val; val2=min(val2,val1);
			  */
			  val2=*(gvals2+i);
			}

		      //free(dosage); free(finaly);
		      

		      // Debug
		      //printf("i=%d, ",i);
		      //for (j=0; j<=(ntypes-1); j++)
		      //printf("%.3e ",*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+j));
		      //printf("%.3e %.3e\n",val1,val2);


		      // Check whether the current node lies on the new frontier.
		      // If the best of the current node is larger than the worst of a frontier, then it is not a frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  //l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			  l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				k=j;
			      /*
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				{
				  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
				    k=j;
				}
			      */
			    }	
			  j++;
			}				  

		      // Debug
		      //printf("i=%d, k=%d\n",i,k);
		      
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      //l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			      l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			      if (l==1)
				k=j;					  
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;	
			    }
			}
		      
		      if (k>=0)
			*(valid+i)=0;

		      /*
		      if (k<0)
			{
			  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
			  *(candidates+ncandidates-1)=i; 
			}
		      */
		    }
		}


	      // Debug
	      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
	      //printf("td6=%.4f\n",val);
	      //free(t1); free(t2);


	      /*
	      // Debug	      
	      printf("nnewcandseqs=%d\n",nnewcandseqs);
	      k=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		if (*(valid+i)==1)
		  k++;
	      printf("nvalid=%d\n",k);
	      */
	      

	      /*
	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      dosage=malloc(sizeof(double)*ndrugs); 
		      finaly=malloc(sizeof(double)*ntypes);
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=1;
		      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val2=val;
		      free(dosage); free(finaly);
		      	      

		      // Check population composition and predicted final population.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=0; k=j;
			  while ((l<=(ntypes-1))&&(k>=0))
			    {
			      if (*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l)<=*(newfrontiers1+j*ntypes+l))
				k=-1;				  
			      l++;
			    }				      
			  //if ((k<0)&&(val2>*(newfrontiers2+j)))
			  //k=j;					  
			  j++;
			}
		      
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      l=0; k=j;
			      while ((l<=(ntypes-1))&&(k>=0))
				{
				  if (fabs(*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l)-*(newfrontiers1+j*ntypes+l))>0.01)
				    k=-1;
				  l++;
				}
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;
			    }
			}

		      if (k>=0)
			*(valid+i)=0;

		    }
		}
	      */

	      free(visited); free(newfrontiers1); free(newfrontiers2); free(newfrontiers3);
	    }     
	  
		  
	 

	  /*
	  // Debug
	  printf("nnewcandseqs=%d\n",nnewcandseqs);	  	  	  	  	  
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      printf("newcandseq %d, depth %d, valid=%d\n",i,(newcandseqs+i)->depth,*(valid+i));
	      printf("dosages: ");
	      for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((newcandseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(newcandseqs+i)->depth; j++)
		{
		  printf("(%d ",j);
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((newcandseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  */


	  // Debug
	  //if (n==5)
	  //printf("nnewcandseqs=%d, valid[0]=%d, valid[148]=%d\n",nnewcandseqs,*(valid+0),*(valid+148));
	
	  // Update the candidate sequences.
	  // Remove the new candidate sequeces that are inferior to others.
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      free((candseqs+i)->dosages); free((candseqs+i)->ys);
	    }	  
	  ncandseqs=0; 

	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  ncandseqs++; candseqs=realloc(candseqs,sizeof(struct sequence)*ncandseqs);
		  (candseqs+ncandseqs-1)->depth=(newcandseqs+i)->depth;
		  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(candseqs+ncandseqs-1)->depth);
		  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*((candseqs+ncandseqs-1)->depth+1));
		  for (j=0; j<=((candseqs+ncandseqs-1)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((candseqs+ncandseqs-1)->dosages+j*ndrugs+k)=*((newcandseqs+i)->dosages+j*ndrugs+k);
		  for (j=0; j<=(candseqs+ncandseqs-1)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((candseqs+ncandseqs-1)->ys+j*ntypes+k)=*((newcandseqs+i)->ys+j*ntypes+k);

		  cvals1=realloc(cvals1,sizeof(double)*ncandseqs);
		  *(cvals1+ncandseqs-1)=*(gvals1+i);
		  cvals2=realloc(cvals2,sizeof(double)*ncandseqs);
		  *(cvals2+ncandseqs-1)=*(gvals2+i);

		}
	    }
	  
	

	  // Clear up memory for new candidates.
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      free((newcandseqs+i)->dosages);
	      free((newcandseqs+i)->ys);
	    }
	  free(newcandseqs); newcandseqs=malloc(sizeof(struct sequence));
	  free(valid); 
	
	  // Debug
	  //printf("ncandseqs=%d\n",ncandseqs);
	  
	  
	  // If ncandseqs exceeds the maximum number, then only keep maxncandseqs of them and report truncation.
	  // Sort candseqs by their final populations with the best static treatments.
	  // Sort candseqs by the geometric means of the two boundary values (best static treatment and full dosage treatments).
	  if (ncandseqs>=maxncandseqs)
	    {
	      struct pair *pairs;	     
	      double val1, val2;
	      //double *dosage, *finaly;

	      // Evalaute the best outcome of static therapies.
	      // Sort the candidates by these outcomes.
	      //dosage=malloc(sizeof(double)*ndrugs); 
	      //finaly=malloc(sizeof(double)*ntypes);
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  /*
		  val1=1e5000;
		  for (j=0; j<=(nstates-1); j++)
		    {
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=*(dstates+j*ndrugs+l);				    
		      recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val1=min(val1,val);
		    }
		  */
		  val1=*(cvals1+i); val2=*(cvals2+i);
		  (pairs+i)->i=i; 
		  //(pairs+i)->s=val1;
		  (pairs+i)->s=sqrt(val1*val2);

		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      m=maxncandseqs; maxncandseqs=m/10;

	      /*
	      // Debug
	      if (n==0)
		{
		  printf("ncandseqs=%d\n",ncandseqs);
		  for (l=0; l<=(ncandseqs-1); l++)
		    {
		      i=(pairs+l)->i;
		      printf("candseq (%d %d), depth %d\n",l,i,(candseqs+i)->depth);
		      
		      printf("dosages: ");
		      for (j=0; j<=((candseqs+i)->depth-1); j++)
			{
			  printf("( ");
			  for (k=0; k<=(ndrugs-1); k++)
			    printf("%.1f ",*((candseqs+i)->dosages+j*ndrugs+k));
			  printf(") ");
			}
		      printf("\n");		  
		      printf("ys: ");
		      for (j=0; j<=(candseqs+i)->depth; j++)
			{
			  val=0;
			  printf("(%d ",j);
			  for (k=0; k<=(ntypes-1); k++)
			    {
			      printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
			      val+=*((candseqs+i)->ys+j*ntypes+k);
			    }
			  printf(") (%.3e) ",val);
			}
		      printf("\n");
		      printf("bds: (%.3e, %.3e)\n",*(cvals1+i),*(cvals2+i));		  
		    }
		}
	      */


	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs; *truncated=1; 
	      free(pairs); //free(dosage); free(finaly);
	      maxncandseqs=m;


	      /*
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  k=rand(); (pairs+i)->i=i; (pairs+i)->s=(double)k;
		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs;
	      free(pairs); *truncated=1;	    
	      */
	    }	  	 


	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      printf("candseq %d, depth %d\n",i,(candseqs+i)->depth);
	      printf("dosages: ");
	      for (j=0; j<=((candseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((candseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(candseqs+i)->depth; j++)
		{
		  printf("( ");
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  */
	  
	  /*
	  // Debug
	  if (n<=5)
	    {
	      printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  double val1, val2;
		  double *dosage, *finaly;
		  
		  // Evaluate worst and best case final population.
		  dosage=malloc(sizeof(double)*ndrugs); 
		  finaly=malloc(sizeof(double)*ntypes);
		  val1=0;
		  for (j=0; j<=(ndrugs-1); j++)
		    {
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=0;
		      *(dosage+j)=1;
		      recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val1=max(val1,val);
		    }
		  for (l=0; l<=(ndrugs-1); l++)
		    *(dosage+l)=1;
		  recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		  val=0;
		  for (l=0; l<=(ntypes-1); l++)
		    val+=*(finaly+l);
		  val2=val;
		  free(dosage); free(finaly);

		  val=0;
		  printf("candseq %d, depth %d, y=( ",i,(candseqs+i)->depth);
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k));
		      val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k);
		    }
		  printf(") %.3e %.3e %.3e\n",val,val1,val2);
		}
	    }
	  */

	  // Debug
	  //printf("n=%d, stopT=%f, extinction=%d\n",n,stopT,extinction);

	  // Stop when some sequences yield zero population.	  
	  // Report the survival time to be max time span + timeinterval.
	  // This number distinguishes the cases where the patients are not cured but can survive for the max time span.	  
	  if (extinction==1)
	    {
	      //stopT=*(t+nintervals-1);
	      stopT=*(t+nintervals-1)+timeinterval;
	      finaltreatmenttime=*(t+maxdepth);
	    }

	  // Stop when the best candidate treatment sequences lead to mortality.
	  minval=1e5000;
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((candseqs+i)->ys+ntypes*(candseqs+i)->depth+j);
	      minval=min(minval,val);
	    }

	  
	  //if (val>=mortal)
	  if (minval>=mortal)
	    {
	      // Debug
	      //printf("minval=%.3e, maxdepth=%d, t=%.0f\n",minval,maxdepth,*(t+n));

	      if ((maxdepth>=1)&&(maxdepth<=nintervals))
		stopT=*(t+maxdepth-1);
	      else if (maxdepth>nintervals)
		stopT=*(t+nintervals-1);
	      else
		stopT=0;
	      finaltreatmenttime=stopT;
	    }	  
	  

	  /*
	  // Debug
	  //if (n==5)
	  if (n==0)
	  //if (n==3)
	    {
	      printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  printf("candseq %d, depth %d, dosages ",i,(candseqs+i)->depth);
		  for (j=0; j<=((candseqs+i)->depth-1); j++)
		    {
		      printf("( ");
		      for (k=0; k<=(ndrugs-1); k++)
			printf("%.2f ",*((candseqs+i)->dosages+j*ndrugs+k));
		      printf(") ");
		    }
		  printf("\n");
		  printf("pops ");
		  for (j=0; j<=((candseqs+i)->depth); j++)
		    {
		      printf("( ");
		      for (k=0; k<=(ntypes-1); k++)
			printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
		      printf(") ");
		    }
		  printf("\n");
		}
	      exit(0);
	    }
	  */	  
	}


      // Debug
      //stopT=-1;
      
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
	  if (*(y+i*nintervals+nintervals-1)>1)
	    k=i;
	  i++;
	}
      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }

  //if (stopT<0)
  //stopT=*(t+nintervals-1);
  */
  
 
  // Extract the best sequence from candidates.
  // In principle, all the final candidates are equivalent.
  // Pick up the one with minimum total population.
  n=0; minind=-1; minval=1e5000;
  for (n=0; n<=(ncandseqs-1); n++)
    {
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*((candseqs+n)->ys+ntypes*(candseqs+n)->depth+i);
      if (val<minval)
	{
	  minind=n; minval=val;
	}
    }
  
  // Write the selected treatment sequence to dosage and population to y.
  for (n=0; n<=((candseqs+minind)->depth-1); n++)    
    for (i=0; i<=(ndrugs-1); i++)	
      *(dosages+i*nintervals+n)=*((candseqs+minind)->dosages+n*ndrugs+i);
  for (n=0; n<=(candseqs+minind)->depth; n++)
    for (i=0; i<=(ntypes-1); i++)
      *(y+i*nintervals+n)=*((candseqs+minind)->ys+n*ntypes+i);
  if ((minval>=mortal)||(minval<1))
    {
      for (n=(candseqs+minind)->depth; n<=(nintervals-1); n++)
	for (i=0; i<=(ndrugs-1); i++)
	  *(dosages+i*nintervals+n)=0;
      for (n=((candseqs+minind)->depth+1); n<=(nintervals-1); n++)
	for (i=0; i<=(ntypes-1); i++)
	  *(y+i*nintervals+n)=*((candseqs+minind)->ys+(candseqs+minind)->depth*ntypes+i);
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


  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *inity;
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
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

      free(inity);
    }

  /*
  // Debug  
  for (n=0; n<=((candseqs+minind)->depth-1); n++)
    {
      printf("gb t(%d)=%.2f, dosage=( ",n,*(t+n));
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*(dosages+i*nintervals+n));
      printf("), y=( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)	
	{
	  printf("%.4e ",*(y+i*nintervals+n));
	  val+=*(y+i*nintervals+n);
	}
      printf(") %.4e\n",val);
    }
  */

  
  // Release memory.
  free(dstates); free(gvals1); free(gvals2);
  free(cvals1); free(cvals2); free(dosage); free(finaly);

  /*
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); 
  */
  free(nfrontiers); free(frontiers1); 
  free(frontiers2); free(frontiers3); 
  free(resvec); free(bases);
  for (n=0; n<=(ncandseqs-1); n++)
    {
      free((candseqs+n)->dosages);
      free((candseqs+n)->ys);
    }
  free(candseqs); free(tp); free(newcandseqs);
  return stopT;
}


// Check whether vector1 is below vector2 in a high-dimensional space.
// vector1 is belowe vector2 if all components of vector1 <= those of vector2, and vector1!=vector2.
// If vec1 is below vec2 by the current criteiron, but the differences of all components < valuppthre, then nullify the label.
// If the max difference between components of the two vectors < valthre, then choose the strict values for comparison.
int vectorbelow(int ndim, double *vec1, double *vec2)
{
  int i=0, flag=1;
  double val, maxval, newthre=valthre;
  
  while ((i<ndim)&&(flag==1))
    {
      if (*(vec1+i)>(*(vec2+i)+newthre))
	flag=0;
      i++;
    }
  
  if (flag==1)
    {
      i=0; flag=0;
      while ((i<ndim)&&(flag==0))
	{
	  if (fabs(*(vec1+i)-*(vec2+i))>newthre)
	    flag=1;
	  i++;
	}
    }

  maxval=0;
  for (i=0; i<=(ndim-1); i++)	
    if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
      maxval=val;
  
  if (maxval<valuppthre)
    {
      i=0; newthre=0; flag=1;
      while ((i<ndim)&&(flag==1))
	{
	  if (*(vec1+i)>(*(vec2+i)+newthre))
	    flag=0;
	  i++;
	}
      
      if (flag==1)
	{
	  i=0; flag=0;
	  while ((i<ndim)&&(flag==0))
	    {
	      if (fabs(*(vec1+i)-*(vec2+i))>newthre)
		flag=1;
	      i++;
	    }
	}
    }


  
  /*
  if (flag==1)
    {
      maxval=0;
      for (i=0; i<=(ndim-1); i++)	
	if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
	  maxval=val;
      if (maxval<valuppthre)
	flag=0;
    }
  */


  return flag;
}


// Check whether vector1 is below or equal to vector2 in a high-dimensional space.
// vector1 is belowe vector2 if all components of vector1 <= those of vector2.
// If vec1 is not below or equal to vec2 (vec2 is above vec1), but the differences of all components < valuppthre, then re-activate the label.
int vectorbeloweq(int ndim, double *vec1, double *vec2)
{
  int i=0, flag=1;
  double val, maxval, newthre=valthre;

  while ((i<ndim)&&(flag==1))
    {
      if (*(vec1+i)>(*(vec2+i)+newthre))
	flag=0;
      i++;
    }
  
  maxval=0;
  for (i=0; i<=(ndim-1); i++)	
    if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
      maxval=val;
  
  if (maxval<valuppthre)
    {
      i=0; newthre=0; flag=1;
      while ((i<ndim)&&(flag==1))
	{
	  if (*(vec1+i)>(*(vec2+i)+newthre))
	    flag=0;
	  i++;
	}
    }


  /*
  if (flag==0)
    {
      maxval=0;
      for (i=0; i<=(ndim-1); i++)	
	if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
	  maxval=val;
      if (maxval<valuppthre)
	flag=1;
    }
  */

  return flag;
}


// Check whether vector1 and vector2 are equal.
int vectoreq(int ndim, double *vec1, double *vec2)
{
  int i=0, flag=1;

  while ((i<ndim)&&(flag==1))
    {
      if (fabs(*(vec1+i)-*(vec2+i))>eqvalthre)
	flag=0;
      i++;
    }
  
  return flag;
}



// Check whether vector1 is below vector2 in a high-dimensional space.
// vector1 is belowe vector2 if all components of vector1 <= those of vector2, and vector1!=vector2.
// If vec1 is below vec2 by the current criteiron, but the differences of all components < valuppthre, then nullify the label.
// If the max difference between components of the two vectors < valthre, then choose the strict values for comparison.
int deb_vectorbelow(int ndim, double *vec1, double *vec2)
{
  int i=0, flag=1;
  double val, maxval, newthre=valthre;
  
  while ((i<ndim)&&(flag==1))
    {
      if (*(vec1+i)>(*(vec2+i)+newthre))
	flag=0;
      i++;
    }
  
  if (flag==1)
    {
      i=0; flag=0;
      while ((i<ndim)&&(flag==0))
	{
	  if (fabs(*(vec1+i)-*(vec2+i))>newthre)
	    flag=1;
	  i++;
	}
    }

  maxval=0;
  for (i=0; i<=(ndim-1); i++)	
    if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
      maxval=val;

  printf("vec1: ");
  for (i=0; i<=(ndim-1); i++)	
    printf("%f ",*(vec1+i));
  printf("\n");
  printf("vec2: ");
  for (i=0; i<=(ndim-1); i++)	
    printf("%f ",*(vec2+i));
  printf("\n");
  printf("flag=%d, maxval=%f\n",flag,maxval);
  
  if (maxval<valuppthre)
    {
      i=0; newthre=0; flag=1;
      while ((i<ndim)&&(flag==1))
	{
	  if (*(vec1+i)>(*(vec2+i)+newthre))
	    flag=0;
	  i++;
	}
      
      if (flag==1)
	{
	  i=0; flag=0;
	  while ((i<ndim)&&(flag==0))
	    {
	      if (fabs(*(vec1+i)-*(vec2+i))>newthre)
		flag=1;
	      i++;
	    }
	}
    }

  printf("newflag=%d\n",flag);

  
  /*
  if (flag==1)
    {
      maxval=0;
      for (i=0; i<=(ndim-1); i++)	
	if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
	  maxval=val;
      if (maxval<valuppthre)
	flag=0;
    }
  */


  return flag;
}


// Check whether vector1 is below or equal to vector2 in a high-dimensional space.
// vector1 is belowe vector2 if all components of vector1 <= those of vector2.
// If vec1 is not below or equal to vec2 (vec2 is above vec1), but the differences of all components < valuppthre, then re-activate the label.
int deb_vectorbeloweq(int ndim, double *vec1, double *vec2)
{
  int i=0, flag=1;
  double val, maxval, newthre=valthre;

  while ((i<ndim)&&(flag==1))
    {
      if (*(vec1+i)>(*(vec2+i)+newthre))
	flag=0;
      i++;
    }
  
  maxval=0;
  for (i=0; i<=(ndim-1); i++)	
    if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
      maxval=val;

  // Debug
  printf("vec1: ");
  for (i=0; i<=(ndim-1); i++)
    printf("%f ",*(vec1+i));
  printf("\n");
  printf("vec2: ");
  for (i=0; i<=(ndim-1); i++)
    printf("%f ",*(vec2+i));
  printf("\n");
  printf("maxval=%f, flag=%d\n",maxval,flag);

  
  //if (maxval<valuppthre)
  if (maxval<1)
    {
      i=0; newthre=0; flag=1;
      while ((i<ndim)&&(flag==1))
	{
	  if (*(vec1+i)>(*(vec2+i)+newthre))
	    flag=0;
	  i++;
	}
    }


  /*
  if (flag==0)
    {
      maxval=0;
      for (i=0; i<=(ndim-1); i++)	
	if ((val=fabs(*(vec1+i)-*(vec2+i)))>maxval)
	  maxval=val;
      if (maxval<valuppthre)
	flag=1;
    }
  */

  return flag;
}



// Find a global optimal treatment sequence using dynamic programming.
// Traverse the decision tree of all treatment sequences and pick the best one.
// Apply branch and bound to trim inferior branches.
// Branch and bound cannot retract previously established branches.
// Stop every nsteps and keep only superior sequences.  Continue traversing.
// The ultimate criterion is the survival time.
// Difference from global_dp_optimize_multidrug_responses1: apply more rigorous criteria of population composition and worst-case final population to set frontiers.
// One criterion is the predicted final population.  Since evaluating population dynamics is time-consuming, make sure that each population composition only evaluates val1 and val2 once.
// Set candseqs to have the following forms.
// At each time step, dosage specifies the treatment dosage through the interval, population specifies the initial population composition.
// Thus the number of population compositions = the number of dosages + 1.
// Difference from global_dp_optimize_multidrug_responses2: do not terminate early when static therapies work.

double global_dp_optimize_multidrug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, combind, *resvec, *bases;
  long curind, ncandseqs, nnewcandseqs, minind;
  double *dstates, val, stopT, *frontiers1, *frontiers2, *frontiers3, pop, incurablepop, finaltreatmenttime=0, minval, maxtime;
  double *gvals1, *gvals2, *cvals1, *cvals2, *dosage, *finaly;  
  struct seqnode *seqnodes; 
  struct sequence *candseqs, *newcandseqs;
  time_t *tp;  
  

  *truncated=0; maxtime=nintervals*timeinterval;
  tp=malloc(sizeof(time_t)); srand((int)(time(tp)));
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
	*(dosages+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }


  n=0; nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers1=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*ntypes);
  frontiers2=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  frontiers3=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  ncandseqs=1; candseqs=malloc(sizeof(struct sequence));
  (candseqs+ncandseqs-1)->depth=0;
  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs);
  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*2);
  for (i=0; i<=(ntypes-1); i++)
    *((candseqs+ncandseqs-1)->ys+i)=*(y+i*nintervals+0);
  nnewcandseqs=0; newcandseqs=malloc(sizeof(struct sequence));
  gvals1=malloc(sizeof(double)); gvals2=malloc(sizeof(double));
  cvals1=malloc(sizeof(double)); cvals2=malloc(sizeof(double));
  dosage=malloc(sizeof(double)*ndrugs);
  finaly=malloc(sizeof(double)*ntypes);


  // Consider all treatment sequences throughout nintervals.
  while ((n<(nintervals-1))&&(stopT<0))  
    {      
      int maxdepth, maxind, minind, candind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1, extinction=0, *valid;

      // (n%nsteps)=0: incur dynamic programming.
      // For each candidate sequence, run the recursive program to generate all descendant nodes of depth nsteps.
      // Discard the nodes by branch-and-bound.
      if ((n%nsteps)==0)
	{		  
  
	  nnewcandseqs=0;
	  valid=malloc(sizeof(int));
	  
	  
	  // Debug
	  //printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  


	  // For each candidate sequence, run recurse_sequence_multidrug_responses to exhaust all valid following sequences up to nsteps.
	  // If the candidate sequence already leads to mortality, then do not proceed.
	  for (candind=0; candind<=(ncandseqs-1); candind++)
	    {
	      // If the prior treatment sequences lead to zero population, then don't consider new ones.
	      if (extinction==0)
		{
		  // Do not proceed recursion if the current candidate already leads to mortality.
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		  
		  if (val>=mortal)
		    {		      
		      nnewcandseqs++;
		      valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
		      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);		      
		      (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth;
		      (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
		      (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
		      for (j=0; j<=((candseqs+candind)->depth-1); j++)
			{
			  for (k=0; k<=(ndrugs-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			}
		      for (j=0; j<=(candseqs+candind)->depth; j++)
			{
			  for (k=0; k<=(ntypes-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			}		      
		    }
				  
		  else
		    {
		      double *lvals1, *lvals2;
		     		      
		      // Initialize memory for seqnodes.
		      nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
		      (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
		      (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
		      (seqnodes+0)->children=malloc(sizeof(long));
		      (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
		      // Also store val1 and val2 under ->y.
		      //(seqnodes+0)->y=malloc(sizeof(double)*ntypes);
		      (seqnodes+0)->y=malloc(sizeof(double)*(ntypes+2));
		      for (i=0; i<=(ntypes-1); i++)
			*((seqnodes+0)->y+i)=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i); 

		      for (i=0; i<=nsteps; i++)
			*(nfrontiers+i)=0;		  
		      curind=0;	  	  	
		      		    		     
		      seqnodes=recurse_sequence_multidrug_responses3(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind);		      

		     
		      lvals1=malloc(sizeof(double)*nseqnodes);
		      lvals2=malloc(sizeof(double)*nseqnodes);
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  *(lvals1+i)=0; *(lvals2+i)=0;
			}


		      maxdepth=0;
		      for (i=0; i<=(nseqnodes-1); i++)
			if ((seqnodes+i)->nchildren==0)
			  maxdepth=max(maxdepth,(seqnodes+i)->depth);

		    
		      // Find the candidate terminal seqnoces that yield almost zero total population.
		      // Each subpopulation < 1, and this trend lasts for maxdepth-1 steps.		      
		      ncandidates=0; candidates=malloc(sizeof(long));
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  // When the population dynamics lasts for maxdepth steps.
			  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			    {	
			      int flag=1, curind=i;

			      // If maxdepth=1, then check the population composition at the current level.
			      // If maxdepth>1, then check the population composition up to (maxdepth-1) levels.
			      if (maxdepth==1)
				m=maxdepth;
			      else
				m=maxdepth-1;				  
			      l=0;
			      while ((l<m)&&(flag==1)&&(curind>=0))
				{
				  j=0;
				  while ((j<=(ntypes-1))&&(flag==1))
				    {				      
				      if (*((seqnodes+curind)->y+j)>=1)				      
					flag=0;
				      j++;
				    }
				  curind=(seqnodes+curind)->pa; l++;
				}

			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }

			  
			  // When the population dynamics ends early.  Add a candidate if each subpopulation < 1.
			  else if ((seqnodes+i)->nchildren==0)
			    {
			      int flag=1; 
			      j=0;
			      while ((j<=(ntypes-1))&&(flag==1))
				{				      
				  if (*((seqnodes+i)->y+j)>=1)
				    flag=0;
				  j++;
				}
			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }			    			  
			}
		    
		      /*
		      // If the no candidates present negligible population, then check whether static therapies can cure the patients.  If so, then label them as curable.
		      if (ncandidates<=0)
			{
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  int flag=0;	
				  l=0; 
				  while ((l<=(nstates-1))&&(flag==0))
				    {
				      for (j=0; j<=(ndrugs-1); j++)
					*(dosage+j)=*(dstates+l*ndrugs+j);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
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
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;				 
				    }
				}
			    }
			}				      				  
		      */


		      // If the candidate sequences lead to zero population, then pick up any instance sequence and stop iteration.  Nullify the sequences derived from other candidates.
		      if (ncandidates>0)
			{
			  extinction=1;
			  for (i=0; i<=(nnewcandseqs-1); i++)
			    *(valid+i)=0;
			  maxdepth=nsteps;
			}
		      
		      // If no seqnodes yield zero total population, then keep the ones not inferior to others.
		      else
			{
			  int nnewfrontiers, nleft, *visited;
			  double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;			  
			  
			  // Establish frontiers of terminal nodes.			  
			  nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
			  newfrontiers2=malloc(sizeof(double));
			  newfrontiers3=malloc(sizeof(double));			  
			  template=malloc(sizeof(double)*ntypes);			  
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome for static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  // val1 and val2 are already evaluated and stored in *((seqnodes+i)->y+ntypes) and *((seqnodes+i)->y+ntypes+1).
				  
				  if (frontier2flag==1)
				    {				      
				      val1=*((seqnodes+i)->y+ntypes);

				      *(lvals1+i)=val1;				      				      
				    }
			
				  if (frontier2flag==1)
				    {				     
				      val2=*((seqnodes+i)->y+ntypes+1);
				      *(lvals2+i)=val2;
				    }

				  // Check population composition and predicted final population.
				  // If the current seqnode lies below every new frontier, then it is a new frontier.
				  // If the best of the current seqnode is larger than the worst of any new frontier, then it is not a new frontier.
				  // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y); 
				      if (l==1)
					k=j;
				      
				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;
					}					  
				      j++;
				    }
							
				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      // An unqualified frontier lies above the current node.
				      // The best of an unqualified frontier is larger than the worst of the current node.
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;

					  l=vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);
					  if (l==1)
					    qualified=0;
					 					 
					  if (frontier2flag==1)
					    {
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						qualified=0; 
					    }
					  					 
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}				      
				    }
				}
			    }				      
			
			  free(template);	  					  			
			
			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome of static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  if (frontier2flag==1)
				    {
				      val1=*(lvals1+i);
				    }
				  
				  if (frontier2flag==1)
				    {
				      val2=*(lvals2+i);
				    }

				  //free(dosage); free(finaly);
				
				  // Check whether the current node lies on the new frontier.
				  // If the best of the current node is larger than the worst of the new frontier, then it is not a nwe frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y); 
				      if (l==1)
					k=j;

				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;
					}	
				      j++;
				    }				  

				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=vectoreq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  if (l==1)
					    k=j;					  
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i; 				      
				    }
				}
			    }
			  			  
			  free(visited); 			  			  
			  free(newfrontiers1); free(newfrontiers2); free(newfrontiers3); 
			 
			}			    

		      // Append the candidate sequences to the existing ones.
		      // Also record the val1 and val2 of each population composition.
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  nnewcandseqs++;

			  // Write the candidates at the prior level.
			  valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
			  newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);
			  //(newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+nsteps;
			  (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+maxdepth;
			  (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
			  (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
			  for (j=0; j<=((candseqs+candind)->depth-1); j++)
			    {
			      for (k=0; k<=(ndrugs-1); k++)
				*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			    }
			  for (j=0; j<=(candseqs+candind)->depth; j++)
			    {
			      for (k=0; k<=(ntypes-1); k++)
				*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			    }

			  // Append the current candidates to the existing ones.
			  curind=*(candidates+i); l=(newcandseqs+nnewcandseqs-1)->depth;
			  while (curind>0)
			    {
			      for (j=0; j<=(ndrugs-1); j++)			    
				*((newcandseqs+nnewcandseqs-1)->dosages+(l-1)*ndrugs+j)=*((seqnodes+curind)->dosage+j);
			      for (j=0; j<=(ntypes-1); j++)
				*((newcandseqs+nnewcandseqs-1)->ys+l*ntypes+j)=*((seqnodes+curind)->y+j); 
			      curind=(seqnodes+curind)->pa; l--;
			    }

			  gvals1=realloc(gvals1,sizeof(double)*nnewcandseqs);			  
			  *(gvals1+nnewcandseqs-1)=*(lvals1+*(candidates+i));
			  gvals2=realloc(gvals2,sizeof(double)*nnewcandseqs);
			  *(gvals2+nnewcandseqs-1)=*(lvals2+*(candidates+i));
		
			}

		      free(candidates);
		      free(lvals1); free(lvals2);

		      // Clear up memory for seqnodes.
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  free((seqnodes+i)->children);
			  free((seqnodes+i)->dosage);
			  free((seqnodes+i)->y);
			}
		      free(seqnodes);	     	      		       
		    }		  			
		}
	    }

	  // Distinct treatment sequences may lead to identical or inferior terminal population composition.  Establish frontiers again on all treatment sequences and filter out inferior/identical ones.
	  if (nnewcandseqs>0)
	    {
	      int nnewfrontiers, nleft, *visited;
	      double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;
	
	      maxdepth=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		maxdepth=max(maxdepth,(newcandseqs+i)->depth);

	      // Establish frontiers of terminal nodes.			  
	      nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
	      newfrontiers2=malloc(sizeof(double));
	      newfrontiers3=malloc(sizeof(double));
	      template=malloc(sizeof(double)*(ntypes+2));
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{		  		  
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      //double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      // val1 and val2 are already calculated and stored in gvals1 and gvals2.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);
		      
		      if (frontier2flag==1)
			{
			  val1=*(gvals1+i);
			}

		      if (frontier2flag==1)
			{
			  val2=*(gvals2+i);
			}


		      // Check population composition and predicted final population.
		      // If the candidate sequence lies below every new frontier, then it is a new frontier.
		      // If the best of the candidate sequence is larger than the worst of a frontier, then it is not a frontier.
		      // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))			      
				k=j;			     
			    }					  
			  j++;
			}
		
		      // The current candidate sequence is a new frontier.
		      if (k<0)
			{
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("before nnewfrontiers=%d\n", nnewfrontiers);


			  // Remove unqualified old frontiers.				      
			  // An unqualified frontier lies above the current candidate sequence.
			  // The best of an unqualified frontier is larger than the worst of the current node.
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      
			      l=vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes);	
			      if (l==1)
				qualified=0;
			      
			      if (frontier2flag==1)
				{
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))				  
				    qualified=0;				  
				}

			      
			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}					      
			    }
		       
			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l);
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;

			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }
			}
		    }
		}				      
	    
	      free(template);


	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))		
		    {
		      double val1, val2;
		      //double *dosage, *finaly;

		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);

		      if (frontier2flag==1)
			{
			  val1=*(gvals1+i);
			}
		     
		      if (frontier2flag==1)
			{
			  val2=*(gvals2+i);
			}

		      // Check whether the current node lies on the new frontier.
		      // If the best of the current node is larger than the worst of a frontier, then it is not a frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				k=j;			      
			    }	
			  j++;
			}				  
  
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			      if (l==1)
				k=j;					  
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;	
			    }
			}
		      
		      if (k>=0)
			*(valid+i)=0;
		    }
		}

	      free(visited); free(newfrontiers1); free(newfrontiers2); free(newfrontiers3);
	    }     
	  		  
	  // Update the candidate sequences.
	  // Remove the new candidate sequeces that are inferior to others.
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      free((candseqs+i)->dosages); free((candseqs+i)->ys);
	    }	  
	  ncandseqs=0; 

	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  ncandseqs++; candseqs=realloc(candseqs,sizeof(struct sequence)*ncandseqs);
		  (candseqs+ncandseqs-1)->depth=(newcandseqs+i)->depth;
		  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(candseqs+ncandseqs-1)->depth);
		  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*((candseqs+ncandseqs-1)->depth+1));
		  for (j=0; j<=((candseqs+ncandseqs-1)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((candseqs+ncandseqs-1)->dosages+j*ndrugs+k)=*((newcandseqs+i)->dosages+j*ndrugs+k);
		  for (j=0; j<=(candseqs+ncandseqs-1)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((candseqs+ncandseqs-1)->ys+j*ntypes+k)=*((newcandseqs+i)->ys+j*ntypes+k);

		  cvals1=realloc(cvals1,sizeof(double)*ncandseqs);
		  *(cvals1+ncandseqs-1)=*(gvals1+i);
		  cvals2=realloc(cvals2,sizeof(double)*ncandseqs);
		  *(cvals2+ncandseqs-1)=*(gvals2+i);

		}
	    }
	  
	
	  // Clear up memory for new candidates.
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      free((newcandseqs+i)->dosages);
	      free((newcandseqs+i)->ys);
	    }
	  free(newcandseqs); newcandseqs=malloc(sizeof(struct sequence));
	  free(valid); 
	
	  
	  // If ncandseqs exceeds the maximum number, then only keep maxncandseqs of them and report truncation.
	  // Sort candseqs by their final populations with the best static treatments.
	  // Sort candseqs by the geometric means of the two boundary values (best static treatment and full dosage treatments).
	  if (ncandseqs>=maxncandseqs)
	    {
	      struct pair *pairs;	     
	      double val1, val2;
	      
	      // Evalaute the best outcome of static therapies.
	      // Sort the candidates by these outcomes.
	      //dosage=malloc(sizeof(double)*ndrugs); 
	      //finaly=malloc(sizeof(double)*ntypes);
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  val1=*(cvals1+i); val2=*(cvals2+i);
		  (pairs+i)->i=i; 
		  (pairs+i)->s=sqrt(val1*val2);

		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      m=maxncandseqs; maxncandseqs=m/10;	     

	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs; *truncated=1; 
	      free(pairs); 
	      maxncandseqs=m;
	    }	  	 


	  // Stop when some sequences yield zero population.	  
	  // Report the survival time to be max time span + timeinterval.
	  // This number distinguishes the cases where the patients are not cured but can survive for the max time span.	  
	  if (extinction==1)
	    {
	      stopT=*(t+nintervals-1)+timeinterval;
	      finaltreatmenttime=*(t+maxdepth);
	    }

	  // Stop when the best candidate treatment sequences lead to mortality.
	  minval=1e5000;
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((candseqs+i)->ys+ntypes*(candseqs+i)->depth+j);
	      minval=min(minval,val);
	    }

	  if (minval>=mortal)
	    {
	      if ((maxdepth>=1)&&(maxdepth<=nintervals))
		stopT=*(t+maxdepth-1);
	      else if (maxdepth>nintervals)
		stopT=*(t+nintervals-1);
	      else
		stopT=0;
	      finaltreatmenttime=stopT;
	    }	  	
	}

      n++;
    }

  // Extract the best sequence from candidates.
  // In principle, all the final candidates are equivalent.
  // Pick up the one with minimum total population.
  n=0; minind=-1; minval=1e5000;
  for (n=0; n<=(ncandseqs-1); n++)
    {
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*((candseqs+n)->ys+ntypes*(candseqs+n)->depth+i);
      if (val<minval)
	{
	  minind=n; minval=val;
	}
    }
  
  // Write the selected treatment sequence to dosage and population to y.
  for (n=0; n<=((candseqs+minind)->depth-1); n++)    
    for (i=0; i<=(ndrugs-1); i++)	
      *(dosages+i*nintervals+n)=*((candseqs+minind)->dosages+n*ndrugs+i);
  for (n=0; n<=(candseqs+minind)->depth; n++)
    for (i=0; i<=(ntypes-1); i++)
      *(y+i*nintervals+n)=*((candseqs+minind)->ys+n*ntypes+i);
  if ((minval>=mortal)||(minval<1))
    {
      for (n=(candseqs+minind)->depth; n<=(nintervals-1); n++)
	for (i=0; i<=(ndrugs-1); i++)
	  *(dosages+i*nintervals+n)=0;
      for (n=((candseqs+minind)->depth+1); n<=(nintervals-1); n++)
	for (i=0; i<=(ntypes-1); i++)
	  *(y+i*nintervals+n)=*((candseqs+minind)->ys+(candseqs+minind)->depth*ntypes+i);
    }


  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *inity;
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
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

      free(inity);
    }

  /*
  // Debug    
  for (n=0; n<=((candseqs+minind)->depth-1); n++)
    {
      printf("gb t(%d)=%.2f, dosage=( ",n,*(t+n));
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*(dosages+i*nintervals+n));
      printf("), y=( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)	
	{
	  printf("%.4e ",*(y+i*nintervals+n));
	  val+=*(y+i*nintervals+n);
	}
      printf(") %.4e\n",val);
    }
  */

  
  // Release memory.
  free(dstates); free(gvals1); free(gvals2);
  free(cvals1); free(cvals2); free(dosage); free(finaly);

  /*
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); 
  */
  free(nfrontiers); free(frontiers1); 
  free(frontiers2); free(frontiers3); 
  free(resvec); free(bases);
  for (n=0; n<=(ncandseqs-1); n++)
    {
      free((candseqs+n)->dosages);
      free((candseqs+n)->ys);
    }
  free(candseqs); free(tp); free(newcandseqs);
  return stopT;
}



// Find a global optimal treatment sequence using dynamic programming.
// Traverse the decision tree of all treatment sequences and pick the best one.
// Apply branch and bound to trim inferior branches.
// Branch and bound cannot retract previously established branches.
// Stop every nsteps and keep only superior sequences.  Continue traversing.
// The ultimate criterion is the survival time.
// Difference from global_dp_optimize_multidrug_responses1: apply more rigorous criteria of population composition and worst-case final population to set frontiers.
// One criterion is the predicted final population.  Since evaluating population dynamics is time-consuming, make sure that each population composition only evaluates val1 and val2 once.
// Set candseqs to have the following forms.
// At each time step, dosage specifies the treatment dosage through the interval, population specifies the initial population composition.
// Thus the number of population compositions = the number of dosages + 1.
// Difference from global_dp_optimize_multidrug_responses2: accelerate matrix exponentiation.

double global_dp_optimize_multidrug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, combind, *resvec, *bases;
  long curind, ncandseqs, nnewcandseqs, minind;
  double *dstates, val, stopT, *frontiers1, *frontiers2, *frontiers3, pop, incurablepop, finaltreatmenttime=0, minval, maxtime;
  double *gvals1, *gvals2, *cvals1, *cvals2, *dosage, *finaly;  
  struct seqnode *seqnodes; 
  struct sequence *candseqs, *newcandseqs;
  time_t *tp, *t1, *t2;  
        

  *truncated=0; maxtime=nintervals*timeinterval;
  tp=malloc(sizeof(time_t)); srand((int)(time(tp)));
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

  /*
  // Debug
  nstates=2; dstates=malloc(sizeof(double)*nstates*ndrugs);
  *(dstates+0)=0; *(dstates+1)=1;
  *(dstates+2)=1; *(dstates+3)=0;
  */

  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosages+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }


  n=0; nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers1=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*ntypes);
  frontiers2=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  frontiers3=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  ncandseqs=1; candseqs=malloc(sizeof(struct sequence));
  (candseqs+ncandseqs-1)->depth=0;
  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs);
  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*2);
  for (i=0; i<=(ntypes-1); i++)
    *((candseqs+ncandseqs-1)->ys+i)=*(y+i*nintervals+0);
  nnewcandseqs=0; newcandseqs=malloc(sizeof(struct sequence));
  gvals1=malloc(sizeof(double)); gvals2=malloc(sizeof(double));
  cvals1=malloc(sizeof(double)); cvals2=malloc(sizeof(double));
  dosage=malloc(sizeof(double)*ndrugs);
  finaly=malloc(sizeof(double)*ntypes);


  // Consider all treatment sequences throughout nintervals.
  while ((n<(nintervals-1))&&(stopT<0))  
    {      
      int maxdepth, maxind, minind, candind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1, extinction=0, *valid;

      // (n%nsteps)=0: incur dynamic programming.
      // For each candidate sequence, run the recursive program to generate all descendant nodes of depth nsteps.
      // Discard the nodes by branch-and-bound.
      if ((n%nsteps)==0)
	{		  
  
	  nnewcandseqs=0;
	  valid=malloc(sizeof(int));
	  
	  
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  /*
	  for (i=0; i<=(ncandseqs-1); i++)
	    {	      
	      minval=1e500;
	      for (l=0; l<=(nstates-1); l++)
		{
		  for (j=0; j<=(ndrugs-1); j++)
		    *(dosage+j)=*(dstates+l*ndrugs+j);
		  recurse_multidrug_response_trimrates2_2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    val+=*(finaly+j);
		  minval=min(minval,val);
		}
	      
	      val=0;
	      printf("candseq %d, ( ",i);
	      for (j=0; j<=(ntypes-1); j++)
		{
		  printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j));
		  val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+j);
		}
	      printf("), %.3e, %.3e\n",val,minval);	      	      	      
	    }
	  */


	  // For each candidate sequence, run recurse_sequence_multidrug_responses to exhaust all valid following sequences up to nsteps.
	  // If the candidate sequence already leads to mortality, then do not proceed.
	  for (candind=0; candind<=(ncandseqs-1); candind++)
	    {
	      // If the prior treatment sequences lead to zero population, then don't consider new ones.
	      if (extinction==0)
		{
		  // Do not proceed recursion if the current candidate already leads to mortality.
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		  
		  if (val>=mortal)
		    {		      
		      nnewcandseqs++;
		      valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
		      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);		      
		      (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth;
		      (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
		      (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
		      for (j=0; j<=((candseqs+candind)->depth-1); j++)
			{
			  for (k=0; k<=(ndrugs-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			}
		      for (j=0; j<=(candseqs+candind)->depth; j++)
			{
			  for (k=0; k<=(ntypes-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			}		      
		    }
				  
		  else
		    {
		      double *lvals1, *lvals2;

		      /*
		      // Clear up and initialize memory for seqnodes.
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
		      */
		      
		      // Initialize memory for seqnodes.
		      nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
		      (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
		      (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
		      (seqnodes+0)->children=malloc(sizeof(long));
		      (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
		      // Also store val1 and val2 under ->y.
		      //(seqnodes+0)->y=malloc(sizeof(double)*ntypes);
		      (seqnodes+0)->y=malloc(sizeof(double)*(ntypes+2));
		      for (i=0; i<=(ntypes-1); i++)
			*((seqnodes+0)->y+i)=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i); 

		      for (i=0; i<=nsteps; i++)
			*(nfrontiers+i)=0;		  
		      curind=0;	  	  	
		      		    
		      //seqnodes=recurse_sequence_multidrug_responses(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);

		      //seqnodes=recurse_sequence_multidrug_responses2(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers,&curind);

		      // Debug
		      //t1=malloc(sizeof(time_t));  
		      //time(t1);

		      seqnodes=recurse_sequence_multidrug_responses4(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind);		      

		      // Debug
		      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
		      //printf("td0=%.4f\n",val);
		      //free(t1); free(t2);


		      // Debug
		      //printf("candind=%d, nseqnodes=%d\n",candind,nseqnodes);

		      /*
		      // Debug		      
		      //if (n==5)
		      if (n==0)
		      //if (n==5)
		      //if (n==10)
		      //if ((n==5)&&(candind==3))
		      //if ((n==10)||(n==15))	
		      //if ((n==5)&&((candind==0)||(candind==5)))
			{
			  printf("n=%d, candind=%d, nseqnodes=%d\n",n,candind,nseqnodes); 		      		      
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      printf("i=%d, pa=%d, nchildren=%d, depth=%d, ( ",i,(seqnodes+i)->pa,(seqnodes+i)->nchildren,(seqnodes+i)->depth);
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)
				{
				  printf("%.4e ",*((seqnodes+i)->y+j));
				  val+=*((seqnodes+i)->y+j);
				}
			      printf(") %.4e, ",val);
			      printf("dosage: ( ");
			      for (j=0; j<=(ndrugs-1); j++)
				printf("%.2f ",*((seqnodes+i)->dosage+j));
			      printf(")\n");
			    }
			  //exit(0);
			}
		      */      
      

		      /*
		      // Debug
		      //if ((n==0)&&(candind==0))
		      //if (n==5)
		      //if ((n==5)&&((candind==0)||(candind==5)))
		      //if (n==0)
			{
			  printf("n=%d, candind=%d\n",n,candind);

			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if ((seqnodes+i)->nchildren==0)
				{
				  int curind;

				  printf("i=%d, depth=%d, pa=%d, dosage sequence is ",i,(seqnodes+i)->depth,(seqnodes+i)->pa);
				  curind=i;
				  while (curind>=0)
				    {
				      printf("( ");
				      for (j=0; j<=(ndrugs-1); j++)
					printf("%.2f ",*((seqnodes+curind)->dosage+j));
				      printf(") ");
				      curind=(seqnodes+curind)->pa;
				    }
				  val=0;
				  printf("pop is ( ");
				  for (j=0; j<=(ntypes-1); j++)
				    {
				      val+=*((seqnodes+i)->y+j);
				      printf("%.3e ",*((seqnodes+i)->y+j));
				    }
				  printf(") %.3e, ",val);
				  printf("bd is (%.3e, %.3e)\n",*((seqnodes+i)->y+ntypes),*((seqnodes+i)->y+ntypes+1));
				}
			    }
			  //exit(0);
			}
		      */



		      lvals1=malloc(sizeof(double)*nseqnodes);
		      lvals2=malloc(sizeof(double)*nseqnodes);
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  *(lvals1+i)=0; *(lvals2+i)=0;
			}


		      maxdepth=0;
		      for (i=0; i<=(nseqnodes-1); i++)
			if ((seqnodes+i)->nchildren==0)
			  maxdepth=max(maxdepth,(seqnodes+i)->depth);

		    
		      // Find the candidate terminal seqnoces that yield almost zero total population.
		      // Each subpopulation < 1, and this trend lasts for maxdepth-1 steps.		      
		      ncandidates=0; candidates=malloc(sizeof(long));
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  // When the population dynamics lasts for maxdepth steps.
			  //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			    {	
			      int flag=1, curind=i;

			      // If maxdepth=1, then check the population composition at the current level.
			      // If maxdepth>1, then check the population composition up to (maxdepth-1) levels.
			      if (maxdepth==1)
				m=maxdepth;
			      else
				m=maxdepth-1;				  
			      l=0;
			      //while ((l<(maxdepth-1))&&(flag==1)&&(curind>=0))
			      while ((l<m)&&(flag==1)&&(curind>=0))
				{
				  j=0;
				  while ((j<=(ntypes-1))&&(flag==1))
				    {				      
				      if (*((seqnodes+curind)->y+j)>=1)
				      //if (*((seqnodes+curind)->y+j)>1)
					flag=0;
				      j++;
				    }
				  //if (*((seqnodes+curind)->y+ntypes-1)>=0.1)
				  //flag=0;
				  curind=(seqnodes+curind)->pa; l++;
				}

			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}

			      /*
			      val=0;
			      for (j=0; j<=(ntypes-1); j++)		
				val+=*((seqnodes+i)->y+j);
			      if (val<1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;
				}
			      */
			    }

			  
			  // When the population dynamics ends early.  Add a candidate if each subpopulation < 1.
			  else if ((seqnodes+i)->nchildren==0)
			    {
			      int flag=1; 
			      j=0;
			      while ((j<=(ntypes-1))&&(flag==1))
				{				      
				  if (*((seqnodes+i)->y+j)>=1)
				    flag=0;
				  j++;
				}
			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }
			    			  
			}

		      
		      // Debug
		      //printf("n=%d, candind=%d, ncandidates=%d\n",n,candind,ncandidates);


		      
		      // If the no candidates present negligible population, then check whether static therapies can cure the patients.  If so, then label them as curable.
		      if (ncandidates<=0)
			{

			  // Debug
			  //t1=malloc(sizeof(time_t));  
			  //time(t1);

			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  int flag=0;	
				  l=0; 
				  while ((l<=(nstates-1))&&(flag==0))
				    {
				      for (j=0; j<=(ndrugs-1); j++)
					*(dosage+j)=*(dstates+l*ndrugs+j);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
				      
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
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;				 

				      // Debug
				      //printf("finaly: ( ");
				      //for (l=0; l<=(ntypes-1); l++)
				      //printf("%.3e ",*(finaly+l));
				      //printf(")\n");

				    }
				}
			    }

			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td1=%.4f\n",val);
			  //free(t1); free(t2);
			}				      				  
		      

		      
		      // Debug
		      //ncandidates=0;


		      // If the candidate sequences lead to zero population, then pick up any instance sequence and stop iteration.  Nullify the sequences derived from other candidates.
		      if (ncandidates>0)
			{
			  extinction=1;
			  for (i=0; i<=(nnewcandseqs-1); i++)
			    *(valid+i)=0;
			  maxdepth=nsteps;
			}
		      
		      // If no seqnodes yield zero total population, then keep the ones not inferior to others.
		      else
			{
			  int nnewfrontiers, nleft, *visited;
			  double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;			  
			  
			  // Debug
			  //t1=malloc(sizeof(time_t)); time(t1);

			  // Establish frontiers of terminal nodes.			  
			  nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
			  newfrontiers2=malloc(sizeof(double));
			  newfrontiers3=malloc(sizeof(double));			  
			  template=malloc(sizeof(double)*ntypes);			  
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome for static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  // val1 and val2 are already evaluated and stored in *((seqnodes+i)->y+ntypes) and *((seqnodes+i)->y+ntypes+1).
				  
				  if (frontier2flag==1)
				    {
				      /*
				      val1=1e5000;
				      for (j=0; j<=(nstates-1); j++)
					{
					  for (l=0; l<=(ndrugs-1); l++)
					    *(dosage+l)=*(dstates+j*ndrugs+l);				    
					  //recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
					  recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
					  val=0;
					  for (l=0; l<=(ntypes-1); l++)
					    val+=*(finaly+l);
					  val1=min(val1,val);
					}
				      */

				      val1=*((seqnodes+i)->y+ntypes);

				      *(lvals1+i)=val1;				      				      
				    }
				  
				  /*
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  */

				  if (frontier2flag==1)
				    {
				      /*
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=1;
				      //recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val2=val; val2=min(val2,val1);
				      */

				      val2=*((seqnodes+i)->y+ntypes+1);

				      *(lvals2+i)=val2;
				    }

				  //free(dosage); free(finaly);
				
				  // Debug
				  //printf("n=%d, candind=%d, i=%d, val1=%.3e, val2=%.3e\n",n,candind,i,val1,val2);
				  

				  // Debug
				  //printf("i=%d, nnewfrontiers=%d\n",i,nnewfrontiers);
				  
				  // Check population composition and predicted final population.
				  // If the current seqnode lies below every new frontier, then it is a new frontier.
				  // If the best of the current seqnode is larger than the worst of any new frontier, then it is not a new frontier.
				  // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
				      
				      /*
				      // Debug
				      //if ((n==10)&&(candind==14)&&(i==79))
				      //if ((n==5)&&(candind==3)&&(i==100)&&(j==2))
				      if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
					{
					  printf("node %d val (%f %f %f %f), ",i,*((seqnodes+i)->y+0),*((seqnodes+i)->y+1),*((seqnodes+i)->y+2),*((seqnodes+i)->y+3));
					  printf("ft %d val (%f %f %f %f), ",j,*(newfrontiers1+j*ntypes+0),*(newfrontiers1+j*ntypes+1),*(newfrontiers1+j*ntypes+2),*(newfrontiers1+j*ntypes+3));
					  m=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  printf("beloweq=%d\n",m);
					}
				      */


				      if (l==1)
					k=j;
				      
				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;

					  /*
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    {
					      //if (*(newfrontiers2+j)>*(newfrontiers3+j))
					      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
						k=j;
					    }
					  */
					}					  
				      j++;
				    }
				
				  // Debug
				  //if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
				  //printf("node %d, k=%d, nnewfrontiers=%d\n",i,k,nnewfrontiers);

				  /*
				  // Check population composition and predicted final population.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=0; k=j;
				      while ((l<=(ntypes-1))&&(k>=0))
					{
					  if (*((seqnodes+i)->y+l)<=*(newfrontiers1+j*ntypes+l))
					    k=-1;
					  l++;
					}
				      if (k<0)
					{
					  l=0; m=-1;
					  while ((l<=(ntypes-1))&&(m<0))
					    {
					      if (fabs(*((seqnodes+i)->y+l)-*(newfrontiers1+j*ntypes+l))>0.01)
						m=l;
					      l++;
					    }
					  if (m<0)
					    k=j;
					}
				      //if ((k<0)&&(val2>*(newfrontiers2+j)))
				      //k=j;					  
				      j++;
				    }
				  */
				
				  
				  // Debug
				  //if ((n==10)&&(candind==14)&&(i==79))
				  //printf("k=%d\n",k);


				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      // An unqualified frontier lies above the current node.
				      // The best of an unqualified frontier is larger than the worst of the current node.
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;

					  l=vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);
					  if (l==1)
					    qualified=0;
					 
					  // Debug
					  //if ((n==10)&&(candind==14)&&(i==79)&&(j==0))
					  //l=deb_vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);

					  /*
					  // Debug
					  if ((n==0)&&(candind==0)&&(i==31))
					    {
					      printf("j=%d, ft=( ",j);
					      for (l=0; l<=(ntypes-1); l++)
						printf("%f ",*(newfrontiers1+j*ntypes+l));
					      printf("), qualified=%d\n",qualified);
					    }
					  */

					  if (frontier2flag==1)
					    {
					      //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						qualified=0;
					      
					      /*
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						{
						  // If the upper and lower bounds are identical, then don't discard the old frontiers.
						  //if (*(newfrontiers2+j)>*(newfrontiers3+j))
						  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(*(newfrontiers3+j)>=10))
						    qualified=0;
						}
					      */
					    }
					  					 
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}

				      /*
				      // Debug
				      if ((n==5)&&(candind==0))
					{
					  printf("n=%d, candind=%d, i=%d, nnewfrontiers=%d, ft1=( ",n,candind,i,nnewfrontiers);
					  for (l=0; l<=(ntypes-1); l++)
					    printf("%.3e ",*((seqnodes+i)->y+l));
					  printf(")\n");
					}
				      */
				      
				      /*
				      // Debug
				      if ((n==10)&&(candind==14)&&(i==79))
					{
					  printf("nnewfrontiers=%d: ",nnewfrontiers);
					  for (j=0; j<=(nnewfrontiers-1); j++)
					    {
					      printf("( ");
					      for (l=0; l<=(ntypes-1); l++)
						printf("%f ",*(newfrontiers1+j*ntypes+l));
					      printf("%.4e %.4e)\n",*(newfrontiers2+j),*(newfrontiers3+j));
					    }
					  //exit(0);
					}
				      */
				    }


				  /*
				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;
					  l=0;
					  while ((l<=(ntypes-1))&&(qualified==1))
					    {
					      if (*(newfrontiers1+j*ntypes+l)>*((seqnodes+i)->y+l))
						qualified=0;
					      l++;
					    }
					  //if ((qualified==1)&&(val2>*(newfrontiers2+j)))
					  //qualified=0;
					  					  
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}
				      
				      //
				      // Debug
				      //printf("nnewfrontiers=%d: ",nnewfrontiers);
				      //for (j=0; j<=(nnewfrontiers-1); j++)
				      //{
				      //  printf("( ");
				      //  for (l=0; l<=(ntypes-1); l++)
				      //    printf("%.3e ",*(newfrontiers1+j*ntypes+l));
				      //  printf("%.3e %.3e)\n",*(newfrontiers2+j),*(newfrontiers3+j));
				      //}
				      //exit(0);
				      //
				    }
				  */
				}
			    }				      
			
			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td2=%.4f\n",val);
			  //free(t1); free(t2);


			  free(template);
	  		
			  /*
			  // Debug
			  if (n==0)
			    {
			      for (i=0; i<=(nseqnodes-1); i++)
				{
				  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				    printf("i=%d, val1=%.3e, val2=%.3e\n",i,*(lvals1+i),*(lvals2+i));
				}
			    }
			  */	  			  
			 
			  /*
			  // Debug
			  printf("nnewfrontiers=%d\n",nnewfrontiers);			  
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    printf("%.3e %.3e\n",*(newfrontiers+i*2+0),*(newfrontiers+i*2+1));
			  exit(0);
			  */

			  // Debug
			  //printf("n=%d, candind=%d, nnewfrontiers=%d\n",n,candind,nnewfrontiers);

			  /*
			  // Debug
			  if ((n==0)&&(candind==0))
			    {
			      printf("nnewfrontiers=%d\n",nnewfrontiers);
			      
			      for (i=0; i<=(nnewfrontiers-1); i++)
				{
				  printf("(%d ",i);
				  for (j=0; j<=(ntypes-1); j++)
				    printf("%f ",*(newfrontiers1+i*ntypes+j));
				  printf(") ");
				  printf("%f %f ",*(newfrontiers2+i),*(newfrontiers3+i));
				  printf("\n");
				}
			    }			  			  
			  */
			
			  // Debug
			  //t1=malloc(sizeof(time_t));  
			  //time(t1);

			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome of static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  if (frontier2flag==1)
				    {
				      /*
				      val1=1e5000;
				      for (j=0; j<=(nstates-1); j++)
					{
					  for (l=0; l<=(ndrugs-1); l++)
					    *(dosage+l)=*(dstates+j*ndrugs+l);				      
					  recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
					  val=0;
					  for (l=0; l<=(ntypes-1); l++)
					    val+=*(finaly+l);
					  val1=min(val1,val);  					  
					}
				      */
				      val1=*(lvals1+i);
				    }

				  /*
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  */

				  if (frontier2flag==1)
				    {
				      /*
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val2=val; val2=min(val2,val1); 
				      */
				      val2=*(lvals2+i);
				    }

				  //free(dosage); free(finaly);
				
				  // Check whether the current node lies on the new frontier.
				  // If the best of the current node is larger than the worst of the new frontier, then it is not a nwe frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y); 
				      if (l==1)
					k=j;

				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;

					  /*
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    {
					      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
						k=j;
					    }
					  */
					}	
				      j++;
				    }				  

				  // Debug
				  //if ((n==0)&&(candind==0)&&((i==30)||(i==31)))
				  //printf("i=%d, k=%d\n",i,k);

				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=vectoreq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  if (l==1)
					    k=j;					  
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i; 				      
				    }
				}
			    }

			  // Debug
			  //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
			  //printf("td3=%.4f\n",val);
			  //free(t1); free(t2);


			  /*
			  // Debug
			  //if ((n==10)&&(candind==14))
			  //if ((n==5)&&(candind==3))
			  //if ((n==0)&&(candind==0))
			  if ((n==5)&&(candind==0))
			    {
			      printf("n=%d, candind=%d, ncandidates=%d\n",n,candind,ncandidates);
			      for (i=0; i<=(ncandidates-1); i++)
				printf("(%d %d) ",i,*(candidates+i));
			      printf("\n");
			    }
			  */


			  /*
			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  dosage=malloc(sizeof(double)*ndrugs); 
				  finaly=malloc(sizeof(double)*ntypes);
				  val1=0;
				  for (j=0; j<=(ndrugs-1); j++)
				    {
				      for (l=0; l<=(ndrugs-1); l++)
					*(dosage+l)=0;
				      *(dosage+j)=1;
				      recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				      val=0;
				      for (l=0; l<=(ntypes-1); l++)
					val+=*(finaly+l);
				      val1=max(val1,val);
				    }
				  for (l=0; l<=(ndrugs-1); l++)
				    *(dosage+l)=1;
				  recurse_multidrug_response_trimrates2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
				  val=0;
				  for (l=0; l<=(ntypes-1); l++)
				    val+=*(finaly+l);
				  val2=val; 
				  free(dosage); free(finaly);
				  
				  // Check population composition and predicted final population.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=0; k=j;
				      while ((l<=(ntypes-1))&&(k>=0))
					{
					  if (*((seqnodes+i)->y+l)<=*(newfrontiers1+j*ntypes+l))
					    k=-1;
					  l++;
					}				      
				      //if ((k<0)&&(val2>*(newfrontiers2+j)))
				      //k=j;					  
				      j++;
				    }

				  
				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=0; k=j;
					  while ((l<=(ntypes-1))&&(k>=0))
					    {
					      if (fabs(*((seqnodes+i)->y+l)-*(newfrontiers1+j*ntypes+l))>0.01)
						k=-1;
					      l++;
					    }
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  // Debug
					  printf("i=%d, visited[%d]=1\n",i,k);

					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;
				      				      
				    }
				}
			    }
			  */
			  			  
			  free(visited); 			  			  
			  free(newfrontiers1); free(newfrontiers2); free(newfrontiers3); 

			  // Debug
			  //printf("n=%d, candind=%d, ncandidates=%d, nnewfrontiers=%d\n",n,candind,ncandidates,nnewfrontiers);

			  
			}			    

		      /*
		      // Debug
		      if ((n==0)&&(candind==0))
			{
			  printf("ncandidates=%d, candidates: ",ncandidates);
			  for (i=0; i<=(ncandidates-1); i++)
			    printf("%d ",*(candidates+i));
			  printf("\n");
			}
		      */
      		    
		      /*
		      // Debug
		      printf("n=%d, candind=%d, ncandidates=%d, nnewcandseqs=%d\n",n,candind,ncandidates,nnewcandseqs);
		      
		      //for (i=0; i<=(ncandidates-1); i++)
		      //{
		      //  printf("%d ",*(candidates+i));
		      //}
		      //printf("\n");

		      for (i=0; i<=(ncandidates-1); i++)
			{
			  double val1, val2, tt;
			  double *dosage, *finaly;

			  // Evaluate worst and best case final population.
			  // Here the worst case is the best outcome of static therapies.
			  dosage=malloc(sizeof(double)*ndrugs); 
			  finaly=malloc(sizeof(double)*ntypes);

			  val1=1e5000;
			  for (j=0; j<=(nstates-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l);			    
			      recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  
			  //
			  //val1=0;
			  //for (j=0; j<=(ndrugs-1); j++)
			  //{
			  //  for (l=0; l<=(ndrugs-1); l++)
			  //	*(dosage+l)=0;
			  //  *(dosage+j)=1;
			  //  recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  //  val=0;
			  //    for (l=0; l<=(ntypes-1); l++)
			  //	val+=*(finaly+l);
			  //  val1=max(val1,val);
			  //}
			  //

			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val;
			  //free(dosage); free(finaly);


			  printf("node %d, ( ",*(candidates+i)); val=0;			  
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      val+=*((seqnodes+*(candidates+i))->y+j);
			      printf("%.3e ",*((seqnodes+*(candidates+i))->y+j));
			    }
			  printf(") %.3e  %.3e %.3e\n",val,val1,val2);

			  //
			  //tt=maxtime;
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0;
			  //*(dosage+0)=1;
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage1: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0;
			  //*(dosage+1)=1;
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage2: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //for (j=0; j<=(ndrugs-1); j++)
			  //*(dosage+j)=0.5;			  
			  //recurse_multidrug_response_trimrates2((seqnodes+*(candidates+i))->y,T,g0,Sg,a0,Sa,0,tt,dosage,finaly);
			  //printf("dosage12: ( "); val=0;
			  //for (j=0; j<=(ntypes-1); j++)
			  //{
			  //  printf("%.3e ",*(finaly+j)); val+=*(finaly+j);
			  //}
			  //printf(") %.3e\n",val);
			  //
			}
			  
		      //exit(0);
		      //m=nnewcandseqs;
		      */

		      /*
		      // Debug
		      if ((n==5)&&(candind==3))
			{
			  printf("ncandidates=%d: ",ncandidates);
			  for (i=0; i<=(ncandidates-1); i++)
			    printf("(%d %d)",*(candidates+i),(seqnodes+*(candidates+i))->depth);
			  printf("\n");
			}
		      */

		      // Debug
		      //t1=malloc(sizeof(time_t)); time(t1);
		      
		      
		      // Append the candidate sequences to the existing ones.
		      // Also record the val1 and val2 of each population composition.
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  nnewcandseqs++;

			  // Debug
			  //if ((n==5)&&(candind==0)&&((i==0)||(i==148)))
			  //printf("nnewcandseqs-1=%d\n",nnewcandseqs-1);
			  
			  // Write the candidates at the prior level.
			  valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
			  newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);
			  //(newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+nsteps;
			  (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+maxdepth;
			  (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
			  (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
			  for (j=0; j<=((candseqs+candind)->depth-1); j++)
			    {
			      for (k=0; k<=(ndrugs-1); k++)
				*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			    }
			  for (j=0; j<=(candseqs+candind)->depth; j++)
			    {
			      for (k=0; k<=(ntypes-1); k++)
				*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			    }

			  // Append the current candidates to the existing ones.
			  curind=*(candidates+i); l=(newcandseqs+nnewcandseqs-1)->depth;
			  while (curind>0)
			    {
			      for (j=0; j<=(ndrugs-1); j++)			    
				*((newcandseqs+nnewcandseqs-1)->dosages+(l-1)*ndrugs+j)=*((seqnodes+curind)->dosage+j);
			      for (j=0; j<=(ntypes-1); j++)
				*((newcandseqs+nnewcandseqs-1)->ys+l*ntypes+j)=*((seqnodes+curind)->y+j);

			      /*
			      // Debug
			      if ((n==5)&&(candind==0)&&(i==0))
				{
				  printf("curind=%d, dosage=( ",curind);
				  for (j=0; j<=(ndrugs-1); j++)
				    printf("%.2f ",*((seqnodes+curind)->dosage+j));
				  printf(")\n");
				}
			      */

			      /*
			      // Debug
			      if (n==0)
				{
				  printf("candind %d, candidate %d, curind=%d, y ",candind,i,curind);
				  for (j=0; j<=(ntypes-1); j++)
				    printf("%.4e ",*((seqnodes+curind)->y+j));
				  printf("\n");
				}
			      */			      			      

			      curind=(seqnodes+curind)->pa; l--;
			    }

			  gvals1=realloc(gvals1,sizeof(double)*nnewcandseqs);			  
			  *(gvals1+nnewcandseqs-1)=*(lvals1+*(candidates+i));
			  gvals2=realloc(gvals2,sizeof(double)*nnewcandseqs);
			  *(gvals2+nnewcandseqs-1)=*(lvals2+*(candidates+i));

			  /*
			  // Debug
			  //if ((n==5)&&((candind==0)||(candind==12)))
			  if (n==5)
			    {
			      printf("n=%d, candind=%d, i=%d, curind=%d, index=%d, depth=%d\n",n,candind,i,*(candidates+i),nnewcandseqs-1,(newcandseqs+nnewcandseqs-1)->depth);
			      //printf("curdosage: ");
			      //for (k=0; k<=(ndrugs-1); k++)
			      //printf("%.2f ",*((seqnodes+*(candidates+i))->dosage+k));
			      //printf("\n");
			      printf("dosage sequence ");
			      for (j=0; j<=((newcandseqs+nnewcandseqs-1)->depth-1); j++)
				{
				  printf("( ");
				  for (k=0; k<=(ndrugs-1); k++)
				    printf("%.2f ",*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k));
				  printf(") ");
				}
			      printf("\n");
			      printf("pop sequence ");
			      for (j=0; j<=((newcandseqs+nnewcandseqs-1)->depth); j++)
				{
				  printf("( ");
				  for (k=0; k<=(ntypes-1); k++)
				    printf("%.3e ",*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k));
				  printf(") ");
				}
			      printf("\n\n");
			    }
			  */
			  

			  /*
			  // Debug
			  printf("n=%d, candind=%d, i=%d, depth=%d, ( ",n,candind,i,(seqnodes+*(candidates+i))->depth);
			  val=0;
			  for (j=0; j<=(ntypes-1); j++)
			    {
			      printf("%.3e ",*((seqnodes+*(candidates+i))->y+j));
			      val+=*((seqnodes+*(candidates+i))->y+j);
			    }
			  printf("), %.3e\n",val);
			  */
			  

			}

		      // Debug
		      //printf("after appending\n");
		      

		      // Debug
		      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
		      //printf("td4=%.4f\n",val);
		      //free(t1); free(t2);


		      free(candidates);
		      free(lvals1); free(lvals2);

		      // Clear up memory for seqnodes.
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  free((seqnodes+i)->children);
			  free((seqnodes+i)->dosage);
			  free((seqnodes+i)->y);
			}
		      free(seqnodes);	     	      
		       
		    }
		  
		  /*
		  // Debug		  
		  printf("n=%d, candind=%d, depth=%d, extinction=%d, nnewcandseqs=%d\n",n,candind,(candseqs+candind)->depth,extinction,nnewcandseqs);
		  val=0; printf("( ");
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      printf("%.3e ",*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i));
		      val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		    }
		  printf(") %.3e\n",val);
		  */		
		}
	    }

	  /*
	  // Debug	  
	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  printf("n=%d, nnewcandseqs=%d, maxdepth=%d, extinction=%d\n",n,nnewcandseqs,maxdepth,extinction);	  	  
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      //if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  printf("newcand %d, valid=%d, depth=%d, ( ",i,*(valid+i),(newcandseqs+i)->depth);
		  val=0;
		  for (j=0; j<=(ntypes-1); j++)
		    {
		      printf("%.3e ",*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j));
		      val+=*((newcandseqs+i)->ys+(newcandseqs+i)->depth*ntypes+j);
		    }
		  printf(") %.3e\n",val);
		}
	    }
	  */

	  // Debug
	  //printf("n=%d, nnewcandseqs=%d\n",n,nnewcandseqs);
	
		
	  // Distinct treatment sequences may lead to identical or inferior terminal population composition.  Establish frontiers again on all treatment sequences and filter out inferior/identical ones.
	  if (nnewcandseqs>0)
	    {
	      int nnewfrontiers, nleft, *visited;
	      double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;
	
	      maxdepth=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		maxdepth=max(maxdepth,(newcandseqs+i)->depth);

	      // Debug
	      //t1=malloc(sizeof(time_t)); time(t1);

	      // Establish frontiers of terminal nodes.			  
	      nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
	      newfrontiers2=malloc(sizeof(double));
	      newfrontiers3=malloc(sizeof(double));
	      template=malloc(sizeof(double)*(ntypes+2));
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{		  		  
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      //double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      // val1 and val2 are already calculated and stored in gvals1 and gvals2.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);
		      
		      if (frontier2flag==1)
			{
			  /*
			  val1=1e5000;
			  for (j=0; j<=(ndrugs-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l);
			      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  */
			  val1=*(gvals1+i);
			}

		      /*
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      */

		      if (frontier2flag==1)
			{
			  /*
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val; val2=min(val2,val1);
			  */
			  val2=*(gvals2+i);
			}

		      //free(dosage); free(finaly);
		    

		      // Check population composition and predicted final population.
		      // If the candidate sequence lies below every new frontier, then it is a new frontier.
		      // If the best of the candidate sequence is larger than the worst of a frontier, then it is not a frontier.
		      // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  //l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			  l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))			      
				k=j;

			      /*
			      if ((k<0)&&(val2>*(newfrontiers2+j)))	
				{
				  //if (*(newfrontiers2+j)>*(newfrontiers3+j))
				  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
				    k=j;
				}
			      */
			    }					  
			  j++;
			}

		      /*
		      // Debug
		      if ((n==5)&&(i==31))
			{
			  j=0;
			  printf("flag1=%d, flag2=%d, flag3=%d\n",vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes),(val2>*(newfrontiers2+j)),*(newfrontiers2+j)>*(newfrontiers3+j));
			}
		      */

		      /*
		      // Debug
		      if ((n==5)&&(i==10))
			{
			  double *dosage, *finaly;

			  dosage=malloc(sizeof(double)*ndrugs);
			  finaly=malloc(sizeof(double)*ntypes);
			  printf("n=%d, i=%d, k=%d\n",n,i,k);
			  printf("ys: ( ");
			  for (l=0; l<=(ntypes-1); l++)
			    //printf("%.3e ",*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l));
			    printf("%.3e ",*((newcandseqs+i)->ys+maxdepth*ntypes+l));
			  printf("), val1=%.3e, val2=%.3e\n",val1,val2);
			  printf("ft1: ( ");
			  for (l=0; l<=(ntypes-1); l++)
			    printf("%.3e ",*(newfrontiers1+k*ntypes+l));
			  printf("), ft2=%.3e, ft3=%.3e\n",*(newfrontiers2+k),*(newfrontiers3+k));
			  //printf("flag1=%d, flag2=%d\n",vectorbeloweq(ntypes,newfrontiers1+k*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes),val2>*(newfrontiers2+k));
			  printf("flag1=%d, flag2=%d\n",vectorbeloweq(ntypes,newfrontiers1+k*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes),val2>*(newfrontiers2+k));
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1.0;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_1=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+0)=1;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_2=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+1)=1;
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_3=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.5;			  
			  //recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+(maxdepth-1)*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  recurse_multidrug_response_trimrates2_2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv1_4=%.3e\n",val);

			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1.0;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_1=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+0)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_2=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.0;
			  *(dosage+1)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_3=%.3e\n",val);
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0.5;
			  *(dosage+0)=1;
			  recurse_multidrug_response_trimrates2_2(newfrontiers1+k*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  printf("vv2_4=%.3e\n",val);
			  
			}
		      */
			
		      
		      /*
		      // Check population composition and predicted final population.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=0; k=j;
			  while ((l<=(ntypes-1))&&(k>=0))
			    {
			      if (*((newcandseqs+i)->ys+maxdepth*ntypes+l)<*(newfrontiers1+j*ntypes+l))			      
				k=-1;
			      l++;
			    }				      
			  //if ((k<0)&&(val2>*(newfrontiers2+j)))
			  //k=j;					  
			  j++;
			}
		      */

		      // The current candidate sequence is a new frontier.
		      if (k<0)
			{
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("before nnewfrontiers=%d\n", nnewfrontiers);


			  // Remove unqualified old frontiers.				      
			  // An unqualified frontier lies above the current candidate sequence.
			  // The best of an unqualified frontier is larger than the worst of the current node.
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      
			      //l=vectorbelow(ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes,newfrontiers1+j*ntypes);
			      l=vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes);	
			      if (l==1)
				qualified=0;
			      
			      if (frontier2flag==1)
				{
				  //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))				  
				    qualified=0;
				  /*
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))
				    {
				      // If the upper and lower bounds are identical, then don't discard the old frontiers.
				      //if (*(newfrontiers2+j)>*(newfrontiers3+j))				      
				      if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(*(newfrontiers3+j)>=10))
					qualified=0;
				    }
				  */
				}


			      // Debug
			      //if ((n==5)&&(i==29)&&(j==(nnewfrontiers-1)))
			      //printf("qualified=%d, flag1=%d, flag2=%d, flag3=%d, diff=%f\n",qualified,vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes),*(newfrontiers3+j)>val1,*(newfrontiers2+j)>*(newfrontiers3+j),*(newfrontiers2+j)-*(newfrontiers3+j));

			      
			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}					      
			    }
		       
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("middle, nleft=%d\n", nleft);
			  

			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l);
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l);
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;

			  // Debug
			  //if (n==5)
			  //{
			  //  printf("n=%d, i=%d, y=( ",n,i);
			  //  for (l=0; l<=(ntypes-1); l++)
			  //	printf("%.3e ",*((newcandseqs+i)->ys+maxdepth*ntypes+l));
			  //  printf("), val1=%.3e, val2=%.3e\n",val1,val2);
			  //}

			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }

			  /*
			  // Debug
			  if (n==5)
			    {
			      //printf("n=%d, i=%d, newfrontier[%d]=%d, maxdepth=%d: ( ",n,i,nnewfrontiers-1,i,maxdepth);
			      printf("n=%d, i=%d, nnewfrontiers=%d, maxdepth=%d\n",n,i,nnewfrontiers,maxdepth);
			      for (j=0; j<=(nnewfrontiers-1); j++)
				{
				  printf("ft[%d]: ",j);
				  for (l=0; l<=(ntypes-1); l++)
				    printf("%f ",*(newfrontiers1+j*ntypes+l));
				  printf("%f %f\n",*(newfrontiers2+j),*(newfrontiers3+j));
				}
			      //for (l=0; l<=(ntypes-1); l++)
			      //printf("%.3e ",*(template+(nnewfrontiers-1)*(ntypes+2)+l));
			      //printf(")\n");
			    }
			  */
			  
			}


		      /*
		      if (k<0)
			{	
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      l=0;
			      while ((l<=(ntypes-1))&&(qualified==1))
				{
				  if (*(newfrontiers1+j*ntypes+l)>*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l))
				    qualified=0;
				  l++;
				}
			      //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
			      //qualified=0;

			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}	
			    }
			  
			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    //*(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l); 
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l); 
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers); 			 
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }					 
			}
		      */
		    }
		}				      
	    
	      free(template);

	      // Debug
	      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
	      //printf("td5=%.4f\n",val);
	      //free(t1); free(t2);

	      /*
	      // Debug
	      printf("nnewfrontiers=%d\n",nnewfrontiers);	      	      
	      for (i=0; i<=(nnewfrontiers-1); i++)
		printf("%d %f %f\n",i,*(newfrontiers+i*2+0),*(newfrontiers+i*2+1));
	      */

	      
	      /*
	      // Debug
	      if (n==5)
		{
		  printf("nnewfrontiers=%d\n",nnewfrontiers);	      
		  for (i=0; i<=(nnewfrontiers-1); i++)
		    {
		      printf("%d ",i);
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.3e ",*(newfrontiers1+i*ntypes+j));
		      printf("%.3e %.3e\n",*(newfrontiers2+i),*(newfrontiers3+i));
		    }
		  printf("nnewcandseqs=%d, maxdepth=%d\n",nnewcandseqs,maxdepth);
		  for (i=0; i<=(nnewcandseqs-1); i++)
		    if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		      printf("%d ",i);
		  printf("\n");

		}
	      */
	      	      
	      // Debug
	      //t1=malloc(sizeof(time_t)); time(t1);


	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))		
		    {
		      double val1, val2;
		      //double *dosage, *finaly;

		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);

		      if (frontier2flag==1)
			{
			  /*
			  val1=1e5000;
			  for (j=0; j<=(ndrugs-1); j++)
			    {
			      for (l=0; l<=(ndrugs-1); l++)
				*(dosage+l)=*(dstates+j*ndrugs+l); 			  
			      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			      val=0;
			      for (l=0; l<=(ntypes-1); l++)
				val+=*(finaly+l);
			      val1=min(val1,val);
			    }
			  */
			  val1=*(gvals1+i);
			}

		      /*
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      */

		      if (frontier2flag==1)
			{
			  /*
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val2=val; val2=min(val2,val1);
			  */
			  val2=*(gvals2+i);
			}

		      //free(dosage); free(finaly);
		      

		      // Debug
		      //printf("i=%d, ",i);
		      //for (j=0; j<=(ntypes-1); j++)
		      //printf("%.3e ",*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+j));
		      //printf("%.3e %.3e\n",val1,val2);


		      // Check whether the current node lies on the new frontier.
		      // If the best of the current node is larger than the worst of a frontier, then it is not a frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  //l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			  l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				k=j;
			      /*
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				{
				  if (((*(newfrontiers2+j)-*(newfrontiers3+j))>0.1)&&((val1-val2)>0.1)&&(val2>=10))
				    k=j;
				}
			      */
			    }	
			  j++;
			}				  

		      // Debug
		      //printf("i=%d, k=%d\n",i,k);
		      
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      //l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			      l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			      if (l==1)
				k=j;					  
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;	
			    }
			}
		      
		      if (k>=0)
			*(valid+i)=0;

		      /*
		      if (k<0)
			{
			  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
			  *(candidates+ncandidates-1)=i; 
			}
		      */
		    }
		}


	      // Debug
	      //t2=malloc(sizeof(time_t)); time(t2); val=difftime(*t2,*t1);
	      //printf("td6=%.4f\n",val);
	      //free(t1); free(t2);


	      /*
	      // Debug	      
	      printf("nnewcandseqs=%d\n",nnewcandseqs);
	      k=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		if (*(valid+i)==1)
		  k++;
	      printf("nvalid=%d\n",k);
	      */
	      

	      /*
	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      dosage=malloc(sizeof(double)*ndrugs); 
		      finaly=malloc(sizeof(double)*ntypes);
		      val1=0;
		      for (j=0; j<=(ndrugs-1); j++)
			{
			  for (l=0; l<=(ndrugs-1); l++)
			    *(dosage+l)=0;
			  *(dosage+j)=1;
			  recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
			  val=0;
			  for (l=0; l<=(ntypes-1); l++)
			    val+=*(finaly+l);
			  val1=max(val1,val);
			}
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=1;
		      recurse_multidrug_response_trimrates2((newcandseqs+i)->ys+maxdepth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val2=val;
		      free(dosage); free(finaly);
		      	      

		      // Check population composition and predicted final population.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=0; k=j;
			  while ((l<=(ntypes-1))&&(k>=0))
			    {
			      if (*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l)<=*(newfrontiers1+j*ntypes+l))
				k=-1;				  
			      l++;
			    }				      
			  //if ((k<0)&&(val2>*(newfrontiers2+j)))
			  //k=j;					  
			  j++;
			}
		      
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      l=0; k=j;
			      while ((l<=(ntypes-1))&&(k>=0))
				{
				  if (fabs(*((newcandseqs+i)->ys+(maxdepth-1)*ntypes+l)-*(newfrontiers1+j*ntypes+l))>0.01)
				    k=-1;
				  l++;
				}
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;
			    }
			}

		      if (k>=0)
			*(valid+i)=0;

		    }
		}
	      */

	      free(visited); free(newfrontiers1); free(newfrontiers2); free(newfrontiers3);
	    }     
	  
		  
	 

	  /*
	  // Debug
	  printf("nnewcandseqs=%d\n",nnewcandseqs);	  	  	  	  	  
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      printf("newcandseq %d, depth %d, valid=%d\n",i,(newcandseqs+i)->depth,*(valid+i));
	      printf("dosages: ");
	      for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((newcandseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(newcandseqs+i)->depth; j++)
		{
		  printf("(%d ",j);
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((newcandseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  */


	  // Debug
	  //if (n==5)
	  //printf("nnewcandseqs=%d, valid[0]=%d, valid[148]=%d\n",nnewcandseqs,*(valid+0),*(valid+148));
	
	  // Update the candidate sequences.
	  // Remove the new candidate sequeces that are inferior to others.
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      free((candseqs+i)->dosages); free((candseqs+i)->ys);
	    }	  
	  ncandseqs=0; 

	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  ncandseqs++; candseqs=realloc(candseqs,sizeof(struct sequence)*ncandseqs);
		  (candseqs+ncandseqs-1)->depth=(newcandseqs+i)->depth;
		  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(candseqs+ncandseqs-1)->depth);
		  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*((candseqs+ncandseqs-1)->depth+1));
		  for (j=0; j<=((candseqs+ncandseqs-1)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((candseqs+ncandseqs-1)->dosages+j*ndrugs+k)=*((newcandseqs+i)->dosages+j*ndrugs+k);
		  for (j=0; j<=(candseqs+ncandseqs-1)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((candseqs+ncandseqs-1)->ys+j*ntypes+k)=*((newcandseqs+i)->ys+j*ntypes+k);

		  cvals1=realloc(cvals1,sizeof(double)*ncandseqs);
		  *(cvals1+ncandseqs-1)=*(gvals1+i);
		  cvals2=realloc(cvals2,sizeof(double)*ncandseqs);
		  *(cvals2+ncandseqs-1)=*(gvals2+i);

		}
	    }
	  
	

	  // Clear up memory for new candidates.
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      free((newcandseqs+i)->dosages);
	      free((newcandseqs+i)->ys);
	    }
	  free(newcandseqs); newcandseqs=malloc(sizeof(struct sequence));
	  free(valid); 
	
	  // Debug
	  //printf("ncandseqs=%d\n",ncandseqs);
	  
	  
	  // If ncandseqs exceeds the maximum number, then only keep maxncandseqs of them and report truncation.
	  // Sort candseqs by their final populations with the best static treatments.
	  // Sort candseqs by the geometric means of the two boundary values (best static treatment and full dosage treatments).
	  if (ncandseqs>=maxncandseqs)
	    {
	      struct pair *pairs;	     
	      double val1, val2;
	      //double *dosage, *finaly;

	      // Evalaute the best outcome of static therapies.
	      // Sort the candidates by these outcomes.
	      //dosage=malloc(sizeof(double)*ndrugs); 
	      //finaly=malloc(sizeof(double)*ntypes);
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  /*
		  val1=1e5000;
		  for (j=0; j<=(nstates-1); j++)
		    {
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=*(dstates+j*ndrugs+l);				    
		      recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val1=min(val1,val);
		    }
		  */
		  val1=*(cvals1+i); val2=*(cvals2+i);
		  (pairs+i)->i=i; 
		  //(pairs+i)->s=val1;
		  (pairs+i)->s=sqrt(val1*val2);

		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      m=maxncandseqs; maxncandseqs=m/10;

	      /*
	      // Debug
	      if (n==0)
		{
		  printf("ncandseqs=%d\n",ncandseqs);
		  for (l=0; l<=(ncandseqs-1); l++)
		    {
		      i=(pairs+l)->i;
		      printf("candseq (%d %d), depth %d\n",l,i,(candseqs+i)->depth);
		      
		      printf("dosages: ");
		      for (j=0; j<=((candseqs+i)->depth-1); j++)
			{
			  printf("( ");
			  for (k=0; k<=(ndrugs-1); k++)
			    printf("%.1f ",*((candseqs+i)->dosages+j*ndrugs+k));
			  printf(") ");
			}
		      printf("\n");		  
		      printf("ys: ");
		      for (j=0; j<=(candseqs+i)->depth; j++)
			{
			  val=0;
			  printf("(%d ",j);
			  for (k=0; k<=(ntypes-1); k++)
			    {
			      printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
			      val+=*((candseqs+i)->ys+j*ntypes+k);
			    }
			  printf(") (%.3e) ",val);
			}
		      printf("\n");
		      printf("bds: (%.3e, %.3e)\n",*(cvals1+i),*(cvals2+i));		  
		    }
		}
	      */


	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs; *truncated=1; 
	      free(pairs); //free(dosage); free(finaly);
	      maxncandseqs=m;


	      /*
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  k=rand(); (pairs+i)->i=i; (pairs+i)->s=(double)k;
		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs;
	      free(pairs); *truncated=1;	    
	      */
	    }	  	 


	  /*
	  // Debug
	  printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      printf("candseq %d, depth %d\n",i,(candseqs+i)->depth);
	      printf("dosages: ");
	      for (j=0; j<=((candseqs+i)->depth-1); j++)
		{
		  printf("( ");
		  for (k=0; k<=(ndrugs-1); k++)
		    printf("%.1f ",*((candseqs+i)->dosages+j*ndrugs+k));
		  printf(") ");
		}
	      printf("\n");
	      printf("ys: ");
	      for (j=0; j<=(candseqs+i)->depth; j++)
		{
		  printf("( ");
		  for (k=0; k<=(ntypes-1); k++)
		    printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
		  printf(") ");
		}
	      printf("\n");
	    }
	  */
	  
	  /*
	  // Debug
	  if (n<=5)
	    {
	      printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  double val1, val2;
		  double *dosage, *finaly;
		  
		  // Evaluate worst and best case final population.
		  dosage=malloc(sizeof(double)*ndrugs); 
		  finaly=malloc(sizeof(double)*ntypes);
		  val1=0;
		  for (j=0; j<=(ndrugs-1); j++)
		    {
		      for (l=0; l<=(ndrugs-1); l++)
			*(dosage+l)=0;
		      *(dosage+j)=1;
		      recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		      val=0;
		      for (l=0; l<=(ntypes-1); l++)
			val+=*(finaly+l);
		      val1=max(val1,val);
		    }
		  for (l=0; l<=(ndrugs-1); l++)
		    *(dosage+l)=1;
		  recurse_multidrug_response_trimrates2((candseqs+i)->ys+(candseqs+i)->depth*ntypes,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
		  val=0;
		  for (l=0; l<=(ntypes-1); l++)
		    val+=*(finaly+l);
		  val2=val;
		  free(dosage); free(finaly);

		  val=0;
		  printf("candseq %d, depth %d, y=( ",i,(candseqs+i)->depth);
		  for (k=0; k<=(ntypes-1); k++)
		    {
		      printf("%.3e ",*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k));
		      val+=*((candseqs+i)->ys+(candseqs+i)->depth*ntypes+k);
		    }
		  printf(") %.3e %.3e %.3e\n",val,val1,val2);
		}
	    }
	  */

	  // Debug
	  //printf("n=%d, stopT=%f, extinction=%d\n",n,stopT,extinction);

	  // Stop when some sequences yield zero population.	  
	  // Report the survival time to be max time span + timeinterval.
	  // This number distinguishes the cases where the patients are not cured but can survive for the max time span.	  
	  if (extinction==1)
	    {
	      //stopT=*(t+nintervals-1);
	      stopT=*(t+nintervals-1)+timeinterval;
	      finaltreatmenttime=*(t+maxdepth);
	    }

	  // Stop when the best candidate treatment sequences lead to mortality.
	  minval=1e5000;
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((candseqs+i)->ys+ntypes*(candseqs+i)->depth+j);
	      minval=min(minval,val);
	    }

	  
	  //if (val>=mortal)
	  if (minval>=mortal)
	    {
	      // Debug
	      //printf("minval=%.3e, maxdepth=%d, t=%.0f\n",minval,maxdepth,*(t+n));

	      if ((maxdepth>=1)&&(maxdepth<=nintervals))
		stopT=*(t+maxdepth-1);
	      else if (maxdepth>nintervals)
		stopT=*(t+nintervals-1);
	      else
		stopT=0;
	      finaltreatmenttime=stopT;
	    }	  
	  

	  /*
	  // Debug
	  //if (n==5)
	  if (n==0)
	  //if (n==3)
	    {
	      printf("n=%d, ncandseqs=%d\n",n,ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  printf("candseq %d, depth %d, dosages ",i,(candseqs+i)->depth);
		  for (j=0; j<=((candseqs+i)->depth-1); j++)
		    {
		      printf("( ");
		      for (k=0; k<=(ndrugs-1); k++)
			printf("%.2f ",*((candseqs+i)->dosages+j*ndrugs+k));
		      printf(") ");
		    }
		  printf("\n");
		  printf("pops ");
		  for (j=0; j<=((candseqs+i)->depth); j++)
		    {
		      printf("( ");
		      for (k=0; k<=(ntypes-1); k++)
			printf("%.3e ",*((candseqs+i)->ys+j*ntypes+k));
		      printf(") ");
		    }
		  printf("\n");
		}
	      exit(0);
	    }
	  */	  
	}


      // Debug
      //stopT=-1;
      
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
	  if (*(y+i*nintervals+nintervals-1)>1)
	    k=i;
	  i++;
	}
      if (k<0)
	stopT=*(t+nintervals-1)+timeinterval;
      else
	stopT=*(t+nintervals-1);
    }

  //if (stopT<0)
  //stopT=*(t+nintervals-1);
  */
  
 
  // Extract the best sequence from candidates.
  // In principle, all the final candidates are equivalent.
  // Pick up the one with minimum total population.
  n=0; minind=-1; minval=1e5000;
  for (n=0; n<=(ncandseqs-1); n++)
    {
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*((candseqs+n)->ys+ntypes*(candseqs+n)->depth+i);
      if (val<minval)
	{
	  minind=n; minval=val;
	}
    }
  
  // Write the selected treatment sequence to dosage and population to y.
  for (n=0; n<=((candseqs+minind)->depth-1); n++)    
    for (i=0; i<=(ndrugs-1); i++)	
      *(dosages+i*nintervals+n)=*((candseqs+minind)->dosages+n*ndrugs+i);
  for (n=0; n<=(candseqs+minind)->depth; n++)
    for (i=0; i<=(ntypes-1); i++)
      *(y+i*nintervals+n)=*((candseqs+minind)->ys+n*ntypes+i);
  if ((minval>=mortal)||(minval<1))
    {
      for (n=(candseqs+minind)->depth; n<=(nintervals-1); n++)
	for (i=0; i<=(ndrugs-1); i++)
	  *(dosages+i*nintervals+n)=0;
      for (n=((candseqs+minind)->depth+1); n<=(nintervals-1); n++)
	for (i=0; i<=(ntypes-1); i++)
	  *(y+i*nintervals+n)=*((candseqs+minind)->ys+(candseqs+minind)->depth*ntypes+i);
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


  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *inity;
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
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

      free(inity);
    }

  /*
  // Debug  
  for (n=0; n<=((candseqs+minind)->depth-1); n++)
    {
      printf("gb t(%d)=%.2f, dosage=( ",n,*(t+n));
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*(dosages+i*nintervals+n));
      printf("), y=( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)	
	{
	  printf("%.4e ",*(y+i*nintervals+n));
	  val+=*(y+i*nintervals+n);
	}
      printf(") %.4e\n",val);
    }
  */

  
  // Release memory.
  free(dstates); free(gvals1); free(gvals2);
  free(cvals1); free(cvals2); free(dosage); free(finaly);

  /*
  for (i=0; i<=(nseqnodes-1); i++)
    {
      free((seqnodes+i)->children);
      free((seqnodes+i)->dosage);
      free((seqnodes+i)->y);
    }
  free(seqnodes); 
  */
  free(nfrontiers); free(frontiers1); 
  free(frontiers2); free(frontiers3); 
  free(resvec); free(bases);
  for (n=0; n<=(ncandseqs-1); n++)
    {
      free((candseqs+n)->dosages);
      free((candseqs+n)->ys);
    }
  free(candseqs); free(tp); free(newcandseqs);
  return stopT;
}


// Recursively exhaust drug responses over all sequences up to depth=nsteps.
// Difference from recurse_sequence_two_drug_responses: apply an improved version of recurse_two_drug_response_trimrates.
// Difference from recurse_sequence_two_drug_responses2: apply branch and bound heuristic to trim unnecessary recursion.
// Report/update the frontier (total number, R12 number) pairs in each time step.  Ignore the pairs beyond the boundaries encircled by the frontiers.
// Difference from recurse_sequence_two_drug_responses3: multi-drug version.
// Difference from recurse_sequence_multidrug_responses: remove the frontiers which are no longer frontiers.
// Difference from recurse_sequence_multidrug_responses2: change the criteria for frontiers.
// Old criteria are total population and multiple-resistant population.
// New criteria are population composition and time to mortality.
// Time to mortality cannot be directly assessed for adaptive treatment sequences.  Thus consider the hypothetically worst and best cases.
// Worst case: stick to one drug throughout the whole cycle.
// Best case: administer max dosage for all drugs throughout the whole cycle. 
// Instead of evaluating time to death, calculate total population at maxt under both cases.
// Since population dynamics here are monotic functions, the one with smaller final population should have a longer survival time.
// If the worst case of parameter A has smaller final population than the best case of parameter B, then A is superior to B.
// frontiers1: population composition.
// frontiers2: worst case final total population.
// frontiers3: best case final total population.
// Also return val1 and val2 of each seqnode, so that there is no need to re-evaluate them.
// Store them in ->y fields after the values of all subpopulations.
// Difference from recurse_sequence_multidrug_responses3: Accelerate matrix exponentiation.
struct seqnode * recurse_sequence_multidrug_responses4(int *nseqnodes_pnt, struct seqnode *seqnodes, int nsteps, int nstates, double *dstates, double *T, double *g0, double *Sg, double *a0, double *Sa, double timeinterval, int *nfrontiers, double *frontiers1, double *frontiers2, double *frontiers3, double maxtime, long *curind_pnt)
{
  int i, j, k, l, m, n, inferior, depth, flag;
  long curind, nseqnodes, curind2;
  double val, *y0, *vec, *vec2, val1, val2;
  double *dosage, *finaly;

  nseqnodes=*nseqnodes_pnt; curind=*curind_pnt;
        
  
  // Evaluate population composition at maxtime.
  // Worst case: the highest outcomes of single drugs: val1.
  // Best case: max dosage for each drug: val2.
  // Here the worst case is actually a lower bound for the final population achievable from the current population composition.
  // Thus choose the best time-invariant therapy (monotherapy or combination) and evaluate its final population.
  dosage=malloc(sizeof(double)*ndrugs); val1=0;
  finaly=malloc(sizeof(double)*ntypes);

  
  if (frontier2flag==1)
    {
      val1=1e5000;
      for (n=0; n<=(nstates-1); n++)
	{
	  for (i=0; i<=(ndrugs-1); i++)
	    *(dosage+i)=*(dstates+n*ndrugs+i);
	  //recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
	  //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);
	  //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
	  //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);
	  recurse_multidrug_response_trimrates3((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);	  
	  val=0;
	  for (i=0; i<=(ntypes-1); i++)
	    val+=*(finaly+i);
	  val1=min(val1,val);      
	}
    }
  

  /*
  val1=0;
  for (n=0; n<=(ndrugs-1); n++)
    {
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i)=0;
      *(dosage+n)=1;
      recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(finaly+i);
      val1=max(val1,val);      
    }
  */

  if (frontier2flag==1)
    {      
      for (i=0; i<=(ndrugs-1); i++)
	*(dosage+i)=1;
      //recurse_multidrug_response_trimrates2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly);
      //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,1);      
      //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
      //recurse_multidrug_response_trimrates2_2((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);
      recurse_multidrug_response_trimrates3((seqnodes+curind)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,curevamode);
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(finaly+i);
      val2=val; val2=min(val2,val1);                 
    }

  free(dosage); free(finaly);

  
  if (frontier2flag==1)
    {
      *((seqnodes+curind)->y+ntypes)=val1;
      *((seqnodes+curind)->y+ntypes+1)=val2;
    }
      

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

  // Check whether the population composition of the current node is inferior to frontiers1.
  // A is inferior to B if all the components of A >= the corresponding components of B, and A B are not identical.
  // If A is not inferior in population composition, then check whether the final total population in the best case of A is larger than that of the worst case of B.  If so then A is inferior to B.
  depth=(seqnodes+curind)->depth; inferior=0; n=0;
  while ((n<*(nfrontiers+depth))&&(inferior==0))
    {
      // The current node is inferior if the frontier is below the current node.
      inferior=vectorbelow(ntypes,frontiers1+depth*maxnfrontiers*ntypes+n*ntypes,(seqnodes+curind)->y);
     
      // Consider the cases when the final population of the best case in the current node is worse than the frontier.
      if (frontier2flag==1)
	{
	  if ((inferior==0)&&(val2>*(frontiers2+depth*maxnfrontiers+n)))
	    inferior=1;
	}         
      n++;
    }   
  

  // Return when the current node is inferior to the frontier.
  if (inferior==1)
    return seqnodes;
  
  // Augment/update the frontier points if the current node is not inferior to any of them.
  // Also remove the frontier points inferior to the current node.
  else
    {
      int nleft;
      double *left1, *left2, *left3;      
      nleft=0; left1=malloc(sizeof(double)*maxnfrontiers*ntypes);
      left2=malloc(sizeof(double)*maxnfrontiers);
      left3=malloc(sizeof(double)*maxnfrontiers);
  
      
      // Remove the frontiers that are inferior or equal to the current node.
      for (n=0; n<*(nfrontiers+depth); n++)
	{	  
	  flag=vectorbeloweq(ntypes,(seqnodes+curind)->y,frontiers1+depth*maxnfrontiers*ntypes+n*ntypes);
	  if (frontier2flag==1)
	    {
	      //if (val2<=*(frontiers2+depth*maxnfrontiers+n))
	      if ((flag==0)&&(*(frontiers3+depth*maxnfrontiers+n)>val1))
		flag=1;
	    }

	  if (flag==0)
	    {
	      nleft++; 
	      for (i=0; i<=(ntypes-1); i++)		
		*(left1+(nleft-1)*ntypes+i)=*(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i);
	      *(left2+(nleft-1))=val1; *(left3+(nleft-1))=val2;	    	      	      
	    }
	}
    
      for (n=0; n<=(nleft-1); n++)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i)=*(left1+n*ntypes+i);
	  *(frontiers2+depth*maxnfrontiers+n)=*(left2+n);
	  *(frontiers3+depth*maxnfrontiers+n)=*(left3+n);	  
	}
      *(nfrontiers+depth)=nleft; 

      if (nleft<maxnfrontiers)
	{
	  *(nfrontiers+depth)+=1;
	  n=*(nfrontiers+depth)-1;
	  for (i=0; i<=(ntypes-1); i++)
	    *(frontiers1+depth*maxnfrontiers*ntypes+n*ntypes+i)=*((seqnodes+curind)->y+i);
	  *(frontiers2+depth*maxnfrontiers+n)=val1;
	  *(frontiers3+depth*maxnfrontiers+n)=val2;	  
	}
      free(left1); free(left2); free(left3);
         
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
      //(seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*ntypes);
      
      // Also store val1 and val2 under ->y.
      (seqnodes+nseqnodes+n)->y=malloc(sizeof(double)*(ntypes+2));
      
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
      seqnodes=recurse_sequence_multidrug_responses4(nseqnodes_pnt,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind2);
    }

  return seqnodes;
}


// Find a global optimal treatment sequence using dynamic programming.
// Traverse the decision tree of all treatment sequences and pick the best one.
// Apply branch and bound to trim inferior branches.
// Branch and bound cannot retract previously established branches.
// Stop every nsteps and keep only superior sequences.  Continue traversing.
// The ultimate criterion is the survival time.
// Difference from global_dp_optimize_multidrug_responses1: apply more rigorous criteria of population composition and worst-case final population to set frontiers.
// One criterion is the predicted final population.  Since evaluating population dynamics is time-consuming, make sure that each population composition only evaluates val1 and val2 once.
// Set candseqs to have the following forms.
// At each time step, dosage specifies the treatment dosage through the interval, population specifies the initial population composition.
// Thus the number of population compositions = the number of dosages + 1.
// Difference from global_dp_optimize_multidrug_responses2: only monotherapies are available.

double global_dp_optimize_multidrug_monotherapy_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int nsteps, double *dosages, int *truncated)
{
  int i, j, k, l, m, n, nseqnodes=0, *nfrontiers, nstates, combind, *resvec, *bases;
  long curind, ncandseqs, nnewcandseqs, minind;
  double *dstates, val, stopT, *frontiers1, *frontiers2, *frontiers3, pop, incurablepop, finaltreatmenttime=0, minval, maxtime;
  double *gvals1, *gvals2, *cvals1, *cvals2, *dosage, *finaly;  
  struct seqnode *seqnodes; 
  struct sequence *candseqs, *newcandseqs;
  time_t *tp, *t1, *t2;  
        

  *truncated=0; maxtime=nintervals*timeinterval;
  tp=malloc(sizeof(time_t)); srand((int)(time(tp)));
  resvec=malloc(sizeof(int)*ndrugs);
  bases=malloc(sizeof(int)*ndrugs);
  for (n=0; n<=(ndrugs-1); n++)
    *(bases+n)=2;

  // Change the number of possible drug combinations to monotherapies.
  nstates=ndrugs; dstates=malloc(sizeof(double)*nstates*ndrugs);
  for (i=0; i<=(nstates*ndrugs-1); i++)
    *(dstates+i)=0;
  for (i=0; i<=(ndrugs-1); i++)
    *(dstates+i*ndrugs+(ndrugs-1-i))=1;
		  
  /*
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
  */

 
  // Initialize the records.
  stopT=-1;

  for (n=0; n<=(nintervals-1); n++)
    {
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n)=0;
      for (i=0; i<=(ndrugs-1); i++)
	*(dosages+i*nintervals+n)=-1;
    }
  pop=0; incurablepop=*(x0+ntypes-1);
  for (i=0; i<=(ntypes-1); i++)
    {
      *(y+i*nintervals+0)=*(x0+i);
      pop+=*(x0+i);
    }


  n=0; nfrontiers=malloc(sizeof(int)*(nsteps+1));
  frontiers1=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers*ntypes);
  frontiers2=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  frontiers3=malloc(sizeof(double)*(nsteps+1)*maxnfrontiers);
  ncandseqs=1; candseqs=malloc(sizeof(struct sequence));
  (candseqs+ncandseqs-1)->depth=0;
  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs);
  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*2);
  for (i=0; i<=(ntypes-1); i++)
    *((candseqs+ncandseqs-1)->ys+i)=*(y+i*nintervals+0);
  nnewcandseqs=0; newcandseqs=malloc(sizeof(struct sequence));
  gvals1=malloc(sizeof(double)); gvals2=malloc(sizeof(double));
  cvals1=malloc(sizeof(double)); cvals2=malloc(sizeof(double));
  dosage=malloc(sizeof(double)*ndrugs);
  finaly=malloc(sizeof(double)*ntypes);


  // Consider all treatment sequences throughout nintervals.
  while ((n<(nintervals-1))&&(stopT<0))  
    {      
      int maxdepth, maxind, minind, candind;
      long ncandidates, *candidates;
      double maxval, minval;
      int opt=-1, extinction=0, *valid;

      // (n%nsteps)=0: incur dynamic programming.
      // For each candidate sequence, run the recursive program to generate all descendant nodes of depth nsteps.
      // Discard the nodes by branch-and-bound.
      if ((n%nsteps)==0)
	{		  
  
	  nnewcandseqs=0;
	  valid=malloc(sizeof(int));
	 

	  // For each candidate sequence, run recurse_sequence_multidrug_responses to exhaust all valid following sequences up to nsteps.
	  // If the candidate sequence already leads to mortality, then do not proceed.
	  for (candind=0; candind<=(ncandseqs-1); candind++)
	    {
	      // If the prior treatment sequences lead to zero population, then don't consider new ones.
	      if (extinction==0)
		{
		  // Do not proceed recursion if the current candidate already leads to mortality.
		  val=0;
		  for (i=0; i<=(ntypes-1); i++)
		    val+=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i);
		  
		  if (val>=mortal)
		    {		      
		      nnewcandseqs++;
		      valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
		      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);		      
		      (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth;
		      (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
		      (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
		      for (j=0; j<=((candseqs+candind)->depth-1); j++)
			{
			  for (k=0; k<=(ndrugs-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			}
		      for (j=0; j<=(candseqs+candind)->depth; j++)
			{
			  for (k=0; k<=(ntypes-1); k++)
			    *((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			}		      
		    }
				  
		  else
		    {
		      double *lvals1, *lvals2;

		      // Initialize memory for seqnodes.
		      nseqnodes=1; seqnodes=malloc(sizeof(struct seqnode));
		      (seqnodes+0)->index=0; (seqnodes+0)->depth=0;
		      (seqnodes+0)->pa=-1; (seqnodes+0)->nchildren=0;
		      (seqnodes+0)->children=malloc(sizeof(long));
		      (seqnodes+0)->dosage=malloc(sizeof(double)*ndrugs);
		      // Also store val1 and val2 under ->y.
		      (seqnodes+0)->y=malloc(sizeof(double)*(ntypes+2));
		      for (i=0; i<=(ntypes-1); i++)
			*((seqnodes+0)->y+i)=*((candseqs+candind)->ys+ntypes*(candseqs+candind)->depth+i); 

		      for (i=0; i<=nsteps; i++)
			*(nfrontiers+i)=0;		  
		      curind=0;	  	  	
		      		    		    
		      seqnodes=recurse_sequence_multidrug_responses3(&nseqnodes,seqnodes,nsteps,nstates,dstates,T,g0,Sg,a0,Sa,timeinterval,nfrontiers,frontiers1,frontiers2,frontiers3,maxtime,&curind);		      

		     
		      lvals1=malloc(sizeof(double)*nseqnodes);
		      lvals2=malloc(sizeof(double)*nseqnodes);
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  *(lvals1+i)=0; *(lvals2+i)=0;
			}


		      maxdepth=0;
		      for (i=0; i<=(nseqnodes-1); i++)
			if ((seqnodes+i)->nchildren==0)
			  maxdepth=max(maxdepth,(seqnodes+i)->depth);

		    
		      // Find the candidate terminal seqnoces that yield almost zero total population.
		      // Each subpopulation < 1, and this trend lasts for maxdepth-1 steps.		      
		      ncandidates=0; candidates=malloc(sizeof(long));
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  // When the population dynamics lasts for maxdepth steps.
			  if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
			    {	
			      int flag=1, curind=i;

			      // If maxdepth=1, then check the population composition at the current level.
			      // If maxdepth>1, then check the population composition up to (maxdepth-1) levels.
			      if (maxdepth==1)
				m=maxdepth;
			      else
				m=maxdepth-1;				  
			      l=0;
			      while ((l<m)&&(flag==1)&&(curind>=0))
				{
				  j=0;
				  while ((j<=(ntypes-1))&&(flag==1))
				    {				      
				      if (*((seqnodes+curind)->y+j)>=1)
					flag=0;
				      j++;
				    }
				  curind=(seqnodes+curind)->pa; l++;
				}

			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }

			  
			  // When the population dynamics ends early.  Add a candidate if each subpopulation < 1.
			  else if ((seqnodes+i)->nchildren==0)
			    {
			      int flag=1; 
			      j=0;
			      while ((j<=(ntypes-1))&&(flag==1))
				{				      
				  if (*((seqnodes+i)->y+j)>=1)
				    flag=0;
				  j++;
				}
			      if (flag==1)
				{
				  ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				  *(candidates+ncandidates-1)=i;				 
				}
			    }
			    			  
			}

		    
		      // If the no candidates present negligible population, then check whether static therapies can cure the patients.  If so, then label them as curable.
		      if (ncandidates<=0)
			{

			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  int flag=0;	
				  l=0; 
				  while ((l<=(nstates-1))&&(flag==0))
				    {
				      for (j=0; j<=(ndrugs-1); j++)
					*(dosage+j)=*(dstates+l*ndrugs+j);
				      recurse_multidrug_response_trimrates2_2((seqnodes+i)->y,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
				      
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
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i;				 
				    }
				}
			    }
			}				      				  
		      

		      // If the candidate sequences lead to zero population, then pick up any instance sequence and stop iteration.  Nullify the sequences derived from other candidates.
		      if (ncandidates>0)
			{
			  extinction=1;
			  for (i=0; i<=(nnewcandseqs-1); i++)
			    *(valid+i)=0;
			  maxdepth=nsteps;
			}
		      
		      // If no seqnodes yield zero total population, then keep the ones not inferior to others.
		      else
			{
			  int nnewfrontiers, nleft, *visited;
			  double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;			  
			 
			  // Establish frontiers of terminal nodes.			  
			  nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
			  newfrontiers2=malloc(sizeof(double));
			  newfrontiers3=malloc(sizeof(double));			  
			  template=malloc(sizeof(double)*ntypes);			  
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome for static therapies.
				  //dosage=malloc(sizeof(double)*ndrugs); 
				  //finaly=malloc(sizeof(double)*ntypes);
				  
				  // val1 and val2 are already evaluated and stored in *((seqnodes+i)->y+ntypes) and *((seqnodes+i)->y+ntypes+1).
				  
				  if (frontier2flag==1)
				    {
				      
				      val1=*((seqnodes+i)->y+ntypes);
				      *(lvals1+i)=val1;				      				      
				    }
				  				 
				  if (frontier2flag==1)
				    {
				     val2=*((seqnodes+i)->y+ntypes+1);
				      *(lvals2+i)=val2;
				    }

				  // Check population composition and predicted final population.
				  // If the current seqnode lies below every new frontier, then it is a new frontier.
				  // If the best of the current seqnode is larger than the worst of any new frontier, then it is not a new frontier.
				  // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
				      				      
				      if (l==1)
					k=j;
				      
				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;					 
					}					  
				      j++;
				    }
				

				  // The current seqnode is a new frontier.				  
				  if (k<0)
				    {				      
				      // Remove unqualified old frontiers.				      
				      // An unqualified frontier lies above the current node.
				      // The best of an unqualified frontier is larger than the worst of the current node.
				      template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
				      nleft=0;
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  int qualified=1;

					  l=vectorbelow(ntypes,(seqnodes+i)->y,newfrontiers1+j*ntypes);
					  if (l==1)
					    qualified=0;
					 					 

					  if (frontier2flag==1)
					    {
					      //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
					      if ((qualified==1)&&(*(newfrontiers3+j)>val1))
						qualified=0; 
					    }
					  					 
					  if (qualified==1)
					    {
					      nleft++; 
					      for (l=0; l<=(ntypes-1); l++)
						*(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
					      *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
					      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
					    }					      
					}
				      nleft++;
				      for (l=0; l<=(ntypes-1); l++)
					*(template+(nleft-1)*(ntypes+2)+l)=*((seqnodes+i)->y+l);
				      *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
				      *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
				      nnewfrontiers=nleft;
				      newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
				      newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
				      newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
				      for (j=0; j<=(nnewfrontiers-1); j++)
					{
					  for (l=0; l<=(ntypes-1); l++)
					    *(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
					  *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
					  *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
					}
				    }
				}
			    }				      
			
			  free(template);
	  		
			  // Find the candidates on the new frontiers.
			  // Each new frontier takes one candidate.
			  visited=malloc(sizeof(int)*nnewfrontiers);
			  for (i=0; i<=(nnewfrontiers-1); i++)
			    *(visited+i)=0;
			  for (i=0; i<=(nseqnodes-1); i++)
			    {
			      //if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==nsteps))
			      if (((seqnodes+i)->nchildren==0)&&((seqnodes+i)->depth==maxdepth))
				{
				  double val1, val2;
				  //double *dosage, *finaly;

				  // Evaluate worst and best case final population.
				  // Here the worst case is the best outcome of static therapies.
				  
				  if (frontier2flag==1)
				    {
				      val1=*(lvals1+i);
				    }

				  if (frontier2flag==1)
				    {
				      val2=*(lvals2+i);
				    }

				  // Check whether the current node lies on the new frontier.
				  // If the best of the current node is larger than the worst of the new frontier, then it is not a nwe frontier.
				  j=0; k=-1;
				  while ((j<=(nnewfrontiers-1))&&(k<0))
				    {
				      l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y); 
				      if (l==1)
					k=j;

				      if (frontier2flag==1)
					{
					  if ((k<0)&&(val2>*(newfrontiers2+j)))
					    k=j;

					}	
				      j++;
				    }				  

				  // Keep only one candidate for each new boundary.
				  if (k<0)
				    {
				      j=0; k=-1;
				      while ((j<=(nnewfrontiers-1))&&(k<0))
					{
					  l=vectoreq(ntypes,newfrontiers1+j*ntypes,(seqnodes+i)->y);
					  if (l==1)
					    k=j;					  
					  j++;
					}
				      if ((k>=0)&&(*(visited+k)==0))
					{
					  *(visited+k)=1; k=-1;	
					}
				    }
				 				  					    
				  if (k<0)
				    {
				      ncandidates++; candidates=realloc(candidates,sizeof(long)*ncandidates);
				      *(candidates+ncandidates-1)=i; 				      
				    }
				}
			    }

			  			  
			  free(visited); 			  			  
			  free(newfrontiers1); free(newfrontiers2); free(newfrontiers3); 

			}			    
 
		      
		      // Append the candidate sequences to the existing ones.
		      // Also record the val1 and val2 of each population composition.
		      for (i=0; i<=(ncandidates-1); i++)
			{
			  nnewcandseqs++;

			  // Write the candidates at the prior level.
			  valid=realloc(valid,sizeof(int)*nnewcandseqs); *(valid+nnewcandseqs-1)=1;
			  newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*nnewcandseqs);
			  //(newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+nsteps;
			  (newcandseqs+nnewcandseqs-1)->depth=(candseqs+candind)->depth+maxdepth;
			  (newcandseqs+nnewcandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+nnewcandseqs-1)->depth);
			  (newcandseqs+nnewcandseqs-1)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+nnewcandseqs-1)->depth+1));
			  for (j=0; j<=((candseqs+candind)->depth-1); j++)
			    {
			      for (k=0; k<=(ndrugs-1); k++)
				*((newcandseqs+nnewcandseqs-1)->dosages+j*ndrugs+k)=*((candseqs+candind)->dosages+j*ndrugs+k);
			    }
			  for (j=0; j<=(candseqs+candind)->depth; j++)
			    {
			      for (k=0; k<=(ntypes-1); k++)
				*((newcandseqs+nnewcandseqs-1)->ys+j*ntypes+k)=*((candseqs+candind)->ys+j*ntypes+k);
			    }

			  // Append the current candidates to the existing ones.
			  curind=*(candidates+i); l=(newcandseqs+nnewcandseqs-1)->depth;
			  while (curind>0)
			    {
			      for (j=0; j<=(ndrugs-1); j++)			    
				*((newcandseqs+nnewcandseqs-1)->dosages+(l-1)*ndrugs+j)=*((seqnodes+curind)->dosage+j);
			      for (j=0; j<=(ntypes-1); j++)
				*((newcandseqs+nnewcandseqs-1)->ys+l*ntypes+j)=*((seqnodes+curind)->y+j);

			      curind=(seqnodes+curind)->pa; l--;
			    }

			  gvals1=realloc(gvals1,sizeof(double)*nnewcandseqs);			  
			  *(gvals1+nnewcandseqs-1)=*(lvals1+*(candidates+i));
			  gvals2=realloc(gvals2,sizeof(double)*nnewcandseqs);
			  *(gvals2+nnewcandseqs-1)=*(lvals2+*(candidates+i));

			}

		      free(candidates);
		      free(lvals1); free(lvals2);

		      // Clear up memory for seqnodes.
		      for (i=0; i<=(nseqnodes-1); i++)
			{
			  free((seqnodes+i)->children);
			  free((seqnodes+i)->dosage);
			  free((seqnodes+i)->y);
			}
		      free(seqnodes);	     	      
		       
		    }
			
		}
	    }

	
		
	  // Distinct treatment sequences may lead to identical or inferior terminal population composition.  Establish frontiers again on all treatment sequences and filter out inferior/identical ones.
	  if (nnewcandseqs>0)
	    {
	      int nnewfrontiers, nleft, *visited;
	      double *newfrontiers1, *newfrontiers2, *newfrontiers3, *template;
	
	      maxdepth=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		maxdepth=max(maxdepth,(newcandseqs+i)->depth);

	      // Establish frontiers of terminal nodes.			  
	      nnewfrontiers=0; newfrontiers1=malloc(sizeof(double)*ntypes);
	      newfrontiers2=malloc(sizeof(double));
	      newfrontiers3=malloc(sizeof(double));
	      template=malloc(sizeof(double)*(ntypes+2));
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{		  		  
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		    {
		      double val1, val2;
		      //double *dosage, *finaly;
		      
		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      // val1 and val2 are already calculated and stored in gvals1 and gvals2.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);
		      
		      if (frontier2flag==1)
			{
			 val1=*(gvals1+i);
			}


		      if (frontier2flag==1)
			{
			 val2=*(gvals2+i);
			}

		      //free(dosage); free(finaly);
		    

		      // Check population composition and predicted final population.
		      // If the candidate sequence lies below every new frontier, then it is a new frontier.
		      // If the best of the candidate sequence is larger than the worst of a frontier, then it is not a frontier.
		      // However, if the upper and lower bounds of frontiers are identical and small, and the difference between val2 and the frontier lower bound is small, then add it as a new frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  l=vectorbeloweq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))			      
				k=j;
			    }					  
			  j++;
			}
		    
		      // The current candidate sequence is a new frontier.
		      if (k<0)
			{
			  // Debug
			  //if ((n==5)&&(i==31))
			  //printf("before nnewfrontiers=%d\n", nnewfrontiers);


			  // Remove unqualified old frontiers.				      
			  // An unqualified frontier lies above the current candidate sequence.
			  // The best of an unqualified frontier is larger than the worst of the current node.
			  template=realloc(template,sizeof(double)*(nnewfrontiers+1)*(ntypes+2));
			  nleft=0;
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      int qualified=1;
			      
			      l=vectorbelow(ntypes,(newcandseqs+i)->ys+maxdepth*ntypes,newfrontiers1+j*ntypes);	
			      if (l==1)
				qualified=0;
			      
			      if (frontier2flag==1)
				{
				  //if ((qualified==1)&&(*(newfrontiers2+j)>val2))
				  if ((qualified==1)&&(*(newfrontiers3+j)>val1))				  
				    qualified=0;				  
				}

			      
			      if (qualified==1)
				{
				  nleft++; 
				  for (l=0; l<=(ntypes-1); l++)
				    *(template+(nleft-1)*(ntypes+2)+l)=*(newfrontiers1+j*ntypes+l);
				  *(template+(nleft-1)*(ntypes+2)+ntypes)=*(newfrontiers2+j);
				  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=*(newfrontiers3+j);
				}					      
			    }		       

			  nleft++;
			  for (l=0; l<=(ntypes-1); l++)
			    *(template+(nleft-1)*(ntypes+2)+l)=*((newcandseqs+i)->ys+maxdepth*ntypes+l);
			  *(template+(nleft-1)*(ntypes+2)+ntypes)=val1;
			  *(template+(nleft-1)*(ntypes+2)+ntypes+1)=val2;
			 
			  nnewfrontiers=nleft;
			  newfrontiers1=realloc(newfrontiers1,sizeof(double)*nnewfrontiers*ntypes);
			  newfrontiers2=realloc(newfrontiers2,sizeof(double)*nnewfrontiers);
			  newfrontiers3=realloc(newfrontiers3,sizeof(double)*nnewfrontiers);
			  for (j=0; j<=(nnewfrontiers-1); j++)
			    {
			      for (l=0; l<=(ntypes-1); l++)
				*(newfrontiers1+j*ntypes+l)=*(template+j*(ntypes+2)+l);
			      *(newfrontiers2+j)=*(template+j*(ntypes+2)+ntypes);
			      *(newfrontiers3+j)=*(template+j*(ntypes+2)+ntypes+1);
			    }
			}		     
		    }
		}				      
	    
	      free(template);
	    

	      // Find the candidates on the new frontiers.
	      // Each new frontier takes one candidate.
	      visited=malloc(sizeof(int)*nnewfrontiers);
	      for (i=0; i<=(nnewfrontiers-1); i++)
		*(visited+i)=0;
	      for (i=0; i<=(nnewcandseqs-1); i++)
		{
		  if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))		
		    {
		      double val1, val2;
		      //double *dosage, *finaly;

		      // Evaluate worst and best case final population.
		      // Here the worst case is the best outcome for static therapies.
		      //dosage=malloc(sizeof(double)*ndrugs); 
		      //finaly=malloc(sizeof(double)*ntypes);

		      if (frontier2flag==1)
			{
			  val1=*(gvals1+i);
			}

		      if (frontier2flag==1)
			{
			  val2=*(gvals2+i);
			}


		      // Check whether the current node lies on the new frontier.
		      // If the best of the current node is larger than the worst of a frontier, then it is not a frontier.
		      j=0; k=-1;
		      while ((j<=(nnewfrontiers-1))&&(k<0))
			{
			  //l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+(maxdepth-1)*ntypes);
			  l=vectorbelow(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			  if (l==1)
			    k=j;
			  if (frontier2flag==1)
			    {
			      if ((k<0)&&(val2>*(newfrontiers2+j)))
				k=j;			      
			    }	
			  j++;
			}				  
    
		      // Keep only one candidate for each new boundary.
		      if (k<0)
			{
			  j=0; k=-1;
			  while ((j<=(nnewfrontiers-1))&&(k<0))
			    {
			      l=vectoreq(ntypes,newfrontiers1+j*ntypes,(newcandseqs+i)->ys+maxdepth*ntypes);
			      if (l==1)
				k=j;					  
			      j++;
			    }
			  if ((k>=0)&&(*(visited+k)==0))
			    {
			      *(visited+k)=1; k=-1;	
			    }
			}
		      
		      if (k>=0)
			*(valid+i)=0;		     
		    }
		}



	      free(visited); free(newfrontiers1); free(newfrontiers2); free(newfrontiers3);
	    }     
	  
	
	  // Update the candidate sequences.
	  // Remove the new candidate sequeces that are inferior to others.
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      free((candseqs+i)->dosages); free((candseqs+i)->ys);
	    }	  
	  ncandseqs=0; 

	  maxdepth=0;
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    maxdepth=max(maxdepth,(newcandseqs+i)->depth);
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      if ((*(valid+i)==1)&&((newcandseqs+i)->depth==maxdepth))
		{
		  ncandseqs++; candseqs=realloc(candseqs,sizeof(struct sequence)*ncandseqs);
		  (candseqs+ncandseqs-1)->depth=(newcandseqs+i)->depth;
		  (candseqs+ncandseqs-1)->dosages=malloc(sizeof(double)*ndrugs*(candseqs+ncandseqs-1)->depth);
		  (candseqs+ncandseqs-1)->ys=malloc(sizeof(double)*ntypes*((candseqs+ncandseqs-1)->depth+1));
		  for (j=0; j<=((candseqs+ncandseqs-1)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((candseqs+ncandseqs-1)->dosages+j*ndrugs+k)=*((newcandseqs+i)->dosages+j*ndrugs+k);
		  for (j=0; j<=(candseqs+ncandseqs-1)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((candseqs+ncandseqs-1)->ys+j*ntypes+k)=*((newcandseqs+i)->ys+j*ntypes+k);

		  cvals1=realloc(cvals1,sizeof(double)*ncandseqs);
		  *(cvals1+ncandseqs-1)=*(gvals1+i);
		  cvals2=realloc(cvals2,sizeof(double)*ncandseqs);
		  *(cvals2+ncandseqs-1)=*(gvals2+i);

		}
	    }
	  
	

	  // Clear up memory for new candidates.
	  for (i=0; i<=(nnewcandseqs-1); i++)
	    {
	      free((newcandseqs+i)->dosages);
	      free((newcandseqs+i)->ys);
	    }
	  free(newcandseqs); newcandseqs=malloc(sizeof(struct sequence));
	  free(valid); 
	  
	  
	  // If ncandseqs exceeds the maximum number, then only keep maxncandseqs of them and report truncation.
	  // Sort candseqs by their final populations with the best static treatments.
	  // Sort candseqs by the geometric means of the two boundary values (best static treatment and full dosage treatments).
	  if (ncandseqs>=maxncandseqs)
	    {
	      struct pair *pairs;	     
	      double val1, val2;
	      //double *dosage, *finaly;

	      // Evalaute the best outcome of static therapies.
	      // Sort the candidates by these outcomes.
	      //dosage=malloc(sizeof(double)*ndrugs); 
	      //finaly=malloc(sizeof(double)*ntypes);
	      pairs=malloc(sizeof(struct pair)*ncandseqs);
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  val1=*(cvals1+i); val2=*(cvals2+i);
		  (pairs+i)->i=i; 
		  (pairs+i)->s=sqrt(val1*val2);

		}
	      qsort(pairs,ncandseqs,sizeof(struct pair),(* pair_cmp_ascend));	      
	      m=maxncandseqs; maxncandseqs=m/10;
	      
	      newcandseqs=realloc(newcandseqs,sizeof(struct sequence)*maxncandseqs);	     
	      for (i=0; i<=(maxncandseqs-1); i++)
		{
		  l=(pairs+i)->i;		  
		  (newcandseqs+i)->depth=(candseqs+l)->depth;
		  (newcandseqs+i)->dosages=malloc(sizeof(double)*ndrugs*(newcandseqs+i)->depth);
		  (newcandseqs+i)->ys=malloc(sizeof(double)*ntypes*((newcandseqs+i)->depth+1));
		  for (j=0; j<=((newcandseqs+i)->depth-1); j++)
		    for (k=0; k<=(ndrugs-1); k++)
		      *((newcandseqs+i)->dosages+j*ndrugs+k)=*((candseqs+l)->dosages+j*ndrugs+k);		  
		  for (j=0; j<=(newcandseqs+i)->depth; j++)
		    for (k=0; k<=(ntypes-1); k++)
		      *((newcandseqs+i)->ys+j*ntypes+k)=*((candseqs+l)->ys+j*ntypes+k);
		}
	      for (i=0; i<=(ncandseqs-1); i++)
		{
		  free((candseqs+i)->dosages); free((candseqs+i)->ys);
		}
	      free(candseqs); candseqs=newcandseqs; newcandseqs=malloc(sizeof(struct sequence));	      
	      ncandseqs=maxncandseqs; *truncated=1; 
	      free(pairs); 
	      maxncandseqs=m;
	    }	  	 


	  // Stop when some sequences yield zero population.	  
	  // Report the survival time to be max time span + timeinterval.
	  // This number distinguishes the cases where the patients are not cured but can survive for the max time span.	  
	  if (extinction==1)
	    {
	      stopT=*(t+nintervals-1)+timeinterval;
	      finaltreatmenttime=*(t+maxdepth);
	    }

	  // Stop when the best candidate treatment sequences lead to mortality.
	  minval=1e5000;
	  for (i=0; i<=(ncandseqs-1); i++)
	    {
	      val=0;
	      for (j=0; j<=(ntypes-1); j++)
		val+=*((candseqs+i)->ys+ntypes*(candseqs+i)->depth+j);
	      minval=min(minval,val);
	    }

	  
	  if (minval>=mortal)
	    {
	      if ((maxdepth>=1)&&(maxdepth<=nintervals))
		stopT=*(t+maxdepth-1);
	      else if (maxdepth>nintervals)
		stopT=*(t+nintervals-1);
	      else
		stopT=0;
	      finaltreatmenttime=stopT;
	    }	   
	}

      n++;
    }

  // Extract the best sequence from candidates.
  // In principle, all the final candidates are equivalent.
  // Pick up the one with minimum total population.
  n=0; minind=-1; minval=1e5000;
  for (n=0; n<=(ncandseqs-1); n++)
    {
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*((candseqs+n)->ys+ntypes*(candseqs+n)->depth+i);
      if (val<minval)
	{
	  minind=n; minval=val;
	}
    }
  
  // Write the selected treatment sequence to dosage and population to y.
  for (n=0; n<=((candseqs+minind)->depth-1); n++)    
    for (i=0; i<=(ndrugs-1); i++)	
      *(dosages+i*nintervals+n)=*((candseqs+minind)->dosages+n*ndrugs+i);
  for (n=0; n<=(candseqs+minind)->depth; n++)
    for (i=0; i<=(ntypes-1); i++)
      *(y+i*nintervals+n)=*((candseqs+minind)->ys+n*ntypes+i);
  if ((minval>=mortal)||(minval<1))
    {
      for (n=(candseqs+minind)->depth; n<=(nintervals-1); n++)
	for (i=0; i<=(ndrugs-1); i++)
	  *(dosages+i*nintervals+n)=0;
      for (n=((candseqs+minind)->depth+1); n<=(nintervals-1); n++)
	for (i=0; i<=(ntypes-1); i++)
	  *(y+i*nintervals+n)=*((candseqs+minind)->ys+(candseqs+minind)->depth*ntypes+i);
    }


  // If stop<0, then decide whether the patient recovers or not.
  // Check if the patient is curable with a static therapy.
  // If the patient recovers, then report survival time to max time span + timeinterval.
  if (stopT<0)
    {
      int flag=0;
      double *inity;
      inity=malloc(sizeof(double)*ntypes);
      for (i=0; i<=(ntypes-1); i++)
	*(inity+i)=*(y+i*nintervals+nintervals-2);
      l=0;
      while ((l<=(nstates-1))&&(flag==0))
	{
	  for (j=0; j<=(ndrugs-1); j++)
	    *(dosage+j)=*(dstates+l*ndrugs+j);
	  recurse_multidrug_response_trimrates2_2(inity,T,g0,Sg,a0,Sa,0,maxtime,dosage,finaly,0);
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

      free(inity);
    }

  /*
  // Debug  
  for (n=0; n<=((candseqs+minind)->depth-1); n++)
    {
      printf("gb t(%d)=%.2f, dosage=( ",n,*(t+n));
      for (i=0; i<=(ndrugs-1); i++)
	printf("%.2f ",*(dosages+i*nintervals+n));
      printf("), y=( ");
      val=0;
      for (i=0; i<=(ntypes-1); i++)	
	{
	  printf("%.4e ",*(y+i*nintervals+n));
	  val+=*(y+i*nintervals+n);
	}
      printf(") %.4e\n",val);
    }
  */

  
  // Release memory.
  free(dstates); free(gvals1); free(gvals2);
  free(cvals1); free(cvals2); free(dosage); free(finaly);  
  free(nfrontiers); free(frontiers1); 
  free(frontiers2); free(frontiers3); 
  free(resvec); free(bases);
  for (n=0; n<=(ncandseqs-1); n++)
    {
      free((candseqs+n)->dosages);
      free((candseqs+n)->ys);
    }
  free(candseqs); free(tp); free(newcandseqs);
  return stopT;
}
