// Utility functions involved in drug responses and optimization.
// Implement five strategies of drug administration.  For each strategy update and fix dosage in each time interval.
// Strategy 1: Select the dosage combination to minimize total population.
// Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Consider the model of truncating growth rate effects from fractional cell numbers.


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# define max(i,j) (i>j ? i : j)
# define min(i,j) (i<j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)


extern int ndrugs, ntypes, initialtime;
extern double detected, mortal;

extern void padm(int ndim, double *A, double t, int p, double *E);
extern char * getitemval3(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
extern void diagonal(int ndim, double *vec, double *mat);
extern void recurse_two_drug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
extern void recurse_two_drug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);

double optimize_two_drug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, double popthre, double *dosage);
double drug_response_trimrate(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, int nintervals, double *t, double *dosage, double *y);
double optimize_two_drug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, double popthre, double *dosage);
double optimize_two_drug_responses5(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double *dosage);
double optimize_two_drug_responses6(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double *dosage);
double optimize_two_drug_responses7(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);
double optimize_two_drug_responses8(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage);


// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
double optimize_two_drug_responses3(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, double popthre, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, minval, minval2, *mask;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF;
  double dd=0.05, d1, d2, w, g, stopT, val, pop, *drugrecords;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); 
  //mask=malloc(sizeof(double)*ntypes*ntypes);

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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*4*k*k);  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);

      /*
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
      */

      // Change the drug combinations to optimize certain criterion.
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100; cnt=0;
      for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
      //for (d1=0; d1<=(1.0); d1+=dd)   
	{
	  for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
	  //for (d2=0; d2<=(1-d1); d2+=dd)
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
	      //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	      
	      // Debug
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
		  *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;

	      /*
	      // Adjust tempF by the mask matrix.
	      for (i=0; i<=(ntypes-1); i++)		
		for (j=0; j<=(ntypes-1); j++)
		  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
	      */

	      padm(ntypes,tempF,1,7,template);	
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);


	      /*
	      // Debug
	      //if ((d1<=0)&&(d2>=1.0))
		{
		  printf("dosage: (%.2f %.2f)\n",d1,d2);
		  printf("rate:\n");
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.3e ",*(template+i*ntypes+j));
		      printf("\n");
		    }		  
		  printf("population: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(vec2+i));
		  printf("\n");
		}
	      */    


	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      w=*(vec2+ntypes-1);

	      
	      cnt++; *(drugrecords+(cnt-1)*4+0)=d1; *(drugrecords+(cnt-1)*4+1)=d2;
	      *(drugrecords+(cnt-1)*4+2)=g; *(drugrecords+(cnt-1)*4+3)=w;
	      

	      // Strategy 1: Select the dosage combination to minimize total population.
	      if (mode==1)	      
		{
		  if ((minval-g)>1e-10)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==2)
		{
		  /*
		  // Debug
		  printf("d1=%.2f d2=%.2f population=",d1,d2);
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(vec2+i));
		  printf(", total pop=%.3e, mode1=%d, mode2=%d\n",pop,((pop<=popthre)&&((minval-w)>1e-10)),((pop>popthre)&&((minval-g)>1e-10)));
		  */

		  // If the total population is below the threshold, then optimize R12 population.
		  //if ((pop<=popthre)&&((minval-w)>1e-10))
		  if ((pop<=popthre)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;		      
		    }
		  // If the total population is above the threshold, then optimize the total population.
		  if ((pop>popthre)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }		  
		}

	      // Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==3)
		{		  
		  // Calculate the growth rate matrix.
		  /*
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
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		  */
		  

		  w=*(vec2+ntypes-1);
	      
		  w=0;
		  for (i=0; i<=(ntypes-1); i++)
		    w=max(w,*(template+i*ntypes+i));	      

		  /*
		  // Debug		  
		  if ((d1<=0)&&(d2>=0.95))
		    {
		      printf("rate: ");
		      for (i=0; i<=(ntypes-1); i++)
			{
			  for (j=0; j<=(ntypes-1); j++)
			    printf("%.3e ",*(template+i*ntypes+j));
			  printf("\n");
			}
		      printf("\n");
		      printf("%.2f %.2f %.3e %.3e %d\n",d1,d2,w,minval,(minval-w)>1e-10);
		    }
		  */

		  if ((minval-w)>1e-10)
		    {
		      // Debug
		      //printf("%.2f %.2f %.3e %.3e\n",d1,d2,w,minval);

		      *(mininds+0)=d1; *(mininds+1)=d2; minval=w; 		      
		    }
		}		  
	      
	      // Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
	      if (mode==4)
		{
		  // Minimize R12 population.
		  //if ((n<=initialtime)&&((minval-w)>1e-10))
		  if ((n<=initialtime)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;
		    }
		  // Minimize the total population.
		  else if ((n>initialtime)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12.
	      else if (mode==5)
		{	
		  if ((minval2-w)>1e-10)
		    minval2=w;
	  
		  // If either drug is effective on R12, then the same as strategy 1.
		  //if ((*(Sa+6)>1e-10)||(*(Sa+7)>1e-10))
		    {
		      if ((minval-g)>1e-10)
			{
			  *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			}
		    }

		  /*
		  // Otherwise check whether R12<1.
		  else
		    {
		      // Debug
		      //if (n==0)
		      //printf("d1=%.2f, d2=%.2f, g=%.3e, w=%.3e\n",d1,d2,g,w);

		      // If so then minimize the total population.
		      //if (w<(1-1e-10))
		      if (w<0.9)
			{
			  if ((minval-g)>1e-10)
			    {
			      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			    }
			}
		    }
		  */
		}  
	    }
	}     
      
      
      // Strategy 5: If neither drug is effective on R12, then do the following.
      // If the minimum of R12 population >=1, then choose the strategy to minimize R12.
      // Otherwise stick to strategy 1.
      if ((mode==5)&&(*(Sa+6)<=1e-10)&&(*(Sa+7)<=1e-10))
	{
	  // Debug
	  //printf("minval2=%.3e\n",minval2);


	  if (minval2>=(1-1e-12))
	    {
	      i=0; k=1;
	      while ((i<=(cnt-1))&&(k==1))
		{
		  if (*(drugrecords+i*4+3)<=minval2)
		    {
		      *(mininds+0)=*(drugrecords+i*4+0);
		      *(mininds+1)=*(drugrecords+i*4+1);
		      k=0;
		    }
		  i++;
		}
	    }
	}
      
            

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
           
      /*
      // Debug
      printf("n=%d, dosage=(%.2f %.2f), population=",n,*(mininds+0),*(mininds+1));      
      for (i=0; i<=(ntypes-1); i++)   	              
	printf("%.3e ",*(vec2+i));
      printf("\n");
      */

      // Debug
      //printf("%d %.1f %.1f\n",n,*(dosage+n),*(dosage+nintervals+n));


      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);

      // Debug
      //printf("here: %d %.2f %.2f %.3e %.3e %.3e %.3e %.3e\n",n,*(dosage+n),1-*(dosage+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val);


      if (val>=mortal)
	stopT=*(t+n);
    

      /*
      // Debug
      printf("t=%f, dosage: %f %f, ",*(t+n)/30,*(dosage+n),*(dosage+nintervals+n));
      printf("y: ");      
      for (i=0; i<=(ntypes-1); i++)
	{
	  printf("%.3e ",*(y+i*nintervals+n+1));	  
	}
      printf("%.3e ",val);
      printf("\n");
      */
   
      n++;
    }

  free(drugrecords); //free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
  // Debug
  //printf("stopT=%f\n",stopT);
  

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}



// Given the dynamic profile of drug dosage, calculate the time response of each population and the entire population.
// Use the revised model of truncating the growth rate when population < 1.
// Also return the stop time.
// Break the drug administration into intervals.
double drug_response_trimrate(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, int nintervals, double *t, double *dosage, double *y)
{
  int i, j, k, l, m, n;
  double *localdosage, *y0, stopT, *localresponse, val, dt;
      
  // For each drug administration interval (e.g., one month), calculate the response.
  y0=malloc(sizeof(double)*ntypes); n=0;
  localdosage=malloc(sizeof(double)*ndrugs);
  localresponse=malloc(sizeof(double)*ntypes);
  dt=*(t+1)-*t; stopT=-1;            


  while ((n<(nintervals-1))&&(stopT<0))
    {      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}
      
      
      // Get the local dosage.
      for (i=0; i<=(ndrugs-1); i++)
	*(localdosage+i)=*(dosage+i*nintervals+n);
            
       
      // Run the recursive function of predicting the responses.
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,dt,localdosage,localresponse);
     
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(localresponse+i);
      
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
      
      if (val>=mortal)	
	stopT=*(t+n);	  	
      
      n++;
    }

  if (stopT<0)
    stopT=*(t+n);


  free(localdosage); free(localresponse); free(y0);

  return stopT;
}



// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
// Difference from optimize_two_drug_responses3: apply the mask to F when optimizing drug combinations.
// Difference from optimize_two_drug_responses3: apply the recursive function to predict the responses.
double optimize_two_drug_responses4(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, double popthre, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, minval, minval2, minval2_2, *mask;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, pop, *drugrecords;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*4*k*k);  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);

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
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100; minval2_2=1e100; cnt=0;
      for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
      //for (d1=0; d1<=(1.0); d1+=dd)   
	{
	  for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
	  //for (d2=0; d2<=(1-d1); d2+=dd)
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
	      //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	      
	      /*
	      // Debug
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
		  *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;

	      //
	      // Adjust tempF by the mask matrix.
	      //for (i=0; i<=(ntypes-1); i++)		
	      //for (j=0; j<=(ntypes-1); j++)
	      //  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
	      //

	      padm(ntypes,tempF,1,7,template);	
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
	      */

	      /*
	      // Debug
	      //if ((d1<=0)&&(d2>=1.0))
		{
		  printf("dosage: (%.2f %.2f)\n",d1,d2);
		  printf("rate:\n");
		  for (i=0; i<=(ntypes-1); i++)
		    {
		      for (j=0; j<=(ntypes-1); j++)
			printf("%.3e ",*(template+i*ntypes+j));
		      printf("\n");
		    }		  
		  printf("population: ");
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(vec2+i));
		  printf("\n");
		}
	      */    


	      // Debug
	      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);

	      /*
	      // Debug
	      if (((fabs(d1-0)<=1e-10)&&(fabs(d2-1)<=1e-10))||((fabs(d1-0.6)<=1e-10)&&(fabs(d2-0.4)<=1e-10)))
		{
		  printf("d1=%.2f, d2=%.2f, vec=",d1,d2);
		  for (i=0; i<=3; i++)
		    printf("%.3e ",*(vec2+i));
		  printf("\n");
		}
	      */
	      	      

	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      w=*(vec2+ntypes-1);
	      	      
	      
	      cnt++; *(drugrecords+(cnt-1)*4+0)=d1; *(drugrecords+(cnt-1)*4+1)=d2;
	      *(drugrecords+(cnt-1)*4+2)=g; *(drugrecords+(cnt-1)*4+3)=w;
	      

	      // Strategy 1: Select the dosage combination to minimize total population.
	      if (mode==1)	      
		{
		  if ((minval-g)>1e-10)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==2)
		{
		  /*
		  // Debug
		  printf("d1=%.2f d2=%.2f population=",d1,d2);
		  for (i=0; i<=(ntypes-1); i++)
		    printf("%.3e ",*(vec2+i));
		  printf(", total pop=%.3e, mode1=%d, mode2=%d\n",pop,((pop<=popthre)&&((minval-w)>1e-10)),((pop>popthre)&&((minval-g)>1e-10)));
		  */

		  // If the total population is below the threshold, then optimize R12 population.
		  //if ((pop<=popthre)&&((minval-w)>1e-10))
		  if ((pop<=popthre)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;		      
		    }
		  // If the total population is above the threshold, then optimize the total population.
		  if ((pop>popthre)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }		  
		}

	      // Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==3)
		{		  
		  // Calculate the growth rate matrix.
		  /*
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
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		  */
		  

		  w=*(vec2+ntypes-1);
	      
		  w=0;
		  for (i=0; i<=(ntypes-1); i++)
		    w=max(w,*(template+i*ntypes+i));	      

		  /*
		  // Debug		  
		  if ((d1<=0)&&(d2>=0.95))
		    {
		      printf("rate: ");
		      for (i=0; i<=(ntypes-1); i++)
			{
			  for (j=0; j<=(ntypes-1); j++)
			    printf("%.3e ",*(template+i*ntypes+j));
			  printf("\n");
			}
		      printf("\n");
		      printf("%.2f %.2f %.3e %.3e %d\n",d1,d2,w,minval,(minval-w)>1e-10);
		    }
		  */

		  if ((minval-w)>1e-10)
		    {
		      // Debug
		      //printf("%.2f %.2f %.3e %.3e\n",d1,d2,w,minval);

		      *(mininds+0)=d1; *(mininds+1)=d2; minval=w; 		      
		    }
		}		  
	      
	      // Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
	      if (mode==4)
		{
		  // Minimize R12 population.
		  //if ((n<=initialtime)&&((minval-w)>1e-10))
		  if ((n<=initialtime)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;
		    }
		  // Minimize the total population.
		  else if ((n>initialtime)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12.
	      else if (mode==5)
		{	
		  if ((minval2-w)>1e-10)
		    minval2=w;
	  
		  // If either drug is effective on R12, then the same as strategy 1.
		  //if ((*(Sa+6)>1e-10)||(*(Sa+7)>1e-10))
		    {
		      if ((minval-g)>1e-10)
			{
			  *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			  minval2_2=w;
			}
		    }

		  /*
		  // Otherwise check whether R12<1.
		  else
		    {
		      // Debug
		      //if (n==0)
		      //printf("d1=%.2f, d2=%.2f, g=%.3e, w=%.3e\n",d1,d2,g,w);

		      // If so then minimize the total population.
		      //if (w<(1-1e-10))
		      if (w<0.9)
			{
			  if ((minval-g)>1e-10)
			    {
			      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			    }
			}
		    }
		  */
		}  
	    }
	}     
      
      
      // Strategy 5: If neither drug is effective on R12, then do the following.
      // If the minimum of R12 population >=1, then choose the strategy to minimize R12.
      // Otherwise stick to strategy 1.
      if ((mode==5)&&(*(Sa+6)<=1e-10)&&(*(Sa+7)<=1e-10))
	{
	  // Debug
	  //printf("minval2=%.3e\n",minval2);


	  //if (minval2>=(1-1e-12))
	  if (minval2_2>=(1-1e-12))
	    {
	      i=0; k=1;
	      while ((i<=(cnt-1))&&(k==1))
		{
		  if (*(drugrecords+i*4+3)<=minval2)
		    {
		      *(mininds+0)=*(drugrecords+i*4+0);
		      *(mininds+1)=*(drugrecords+i*4+1);
		      k=0;
		    }
		  i++;
		}
	    }
	}
      
            

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
           
      
      // Debug
      printf("n=%d, dosage=(%.2f %.2f), population=",n,*(mininds+0),*(mininds+1));      
      for (i=0; i<=(ntypes-1); i++)   	              
	printf("%.3e ",*(vec2+i));
      printf("\n");
      

      // Debug
      //printf("%d %.1f %.1f\n",n,*(dosage+n),*(dosage+nintervals+n));


      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);

      // Debug
      //printf("here: %d %.2f %.2f %.3e %.3e %.3e %.3e %.3e\n",n,*(dosage+n),1-*(dosage+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val);


      if (val>=mortal)
	stopT=*(t+n);
    

      /*
      // Debug
      printf("t=%f, dosage: %f %f, ",*(t+n)/30,*(dosage+n),*(dosage+nintervals+n));
      printf("y: ");      
      for (i=0; i<=(ntypes-1); i++)
	{
	  printf("%.3e ",*(y+i*nintervals+n+1));	  
	}
      printf("%.3e ",val);
      printf("\n");
      */
   
      n++;
    }

  free(drugrecords); free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
  // Debug
  //printf("stopT=%f\n",stopT);
  

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}



// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
// Difference from optimize_two_drug_responses3: apply the mask to F when optimizing drug combinations.
// Difference from optimize_two_drug_responses3: apply the recursive function to predict the responses.
// Difference from optimize_two_drug_responses4: can choose either to use F to estimate the population or incur the recursive function.  The recursive function is time consuming and leads to early drop on cluster jobs.  So do not incur the recursive functions if the stop time >= upper bound using the simple method.
double optimize_two_drug_responses5(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, minval, minval2, minval2_2, *mask;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, pop, *drugrecords;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*4*k*k);  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);

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
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100; minval2_2=1e100; cnt=0;
      for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
      //for (d1=0; d1<=(1.0); d1+=dd)   
	{
	  for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
	  //for (d2=0; d2<=(1-d1); d2+=dd)
	    {
	      *(vec+0)=d1; *(vec+1)=d2;

	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(y0+i);

	      //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	      
	      
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
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
		  
		  //
		  // Adjust tempF by the mask matrix.
		  //for (i=0; i<=(ntypes-1); i++)		
		  //for (j=0; j<=(ntypes-1); j++)
		  //  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
		  //
		  
		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		}


	      // compmode=1: incur the recursive function to estimate the population.
	      else if (compmode==1)
		recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	     

	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      w=*(vec2+ntypes-1);
	      	      
	      
	      cnt++; *(drugrecords+(cnt-1)*4+0)=d1; *(drugrecords+(cnt-1)*4+1)=d2;
	      *(drugrecords+(cnt-1)*4+2)=g; *(drugrecords+(cnt-1)*4+3)=w;
	      

	      // Strategy 1: Select the dosage combination to minimize total population.
	      if (mode==1)	      
		{
		  if ((minval-g)>1e-10)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==2)
		{		  
		  // If the total population is below the threshold, then optimize R12 population.
		  //if ((pop<=popthre)&&((minval-w)>1e-10))
		  if ((pop<=popthre)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;		      
		    }
		  // If the total population is above the threshold, then optimize the total population.
		  if ((pop>popthre)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }		  
		}

	      // Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	      else if (mode==3)
		{		  		  

		  w=*(vec2+ntypes-1);
	      
		  w=0;
		  for (i=0; i<=(ntypes-1); i++)
		    w=max(w,*(template+i*ntypes+i));	      

		  if ((minval-w)>1e-10)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=w; 		      
		    }
		}		  
	      
	      // Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
	      if (mode==4)
		{
		  // Minimize R12 population.
		  //if ((n<=initialtime)&&((minval-w)>1e-10))
		  if ((n<=initialtime)&&((minval2-w)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=w;
		    }
		  // Minimize the total population.
		  else if ((n>initialtime)&&((minval-g)>1e-10))
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
		    }
		}

	      // Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12.
	      else if (mode==5)
		{	
		  if ((minval2-w)>1e-10)
		    minval2=w;
	  
		  // If either drug is effective on R12, then the same as strategy 1.
		  //if ((*(Sa+6)>1e-10)||(*(Sa+7)>1e-10))
		    {
		      if ((minval-g)>1e-10)
			{
			  *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			  minval2_2=w;
			}
		    }

		  /*
		  // Otherwise check whether R12<1.
		  else
		    {
		      // Debug
		      //if (n==0)
		      //printf("d1=%.2f, d2=%.2f, g=%.3e, w=%.3e\n",d1,d2,g,w);

		      // If so then minimize the total population.
		      //if (w<(1-1e-10))
		      if (w<0.9)
			{
			  if ((minval-g)>1e-10)
			    {
			      *(mininds+0)=d1; *(mininds+1)=d2; minval=g;
			    }
			}
		    }
		  */
		}  
	    }
	}     
      
      
      // Strategy 5: If neither drug is effective on R12, then do the following.
      // If the minimum of R12 population >=1, then choose the strategy to minimize R12.
      // Otherwise stick to strategy 1.
      if ((mode==5)&&(*(Sa+6)<=1e-10)&&(*(Sa+7)<=1e-10))
	{	  
	  //if (minval2>=(1-1e-12))
	  if (minval2_2>=(1-1e-12))
	    {
	      i=0; k=1;
	      while ((i<=(cnt-1))&&(k==1))
		{
		  if (*(drugrecords+i*4+3)<=minval2)
		    {
		      *(mininds+0)=*(drugrecords+i*4+0);
		      *(mininds+1)=*(drugrecords+i*4+1);
		      k=0;
		    }
		  i++;
		}
	    }
	}
      
            

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
               
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
     

      if (val>=mortal)
	stopT=*(t+n);
    
   
      n++;
    }

  free(drugrecords); free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
 
  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}



// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
// Difference from optimize_two_drug_responses3: apply the mask to F when optimizing drug combinations.
// Difference from optimize_two_drug_responses3: apply the recursive function to predict the responses.
// Difference from optimize_two_drug_responses4: can choose either to use F to estimate the population or incur the recursive function.  The recursive function is time consuming and leads to early drop on cluster jobs.  So do not incur the recursive functions if the stop time >= upper bound using the simple method.
// Difference from optimize_two_drug_responses5: uses 5 strategies:
/* Strategy 1: Minimize total population.
   Strategy 2: Minimize total population if it is above a threshold.  Otherwise minimize the R12 population.   
   Strategy 3: Minimize total population with constraints R12<=1.
   Strategy 4: Maximize expected mortality/incurability time.
*/

double optimize_two_drug_responses6(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, minval, minval2, minval2_2, *mask, maxval;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF, *vec3;
  double dd=0.05, d1, d2, w, g, stopT, val, val2, pop, *drugrecords;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes); vec3=malloc(sizeof(double)*ntypes);
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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*9*k*k);  
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
	}
      else
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(y+i*nintervals+n);
	}

      pop=0;
      for (i=0; i<=(ntypes-1); i++)
	pop+=*(y0+i);

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
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100; minval2_2=1e100; cnt=0;      
      for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
      //for (d1=0; d1<=(1.0); d1+=dd)   
	{
	  for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
	  //for (d2=0; d2<=(1-d1); d2+=dd)
	    {
	      double tauinc, tauS, tauR1, tauR2, tauR12;

	      *(vec+0)=d1; *(vec+1)=d2;

	      for (i=0; i<=(ntypes-1); i++)
		*(vec2+i)=*(y0+i);

	      //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	      
	      
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
		      *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
		 		  
		  //
		  // Adjust tempF by the mask matrix.
		  //for (i=0; i<=(ntypes-1); i++)		
		  //for (j=0; j<=(ntypes-1); j++)
		  //  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
		  //		  		  		  

		  padm(ntypes,tempF,1,7,template);	
		  mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		  
		}


	      // compmode=1: incur the recursive function to estimate the population.
	      else if (compmode==1)
		recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
	     
	      
	      // vec2 is the estimated population vector at the end of the time interval.
	      g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      w=*(vec2+ntypes-1);
	      	      
	      
	      // Record the following quantities for each drug combination.
	      // (1)dosage of drug 1, (2)dosage of drug 2, (3)total population, (4)population of R12 cells, (5)time to incurability, (6)time to death from S cells, (7)time to death from R1 cells, (8)time to death from R2 cells.
	      // If R12 population >1, then tauinc=0.
	      // If a tau value <0, then set the tau value to positive infinity.
	      
	      cnt++; *(drugrecords+(cnt-1)*9+0)=d1; *(drugrecords+(cnt-1)*9+1)=d2;
	      *(drugrecords+(cnt-1)*9+2)=g; *(drugrecords+(cnt-1)*9+3)=w;

	      
	      // Calculate the time to incurability and mortality.
	      tauinc=*(vec2+1)*(*(T+13))+*(vec2+2)*(*(T+14));
	      tauinc*=*(g0+0); tauinc=(double)1/tauinc;

	      if (*(vec2+3)>1)
		tauinc=0;
	      if (tauinc<0)
		tauinc=1e20;

	      tauS=1e20;
	      if (*(vec2+0)>0)
		{
		  tauS=*(g0+0)-*(Sa+0)*d1-*(Sa+1)*d2;
		  tauS=(log(mortal)-log(*(vec2+0)))/tauS;
		}
	      if (tauS<0)
		tauS=1e20;	     
	      tauR1=1e20;
	      if (*(vec2+1)>0)
		{
		  tauR1=*(g0+0)-*(Sa+2)*d1-*(Sa+3)*d2;
		  tauR1=(log(mortal)-log(*(vec2+1)))/tauR1;
		}
	      if (tauR1<0)
		tauR1=1e20;
	      tauR2=1e20;
	      if (*(vec2+2)>0)
		{
		  tauR2=*(g0+0)-*(Sa+4)*d1-*(Sa+5)*d2;
		  tauR2=(log(mortal)-log(*(vec2+2)))/tauR2;
		}
	      if (tauR2<0)
		tauR2=1e20;
	      tauR12=1e20;
	      if (*(vec2+3)>0)
		{
		  tauR12=*(g0+0)-*(Sa+6)*d1-*(Sa+7)*d2;
		  tauR12=(log(mortal)-log(*(vec2+3)))/tauR12;
		}
	      if (tauR12<0)
		tauR12=1e20;

	      *(drugrecords+(cnt-1)*9+4)=tauinc;
	      *(drugrecords+(cnt-1)*9+5)=tauS;
	      *(drugrecords+(cnt-1)*9+6)=tauR1;
	      *(drugrecords+(cnt-1)*9+7)=tauR2;
	      *(drugrecords+(cnt-1)*9+8)=tauR12;
	    }
	}

      // Select the optimal dosage combination.
      *(mininds+0)=0; *(mininds+1)=0;

      
      // Strategy 1: Select the dosage combination to minimize total population.
      if (mode==1)
	{
	  minval=1e20;
	  for (i=0; i<=(cnt-1); i++)
	    {
	      if ((val=*(drugrecords+i*9+2))<minval)
		{
		  *(mininds+0)=*(drugrecords+i*9+0);
		  *(mininds+1)=*(drugrecords+i*9+1);
		  minval=val;
		}
	    }
	}

      // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
      else if (mode==2)
	{
	  minval=1e20;

	  // Total population is above the threshold: minimize total population.
	  if (pop>=popthre)
	    {
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+2))<minval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }
		}
	    }

	  // Total population is below the threshold: minimize R12 population.
	  else
	    {
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+3))<minval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }
		}
	    }
	}

      // Strategy 3: Minimize total population unless estimated R12 population>=1.  In that case minimize R12.  
      // But if the current population > 1 and R12 is not curable, then minimize the total population.
      else if (mode==3)
	{
	  // Find the dosage combination that minimize total population with constraint R12<1.
	  minval=1e20;
	  for (i=0; i<=(cnt-1); i++)
	    {
	      //if (((val=*(drugrecords+i*8+2))<minval)&&((val2=*(drugrecords+i*8+3))<(1+1e-10)))
	      
	      if (((val=*(drugrecords+i*9+2))<minval)&&(((val2=*(drugrecords+i*9+3))<(1+1e-10))||((*(y0+ntypes-1)>(1+1e-10))&&((*(g0+ntypes-1)-*(Sa+6)-*(Sa+7))>0))))
	      
		{
		  *(mininds+0)=*(drugrecords+i*9+0);
		  *(mininds+1)=*(drugrecords+i*9+1);
		  minval=val;
		}		
	    }
	  
	  // If no such dosage combination exists, then find the dosage combination that minimizes R12.
	  if (minval>=1e19)
	    {
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+3))<minval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }		
		}
	    }
	}

      // Strategy 4: If tauinc < tauS, tauR1, tauR2, and tauS, tauR1, tauR2 > 1 month, then maximize tauinc.  Otherwise maximize min(tauS,tauR1,tauR2).
      // But if the current population > 1 and R12 is not curable, then maximize min(tauS,tauR1,tauR2,tauR12).
      else if (mode==4)
	{
	  //
	  // Debug
	  //for (i=0; i<=(cnt-1); i++)
	  //{
	  //  if (((fabs(*(drugrecords+i*9+0)-0.55)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-0.45)<=1e-6))||((fabs(*(drugrecords+i*9+0)-0)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-1)<=1e-6)))
	  //{
	  //  printf("t=%f, dosage=(%.2f,%.2f), pop=(%.3e,%.3e), taus=(%.3e,%.3e,%.3e,%.3e)\n",*(t+n),*(drugrecords+i*9+0),*(drugrecords+i*9+1),*(drugrecords+i*9+2),*(drugrecords+i*9+3),*(drugrecords+i*9+4),*(drugrecords+i*9+5),*(drugrecords+i*9+6),*(drugrecords+i*9+7));
	  //}
	  //}
	  //
	  	  
	  
	  // Current R12 population <=1 or R12 is curable.
	  if ((*(y0+ntypes-1)<(1-1e-10))||((*(g0+ntypes-1)-*(Sa+6))<=0)||((*(g0+ntypes-1)-*(Sa+7))<=0))
	    {
	      // Find the dosage combination to maximize min(tauinc,tauS,tauR1,tauR2,tauR12) with constraint min(tauS,tauR1,tauR2,tauR12)>1 month.
	      maxval=0;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if (((val=min(*(drugrecords+i*9+4),min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))))>maxval)&&(min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))>30))
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      maxval=val;
		    }
		}

	      // If such dosage combination does not exist, then maximize min(tauS,tauR1,tauR2,tauR12).
	      if (maxval<=1e-5)
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  maxval=val;
			}
		    }
		}
	    }
	      
	  // Current R12 population >1 and R12 is incurable.
	  else
	    {
	      // Find the dosage combination to maximize min(tauS,tauR1,tauR2,tauR12).
	      maxval=0;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      maxval=val;
		    }
		}
	    }
	  

	  
	  /*
	  // Find the dosage combination to maximize min(tauinc,tauS,tauR1,tauR2) with constraint min(tauS,tauR1,tauR2)>1 month.
	  maxval=0;
	  for (i=0; i<=(cnt-1); i++)
	    {
	      if (((val=min(*(drugrecords+i*9+4),min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),*(drugrecords+i*9+7)))))>maxval)&&(min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),*(drugrecords+i*9+7)))>30))
		{
		  *(mininds+0)=*(drugrecords+i*9+0);
		  *(mininds+1)=*(drugrecords+i*9+1);
		  maxval=val;
		}
	    }

	  // If such dosage combination does not exist, then maximize min(tauS,tauR1,tauR2).
	  if (maxval<=1e-5)
	    {
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),*(drugrecords+i*9+7))))>maxval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      maxval=val;
		    }
		}
	    }
	  */

	  
	  //
	  // Find the dosage combination to maximize min(tauS,tauR1,tauR2) with constraint min(tauS,tauR1,tauR2)<1 month or tauinc>=min(tauS,tauR1,tauR2).
	  //maxval=0;
	  //for (i=0; i<=(cnt-1); i++)
	  //{
	  //  if (((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),*(drugrecords+i*9+7))))>maxval)&&((*(drugrecords+i*9+4)>=val)||(val<=30)))
	  //{
	  //  *(mininds+0)=*(drugrecords+i*9+0);
	  //  *(mininds+1)=*(drugrecords+i*9+1);
	  //  maxval=val;
	  //}
	  //}

	  // If such dosage combination does not exist, then maximize tauinc.
	  //if (maxval<=1e-5)
	  //{
	  //  for (i=0; i<=(cnt-1); i++)
	  //{
	  //  if ((val=*(drugrecords+i*9+4))>maxval)
	  //    {
	  //      *(mininds+0)=*(drugrecords+i*9+0);
	  //      *(mininds+1)=*(drugrecords+i*9+1);
	  //      maxval=val;
	  //    }
	  //}
	  //}
	  //
	}
            

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
               
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
     
      // Debug
      printf("t=%f, dosage=(%f,%f), y=(%.3e,%.3e,%.3e,%.3e), val=%.3e, mortal=%.3e\n",*(t+n),*(dosage+n),*(dosage+nintervals+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val,mortal);


      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
 
  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}




// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
// Difference from optimize_two_drug_responses3: apply the mask to F when optimizing drug combinations.
// Difference from optimize_two_drug_responses3: apply the recursive function to predict the responses.
// Difference from optimize_two_drug_responses4: can choose either to use F to estimate the population or incur the recursive function.  The recursive function is time consuming and leads to early drop on cluster jobs.  So do not incur the recursive functions if the stop time >= upper bound using the simple method.
// Difference from optimize_two_drug_responses5: uses 5 strategies:
/* Strategy 1: Minimize total population.
   Strategy 2: Minimize total population if it is above a threshold.  Otherwise minimize the R12 population.   
   Strategy 3: Minimize total population with constraints R12<=1.
   Strategy 4: Maximize expected mortality/incurability time.
*/
// Difference from optimize_two_drug_responses6: introduce a strategy 0.  Initial treat the patient with d1 if R1/N<=0.5.  Otherwise treat the patient with d2.  Switch between the two regimens when (a)the total population reaches double the nadir for the current regimen, or (b)the total population disappears with the current regimen and then re-emerges.
// Also the prediction period is imported instead of being fixed.

double optimize_two_drug_responses7(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, *mininds2, *mininds3, minval, minval2, minval3, *mask, maxval;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF, *vec3;
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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*9*k*k);  
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
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

	  // Initially if R1/N <= 0.5, then apply d1.  Otherwise apply d2.
	  if (n==0)
	    {
	      if (*(y0+1)<=(0.5*pop))
		{
		  *(mininds+0)=1; *(mininds+1)=0;
		}
	      else
		{
		  *(mininds+0)=0; *(mininds+1)=1;
		}
	    }
	  
	  // Switch the drug when the total population reaches the double the nadir or the population reemerges.
	  else
	    {
	      int switching=0;
	      double nadir=1e50;

	      // Find the nadir of the current regimen.
	      // A nadir is the local minimum of the population under the current regimen.
	      i=n-1; k=1;
	      while ((i>=0)&&(k==1))
		{
		  if ((fabs(*(dosage+i)-*(dosage+n-1))<=1e-12)&&(fabs(*(dosage+nintervals+i)-*(dosage+nintervals+n-1))<=1e-12))
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

	      // Debug
	      //printf("nadir=%.3e, pop=%.3e, val=%.3e, detected=%.3e, switch=%d\n",nadir,pop,val,detected,switching);


	      if (switching==1)
		{
		  if (*(dosage+n-1)>=0.01)
		    {
		      *(mininds+0)=0; *(mininds+1)=1;
		    }
		  else
		    {
		      *(mininds+0)=1; *(mininds+1)=0;
		    }
		}
	      else
		{
		  *(mininds+0)=*(dosage+n-1);
		  *(mininds+1)=*(dosage+nintervals+n-1);
		}
	    }
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
	  
	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;

	  for (i=0; i<=(ndrugs-1); i++)
	    {
	      *(mininds2+i)=0; *(mininds3+i)=0;
	    }

	  minval=1e100; minval2=1e100; minval3=1e100; cnt=0;      
	  for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
	    //for (d1=0; d1<=(1.0); d1+=dd)   
	    {
	      for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
		//for (d2=0; d2<=(1-d1); d2+=dd)
		{
		  double tauinc, tauS, tauR1, tauR2, tauR12;
		  
		  *(vec+0)=d1; *(vec+1)=d2;
		  
		  for (i=0; i<=(ntypes-1); i++)
		    *(vec2+i)=*(y0+i);
		  
		  //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		  
		  
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
			  //*(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
			  *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*predictionperiod;
		      
		      //
		      // Adjust tempF by the mask matrix.
		      //for (i=0; i<=(ntypes-1); i++)		
		      //for (j=0; j<=(ntypes-1); j++)
		      //  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
		      //		  		  		  
		      
		      padm(ntypes,tempF,1,7,template);	
		      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		      
		    }


		  // compmode=1: incur the recursive function to estimate the population.
		  else if (compmode==1)
		    //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		    recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,predictionperiod,vec,vec2);
	     
	      
		  // vec2 is the estimated population vector at the end of the time interval.
		  g=0;
		  for (i=0; i<=(ntypes-1); i++)
		    g+=*(vec2+i);
		  w=*(vec2+ntypes-1);


		  // Debug
		  //if (((mode==1)||(mode==3))&&(n==13))
		  //if (((mode==1)||(mode==3))&&(n==17))
		  //printf("dosage (%.3f, %.3f), population (%.3e, %.3e, %.3e, %.3e)\n",*(vec+0),*(vec+1),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3));
		  
		  
		  // Record the following quantities for each drug combination.
		  // (1)dosage of drug 1, (2)dosage of drug 2, (3)total population, (4)population of R12 cells, (5)time to incurability, (6)time to death from S cells, (7)time to death from R1 cells, (8)time to death from R2 cells.
		  // If R12 population >1, then tauinc=0.
		  // If a tau value <0, then set the tau value to positive infinity.
		  
		  cnt++; *(drugrecords+(cnt-1)*9+0)=d1; *(drugrecords+(cnt-1)*9+1)=d2;
		  *(drugrecords+(cnt-1)*9+2)=g; *(drugrecords+(cnt-1)*9+3)=w;
		  
		  
		  // Calculate the time to incurability and mortality.
		  tauinc=*(vec2+1)*(*(T+13))+*(vec2+2)*(*(T+14));
		  tauinc*=*(g0+0); tauinc=(double)1/tauinc;
		  
		  if (*(vec2+3)>1)
		    tauinc=0;
		  if (tauinc<0)
		    tauinc=1e20;
		  
		  tauS=1e20;
		  if (*(vec2+0)>0)
		    {
		      tauS=*(g0+0)-*(Sa+0)*d1-*(Sa+1)*d2;
		      tauS=(log(mortal)-log(*(vec2+0)))/tauS;
		    }
		  if (tauS<0)
		    tauS=1e20;	     
		  tauR1=1e20;
		  if (*(vec2+1)>0)
		    {
		      tauR1=*(g0+0)-*(Sa+2)*d1-*(Sa+3)*d2;
		      tauR1=(log(mortal)-log(*(vec2+1)))/tauR1;
		    }
		  if (tauR1<0)
		    tauR1=1e20;
		  tauR2=1e20;
		  if (*(vec2+2)>0)
		    {
		      tauR2=*(g0+0)-*(Sa+4)*d1-*(Sa+5)*d2;
		      tauR2=(log(mortal)-log(*(vec2+2)))/tauR2;
		    }
		  if (tauR2<0)
		    tauR2=1e20;
		  tauR12=1e20;
		  if (*(vec2+3)>0)
		    {
		      tauR12=*(g0+0)-*(Sa+6)*d1-*(Sa+7)*d2;
		      tauR12=(log(mortal)-log(*(vec2+3)))/tauR12;
		    }
		  if (tauR12<0)
		    tauR12=1e20;
		  
		  *(drugrecords+(cnt-1)*9+4)=tauinc;
		  *(drugrecords+(cnt-1)*9+5)=tauS;
		  *(drugrecords+(cnt-1)*9+6)=tauR1;
		  *(drugrecords+(cnt-1)*9+7)=tauR2;
		  *(drugrecords+(cnt-1)*9+8)=tauR12;
		}
	    }
	  
	  // Select the optimal dosage combination.
	  *(mininds+0)=0; *(mininds+1)=0;
	  
	  
	  // Strategy 1: Select the dosage combination to minimize total population.
	  if (mode==1)
	    {
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+2))<minval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }
		}
	    }
	  
	  // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	  else if (mode==2)
	    {
	      minval=1e20;
	      
	      // Total population is above the threshold: minimize total population.
	      if (pop>=popthre)
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+2))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}
		    }
		}
	      
	      // Total population is below the threshold: minimize R12 population.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+3))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}
		    }
		}
	    }
	  
	  // Strategy 3: Minimize total population unless estimated R12 population>=1.  In that case minimize R12.  
	  // But if the current R12 population > 1 and R12 is not curable, then minimize the total population.
	  else if (mode==3)
	    {
	      // Find the dosage combination that minimize total population with constraint R12<1.
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  //if (((val=*(drugrecords+i*8+2))<minval)&&((val2=*(drugrecords+i*8+3))<(1+1e-10)))
		  
		  if (((val=*(drugrecords+i*9+2))<minval)&&(((val2=*(drugrecords+i*9+3))<(1+1e-10))||((*(y0+ntypes-1)>(1+1e-10))&&((*(g0+ntypes-1)-*(Sa+6)-*(Sa+7))>0))))
		    
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }		
		}
	      
	      // If no such dosage combination exists, then find the dosage combination that minimizes R12.
	      if (minval>=1e19)
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+3))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}		
		    }
		}


	      
	      // A revision of strategy 3.
	      // If R12 is curable (g0<max(Sa(d1,R12),Sa(d2,R12))*0.1), then minimize the total population.
	      // If R1 is not curable and R12 population according to the optimal strategy >= 1, then minimize the total population.
	      // If R1 is not curable, R12 population according to the optimal strategy < 1, and the total population minimizer can also control R12 (R12 < 1), then minimize the total population.
	      // Otherwise minimize R12 population.
	      	      
	      minval2=1e20; l=-1;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+2))<minval2)
		    {
		      *(mininds2+0)=*(drugrecords+i*9+0);
		      *(mininds2+1)=*(drugrecords+i*9+1);
		      minval2=val; l=i;
		    }		
		}
	      minval3=1e20; m=-1;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+3))<minval3)
		    {
		      *(mininds3+0)=*(drugrecords+i*9+0);
		      *(mininds3+1)=*(drugrecords+i*9+1);
		      minval3=val; m=i;
		    }		
		}

	      // R12 is curable.  Minimize total population.
	      if ((*(g0+ntypes-1)<=(*(Sa+6)*0.1))||(*(g0+ntypes-1)<=(*(Sa+7)*0.1)))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // R1 is not curable and R12 population according to the optimal strategy >= 1.  Minimize the total population.
	      else if (minval3>(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // R1 is not curable, R12 population according to the optimal strategy < 1, and the total population minimizer can also control R12 (R12 < 1).  Minimize the total population.
	      else if (*(drugrecords+l*9+3)<(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // Otherwise minimize R12 population.
	      else
		{
		  *(mininds+0)=*(mininds3+0);
		  *(mininds+1)=*(mininds3+1);
		}
	      

	      /*
	      // Debug	      
	      if (minval3>(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}
	      else
		{
		  if (*(drugrecords+l*9+3)<(1-1e-10))
		    {
		      *(mininds+0)=*(mininds2+0);
		      *(mininds+1)=*(mininds2+1);
		    }
		  else
		    {
		      *(mininds+0)=*(mininds3+0);
		      *(mininds+1)=*(mininds3+1);
		    }
		}
	      */

	    }
	  
	  // Strategy 4: If tauinc < tauS, tauR1, tauR2, and tauS, tauR1, tauR2 > 1 month, then maximize tauinc.  Otherwise maximize min(tauS,tauR1,tauR2).
	  // But if the current population > 1 and R12 is not curable, then maximize min(tauS,tauR1,tauR2,tauR12).
	  else if (mode==4)
	    {
	      //
	      // Debug
	      //for (i=0; i<=(cnt-1); i++)
	      //{
	      //  if (((fabs(*(drugrecords+i*9+0)-0.55)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-0.45)<=1e-6))||((fabs(*(drugrecords+i*9+0)-0)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-1)<=1e-6)))
	      //{
	      //  printf("t=%f, dosage=(%.2f,%.2f), pop=(%.3e,%.3e), taus=(%.3e,%.3e,%.3e,%.3e)\n",*(t+n),*(drugrecords+i*9+0),*(drugrecords+i*9+1),*(drugrecords+i*9+2),*(drugrecords+i*9+3),*(drugrecords+i*9+4),*(drugrecords+i*9+5),*(drugrecords+i*9+6),*(drugrecords+i*9+7));
	      //}
	      //}
	      //
	      
	      
	      // Current R12 population <=1 or R12 is curable.
	      if ((*(y0+ntypes-1)<(1-1e-10))||((*(g0+ntypes-1)-*(Sa+6))<=0)||((*(g0+ntypes-1)-*(Sa+7))<=0))
		{
		  // Find the dosage combination to maximize min(tauinc,tauS,tauR1,tauR2,tauR12) with constraint min(tauS,tauR1,tauR2,tauR12)>1 month.
		  maxval=0;
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if (((val=min(*(drugrecords+i*9+4),min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))))>maxval)&&(min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))>30))
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  maxval=val;
			}
		    }
		  
		  // If such dosage combination does not exist, then maximize min(tauS,tauR1,tauR2,tauR12).
		  if (maxval<=1e-5)
		    {
		      for (i=0; i<=(cnt-1); i++)
			{
			  if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
			    {
			      *(mininds+0)=*(drugrecords+i*9+0);
			      *(mininds+1)=*(drugrecords+i*9+1);
			      maxval=val;
			    }
			}
		    }
		}
	      
	      // Current R12 population >1 and R12 is incurable.
	      else
		{
		  // Find the dosage combination to maximize min(tauS,tauR1,tauR2,tauR12).
		  maxval=0;
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  maxval=val;
			}
		    }
		}	      	      
	    }	 
	}
                  

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
               
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
     
      // Debug
      printf("t=%f, dosage=(%f,%f), y=(%.3e,%.3e,%.3e,%.3e), val=%.3e, mortal=%.3e\n",*(t+n),*(dosage+n),*(dosage+nintervals+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val,mortal);

      //if (compmode==1)
      //printf("%f%c%f%c%f%c%.3e%c%.3e%c%.3e%c%.3e%c%.3e\n",*(t+n),0x09,*(dosage+n),0x09,*(dosage+nintervals+n),0x09,*(vec2+0),0x09,*(vec2+1),0x09,*(vec2+2),0x09,*(vec2+3),0x09,val);
      //printf("t=%f, dosage=(%f,%f), y=(%.3e,%.3e,%.3e,%.3e), val=%.3e, mortal=%.3e\n",*(t+n),*(dosage+n),*(dosage+nintervals+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val,mortal);



      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
 
  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}



// Find the optimal drug dosage and responses by applying 4 possible strategies.
// Strategy 1: Minimize total population.
// Strategy 2: Minimize R12 population unless the total population is above a threshold.  If the total population is below the threshold, then switch back to minimizing R12 population again.
// Strategy 3: Minimize the max growth rate for 4 cell types.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
// Strategy 4: Initially minimize R12 population.  Then switch to minimize total population.
// Strategy 5: If the drugs are effective on R12, then minimize total population.  If the drugs are ineffective on R12, then minimize total population unless R12 population >=1.  In that case minimize R12 population.
// Difference from optimize_two_drug_responses3: apply the mask to F when optimizing drug combinations.
// Difference from optimize_two_drug_responses3: apply the recursive function to predict the responses.
// Difference from optimize_two_drug_responses4: can choose either to use F to estimate the population or incur the recursive function.  The recursive function is time consuming and leads to early drop on cluster jobs.  So do not incur the recursive functions if the stop time >= upper bound using the simple method.
// Difference from optimize_two_drug_responses5: uses 5 strategies:
/* Strategy 1: Minimize total population.
   Strategy 2: Minimize total population if it is above a threshold.  Otherwise minimize the R12 population.   
   Strategy 3: Minimize total population with constraints R12<=1.
   Strategy 4: Maximize expected mortality/incurability time.
*/
// Difference from optimize_two_drug_responses6: introduce a strategy 0.  Initial treat the patient with d1 if R1/N<=0.5.  Otherwise treat the patient with d2.  Switch between the two regimens when (a)the total population reaches double the nadir for the current regimen, or (b)the total population disappears with the current regimen and then re-emerges.
// Also the prediction period is imported instead of being fixed.
// Difference from optimize_two_drug_responses7: apply an improved recurse_two_drug_response_trimrates function.

double optimize_two_drug_responses8(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, int mode, int compmode, double popthre, double predictionperiod, double *dosage)
{
  int i, j, k, l, m, n, cnt;
  double *d, *mininds, *mininds2, *mininds3, minval, minval2, minval3, *mask, maxval;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF, *vec3;
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

  val=(double)1/dd; k=(int)(val+1);  
  drugrecords=malloc(sizeof(double)*9*k*k);  
  
  
  // Iterate at each time step.
  n=0; y0=malloc(sizeof(double)*ntypes);
  vec=malloc(sizeof(double)*ndrugs); stopT=-1;
  while ((n<(nintervals-1))&&(stopT<0))
    {
      
      // Retrieve the population from the previous step.
      if (n==0)
	{
	  for (i=0; i<=(ntypes-1); i++)
	    *(y0+i)=*(x0+i);
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

	  // Initially if R1/N <= 0.5, then apply d1.  Otherwise apply d2.
	  if (n==0)
	    {
	      if (*(y0+1)<=(0.5*pop))
		{
		  *(mininds+0)=1; *(mininds+1)=0;
		}
	      else
		{
		  *(mininds+0)=0; *(mininds+1)=1;
		}
	    }
	  
	  // Switch the drug when the total population reaches the double the nadir or the population reemerges.
	  else
	    {
	      int switching=0;
	      double nadir=1e50;

	      // Find the nadir of the current regimen.
	      // A nadir is the local minimum of the population under the current regimen.
	      i=n-1; k=1;
	      while ((i>=0)&&(k==1))
		{
		  if ((fabs(*(dosage+i)-*(dosage+n-1))<=1e-12)&&(fabs(*(dosage+nintervals+i)-*(dosage+nintervals+n-1))<=1e-12))
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

	      // Debug
	      //printf("nadir=%.3e, pop=%.3e, val=%.3e, detected=%.3e, switch=%d\n",nadir,pop,val,detected,switching);


	      if (switching==1)
		{
		  if (*(dosage+n-1)>=0.01)
		    {
		      *(mininds+0)=0; *(mininds+1)=1;
		    }
		  else
		    {
		      *(mininds+0)=1; *(mininds+1)=0;
		    }
		}
	      else
		{
		  *(mininds+0)=*(dosage+n-1);
		  *(mininds+1)=*(dosage+nintervals+n-1);
		}
	    }
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
	  
	  // Consider the drug combinations (d1,d2), d1+d2=1.

	  for (i=0; i<=(ndrugs-1); i++)
	    *(mininds+i)=0;

	  for (i=0; i<=(ndrugs-1); i++)
	    {
	      *(mininds2+i)=0; *(mininds3+i)=0;
	    }

	  minval=1e100; minval2=1e100; minval3=1e100; cnt=0;      
	  for (d1=0; d1<=(1.0+1e-10); d1+=dd)      
	    //for (d1=0; d1<=(1.0); d1+=dd)   
	    {
	      //for (d2=0; d2<=(1-d1+1e-10); d2+=dd)	 
		//for (d2=0; d2<=(1-d1); d2+=dd)
		{
		  double tauinc, tauS, tauR1, tauR2, tauR12;
		
		  // Debug: only consider d1+d2=1.
		  d2=1-d1;
  
		  *(vec+0)=d1; *(vec+1)=d2;
		  
		  for (i=0; i<=(ntypes-1); i++)
		    *(vec2+i)=*(y0+i);
		  
		  //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		  
		  
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
			  //*(tempF+i*ntypes+j)=*(F+i*ntypes+j)*timeinterval;
			  *(tempF+i*ntypes+j)=*(F+i*ntypes+j)*predictionperiod;
		      
		      //
		      // Adjust tempF by the mask matrix.
		      //for (i=0; i<=(ntypes-1); i++)		
		      //for (j=0; j<=(ntypes-1); j++)
		      //  *(tempF+i*ntypes+j)*=*(mask+i*ntypes+j);
		      //		  		  		  
		      
		      padm(ntypes,tempF,1,7,template);	
		      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
		      
		    }


		  // compmode=1: incur the recursive function to estimate the population.
		  else if (compmode==1)
		    //recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,vec,vec2);
		    recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,predictionperiod,vec,vec2);
	     
	      
		  // vec2 is the estimated population vector at the end of the time interval.
		  g=0;
		  for (i=0; i<=(ntypes-1); i++)
		    g+=*(vec2+i);
		  w=*(vec2+ntypes-1);


		  // Debug
		  //if (((mode==1)||(mode==3))&&(n==13))
		  //if (((mode==1)||(mode==3))&&(n==17))
		  //printf("dosage (%.3f, %.3f), population (%.3e, %.3e, %.3e, %.3e)\n",*(vec+0),*(vec+1),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3));
		  
		  
		  // Record the following quantities for each drug combination.
		  // (1)dosage of drug 1, (2)dosage of drug 2, (3)total population, (4)population of R12 cells, (5)time to incurability, (6)time to death from S cells, (7)time to death from R1 cells, (8)time to death from R2 cells.
		  // If R12 population >1, then tauinc=0.
		  // If a tau value <0, then set the tau value to positive infinity.
		  
		  cnt++; *(drugrecords+(cnt-1)*9+0)=d1; *(drugrecords+(cnt-1)*9+1)=d2;
		  *(drugrecords+(cnt-1)*9+2)=g; *(drugrecords+(cnt-1)*9+3)=w;
		  
		  
		  // Calculate the time to incurability and mortality.
		  tauinc=*(vec2+1)*(*(T+13))+*(vec2+2)*(*(T+14));
		  tauinc*=*(g0+0); tauinc=(double)1/tauinc;
		  
		  if (*(vec2+3)>1)
		    tauinc=0;
		  if (tauinc<0)
		    tauinc=1e20;
		  
		  tauS=1e20;
		  if (*(vec2+0)>0)
		    {
		      tauS=*(g0+0)-*(Sa+0)*d1-*(Sa+1)*d2;
		      tauS=(log(mortal)-log(*(vec2+0)))/tauS;
		    }
		  if (tauS<0)
		    tauS=1e20;	     
		  tauR1=1e20;
		  if (*(vec2+1)>0)
		    {
		      tauR1=*(g0+0)-*(Sa+2)*d1-*(Sa+3)*d2;
		      tauR1=(log(mortal)-log(*(vec2+1)))/tauR1;
		    }
		  if (tauR1<0)
		    tauR1=1e20;
		  tauR2=1e20;
		  if (*(vec2+2)>0)
		    {
		      tauR2=*(g0+0)-*(Sa+4)*d1-*(Sa+5)*d2;
		      tauR2=(log(mortal)-log(*(vec2+2)))/tauR2;
		    }
		  if (tauR2<0)
		    tauR2=1e20;
		  tauR12=1e20;
		  if (*(vec2+3)>0)
		    {
		      tauR12=*(g0+0)-*(Sa+6)*d1-*(Sa+7)*d2;
		      tauR12=(log(mortal)-log(*(vec2+3)))/tauR12;
		    }
		  if (tauR12<0)
		    tauR12=1e20;
		  
		  *(drugrecords+(cnt-1)*9+4)=tauinc;
		  *(drugrecords+(cnt-1)*9+5)=tauS;
		  *(drugrecords+(cnt-1)*9+6)=tauR1;
		  *(drugrecords+(cnt-1)*9+7)=tauR2;
		  *(drugrecords+(cnt-1)*9+8)=tauR12;
		}
	    }
	  
	  // Select the optimal dosage combination.
	  *(mininds+0)=0; *(mininds+1)=0;
	  
	  
	  // Strategy 1: Select the dosage combination to minimize total population.
	  if (mode==1)
	    {
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+2))<minval)
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }
		}
	    }
	  
	  // Strategy 2: Select the dosage combination to minimize R12 population.  If total population reaches a threshold, then switch to strategy 1.  If total population is below the threshold, then switch back to the original strategy.
	  else if (mode==2)
	    {
	      minval=1e20;
	      
	      // Total population is above the threshold: minimize total population.
	      if (pop>=popthre)
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+2))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}
		    }
		}
	      
	      // Total population is below the threshold: minimize R12 population.
	      else
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+3))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}
		    }
		}
	    }
	  
	  // Strategy 3: Minimize total population unless estimated R12 population>=1.  In that case minimize R12.  
	  // But if the current R12 population > 1 and R12 is not curable, then minimize the total population.
	  else if (mode==3)
	    {
	      // Find the dosage combination that minimize total population with constraint R12<1.
	      minval=1e20;
	      for (i=0; i<=(cnt-1); i++)
		{
		  //if (((val=*(drugrecords+i*8+2))<minval)&&((val2=*(drugrecords+i*8+3))<(1+1e-10)))
		  
		  if (((val=*(drugrecords+i*9+2))<minval)&&(((val2=*(drugrecords+i*9+3))<(1+1e-10))||((*(y0+ntypes-1)>(1+1e-10))&&((*(g0+ntypes-1)-*(Sa+6)-*(Sa+7))>0))))
		    
		    {
		      *(mininds+0)=*(drugrecords+i*9+0);
		      *(mininds+1)=*(drugrecords+i*9+1);
		      minval=val;
		    }		
		}
	      
	      // If no such dosage combination exists, then find the dosage combination that minimizes R12.
	      if (minval>=1e19)
		{
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=*(drugrecords+i*9+3))<minval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  minval=val;
			}		
		    }
		}


	      
	      // A revision of strategy 3.
	      // If R12 is curable (g0<max(Sa(d1,R12),Sa(d2,R12))*0.1), then minimize the total population.
	      // If R1 is not curable and R12 population according to the optimal strategy >= 1, then minimize the total population.
	      // If R1 is not curable, R12 population according to the optimal strategy < 1, and the total population minimizer can also control R12 (R12 < 1), then minimize the total population.
	      // Otherwise minimize R12 population.
	      	      
	      minval2=1e20; l=-1;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+2))<minval2)
		    {
		      *(mininds2+0)=*(drugrecords+i*9+0);
		      *(mininds2+1)=*(drugrecords+i*9+1);
		      minval2=val; l=i;
		    }		
		}
	      minval3=1e20; m=-1;
	      for (i=0; i<=(cnt-1); i++)
		{
		  if ((val=*(drugrecords+i*9+3))<minval3)
		    {
		      *(mininds3+0)=*(drugrecords+i*9+0);
		      *(mininds3+1)=*(drugrecords+i*9+1);
		      minval3=val; m=i;
		    }		
		}

	      // R12 is curable.  Minimize total population.
	      if ((*(g0+ntypes-1)<=(*(Sa+6)*0.1))||(*(g0+ntypes-1)<=(*(Sa+7)*0.1)))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // R1 is not curable and R12 population according to the optimal strategy >= 1.  Minimize the total population.
	      else if (minval3>(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // R1 is not curable, R12 population according to the optimal strategy < 1, and the total population minimizer can also control R12 (R12 < 1).  Minimize the total population.
	      else if (*(drugrecords+l*9+3)<(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}

	      // Otherwise minimize R12 population.
	      else
		{
		  *(mininds+0)=*(mininds3+0);
		  *(mininds+1)=*(mininds3+1);
		}
	      

	      /*
	      // Debug	      
	      if (minval3>(1-1e-10))
		{
		  *(mininds+0)=*(mininds2+0);
		  *(mininds+1)=*(mininds2+1);
		}
	      else
		{
		  if (*(drugrecords+l*9+3)<(1-1e-10))
		    {
		      *(mininds+0)=*(mininds2+0);
		      *(mininds+1)=*(mininds2+1);
		    }
		  else
		    {
		      *(mininds+0)=*(mininds3+0);
		      *(mininds+1)=*(mininds3+1);
		    }
		}
	      */

	    }
	  
	  // Strategy 4: If tauinc < tauS, tauR1, tauR2, and tauS, tauR1, tauR2 > 1 month, then maximize tauinc.  Otherwise maximize min(tauS,tauR1,tauR2).
	  // But if the current population > 1 and R12 is not curable, then maximize min(tauS,tauR1,tauR2,tauR12).
	  else if (mode==4)
	    {
	      //
	      // Debug
	      //for (i=0; i<=(cnt-1); i++)
	      //{
	      //  if (((fabs(*(drugrecords+i*9+0)-0.55)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-0.45)<=1e-6))||((fabs(*(drugrecords+i*9+0)-0)<=1e-6)&&(fabs(*(drugrecords+i*9+1)-1)<=1e-6)))
	      //{
	      //  printf("t=%f, dosage=(%.2f,%.2f), pop=(%.3e,%.3e), taus=(%.3e,%.3e,%.3e,%.3e)\n",*(t+n),*(drugrecords+i*9+0),*(drugrecords+i*9+1),*(drugrecords+i*9+2),*(drugrecords+i*9+3),*(drugrecords+i*9+4),*(drugrecords+i*9+5),*(drugrecords+i*9+6),*(drugrecords+i*9+7));
	      //}
	      //}
	      //
	      
	      
	      // Current R12 population <=1 or R12 is curable.
	      if ((*(y0+ntypes-1)<(1-1e-10))||((*(g0+ntypes-1)-*(Sa+6))<=0)||((*(g0+ntypes-1)-*(Sa+7))<=0))
		{
		  // Find the dosage combination to maximize min(tauinc,tauS,tauR1,tauR2,tauR12) with constraint min(tauS,tauR1,tauR2,tauR12)>1 month.
		  maxval=0;
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if (((val=min(*(drugrecords+i*9+4),min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))))>maxval)&&(min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8))))>30))
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  maxval=val;
			}
		    }
		  
		  // If such dosage combination does not exist, then maximize min(tauS,tauR1,tauR2,tauR12).
		  if (maxval<=1e-5)
		    {
		      for (i=0; i<=(cnt-1); i++)
			{
			  if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
			    {
			      *(mininds+0)=*(drugrecords+i*9+0);
			      *(mininds+1)=*(drugrecords+i*9+1);
			      maxval=val;
			    }
			}
		    }
		}
	      
	      // Current R12 population >1 and R12 is incurable.
	      else
		{
		  // Find the dosage combination to maximize min(tauS,tauR1,tauR2,tauR12).
		  maxval=0;
		  for (i=0; i<=(cnt-1); i++)
		    {
		      if ((val=min(*(drugrecords+i*9+5),min(*(drugrecords+i*9+6),min(*(drugrecords+i*9+7),*(drugrecords+i*9+8)))))>maxval)
			{
			  *(mininds+0)=*(drugrecords+i*9+0);
			  *(mininds+1)=*(drugrecords+i*9+1);
			  maxval=val;
			}
		    }
		}	      	      
	    }	 
	}
                  

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);            
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);
               
      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);
     
      // Debug
      //printf("t=%f, dosage=(%f,%f), y=(%.3e,%.3e,%.3e,%.3e), val=%.3e, mortal=%.3e\n",*(t+n),*(dosage+n),*(dosage+nintervals+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val,mortal);

      //if (compmode==1)
      //printf("%f%c%f%c%f%c%.3e%c%.3e%c%.3e%c%.3e%c%.3e\n",*(t+n),0x09,*(dosage+n),0x09,*(dosage+nintervals+n),0x09,*(vec2+0),0x09,*(vec2+1),0x09,*(vec2+2),0x09,*(vec2+3),0x09,val);
      //printf("t=%f, dosage=(%f,%f), y=(%.3e,%.3e,%.3e,%.3e), val=%.3e, mortal=%.3e\n",*(t+n),*(dosage+n),*(dosage+nintervals+n),*(vec2+0),*(vec2+1),*(vec2+2),*(vec2+3),val,mortal);



      if (val>=mortal)
	stopT=*(t+n);
          
      n++;
    }

  
  free(drugrecords); free(mask);
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);
  
 
  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}
