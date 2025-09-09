// Utility functions involved in drug responses and optimization.

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# define max(i,j) (i>j ? i : j)
# define min(i,j) (i<j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)


extern void padm(int ndim, double *A, double t, int p, double *E);
extern char * getitemval3(char *s, int ind, int *length, char sepch);
extern double atof2(char *s, int length);
extern void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
extern void diagonal(int ndim, double *vec, double *mat);

extern int ndrugs, ntypes;
extern double detected, mortal;


double optimize_two_drug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, double *dosage);
double optimize_two_drug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, double *dosage);
void recurse_two_drug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);
void recurse_two_drug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses);


// Find an optimal treatment strategy and stopping time for a specific paramater setting.
// Divide the time into small intervals.  Within each interval vary the drug dosage to minimize the cell population.
// Do nothing to address the fractional cell number problem.
double optimize_two_drug_responses1(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, double *dosage)
{
  int i, j, k, l, m, n;
  double *d, *mininds, minval, minval2;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF;
  double dd=0.1, d1, d2, w, g, stopT, val;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);

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
	            
      // Change the drug combinations to minimize the total population at the end of each interval.
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100;
      for (d1=0; d1<=1; d1+=dd)      
	{
	  for (d2=0; d2<=(1-d1); d2+=dd)	 
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
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
	      //padm(ntypes,F,timeinterval,7,template);	      
	      padm(ntypes,tempF,1,7,template);	
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
	      w=*(vec2+ntypes-1);

	      // Debug
	      w=0;
	      for (i=0; i<=(ntypes-1); i++)
		w=max(w,*(template+i*ntypes+i));
	      
	       	      
	      if ((minval-w)>1e-10)
		{
		  *(mininds+0)=d1; *(mininds+1)=d2; minval=w; minval2=0;
		  for (i=0; i<=(ntypes-1); i++)
		    minval2+=*(vec2+i);
		}
	    }
	}

      
      // If other drug responses yield lower entire population, then choose it.
      for (d1=0; d1<=1; d1+=dd)
	{
	  for (d2=0; d2<=1; d2+=dd)
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
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
	      //padm(ntypes,F,timeinterval,7,template);
	      padm(ntypes,tempF,1,7,template);
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
	      w=*(vec2+ntypes-1); g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      if (fabs(minval-w)<1e-30)
		{
		  if ((minval2-g)>1e-6)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=g;
		    }
		}	     
	    }
	}           

      // Debug
      //printf("mininds=[%f %f]\n",*(mininds+0),*(mininds+1));

      // Fix dosage at each time interval.
      *(vec+0)=*(mininds+0); *(vec+1)=*(mininds+1);
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
      //padm(ntypes,F,timeinterval,7,template);
      padm(ntypes,tempF,1,7,template);
      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
            
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);

      
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);

      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);

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

      printf("template=\n");
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ntypes-1); j++)
	    printf("%.3e ",*(template+i*ntypes+j));
	  printf("\n");
	}
      printf("tempF=\n");
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ntypes-1); j++)
	    printf("%.3e ",*(tempF+i*ntypes+j));
	  printf("\n");
	}
      */
   
      n++;
    }
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);

  /*
  // Find the time where the total population reaches a mortal threshold.
  stopT=-1; n=0;
  while ((n<nintervals)&&(stopT<0))
    {
      double val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n);
      if (val>=mortal)
	stopT=*(t+n);
      n++;
    }
  */

  // Debug
  //printf("stopT=%f\n",stopT);

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}


// Find an optimal treatment strategy and stopping time for a specific paramater setting.
// Divide the time into small intervals.  Within each interval vary the drug dosage to minimize the cell population.
// Apply different parameters to the equations when fractional cell numbers are encountered.
double optimize_two_drug_responses2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double *t, double timeinterval, int nintervals, double *y, double *dosage)
{
  int i, j, k, l, m, n;
  double *d, *mininds, minval, minval2;
  double *F, *vec, *template, *constF, *y0, *vec2, *template2, *template3, *tempF;
  double dd=0.1, d1, d2, w, g, stopT, val;

  d=malloc(sizeof(double)*ndrugs); mininds=malloc(sizeof(double)*ndrugs);
  F=malloc(sizeof(double)*ntypes*ntypes);
  tempF=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);
  template2=malloc(sizeof(double)*ntypes*ntypes);
  template3=malloc(sizeof(double)*ntypes*ntypes);
  constF=malloc(sizeof(double)*ntypes*ntypes);
  vec2=malloc(sizeof(double)*ntypes);

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
	      
  
      // Change the drug combinations to minimize the total population at the end of each interval.
      for (i=0; i<=(ndrugs-1); i++)
	*(mininds+i)=0;
      minval=1e100; minval2=1e100;
      for (d1=0; d1<=1; d1+=dd)      
	{
	  for (d2=0; d2<=(1-d1); d2+=dd)	 
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
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
	      //padm(ntypes,F,timeinterval,7,template);	      
	      padm(ntypes,tempF,1,7,template);	
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
	      w=*(vec2+ntypes-1);

	      // Debug
	      w=0;
	      for (i=0; i<=(ntypes-1); i++)
		w=max(w,*(template+i*ntypes+i));	      
	      
	      if ((minval-w)>1e-10)
		{
		  *(mininds+0)=d1; *(mininds+1)=d2; minval=w; minval2=0;
		  for (i=0; i<=(ntypes-1); i++)
		    minval2+=*(vec2+i);
		}
	    }
	}

      
      // If other drug responses yield lower entire population, then choose it.
      for (d1=0; d1<=1; d1+=dd)
	{
	  for (d2=0; d2<=1; d2+=dd)
	    {
	      *(vec+0)=d1; *(vec+1)=d2;
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
	      //padm(ntypes,F,timeinterval,7,template);
	      padm(ntypes,tempF,1,7,template);
	      mat_multiply(ntypes,ntypes,1,template,y0,vec2);
	      w=*(vec2+ntypes-1); g=0;
	      for (i=0; i<=(ntypes-1); i++)
		g+=*(vec2+i);
	      if (fabs(minval-w)<1e-30)
		{
		  if ((minval2-g)>1e-6)
		    {
		      *(mininds+0)=d1; *(mininds+1)=d2; minval2=g;
		    }
		}	     
	    }
	}           

      // Fix dosage at each time interval.
      // Predict the drug response by running the modified differential equations.      
      //recurse_two_drug_response_trimrates(x0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,0,timeinterval,mininds,vec2);
      
      *(dosage+n)=*(mininds+0); *(dosage+nintervals+n)=*(mininds+1);

      
      for (i=0; i<=(ntypes-1); i++)
	*(y+i*nintervals+n+1)=*(vec2+i);

      val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n+1);

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
  
  *(dosage+n)=*(dosage+n-1);
  *(dosage+nintervals+n)=*(dosage+nintervals+n-1);

  
  // If stopT<0 then the patient recovers.
  if (stopT<0)
    stopT=*(t+nintervals-1);

  /*
  // Find the time where the total population reaches a mortal threshold.
  stopT=-1; n=0;
  while ((n<nintervals)&&(stopT<0))
    {
      double val=0;
      for (i=0; i<=(ntypes-1); i++)
	val+=*(y+i*nintervals+n);
      if (val>=mortal)
	stopT=*(t+n);
      n++;
    }
  */

  // Debug
  //printf("stopT=%f\n",stopT);

  // Release memory.
  free(d); free(mininds); free(y0); free(tempF);
  free(F); free(vec); free(vec2); free(template); free(template2); free(template3);
  
  
  return stopT;
}




// Predict the drug response of populations using the modified differential equations.
// Trim the rates of the terms if the population size drops below 1.
// Only report the responses at end of the time interval.
// Recursively run the differential equations with the masked terms till timeinterval.  If there are populations crossing the boundary of 1, then step back and find the first time when the boundary crossing occurs.
void recurse_two_drug_response_trimrates(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;

  // Debug
  //printf("startt=%f, endt=%f\n",startt,endt);
    
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

  /*
  // Debug
  printf("initial y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  */

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
  
  /*
  // Debug
  printf("F:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(F+i*ntypes+j));
      printf("\n");
    }
  
  // Debug
  printf("mask:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(mask+i*ntypes+j));
      printf("\n");
    }
  */  


  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));

  /*
  // Debug
  printf("startt=%f, endt=%f\n",startt,endt);
  printf("F2:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(F2+i*ntypes+j));
      printf("\n");
    }
  */

  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template);

  //padm(ntypes,F2,endt-startt,7,template); 

  /*
  // Debug
  printf("expm(F2t):\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(template+i*ntypes+j));
      printf("\n");
    }
  exit(0);
  */

  mat_multiply(ntypes,ntypes,1,template,y0,yn);
      
  /*
  // Debug
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(yn+i));
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

  /*
  // Debug
  printf("startt=%f, endt=%f\n",startt,endt);
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(yn+i));
  printf("\n");  
  printf("alter: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%d ",*(alters+i));
  printf("\n");  
  */

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
		
	      /*
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>1e-5)&&((*(vec3+i)-1)>1e-5))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>1e-5)&&((1-*(vec3+i))>1e-5))
		*(alters+i)=1;
	      */

	      /*
	      if ((*(vec2+i)<1)&&(*(vec3+i)>1))
		*(alters+i)=1;
	      else if ((*(vec2+i)>1)&&(*(vec3+i)<1))
		*(alters+i)=-1;
		*/
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;

	  /*
	  // Debug
	  printf("t1=%f, t2=%f, curt=%f\n",t1,t2,curt);
	  printf("y0: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(vec2+i));
	  printf("\n");
	  printf("yn: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(vec3+i));
	  printf("\n");
	  printf("alters: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%d ",*(alters+i));
	  printf("\n");	  
	  */

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

      
      // Debug
      //printf("final t1=%f, t2=%f, curt=%f\n",t1,t2,curt);
      //exit(0);
      
      
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
	  recurse_two_drug_response_trimrates(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
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
void recurse_two_drug_response_trimrates2(double *x0, double *T, double *g0, double *Sg, double *a0, double *Sa, double startt, double endt, double *dosage, double *responses)
{
  int i, j, k, l, m, n, *alters;
  double *F, *F2, *template, *constF, *y0, *vec2, *vec3, *template2, *template3, *mask, *yn, *tempF;
  double val;

  // Debug
  //printf("startt=%f, endt=%f\n",startt,endt);
    
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

  /*
  // Debug
  printf("initial y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  */

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
  
  /*
  // Debug
  printf("F:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(F+i*ntypes+j));
      printf("\n");
    }
  
  // Debug
  printf("mask:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(mask+i*ntypes+j));
      printf("\n");
    }
  */  


  for (i=0; i<=(ntypes-1); i++)		
    for (j=0; j<=(ntypes-1); j++)
      *(F2+i*ntypes+j)=*(F+i*ntypes+j)*(*(mask+i*ntypes+j));

  /*
  // Debug
  printf("startt=%f, endt=%f\n",startt,endt);
  printf("F2:\n");
  for (i=0; i<=(ntypes-1); i++)		
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(F2+i*ntypes+j));
      printf("\n");
    }
  */

  for (i=0; i<=(ntypes-1); i++)
    for (j=0; j<=(ntypes-1); j++)
      *(tempF+i*ntypes+j)=*(F2+i*ntypes+j)*(endt-startt);
  padm(ntypes,tempF,1,7,template);

  //padm(ntypes,F2,endt-startt,7,template); 

  /*
  // Debug
  printf("expm(F2t):\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(template+i*ntypes+j));
      printf("\n");
    }
  //exit(0);
  */

  mat_multiply(ntypes,ntypes,1,template,y0,yn);
      
  /*
  // Debug
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(yn+i));
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

  /*
  // Debug
  printf("startt=%f, endt=%f\n",startt,endt);
  printf("y0: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(y0+i));
  printf("\n");
  printf("yn: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.3e ",*(yn+i));
  printf("\n");  
  printf("alter: ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%d ",*(alters+i));
  printf("\n");  
  */

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
				

	      /*
	      // Debug
	      // Avoid the deadlock at boundary 1.
	      if (((1-*(vec2+i))>1e-6)&&((*(vec3+i)-1)>1e-6))
		*(alters+i)=1;
	      else if (((*(vec2+i)-1)>1e-6)&&((1-*(vec3+i))>1e-6))
		*(alters+i)=1;
		*/

	      /*
	      if ((*(vec2+i)<1)&&(*(vec3+i)>1))
		*(alters+i)=1;
	      else if ((*(vec2+i)>1)&&(*(vec3+i)<1))
		*(alters+i)=-1;
		*/
	    }
	  k=0;
	  for (i=0; i<=(ntypes-1); i++)
	    if (*(alters+i)!=0)
	      k++;

	  /*
	  // Debug
	  //if ((t1>=1)&&(t1<2)&&(t2>=44))
	    {
	      printf("int start\n");
	      printf("t1=%f, t2=%f, curt=%f\n",t1,t2,curt);
	      printf("y0: ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.3e ",*(vec2+i));
	      printf("\n");
	      printf("yn: ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%.3e ",*(vec3+i));
	      printf("\n");
	      printf("alters: ");
	      for (i=0; i<=(ntypes-1); i++)
		printf("%d ",*(alters+i));
	      printf("\n");	  
	      printf("int end\n");
	    }
	  */

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

      
      // Debug
      //printf("final t1=%f, t2=%f, curt=%f\n",t1,t2,curt);
      //exit(0);
      
      
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

	  /*
	  // Debug
	  printf("vec3: ");
	  for (i=0; i<=(ntypes-1); i++)
	    printf("%.3e ",*(vec3+i));
	  printf("\n");
	  */

	  recurse_two_drug_response_trimrates2(y0,T,g0,Sg,a0,Sa,curt,endt,dosage,responses);
	}      

      // Release memory.
      free(y0); free(vec3);
    }
	    

  return;
}


