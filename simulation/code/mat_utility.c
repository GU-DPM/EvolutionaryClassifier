// Utility functions involved in matrix and rate evaluations.


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include "phylogeny.h"
# define MAXLINEWIDTH 100000
# define NAMELENGTH 100
# define LONGNAMELENGTH 100
# define max(i,j) (i>j ? i : j)
# define min(i,j) (i<j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)

extern int nu2ind(char ch);
extern int aa2ind(char ch);
double * kron(int nrows1, int ncols1, double *mat1, int nrows2, int ncols2, double *mat2);
void calculate_joint_rate(double *q1, double *q2, double *q2d);
void padm(int ndim, double *A, double t, int p, double *E);
void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3);
void mat_inv(int ndim, double *A1, double *A2);
void ludcmp(int ndim, double *A1, int *indx, int *d, double *A2);
void lubksb(int ndim, double *A, int *indx, double *b);
double llratio(int nnodes, int **child, double *t, int *order, int nseq, int *seqs, double *prior2, double *nulltrans, double *activetrans);
struct exp_real felsenstein(int nseq, int *x, int nnodes, int **child, int *order, double *prior2, double *transitions);
struct exp_real add(struct exp_real v1, struct exp_real v2);
struct exp_real multiply1(struct exp_real v1, struct exp_real v2);
struct exp_real multiply2(struct exp_real v1, double v2);
double llratio_quantized_rates(int nnodes, int **child, double *t, int *order, int nseq, int *seqs, int Nbllevels, double *bllevels, double *prior2, double *nulltransitions, double *activetransitions);
struct exp_real felsenstein_quantized_rates(int nseq, int *x, int nnodes, int **child, double *t, int *order, int Nbllevels, double *bllevels, double *prior2, double *transitions);
void diagonal(int ndim, double *vec, double *mat);

extern int nchar, njointchar, nnodes;
extern double epsilon;

// Calculate the Kronecker tensor product of two matrices.
double * kron(int nrows1, int ncols1, double *mat1, int nrows2, int ncols2, double *mat2)
{
  int i, j, k, l, m, n, nrows, ncols;
  double *mat;

  nrows=nrows1*nrows2; ncols=ncols1*ncols2;
  mat=malloc(sizeof(double)*nrows*ncols);
  for (i=0; i<=(nrows1-1); i++)
    {
      for (j=0; j<=(ncols1-1); j++)
	{
	  for (k=0; k<=(nrows2-1); k++)
	    {
	      for (l=0; l<=(ncols2-1); l++)
		{
		  m=i*nrows2+k; n=j*ncols2+l;
		  *(mat+m*ncols+n)=*(mat1+i*ncols1+j)*(*(mat2+k*ncols2+l));
		}
	    }
	}
    }
 
  return mat;
}

// Calculate the independent and dependent rate matrices from the rate matrix of single states.
void calculate_joint_rate(double *q1, double *q2, double *q2d)
{
  int i, j, k, a, b, c, d, rowind, colind, *phi2d;
  double val, perval;

  // Get the independent rate matrix.  Simultaneous transition rates=0, non-simultaneous transition rates=single aa transition rates.
  for (i=0; i<=(njointchar*njointchar-1); i++)
    *(q2+i)=0;

  for (a=0; a<=(nchar-1); a++)
    {
      for (b=0; b<=(nchar-1); b++)
	{
	  for (c=0; c<=(nchar-1); c++)
	    {
	      for (d=0; d<=(nchar-1); d++)
		{
		  rowind=a*nchar+b; colind=c*nchar+d;
		  if ((a!=c)&&(b==d))
		    *(q2+rowind*njointchar+colind)=*(q1+a*nchar+c);
		  else if ((a==c)&&(b!=d))
		    *(q2+rowind*njointchar+colind)=*(q1+b*nchar+d);
		}
	    }
	}
    }
    
  for (a=0; a<=(njointchar-1); a++)
    {
      val=0;
      for (b=0; b<=(njointchar-1); b++)
	if (b!=a)
	  val-=*(q2+a*njointchar+b);
      *(q2+a*njointchar+a)=val;
    }

  // Get the dependent matrix.  Reweight the independent matrix using simultaneous changes.
  phi2d=malloc(sizeof(int)*njointchar*njointchar);
  for (i=0; i<=(njointchar*njointchar-1); i++)
    {
      *(q2d+i)=0; *(phi2d+i)=0;
    }

  for (a=0; a<=(nchar-1); a++)
    {
      for (b=0; b<=(nchar-1); b++)
	{
	  for (c=0; c<=(nchar-1); c++)
	    {
	      for (d=0; d<=(nchar-1); d++)
		{
		  rowind=a*nchar+b; colind=c*nchar+d;
		  
		  // Simultaneous changes.
		  if ((a!=c)&&(b!=d))
		    *(phi2d+rowind*njointchar+colind)=2;
		  
		  // Unilateral changes.
		  else if (((a!=c)&&(b==d))||((a==c)&&(b!=d)))
		    *(phi2d+rowind*njointchar+colind)=1;
		}
	    }
	}
    }

  for (a=0; a<=(njointchar-1); a++)
    {
      val=0; k=0;
      for (b=0; b<=(njointchar-1); b++)
	{
	  if (*(phi2d+a*njointchar+b)==1)
	    val+=*(q2+a*njointchar+b);
	  else if (*(phi2d+a*njointchar+b)==2)
	    k++;
	}
      perval=val*(1-epsilon)/(double)k;      
      for (b=0; b<=(njointchar-1); b++)
	{
	  if (*(phi2d+a*njointchar+b)==1)
	    *(q2d+a*njointchar+b)=*(q2+a*njointchar+b)*epsilon;
	  else if (*(phi2d+a*njointchar+b)==2)
	    *(q2d+a*njointchar+b)=perval;
	}
      val=0;
      for (b=0; b<=(njointchar-1); b++)
	if (b!=a)
	  val-=*(q2d+a*njointchar+b);
      *(q2d+a*njointchar+a)=val;
    }

  free(phi2d);
}

// Pade approximation of matrix exponential.
// Inputs: a square matrix, scaling constant t, approximation degree.
// Output: a square matrix.
void padm(int ndim, double *A, double t, int p, double *E)
{
  int i, j, k, l, m, n, odd;
  double *temp, *temp2, *tempA, *c, s, val, val2, *A2, *P, *Q;

  
  // If t=0, then return the identifty matrix.
  if (fabs(t)<1e-50)
    {
      for (i=0; i<=(ndim*ndim-1); i++)
	*(E+i)=0;
      for (i=1; i<=ndim; i++)
	*(E+(i-1)*ndim+(i-1))=1.0;
      return;
    }
  
  // Pade coefficients.
  c=malloc(sizeof(double)*(p+2));
  *c=0; *(c+1)=1;
  for (k=1; k<=p; k++)
    *(c+k+1)=*(c+k)*((double)(p+1-k)/(double)(k*(2*p+1-k)));
  
  // Scaling.
  tempA=malloc(sizeof(double)*ndim*ndim);
    
  s=0;
  for (i=1; i<=ndim; i++)
    {
      val=0;
      for (j=1; j<=ndim; j++)
	{
	  val2=*(A+(i-1)*ndim+(j-1)); val2=fabs(val2);
	  val+=val2;
	}
      if (val>s)
	s=val;
    }
  
  if (s>0.5)
    {
      val=log(s)/log(2);
      if (val>=0)
	k=(int)val;
      else
	k=(-1)*((int)((-1)*val));
      k=max(0,k+2);
      val=(-1)*k*log(2);
      val=exp(val);
      for (i=0; i<=(ndim*ndim-1); i++)
	*(tempA+i)=*(A+i)*val*t;
      s=k;
    }
  else
    for (i=0; i<=(ndim*ndim-1); i++)
      *(tempA+i)=*(A+i)*t;
  
  
  // Horner evaluation of the irreducible fraction (see ref. above).  
  A2=malloc(sizeof(double)*ndim*ndim);
  mat_multiply(ndim,ndim,ndim,tempA,tempA,A2);
  Q=malloc(sizeof(double)*ndim*ndim);
  P=malloc(sizeof(double)*ndim*ndim);
  for (i=0; i<=(ndim*ndim-1); i++)
    {
      *(Q+i)=0; *(P+i)=0;
    }
  for (i=1; i<=ndim; i++)
    {
      *(Q+(i-1)*ndim+(i-1))=*(c+p+1);
      *(P+(i-1)*ndim+(i-1))=*(c+p);
    }
  odd=1; temp=malloc(sizeof(double)*ndim*ndim);
  temp2=malloc(sizeof(double)*ndim*ndim);   
  for (k=(p-1); k>=1; k--)
    {
      if (odd==1)
	{
	  mat_multiply(ndim,ndim,ndim,Q,A2,temp);
	  for (i=1; i<=ndim; i++)
	    *(temp+(i-1)*ndim+(i-1))+=*(c+k);
	  for (i=0; i<=(ndim*ndim-1); i++)
	    *(Q+i)=*(temp+i);
	}
      else
	{
	  mat_multiply(ndim,ndim,ndim,P,A2,temp);
	  for (i=1; i<=ndim; i++)
	    *(temp+(i-1)*ndim+(i-1))+=*(c+k);
	  for (i=0; i<=(ndim*ndim-1); i++)
	    *(P+i)=*(temp+i);
	}
      odd=1-odd;
    }  

  if (odd==1)
    {
      mat_multiply(ndim,ndim,ndim,Q,tempA,temp);            
      for (i=0; i<=(ndim*ndim-1); i++)
	*(Q+i)=*(temp+i)-*(P+i);            

      mat_inv(ndim,Q,temp2);      
   
      mat_multiply(ndim,ndim,ndim,temp2,P,temp);
      for (i=0; i<=(ndim*ndim-1); i++)
	*(E+i)=*(temp+i)*(-2);
      for (i=1; i<=ndim; i++)
	*(E+(i-1)*ndim+(i-1))-=1;	
    }  
  else
    {
      mat_multiply(ndim,ndim,ndim,P,tempA,temp);
      for (i=0; i<=(ndim*ndim-1); i++)
	{
	  *(P+i)=*(temp+i);
	  *(Q+i)=*(Q+i)-*(P+i);
	}
      mat_inv(ndim,Q,temp2);
      mat_multiply(ndim,ndim,ndim,temp2,P,temp);
      for (i=0; i<=(ndim*ndim-1); i++)
	*(E+i)=*(temp+i)*2;
      for (i=1; i<=ndim; i++)
	*(E+(i-1)*ndim+(i-1))+=1;
    }
  
  
  // Squaring.
  i=(int)s;
  for (k=1; k<=i; k++)
    {
      mat_multiply(ndim,ndim,ndim,E,E,temp);
      for (j=0; j<=(ndim*ndim-1); j++)
	*(E+j)=*(temp+j);
    }  

  free(A2); free(Q); free(P); free(temp); free(temp2);     
  free(tempA); free(c); 
 
}


// Matrix multiplication.  A3=A1*A2.  
// A1: ndim1 x ndim2, A2: ndim2 x ndim3, A3: ndim1 x ndim3.
void mat_multiply(int ndim1, int ndim2, int ndim3, double *A1, double *A2, double *A3)
{
  int i, j, k;
  double val;
  for (i=1; i<=ndim1; i++)
    {
      for (j=1; j<=ndim3; j++)
	{
	  val=0;
	  for (k=1; k<=ndim2; k++)
	    val+=*(A1+(i-1)*ndim2+(k-1))*(*(A2+(k-1)*ndim3+(j-1)));
	  *(A3+(i-1)*ndim3+(j-1))=val;
	}
    }
}

// Square matrix inversion.  Use the C codes in numerical recipe.
// Inputs: matrix dimension ndim, matrix A1.
// Output: inverted matrix A2.
void mat_inv(int ndim, double *A1, double *A2)
{
  int i, j, k, *indx, d;
  double *temp, *cols;
    
  indx=malloc(sizeof(int)*ndim);
  temp=malloc(sizeof(double)*ndim*ndim);
  cols=malloc(sizeof(double)*ndim);
  
  // Perform LU decomposition.
  ludcmp(ndim,A1,indx,&d,temp);    

  // Solve the linear system of equations.
  for (j=1; j<=ndim; j++)
    {
      for (i=1; i<=ndim; i++)
	*(cols+i-1)=0;
      *(cols+j-1)=1.0;
      lubksb(ndim,temp,indx,cols);
      for (i=1; i<=ndim; i++)
	*(A2+(i-1)*ndim+(j-1))=*(cols+i-1);
    }  
  free(indx); free(temp); free(cols);
}


// LU decomposition of a square matrix.  Use the C codes in numerical recipe.
// Inputs: matrix dimension ndim, matrix A1, row permutation affected by partial pivoting indx, eveness/oddnedd of row interchange.
void ludcmp(int ndim, double *A1, int *indx, int *d, double *A2)
{
  int i, j, k, imax=0;
  double big=0.0, dum, sum, temp, *vv;
  
  for (i=0; i<=(ndim*ndim-1); i++)
    *(A2+i)=*(A1+i); 

  vv=malloc(sizeof(double)*ndim); *d=1;
  for (i=1; i<=ndim; i++)
    *(vv+i-1)=i;

  
  for (i=1; i<=ndim; i++)
    {
      big=0.0;
      for (j=1; j<=ndim; j++)
	if ((temp=fabs(*(A2+(i-1)*ndim+(j-1))))>big)
	  big=temp;
      *(vv+i-1)=1.0/big;
    }   
  

  
  for (j=1; j<=ndim; j++)
    {
      for (i=1; i<j; i++)
	{
	  sum=*(A2+(i-1)*ndim+(j-1));
	  for (k=1; k<i; k++)
	    sum-=*(A2+(i-1)*ndim+(k-1))*(*(A2+(k-1)*ndim+(j-1)));
	  *(A2+(i-1)*ndim+(j-1))=sum;
	}
      big=0.0;
      for (i=j; i<=ndim; i++)
	{
	  sum=*(A2+(i-1)*ndim+(j-1));
	  for (k=1; k<j; k++)
	    sum-=*(A2+(i-1)*ndim+(k-1))*(*(A2+(k-1)*ndim+(j-1)));
	  *(A2+(i-1)*ndim+(j-1))=sum;
	  if ((dum=*(vv+i-1)*fabs(sum))>=big)
	    {
	      big=dum;
	      imax=i;
	    }
	}
      
      if (j!=imax)
	{
	  for (k=1; k<=ndim; k++)
	    {
	      dum=*(A2+(imax-1)*ndim+(k-1));
	      *(A2+(imax-1)*ndim+(k-1))=*(A2+(j-1)*ndim+(k-1));
	      *(A2+(j-1)*ndim+(k-1))=dum;
	    }	  
	  *d=-(*d);
	  *(vv+imax-1)=*(vv+j-1);	  

	  // Debug
	  //printf("ndim=%d, vv(%d)=vv(%d)=%f\n",ndim,imax-1,j-1,*(vv+j-1));

	}
      
      *(indx+j-1)=imax;
      if (*(A2+(j-1)*ndim+(j-1))==0.0)
	*(A2+(j-1)*ndim+(j-1))=1e-15;

      if (j!=ndim)
	{	  
	  dum=*(A2+(j-1)*ndim+(j-1)); dum=1.0/dum;
	  for (i=(j+1); i<=ndim; i++)
	    *(A2+(i-1)*ndim+(j-1))*=dum;
	}      
    }
  

  // Debug
  //exit(0);
  
  free(vv);

  // Debug
  //exit(0);

}


// Perform forward and backward substitutions to solve the equation Ax=b.
// Inputs: A, b, ndim, indx.
// Output: x replaces b.
void lubksb(int ndim, double *A, int *indx, double *b)
{
  int i, ii=0, ip, j;
  double sum, val;

  for (i=1; i<=ndim; i++)
    {
      ip=*(indx+i-1);
      sum=*(b+ip-1);
      *(b+ip-1)=*(b+i-1);
      if (ii)
	for (j=ii; j<=(i-1); j++)
	  sum-=*(A+(i-1)*ndim+(j-1))*(*(b+j-1));
      else if (sum)
	ii=i;
      *(b+i-1)=sum;
    }

  for (i=ndim; i>=1; i--)
    {
      sum=*(b+i-1);
      for (j=(i+1); j<=ndim; j++)
	sum-=*(A+(i-1)*ndim+(j-1))*(*(b+j-1));
      val=*(A+(i-1)*ndim+(i-1));
      *(b+i-1)=sum/val;
    }
}

// Use Felsenstein's algorithm to calculate the log likelihood ratio of one pair.
// Use pre-calculated transition probability on each branch.
double llratio(int nnodes, int **child, double *t, int *order, int nseq, int *seqs, double *prior2, double *nulltransitions, double *activetransitions)
{
  int i, j, k, l, m, n, *x;
  double a, b;
  struct exp_real *jLtree, tmpval;
  double retval, val;

  // Convert the sequences into integers.
  x=malloc(sizeof(int)*nseq*2);
  for (i=0; i<=(nseq-1); i++)
    {
      *(x+i)=*(seqs+i*2+0);
      *(x+i+nseq)=*(seqs+i*2+1);     
    }

  // Calculate the likelihood scores.
  jLtree=malloc(sizeof(struct exp_real)*2);
  tmpval=felsenstein(nseq,x,nnodes,child,order,prior2,nulltransitions);
  (jLtree+0)->mantissa=tmpval.mantissa; (jLtree+0)->exponent=tmpval.exponent;
  tmpval=felsenstein(nseq,x,nnodes,child,order,prior2,activetransitions);
  (jLtree+1)->mantissa=tmpval.mantissa; (jLtree+1)->exponent=tmpval.exponent;
  val=(jLtree+1)->mantissa/(jLtree+0)->mantissa; k=(jLtree+1)->exponent-(jLtree+0)->exponent;
  retval=(double)k+log10(val);
  free(jLtree); free(x);
  return retval;
}

// Felsenstein's algorithm of calculating the marginal likelihood of sequences on the leaves of a phylogenetic tree.
struct exp_real felsenstein(int nseq, int *x, int nnodes, int **child, int *order, double *prior2, double *transitions)
{
  int c, i, j, k, l, m, n;
  double *temp, a, b;
  double val, val2, val3;
  struct exp_real *jL, oneoverjoint, oneoversingle, sumval, tmpval, maxval;

  jL=malloc(sizeof(struct exp_real)*nnodes*(njointchar+1));
  temp=malloc(sizeof(double)*njointchar*njointchar);
  val=(double)1/(double)(njointchar);
  val2=log10(val); val2=floor(val2); 
  k=(int)val2; val2=(double)pow(10,k);
  val/=(double)val2;
  oneoverjoint.mantissa=val; oneoverjoint.exponent=k;
  val=(double)1/(double)(nchar);
  val2=log10(val); val2=floor(val2); 
  k=(int)val2; val2=(double)pow(10,k);
  val/=(double)val2;
  oneoversingle.mantissa=val; oneoversingle.exponent=k;

  // Follow the order to calculate the likelihood scores.
  for (i=0; i<=(nnodes-1); i++)
    {
      char ch1, ch2;

      n=*(order+i);
      
      // Leaves: assign the entire mass on the observed state.  Do not consider gaps.
      if ((**(child+n)==0)&&(*(x+nseq+n)!=0)&&(*(x+n)!=0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }

	  j=*(x+nseq+n)+(*(x+n)-1)*nchar;
	  (jL+n*(njointchar+1)+j)->mantissa=1.0;
	  (jL+n*(njointchar+1)+j)->exponent=0;	  
	}

      // Leaves with single gaps.  Apply uniform probability to each consistent state.
      else if ((**(child+n)==0)&&(*(x+nseq+n)==0)&&(*(x+n)!=0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }
	  for (k=1; k<=nchar; k++)
	    {
	      j=k+(*(x+n)-1)*nchar;
	      (jL+n*(njointchar+1)+j)->mantissa=oneoversingle.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoversingle.exponent;
	    }
	}

      else if ((**(child+n)==0)&&(*(x+nseq+n)!=0)&&(*(x+n)==0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }
	  for (k=1; k<=nchar; k++)
	    {
	      j=*(x+nseq+n)+(k-1)*nchar;
	      (jL+n*(njointchar+1)+j)->mantissa=oneoversingle.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoversingle.exponent;
	    }
	}

      // Leaves with double gaps.  Apply uniform probability to each state.
      else if ((**(child+n)==0)&&(*(x+nseq+n)==0)&&(*(x+n)==0))
	{
	  for (j=1; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=oneoverjoint.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoverjoint.exponent;
	    }	  
	}
  
      // Internal nodes: the likelihood score is the marginalized product of the likelihood scores of children weighted by branch transition probabilities.
      else if (**(child+n)>0)
	{
	  for (j=1; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=1.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }

	  for (j=1; j<=**(child+n); j++)
	    {
	      int c=*(*(child+n)+j);
	      struct exp_real *vals;

	      vals=malloc(sizeof(struct exp_real)*(njointchar+1));
	      
	      // Copy the transition probability into temp.
	      for (l=0; l<=(njointchar*njointchar-1); l++)
		*(temp+l)=*(transitions+c*njointchar*njointchar+l);

	      // Calculate the posterior probability.	      
	      for (l=1; l<=njointchar; l++)
		{
		  (vals+l)->mantissa=0.0; (vals+l)->exponent=0.0;
		}
	      for (l=1; l<=njointchar; l++)
		{
		  for (m=1; m<=njointchar; m++)
		    {
		      tmpval=multiply2(*(jL+c*(njointchar+1)+m),*(temp+(l-1)*njointchar+(m-1)));
		      tmpval=add(*(vals+l),tmpval);
		      (vals+l)->mantissa=tmpval.mantissa;
		      (vals+l)->exponent=tmpval.exponent;
		    }
		}
	      for (l=1; l<=njointchar; l++)
		{
		  tmpval=multiply1(*(jL+n*(njointchar+1)+l),*(vals+l));
		  (jL+n*(njointchar+1)+l)->mantissa=tmpval.mantissa;
		  (jL+n*(njointchar+1)+l)->exponent=tmpval.exponent;
		}
	      free(vals);
	    }
	}            
    }

  // Likelihood of the tree.
  sumval.mantissa=0.0; sumval.exponent=0.0;
  for (j=1; j<=njointchar; j++)
    {
      tmpval=multiply2(*(jL+(nnodes-1)*(njointchar+1)+j),*(prior2+j-1));
      tmpval=add(sumval,tmpval);
      sumval.mantissa=tmpval.mantissa; sumval.exponent=tmpval.exponent;
    } 
  free(jL); free(temp);  
  return sumval;
}
  
// Add two numbers in exp_real formats.
struct exp_real add(struct exp_real v1, struct exp_real v2)
{
  int i, j, scale;
  double val;
  struct exp_real retval;

  val=v1.mantissa; val=fabs(val);
  if (val<1e-50)
    {
      retval.mantissa=v2.mantissa; retval.exponent=v2.exponent;
      return retval;
    } 
  val=v2.mantissa; val=fabs(val);
  if (val<1e-50)
    {
      retval.mantissa=v1.mantissa; retval.exponent=v1.exponent;
      return retval;
    } 
  
  i=argmax(v1.exponent,v2.exponent);
  if (i==1)
    {
      j=v2.exponent-v1.exponent;
      val=(double)pow(10,j);
      val*=v2.mantissa;
      retval.mantissa=v1.mantissa+val;
      retval.exponent=v1.exponent;
    }
  else
    {
      j=v1.exponent-v2.exponent;
      val=(double)pow(10,j);
      val*=v1.mantissa;
      retval.mantissa=v2.mantissa+val;
      retval.exponent=v2.exponent;
    }

  if (retval.mantissa>=10.0)
    {
      retval.mantissa/=10.0; retval.exponent+=1;
    }

  return retval;
}


// Multiply two numbers in exp_real formats.
struct exp_real multiply1(struct exp_real v1, struct exp_real v2)
{
  int i;
  double val;
  struct exp_real retval;

  val=v1.mantissa*v2.mantissa;
  i=v1.exponent+v2.exponent;
  
  if (val>=10.0)
    {
      val/=10.0; i+=1;
    }

  if (fabs(val)<1e-50)
    i=0;

  retval.mantissa=val; retval.exponent=i;
  return retval;
}


// Multiply a exp_real number and a double number.
struct exp_real multiply2(struct exp_real v1, double v2)
{
  int i, inv=0;
  double val;
  struct exp_real retval;
 
  retval.mantissa=v1.mantissa*v2; retval.exponent=v1.exponent;
  
  if (retval.mantissa<0)
    {
      inv=1; retval.mantissa*=(-1);
    } 

  if (fabs(retval.mantissa)<1e-50)
    {
      retval.mantissa=0; retval.exponent=0;
    }
  else
    {
      // Obtain the scientific representation of retval.mantissa.
      val=retval.mantissa; i=0;
      if (val>1.0)
	{
	  while (val>10.0)
	    {
	      val/=10; i++;
	    }
	}
      else if (val<1.0)
	{
	  while (val<1.0)
	    {
	      val*=10; i--;
	    }
	}
      retval.mantissa=val; retval.exponent+=i;

      if (inv==1)
	retval.mantissa=(-1)*retval.mantissa;
    } 

  return retval;
}

// Use Felsenstein's algorithm to calculate the log likelihood ratio of one pair.
// Use pre-calculated transition probability on each branch.
// Difference from llratio: use quantized conditional probability.
double llratio_quantized_rates(int nnodes, int **child, double *t, int *order, int nseq, int *seqs, int Nbllevels, double *bllevels, double *prior2, double *nulltransitions, double *activetransitions)
{
  int i, j, k, l, m, n, *x;
  double a, b;
  struct exp_real *jLtree, tmpval;
  double retval, val;

  // Convert the sequences into integers.
  x=malloc(sizeof(int)*nseq*2);
  for (i=0; i<=(nseq-1); i++)
    {
      *(x+i)=*(seqs+i*2+0);
      *(x+i+nseq)=*(seqs+i*2+1);     
    }

  // Calculate the likelihood scores.
  jLtree=malloc(sizeof(struct exp_real)*2);
  tmpval=felsenstein_quantized_rates(nseq,x,nnodes,child,t,order,Nbllevels,bllevels,prior2,nulltransitions);
  (jLtree+0)->mantissa=tmpval.mantissa; (jLtree+0)->exponent=tmpval.exponent;
  tmpval=felsenstein_quantized_rates(nseq,x,nnodes,child,t,order,Nbllevels,bllevels,prior2,activetransitions);
  (jLtree+1)->mantissa=tmpval.mantissa; (jLtree+1)->exponent=tmpval.exponent;
  val=(jLtree+1)->mantissa/(jLtree+0)->mantissa; k=(jLtree+1)->exponent-(jLtree+0)->exponent;
  retval=(double)k+log10(val);
  free(jLtree); free(x);
  return retval;
}

// Felsenstein's algorithm of calculating the marginal likelihood of sequences on the leaves of a phylogenetic tree.
// Difference from felsenstein: use quantized conditional probability.
struct exp_real felsenstein_quantized_rates(int nseq, int *x, int nnodes, int **child, double *t, int *order, int Nbllevels, double *bllevels, double *prior2, double *transitions)
{
  int c, i, j, k, l, m, n;
  double *temp, a, b;
  double val, val2, val3;
  struct exp_real *jL, oneoverjoint, oneoversingle, sumval, tmpval, maxval;

  jL=malloc(sizeof(struct exp_real)*nnodes*(njointchar+1));
  temp=malloc(sizeof(double)*njointchar*njointchar);
  val=(double)1/(double)(njointchar);
  val2=log10(val); val2=floor(val2); 
  k=(int)val2; val2=(double)pow(10,k);
  val/=(double)val2;
  oneoverjoint.mantissa=val; oneoverjoint.exponent=k;
  val=(double)1/(double)(nchar);
  val2=log10(val); val2=floor(val2); 
  k=(int)val2; val2=(double)pow(10,k);
  val/=(double)val2;
  oneoversingle.mantissa=val; oneoversingle.exponent=k;

  // Follow the order to calculate the likelihood scores.
  for (i=0; i<=(nnodes-1); i++)
    {
      char ch1, ch2;

      n=*(order+i);
      
      // Leaves: assign the entire mass on the observed state.  Do not consider gaps.
      if ((**(child+n)==0)&&(*(x+nseq+n)!=0)&&(*(x+n)!=0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }

	  j=*(x+nseq+n)+(*(x+n)-1)*nchar;
	  (jL+n*(njointchar+1)+j)->mantissa=1.0;
	  (jL+n*(njointchar+1)+j)->exponent=0;	  
	}

      // Leaves with single gaps.  Apply uniform probability to each consistent state.
      else if ((**(child+n)==0)&&(*(x+nseq+n)==0)&&(*(x+n)!=0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }
	  for (k=1; k<=nchar; k++)
	    {
	      j=k+(*(x+n)-1)*nchar;
	      (jL+n*(njointchar+1)+j)->mantissa=oneoversingle.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoversingle.exponent;
	    }
	}

      else if ((**(child+n)==0)&&(*(x+nseq+n)!=0)&&(*(x+n)==0))
	{
	  for (j=0; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=0.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }
	  for (k=1; k<=nchar; k++)
	    {
	      j=*(x+nseq+n)+(k-1)*nchar;
	      (jL+n*(njointchar+1)+j)->mantissa=oneoversingle.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoversingle.exponent;
	    }
	}

      // Leaves with double gaps.  Apply uniform probability to each state.
      else if ((**(child+n)==0)&&(*(x+nseq+n)==0)&&(*(x+n)==0))
	{
	  for (j=1; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=oneoverjoint.mantissa;
	      (jL+n*(njointchar+1)+j)->exponent=oneoverjoint.exponent;
	    }	  
	}
  
      // Internal nodes: the likelihood score is the marginalized product of the likelihood scores of children weighted by branch transition probabilities.
      else if (**(child+n)>0)
	{
	  for (j=1; j<=njointchar; j++)
	    {
	      (jL+n*(njointchar+1)+j)->mantissa=1.0;
	      (jL+n*(njointchar+1)+j)->exponent=0;
	    }

	  for (j=1; j<=**(child+n); j++)
	    {
	      int c=*(*(child+n)+j);
	      struct exp_real *vals;

	      vals=malloc(sizeof(struct exp_real)*(njointchar+1));
	      
	      // Evaluate the probability transition matrix using the precomputed quantized matrices.
	      k=0;
	      while ((k<=(Nbllevels-1))&&(*(t+c)>=*(bllevels+k)))
		k++;
	      if (k>(Nbllevels-1))
		k=Nbllevels-1;
	      for (l=0; l<=(njointchar*njointchar-1); l++)
		*(temp+l)=*(transitions+k*njointchar*njointchar+l);		

	      // Copy the transition probability into temp.
	      //for (l=0; l<=(njointchar*njointchar-1); l++)
	      //*(temp+l)=*(transitions+c*njointchar*njointchar+l);

	      // Calculate the posterior probability.	      
	      for (l=1; l<=njointchar; l++)
		{
		  (vals+l)->mantissa=0.0; (vals+l)->exponent=0.0;
		}
	      for (l=1; l<=njointchar; l++)
		{
		  for (m=1; m<=njointchar; m++)
		    {
		      tmpval=multiply2(*(jL+c*(njointchar+1)+m),*(temp+(l-1)*njointchar+(m-1)));
		      tmpval=add(*(vals+l),tmpval);
		      (vals+l)->mantissa=tmpval.mantissa;
		      (vals+l)->exponent=tmpval.exponent;
		    }
		}
	      for (l=1; l<=njointchar; l++)
		{
		  tmpval=multiply1(*(jL+n*(njointchar+1)+l),*(vals+l));
		  (jL+n*(njointchar+1)+l)->mantissa=tmpval.mantissa;
		  (jL+n*(njointchar+1)+l)->exponent=tmpval.exponent;
		}
	      free(vals);
	    }
	}            
    }

  // Likelihood of the tree.
  sumval.mantissa=0.0; sumval.exponent=0.0;
  for (j=1; j<=njointchar; j++)
    {
      tmpval=multiply2(*(jL+(nnodes-1)*(njointchar+1)+j),*(prior2+j-1));
      tmpval=add(sumval,tmpval);
      sumval.mantissa=tmpval.mantissa; sumval.exponent=tmpval.exponent;
    } 
  free(jL); free(temp);  
  return sumval;
}


// Make a diagonal matrix from a vector.
void diagonal(int ndim, double *vec, double *mat)
{
  int i, j;
  for (i=0; i<=(ndim*ndim-1); i++)
    *(mat+i)=0;
  for (i=0; i<=(ndim-1); i++)
    *(mat+i*ndim+i)=*(vec+i);
  return;
}



// Calculate exponential of the growth rate matrix.
// It has a specific format thus can be directly evaluated instead of employing padm.
void drugmat_expm(int ntypes, double *A, double *expA)
{
  int i, j, k, l, m, n, ndistinct, ndrugs, *memberships;
  double *as, *deltas, *V, *iV, *D, val, val0, shift=1e-10, *distinct;
  double *template;
  
  val=log2((double)ntypes); ndrugs=(int)val;
  as=malloc(sizeof(double)*ntypes);
  deltas=malloc(sizeof(double)*ndrugs);
  V=malloc(sizeof(double)*ntypes*ntypes);
  iV=malloc(sizeof(double)*ntypes*ntypes);
  D=malloc(sizeof(double)*ntypes*ntypes);
  template=malloc(sizeof(double)*ntypes*ntypes);

  // Get diagonal entries of A.
  for (n=0; n<=(ntypes-1); n++)
    *(as+n)=*(A+n*ntypes+n);

  // Find distinct diagonal entries.
  ndistinct=0; distinct=malloc(sizeof(double));
  memberships=malloc(sizeof(int)*ntypes);
  for (i=0; i<=(ntypes-1); i++)
    {
      j=0; k=-1; *(memberships+i)=-1;
      while ((j<=(ndistinct-1))&&(k<0))
	{
	  if (fabs(*(as+i)-*(distinct+j))<1e-10)
	    k=j;
	  j++;
	}
      if (k<0)
	{
	  ndistinct++; distinct=realloc(distinct,sizeof(double)*ndistinct);
	  *(distinct+ndistinct-1)=*(as+i); *(memberships+i)=ndistinct-1;
	}
      else
	*(memberships+i)=k;
    }

  /*
  // Debug
  printf("ndistinct=%d, distinct= ",ndistinct);
  for (i=0; i<=(ndistinct-1); i++)
    printf("%.3e ",*(distinct+i));
  printf("\n");
  */

  // Shift the diagonal entries to make them distinct.
  for (i=0; i<=(ndistinct-1); i++)
    {
      k=0;
      for (j=0; j<=(ntypes-1); j++)
	{
	  if (*(memberships+j)==i)
	    {
	      *(as+j)+=(shift*k); 

	      // Debug
	      //printf("i=%d, as[%d]=%.21e, k=%d\n",i,j,*(as+j),k);

	      k++;
	    }
	}
    }

  /*
  // Debug
  printf("as= ");
  for (i=0; i<=(ntypes-1); i++)
    printf("%.21e ",*(as+i));
  printf("\n");
  printf("a1-a3=%.3e\n",*(as+0)-*(as+2));
  */

  free(distinct); free(memberships);

  // Two drug case.
  if (ndrugs==2)
    {
      *(deltas+0)=*(A+1*ntypes+0);
      *(deltas+1)=*(A+2*ntypes+0);      
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ntypes-1); j++)
	    {
	      if (i<j)
		{
		  *(V+i*ntypes+j)=0;
		  *(iV+i*ntypes+j)=0;
		  *(D+i*ntypes+j)=0;
		}
	      else if (i==j)
		{
		  *(V+i*ntypes+j)=1;
		  *(iV+i*ntypes+j)=1;
		  *(D+i*ntypes+j)=*(as+i);
		}
	      else
		{
		  *(V+i*ntypes+j)=0;
		  *(iV+i*ntypes+j)=0;
		  *(D+i*ntypes+j)=0;
		}
	    }
	}
      *(V+1*ntypes+0)=*(deltas+0)/(*(as+0)-*(as+1));
      *(iV+1*ntypes+0)=*(V+1*ntypes+0)*(-1);
      *(V+2*ntypes+0)=*(deltas+1)/(*(as+0)-*(as+2));
      *(iV+2*ntypes+0)=*(V+2*ntypes+0)*(-1);
      *(V+3*ntypes+1)=*(deltas+1)/(*(as+1)-*(as+3));
      *(iV+3*ntypes+1)=*(V+3*ntypes+1)*(-1);
      *(V+3*ntypes+2)=*(deltas+0)/(*(as+2)-*(as+3));
      *(iV+3*ntypes+2)=*(V+3*ntypes+2)*(-1);
      val=*(deltas+0)*(*(deltas+1))/(*(as+0)-*(as+3));
      val*=(1/(*(as+0)-*(as+1))+1/(*(as+0)-*(as+2)));
      *(V+3*ntypes+0)=val;
      val=0;
      for (i=0; i<=(ntypes-2); i++)
	val+=*(V+(ntypes-1)*ntypes+i)*(*(iV+i*ntypes+0));
      *(iV+(ntypes-1)*ntypes+0)=-1*val;
    }

  // Three drug case.
  else if (ndrugs==3)
    {
      int vec1[ntypes], vec2[ntypes], vec3[ntypes];
      *(deltas+0)=*(A+1*ntypes+0);
      *(deltas+1)=*(A+2*ntypes+0);      
      *(deltas+2)=*(A+4*ntypes+0); 
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (j=0; j<=(ntypes-1); j++)
	    {
	      if (i<j)
		{
		  *(V+i*ntypes+j)=0;
		  *(iV+i*ntypes+j)=0;
		  *(D+i*ntypes+j)=0;
		}
	      else if (i==j)
		{
		  *(V+i*ntypes+j)=1;
		  *(iV+i*ntypes+j)=1;
		  *(D+i*ntypes+j)=*(as+i);
		}
	      else
		{
		  *(V+i*ntypes+j)=0;
		  *(iV+i*ntypes+j)=0;
		  *(D+i*ntypes+j)=0;
		}
	    }
	}
      for (i=0; i<=(ntypes-1); i++)
	{
	  for (k=0; k<=(ndrugs-1); k++)
	    vec1[k]=0;
	  k=i; l=ndrugs-1;
	  while (l>=0)
	    {
	      m=k%2; vec1[l]=m; k=(k-m)/2; l--;
	    }
	  for (j=0; j<=(ntypes-1); j++)
	    {
	      int ndiff, diff[ntypes];
	      for (k=0; k<=(ndrugs-1); k++)
		*(vec2+k)=0;
	      k=j; l=ndrugs-1;
	      while (l>=0)
		{
		  m=k%2; *(vec2+l)=m; k=(k-m)/2; l--;
		}
	      ndiff=0;
	      for (k=0; k<=(ndrugs-1); k++)
		{
		  if (vec1[k]!=vec2[k])
		    {
		      ndiff++; diff[ndiff-1]=k;
		    }
		}
	      if ((ndiff==1)&&(vec1[diff[0]]==0)&&(vec2[diff[0]]==1))
		{
		  k=ndrugs-1-diff[0]; *(V+j*ntypes+i)=*(deltas+k)/(*(as+i)-*(as+j));
		  *(iV+j*ntypes+i)=*(V+j*ntypes+i)*(-1);
		  if ((i!=0)&&(j!=7))
		    {
		      val=1;
		      for (l=0; l<=(ndrugs-1); l++)
			val*=*(deltas+l);
		      val/=(*(as+0)-*(as+i));
		      val/=(*(as+0)-*(as+j));
		      val/=(*(as+0)-*(as+7));
		      *(V+(ntypes-1)*ntypes+0)+=val;
		    }
		}
	      else if ((ndiff==2)&&(vec1[diff[0]]==0)&&(vec1[diff[1]]==0)&&(vec2[diff[0]]==1)&&(vec2[diff[1]]==1))
		{
		  k=ndrugs-1-diff[0]; l=ndrugs-1-diff[1];
		  val0=*(deltas+k)*(*(deltas+l));
		  for (m=0; m<=(ntypes-1); m++)
		    vec3[m]=vec1[m];
		  vec3[diff[0]]=1; n=0;
		  for (m=0; m<=(ndrugs-1); m++)
		    n=n*2+vec3[m];
		  val=val0;
		  val/=(*(as+i)-*(as+n));
		  val/=(*(as+i)-*(as+j));
		  *(V+j*ntypes+i)+=val;
		  *(iV+j*ntypes+i)-=val;
		  val=val0;
		  val/=(*(as+i)-*(as+n));
		  val/=(*(as+n)-*(as+j));
		  *(iV+j*ntypes+i)+=val;
		  for (m=0; m<=(ntypes-1); m++)
		    vec3[m]=vec1[m];
		  vec3[diff[1]]=1; n=0;
		  for (m=0; m<=(ndrugs-1); m++)
		    n=n*2+vec3[m];
		  val=val0;
		  val/=(*(as+i)-*(as+n));
		  val/=(*(as+i)-*(as+j));
		  *(V+j*ntypes+i)+=val;
		  *(iV+j*ntypes+i)-=val;
		  val=val0;
		  val/=(*(as+i)-*(as+n));
		  val/=(*(as+n)-*(as+j));
		  *(iV+j*ntypes+i)+=val;
		}
	    }
	}
    
      val=0;
      for (i=0; i<=(ntypes-2); i++)
	val+=*(V+(ntypes-1)*ntypes+i)*(*(iV+i*ntypes+0));
      *(iV+(ntypes-1)*ntypes+0)=-1*val;
    }

  /*
  // Debug
  printf("V=\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(V+i*ntypes+j));
      printf("\n");
    }
  printf("D=\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(D+i*ntypes+j));
      printf("\n");
    }
  printf("iV=\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(iV+i*ntypes+j));
      printf("\n");
    }
  */


  // Apply matrix exponentiation.
  // V*expm(D)*iV.
  for (i=0; i<=(ntypes-1); i++)
    *(D+i*ntypes+i)=exp(*(D+i*ntypes+i));
  mat_multiply(ntypes,ntypes,ntypes,V,D,template);
  mat_multiply(ntypes,ntypes,ntypes,template,iV,expA);

  // Release memory.
  free(as); free(deltas); free(V); free(iV); free(D); free(template);

  /*
  // Debug
  printf("A=\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(A+i*ntypes+j));
      printf("\n");
    }
  printf("expA=\n");
  for (i=0; i<=(ntypes-1); i++)
    {
      for (j=0; j<=(ntypes-1); j++)
	printf("%.3e ",*(expA+i*ntypes+j));
      printf("\n");
    }
  exit(0);
  */

  return;
}
  
