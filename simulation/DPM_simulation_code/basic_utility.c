// Basic utility functions.

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include "db.h"
# define MAXLINEWIDTH 100000
# define ENTRYWIDTH 50
# define PHYLIPWIDTH 10
# define SPECIESNAMELENGTH 50
# define max(i,j) (i>j ? i : j)
# define argmax(i,j) (j>i ? 2 : 1)

int inlist(int Nv, int *v, int target);
char * getitemval2(char *s, int ind, int *length);
char * getitemval3(char *s, int ind, int *length, char sepch);
char * getitemval4(char *s, int ind, int *length, char sepch);
double atof2(char *s, int length);
int * setdiff(int NA, int *A, int NB, int *B, int *NC);
int nu2ind(char ch);
int aa2ind(char ch);
int pair_cmp_ascend(struct pair *p1, struct pair *p2);
int pair_cmp_descend(struct pair *p1, struct pair *p2);
int triplet_cmp_ascend(struct triplet *p1, struct triplet *p2);
int triplet_cmp_descend(struct triplet *p1, struct triplet *p2);
double logfactorial(int n);
double logchoose(int n, int k);
int label_conn_comps(int nnodes1, int nnodes2, int *graph, int *conn1, int *conn2);
void recurse_label_comps(int nnodes1, int nnodes2, int *graph, int *conn1, int *conn2, int type, int target);
double mean(int n, double *v);
double variance(int n, double *v);
double corr_coef(int n, double *v1, double *v2);

extern int nchar, njointchar;
extern char *numap, *aamap;

// Check if a number is on a list.
int inlist(int Nv, int *v, int target)
{
  int i=0, indicator=0;

  while ((i<Nv)&&(indicator==0))
    {
      if (*(v+i)==target)
	indicator=1;
      i++;
    }

  return indicator;
}

// Get the designated item in a string.
// Difference from getitemval: items can be separated by tab or space.
char * getitemval2(char *s, int ind, int *length)
{
  int i, tabcnt, curind, cnt;
  char ch, *item, *template;

  template=malloc(sizeof(char)*MAXLINEWIDTH); tabcnt=0; curind=0;

  while ((tabcnt<ind)&&(curind<strlen(s)))
    {
      ch=*(s+curind);

      if ((ch==0x09)||(ch==' '))
	tabcnt++;
      // If the line ends before tabcnt reaching ind, then return null.
      else if ((ch==0x0a)&&(tabcnt<ind))
	tabcnt=ind+1;
      curind++;
    }

  if (tabcnt==ind)
    {      
      cnt=0;
      while (((ch=*(s+curind))!=0x09)&&(ch!=0x0a)&&(ch!=0x0d)&&(ch!=' '))
	{
	  *(template+cnt)=ch;
	  curind++; cnt++;
	}
    
      if (cnt>0)
	{
	  //item=malloc(sizeof(char)*cnt);
	  item=malloc(sizeof(char)*(cnt+1));
	  for (i=0; i<=(cnt-1); i++)
	    *(item+i)=*(template+i);

	  *(item+cnt)='\0';

	  *length=cnt;
	}
      else
	*length=cnt;
           
    }
  else
    {
      item=NULL; *length=0;
    }

  free(template);
  return item;
}

// Get the designated item in a string.
// Difference from getitemval2: specify the separator character.
char * getitemval3(char *s, int ind, int *length, char sepch)
{
  int i, tabcnt, curind, cnt;
  char ch, *item, *template;

  template=malloc(sizeof(char)*MAXLINEWIDTH); tabcnt=0; curind=0;

  while ((tabcnt<ind)&&(curind<strlen(s)))
    {
      ch=*(s+curind);

      if (ch==sepch)
	tabcnt++;
      // If the line ends before tabcnt reaching ind, then return null.
      else if ((ch==0x0a)&&(tabcnt<ind))
	tabcnt=ind+1;
      curind++;
    }

  if (tabcnt==ind)
    {
      
      cnt=0;
      while (((ch=*(s+curind))!=sepch)&&(ch!=0x0a)&&(ch!=0x0d))
	{
	  *(template+cnt)=ch;
	  curind++; cnt++;
	}
    
      if (cnt>0)
	{
	  item=malloc(sizeof(char)*cnt);
	  for (i=0; i<=(cnt-1); i++)
	    *(item+i)=*(template+i);
	  *length=cnt;
	}
      else
	*length=cnt;
           
    }
  else
    {
      item=NULL; *length=0;
    }

  free(template);
  
  return item;
}

// Convert a string into a real number.
// Also works for exponential representation.
double atof2(char *s, int length)
{
  double val, expval, power;
  int i, j, cnt, sign, afterpoint, afterexp, expsign;
  char ch;

  val=0; sign=1; cnt=0; afterpoint=0; afterexp=0; 
  expval=0, expsign=1; power=1.0;
  //while (cnt<strlen(s))
  while (cnt<length)
    {
      ch=*(s+cnt);
      if ((cnt==0)&&(ch=='-'))
	sign=-1;
      else if ((ch=='-')&&(afterexp==1))
	expsign=-1;
      else if (ch=='.')
	{
	  afterpoint=1;
	  power=0.1;
	}
      else if ((ch=='e')||(ch=='E'))
	{	
	  afterexp=1;
	  afterpoint=0;
	}
      else if ((ch>='0')&&(ch<='9')&&(afterpoint==0))
	{
	  if (afterexp==0)
	    val=val*10+(ch-'0');
	  else
	    expval=expval*10+(ch-'0');
	}
      else if ((ch>='0')&&(ch<='9')&&(afterpoint==1))
	{
	  if (afterexp==0)
	    val=val+(ch-'0')*power;
	  else
	    expval=expval+(ch-'0')*power;
	  power=power*0.1;
	}
  
      cnt++;
    }

  
  if (sign==-1)
    val=val*(-1);
  if (expsign==-1)
    expval=expval*(-1);

  val=val*pow(10,expval);

  return val;
}

// Obtain a set difference A\B.
int * setdiff(int NA, int *A, int NB, int *B, int *NC)
{
  int i, j, k, *C;
 
  *NC=0; C=malloc(sizeof(int)*max(NA,NB));
  for (i=0; i<=(NA-1); i++)
    {
      j=0; k=1;
      while ((j<=(NB-1))&&(k==1))
	{
	  if (*(A+i)==*(B+j))
	    k=0;
	  j++;
	}
      if (k==1)
	{
	  *NC+=1; *(C+*NC-1)=*(A+i);
	}
    }

  return C;
}

// Convert a nucleic acid character into an index number.
int nu2ind(char ch)
{
  int ind;
  switch (ch) {
    case 'A': ind=1; break;
    case 'C': ind=2; break;
    case 'G': ind=3; break;
    case 'U': ind=4; break;
    case 'T': ind=4; break;
    default: ind=0; break;
  }
  return ind;
}

// Convert an amino acid character into an index number.
int aa2ind(char ch)
{
  int ind;
  switch (ch) {
    case 'A': ind=1; break;
    case 'R': ind=2; break;
    case 'N': ind=3; break;
    case 'D': ind=4; break;
    case 'C': ind=5; break;
    case 'Q': ind=6; break;
    case 'E': ind=7; break;
    case 'G': ind=8; break;
    case 'H': ind=9; break;
    case 'I': ind=10; break;
    case 'L': ind=11; break;
    case 'K': ind=12; break;
    case 'M': ind=13; break;
    case 'F': ind=14; break;
    case 'P': ind=15; break;
    case 'S': ind=16; break;
    case 'T': ind=17; break;
    case 'W': ind=18; break;
    case 'Y': ind=19; break;
    case 'V': ind=20; break;
    default: ind=0; break;
  }

  return ind;
}


// Compare two pair structures in an ascending order.
int pair_cmp_ascend(struct pair *p1, struct pair *p2)
{
  int retval=0;

  if (p1->s<p2->s)
    retval=-1;
  else if (p1->s>p2->s)
    retval=1;
  else
    retval=0;
  return retval;
}

// Compare two pair structures in a descending order.
int pair_cmp_descend(struct pair *p1, struct pair *p2)
{
  int retval=0;

  if (p1->s<p2->s)
    retval=1;
  else if (p1->s>p2->s)
    retval=-1;
  else
    retval=0;
  return retval;
}

// Compare two triplet structures in an ascending order.
int triplet_cmp_ascend(struct triplet *p1, struct triplet *p2)
{
  int retval=0;

  if (p1->s<p2->s)
    retval=-1;
  else if (p1->s>p2->s)
    retval=1;
  else
    retval=0;
  return retval;
}

// Compare two triplet structures in a descending order.
int triplet_cmp_descend(struct triplet *p1, struct triplet *p2)
{
  int retval=0;

  if (p1->s<p2->s)
    retval=1;
  else if (p1->s>p2->s)
    retval=-1;
  else
    retval=0;
  return retval;
}


// log of a factorial.
double logfactorial(int n)
{
  int i;
  double val;

  if (n<=1)
    return 0.0;

  val=0;
  for (i=1; i<=n; i++)
    val+=log(i);
  return val;
}

// log of a choose function C(n,k).
double logchoose(int n, int k)
{
  double val;

  if (k>=n)
    return 0;

  val=logfactorial(n)-logfactorial(k)-logfactorial(n-k);
  return val;
}


// Find connected components of a graph.
// Label nodes in each component with a different number.
// Only consider the non-orphan components.
int label_conn_comps(int nnodes1, int nnodes2, int *graph, int *conn1, int *conn2)
{
  int i, j, k, l, n, flag, nconn;
  for (i=0; i<=(nnodes1-1); i++)
    *(conn1+i)=-1;
  for (i=0; i<=(nnodes2-1); i++)
    *(conn2+i)=-1;
  flag=1; nconn=0;
  while (flag==1)
    {
      // Pick up an unlabeled node.
      i=0; j=-1;
      while ((i<nnodes1)&&(j==-1))
	{
	  if (*(conn1+i)==-1)
	    j=i;
	  i++;
	}

      // If such node exists then recurivsely label the component.
      if (j>=0)
	{
	  nconn++; *(conn1+j)=nconn-1;
	  recurse_label_comps(nnodes1,nnodes2,graph,conn1,conn2,0,j);
	}
      else
	flag=0;
    }

  // Count the number of non-orphan components.
  n=0;
  for (i=0; i<=nconn; i++)
    {
      k=0;
      for (j=0; j<=(nnodes1-1); j++)
	if (*(conn1+j)==i)
	  k++;
      for (j=0; j<=(nnodes2-1); j++)
	if (*(conn2+j)==i)
	  k++;
      if (k>1)
	n++;
    }


  //return nconn;
  return n;
}


// Recursively label the nodes in a component.
void recurse_label_comps(int nnodes1, int nnodes2, int *graph, int *conn1, int *conn2, int type, int target)
{
  int i, j, k, l;
  
  // Label neighbors which are unlabeled.
  // Target node is type 0.
  if (type==0)
    {
      for (i=0; i<=(nnodes2-1); i++)
	{
	  //if ((i!=target)&&((*(graph+i*nnodes+target)>0)||(*(graph+target*nnodes+i)>0))&&(*(conn+i)==-1))
	  //if ((i!=target)&&(*(graph+i*nnodes+target)>0)&&(*(conn+i)==-1))
	  if ((*(graph+target*nnodes2+i)>0)&&(*(conn2+i)==-1))
	    {
	      *(conn2+i)=*(conn1+target);
	      recurse_label_comps(nnodes1,nnodes2,graph,conn1,conn2,1,i);
	    }
	}            
    }
  // Target node is type 1.
  else if (type==1)
    {
      for (i=0; i<=(nnodes1-1); i++)
	{
	  //if ((i!=target)&&((*(graph+i*nnodes+target)>0)||(*(graph+target*nnodes+i)>0))&&(*(conn+i)==-1))
	  //if ((i!=target)&&(*(graph+i*nnodes+target)>0)&&(*(conn+i)==-1))
	  if ((*(graph+i*nnodes2+target)>0)&&(*(conn1+i)==-1))
	    {
	      *(conn1+i)=*(conn2+target);
	      recurse_label_comps(nnodes1,nnodes2,graph,conn1,conn2,0,i);
	    }
	}            
    }
  return;
}

// Mean of a vector.
double mean(int n, double *v)
{
  int i;
  double retval;
  
  retval=0;
  for (i=0; i<=(n-1); i++)
    retval+=*(v+i);
  retval/=(double)n;
  return retval;
}

// Variance of a vector.
double variance(int n, double *v)
{
  int i;
  double retval, m, val;

  retval=0; m=mean(n,v);
  for (i=0; i<=(n-1); i++)
    {
      val=(*(v+i)-m)*(*(v+i)-m);
      retval+=val;
    }
  retval/=(double)n;
  return retval;
}


// Calculate the correlation coefficient of two vectors.
double corr_coef(int n, double *v1, double *v2)
{
  int i, j, k;
  double val, retval, m1, m2, var1, var2;
  
  retval=0; m1=mean(n,v1); m2=mean(n,v2);
  var1=variance(n,v1); var2=variance(n,v2);
  
  for (i=0; i<=(n-1); i++)
    {
      val=(*(v1+i)-m1)*(*(v2+i)-m2);
      retval+=val;
    }
  retval/=(double)n;

  retval/=sqrt(var1*var2);
  return retval;
}

// Get the designated item in a string.
// Difference from getitemval2: specify the separator character.
// Difference from getitemval3: can handle the last entry in a line.
char * getitemval4(char *s, int ind, int *length, char sepch)
{
  int i, tabcnt, curind, cnt;
  char ch, *item, *template;

  template=malloc(sizeof(char)*MAXLINEWIDTH); tabcnt=0; curind=0;

  while ((tabcnt<ind)&&(curind<strlen(s)))
    {
      ch=*(s+curind);

      if (ch==sepch)
	tabcnt++;
      // If the line ends before tabcnt reaching ind, then return null.
      else if ((ch==0x0a)&&(tabcnt<ind))
	tabcnt=ind+1;
      curind++;
    }

  if (tabcnt==ind)
    {
      
      cnt=0;
      while (((ch=*(s+curind))!=sepch)&&(ch!=0x0a)&&(ch!=0x0d))
	{
	  *(template+cnt)=ch;
	  curind++; cnt++;
	}
    
      if (cnt>0)
	{
	  item=malloc(sizeof(char)*cnt);
	  for (i=0; i<=(cnt-1); i++)
	    *(item+i)=*(template+i);
	  *length=cnt;
	}
      else
	*length=cnt;
           
    }
  else
    {
      item=NULL; *length=0;
    }

  if (*length>0)
    {
      if (((ch=*(item+*length-1))==0x0a)||(ch=='\n'))
	*length-=1;
    }      

  free(template);
  
  return item;
}
