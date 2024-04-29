// Include data structures for db_utility.c

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# define MAXLINEWIDTH 100000
# define NAMELENGTH 100
# define LONGNAMELENGTH 100
# define LINEWIDTH 50

struct protein_species_info {
  int ind;
  char *proteinaccs;
  char *speciesnames;
};

struct protein_position_info {
  int ind;
  int fileind;
  char *proteinaccs;
  long position;
};

struct index_node {
  int ind;
  int parent;
  int Nchildren;
  int *children;
  char indchar;
  int start;
};

struct hierarchy {
  char *speciesname;
  int nmembers;
  int *inds;
  int *ranges;
  char *membernames;
  char *memberswissprotaccs;
  char *memberseqs;
};

struct pair {
  int i;
  int j;
  double s;
};

struct triplet {
  int i;
  int j;
  int k;
  double s;
};

struct seqnode {
  int index;
  int depth;
  int pa;
  int nchildren;
  long *children;
  double *dosage;
  double *y;
};

struct variable {
  char *name;
  int nvalues;
  double *values;
};

struct pair2 {
  int i;
  double s;
};

struct sequence {
  int depth;
  double *dosages;
  double *ys;
};

