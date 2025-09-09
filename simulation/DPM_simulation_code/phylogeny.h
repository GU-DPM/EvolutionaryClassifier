// Include data structures for phylogeny_utility.c

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# define MAXLINEWIDTH 100000
# define ENTRYWIDTH 50
# define PHYLIPWIDTH 10
# define SPECIESNAMELENGTH 20


// A cellstring structure.
struct cellstr {
  int length;  // string length.
  char *s;     // string.
};

// A phylogenetic node structure.
struct phylogenetic_node {
  int ind;     // node index.
  int parent;  // parent node index.
  double bl;   // length of the branch from the current node's parent to the current node.
  struct cellstr *name;  // species name of the leaf node. 
  int nchild;  // number of children.
  int *child;  // children indices.
};

// A data structure of a real number in exponential (10-base) form.
struct exp_real {
  double mantissa;  // the leading term of the number.
  int exponent;     // the 10-base exponential term.
};

