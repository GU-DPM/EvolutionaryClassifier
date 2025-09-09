// An integer set.
struct intset {
  int nmembers; // number of members.
  int *members; // members.
};
  
// A data structure of phylogenetic tree information.
// Replace double pointers with the intset structure.
struct treeinfo2 {
  int nnodes;  // number of nodes.
  struct intset *child; // children.  
  double *t;   // branch length.
  int *order;  // order.
  int *pa;     // parents.
};

