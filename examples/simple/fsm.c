#include <igraph.h>
#include <stdio.h>

int match_rings() {

  igraph_t ro, rc;
  igraph_bool_t iso;
  igraph_integer_t supp;
  igraph_ring(&ro, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 0);
  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  igraph_shallow_support(&ro, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh %d\n", supp);

  igraph_shallow_support(&rc, &rc, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh %d\n", supp);

  igraph_shallow_support(&rc, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh %d\n\n", supp);

  igraph_mib_support(&ro, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib %d\n", supp);

  igraph_mib_support(&rc, &rc, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib %d\n", supp);

  igraph_mib_support(&rc, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib %d\n\n", supp);

  igraph_mib_support_fast(&ro, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibF %d\n", supp);

  igraph_mib_support_fast(&rc, &rc, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibF %d\n", supp);

  igraph_mib_support_fast(&rc, &ro, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibF %d\n\n", supp);

  igraph_destroy(&ro);
  igraph_destroy(&rc);
  return 0;
}

/* ----------------------------------------------------------- */

int main() {
  match_rings();
  return 0;
}

