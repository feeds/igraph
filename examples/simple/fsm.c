#include <igraph.h>
#include <stdio.h>

int match_rings_induced() {
  igraph_t ro, rc;
  igraph_bool_t iso;
  igraph_integer_t supp;
  igraph_ring(&ro, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 0);
  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  printf("INDUCED\n");

  igraph_shallow_support(&ro, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("sh ro-ro %d\n", supp);

  igraph_shallow_support(&rc, &rc, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("sh rc-rc %d\n", supp);

  igraph_shallow_support(&rc, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("sh rc-ro %d\n\n", supp);

  igraph_mib_support(&ro, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mib ro-ro %d\n", supp);

  igraph_mib_support(&rc, &rc, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mib rc-rc %d\n", supp);

  igraph_mib_support(&rc, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mib rc-ro %d\n\n", supp);

  igraph_mib_support_slow(&ro, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mibS ro-ro %d\n", supp);

  igraph_mib_support_slow(&rc, &rc, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mibS rc-rc %d\n", supp);

  igraph_mib_support_slow(&rc, &ro, 0, 0, 0, 0, 0, 0, 1, &supp);
  printf("mibS rc-ro %d\n\n", supp);

  igraph_destroy(&ro);
  igraph_destroy(&rc);
  return 0;
}

int match_rings_noninduced() {
  igraph_t ro, rc;
  igraph_bool_t iso;
  igraph_integer_t supp;
  igraph_ring(&ro, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 0);
  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  printf("NON-INDUCED\n");

  igraph_shallow_support(&ro, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh ro-ro %d\n", supp);

  igraph_shallow_support(&rc, &rc, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh rc-rc %d\n", supp);

  igraph_shallow_support(&rc, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("sh rc-ro %d\n\n", supp);

  igraph_mib_support(&ro, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib ro-ro %d\n", supp);

  igraph_mib_support(&rc, &rc, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib rc-rc %d\n", supp);

  igraph_mib_support(&rc, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mib rc-ro %d\n\n", supp);

  igraph_mib_support_slow(&ro, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibS ro-ro %d\n", supp);

  igraph_mib_support_slow(&rc, &rc, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibS rc-rc %d\n", supp);

  igraph_mib_support_slow(&rc, &ro, 0, 0, 0, 0, 0, 0, 0, &supp);
  printf("mibS rc-ro %d\n\n", supp);

  igraph_destroy(&ro);
  igraph_destroy(&rc);
  return 0;
}

/* ----------------------------------------------------------- */

int main() {
  match_rings_induced();
  match_rings_noninduced();
  return 0;
}

