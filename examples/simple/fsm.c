#include <igraph.h>
#include <stdio.h>

int match_rings_induced() {
  igraph_t ro, rc;
  igraph_bool_t iso;
  igraph_integer_t supp;
  igraph_ring(&ro, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 0);
  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  printf("INDUCED\n");

  igraph_shallow_support(&ro, &ro, 0, 0, 0, 0, 1, &supp, -1);
  printf("sh ro-ro %d\n", supp);

  igraph_shallow_support(&rc, &rc, 0, 0, 0, 0, 1, &supp, -1);
  printf("sh rc-rc %d\n", supp);

  igraph_shallow_support(&rc, &ro, 0, 0, 0, 0, 1, &supp, -1);
  printf("sh rc-ro %d\n\n", supp);

  igraph_mib_support(&ro, &ro, 0, 0, 0, 0, 1, &supp, -1);
  printf("mib ro-ro %d\n", supp);

  igraph_mib_support(&rc, &rc, 0, 0, 0, 0, 1, &supp, -1);
  printf("mib rc-rc %d\n", supp);

  igraph_mib_support(&rc, &ro, 0, 0, 0, 0, 1, &supp, -1);
  printf("mib rc-ro %d\n\n", supp);

  igraph_mib_support_slow(&ro, &ro, 0, 0, 0, 0, 1, &supp, -1);
  printf("mibS ro-ro %d\n", supp);

  igraph_mib_support_slow(&rc, &rc, 0, 0, 0, 0, 1, &supp, -1);
  printf("mibS rc-rc %d\n", supp);

  igraph_mib_support_slow(&rc, &ro, 0, 0, 0, 0, 1, &supp, -1);
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

  igraph_shallow_support(&ro, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("sh ro-ro %d\n", supp);

  igraph_shallow_support(&rc, &rc, 0, 0, 0, 0, 0, &supp, -1);
  printf("sh rc-rc %d\n", supp);

  igraph_shallow_support(&rc, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("sh rc-ro %d\n\n", supp);

  igraph_mib_support(&ro, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("mib ro-ro %d\n", supp);

  igraph_mib_support(&rc, &rc, 0, 0, 0, 0, 0, &supp, -1);
  printf("mib rc-rc %d\n", supp);

  igraph_mib_support(&rc, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("mib rc-ro %d\n\n", supp);

  igraph_mib_support_slow(&ro, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("mibS ro-ro %d\n", supp);

  igraph_mib_support_slow(&rc, &rc, 0, 0, 0, 0, 0, &supp, -1);
  printf("mibS rc-rc %d\n", supp);

  igraph_mib_support_slow(&rc, &ro, 0, 0, 0, 0, 0, &supp, -1);
  printf("mibS rc-ro %d\n\n", supp);

  igraph_destroy(&ro);
  igraph_destroy(&rc);
  return 0;
}

int gspan() {
  igraph_t rc;
  igraph_vector_int_t vcolor, ecolor;

  igraph_vector_ptr_t graphs;
  igraph_vector_ptr_t vertex_colors;
  igraph_vector_ptr_t edge_colors;

  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  igraph_vector_int_init(&vcolor, 10);
  igraph_vector_int_fill(&vcolor, 5);
  VECTOR(vcolor)[0] = 1;
  VECTOR(vcolor)[1] = 1;
  VECTOR(vcolor)[2] = 3;
  VECTOR(vcolor)[3] = 4;

  igraph_vector_int_init(&ecolor, 10);
  igraph_vector_int_fill(&ecolor, 2);
  VECTOR(ecolor)[6] = 7;
  VECTOR(ecolor)[7] = 7;
  VECTOR(ecolor)[8] = 8;

  igraph_vector_ptr_init(&graphs, 1);
  igraph_vector_ptr_init(&vertex_colors, 1);
  igraph_vector_ptr_init(&edge_colors, 1);

  VECTOR(graphs)[0] = &rc;
  VECTOR(vertex_colors)[0] = &vcolor;
  VECTOR(edge_colors)[0] = &ecolor;

  printf("GSPAN\n\n");

  printf("no labels\n");
  igraph_gspan(&graphs, NULL, NULL, &igraph_db_mib_support, /*min_supp=*/ 1, NULL, NULL);

  printf("\nvertex labels\n");
  igraph_gspan(&graphs, &vertex_colors, NULL,
               &igraph_db_mib_support, /*min_supp=*/ 1, NULL, NULL);

  printf("\nedge labels\n");
  igraph_gspan(&graphs, NULL, &edge_colors,
               &igraph_db_mib_support, /*min_supp=*/ 1, NULL, NULL);

  printf("\nvertex and edge labels\n");
  igraph_gspan(&graphs, &vertex_colors, &edge_colors,
               &igraph_db_mib_support, /*min_supp=*/ 1, NULL, NULL);

  igraph_destroy(&rc);
  igraph_vector_ptr_destroy(&graphs);
  igraph_vector_ptr_destroy(&vertex_colors);
  igraph_vector_ptr_destroy(&edge_colors);

  return 0;
}

/* ----------------------------------------------------------- */

int main() {
  //match_rings_induced();
  //match_rings_noninduced();
  gspan();
  return 0;
}

