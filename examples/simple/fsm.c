#include <igraph.h>
#include <stdio.h>

void igraph_i_print(const igraph_t *g, const igraph_vector_int_t *vcolors,
                    const igraph_vector_int_t *ecolors) {
  long int i;
  if (vcolors != NULL) {
    if (ecolors != NULL) {
      for (i = 0; i < igraph_ecount(g); i++) {
        printf("%ld(%d) --%d-- %ld(%d)\n", (long int) VECTOR(g->from)[i],
                                           VECTOR(*vcolors)[(long int) VECTOR(g->from)[i]],
                                           VECTOR(*ecolors)[i],
                                           (long int) VECTOR(g->to)[i],
                                           VECTOR(*vcolors)[(long int) VECTOR(g->to)[i]]);
      }
    } else {
      for (i = 0; i < igraph_ecount(g); i++) {
        printf("%ld(%d) -- %ld(%d)\n", (long int) VECTOR(g->from)[i],
                                       VECTOR(*vcolors)[(long int) VECTOR(g->from)[i]],
                                       (long int) VECTOR(g->to)[i],
                                       VECTOR(*vcolors)[(long int) VECTOR(g->to)[i]]);
      }
    }
  } else {
    if (ecolors != NULL) {
      for (i = 0; i < igraph_ecount(g); i++) {
        printf("%ld --%d-- %ld\n", (long int) VECTOR(g->from)[i],
                                   VECTOR(*ecolors)[i],
                                   (long int) VECTOR(g->to)[i]);
      }
    } else {
      for (i = 0; i < igraph_ecount(g); i++) {
        printf("%ld -- %ld\n", (long int) VECTOR(g->from)[i], (long int) VECTOR(g->to)[i]);
      }
    }
  }
}

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
  igraph_vector_int_t vcolor, ecolor, result_supps;
  igraph_vector_ptr_t graphs, result_graphs;
  igraph_vector_ptr_t vertex_colors, result_vertex_colors;
  igraph_vector_ptr_t edge_colors, result_edge_colors;
  long int i;

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

  igraph_vector_ptr_init(&result_graphs, 0);
  igraph_vector_ptr_init(&result_vertex_colors, 0);
  igraph_vector_ptr_init(&result_edge_colors, 0);
  igraph_vector_int_init(&result_supps, 0);

  printf("no labels\n\n");
  igraph_gspan(&graphs, NULL, NULL, &igraph_db_mib_support, /*min_supp=*/ 1, /*max_edges=*/ 5,
      &result_graphs, NULL, NULL, &result_supps);
  for (i = 0; i < igraph_vector_ptr_size(&result_graphs); i++) {
    printf("supp=%ld\n", (long int) VECTOR(result_supps)[i]);
    igraph_i_print((igraph_t *) VECTOR(result_graphs)[i], NULL, NULL);
    printf("\n");
  }

  printf("\nvertex labels\n\n");
  igraph_gspan(&graphs, &vertex_colors, NULL,
               &igraph_db_mib_support, /*min_supp=*/ 1, /*max_edges=*/ 5,
               &result_graphs, &result_vertex_colors, NULL, &result_supps);
  for (i = 0; i < igraph_vector_ptr_size(&result_graphs); i++) {
    printf("supp=%ld\n", (long int) VECTOR(result_supps)[i]);
    igraph_i_print((igraph_t *) VECTOR(result_graphs)[i],
                   (igraph_vector_int_t *) VECTOR(result_vertex_colors)[i], NULL);
    printf("\n");
  }

  printf("\nedge labels\n\n");
  igraph_gspan(&graphs, NULL, &edge_colors,
               &igraph_db_mib_support, /*min_supp=*/ 1, /*max_edges=*/ 5,
               &result_graphs, NULL, &result_edge_colors, &result_supps);
  for (i = 0; i < igraph_vector_ptr_size(&result_graphs); i++) {
    printf("supp=%ld\n", (long int) VECTOR(result_supps)[i]);
    igraph_i_print((igraph_t *) VECTOR(result_graphs)[i], NULL,
                   (igraph_vector_int_t *) VECTOR(result_edge_colors)[i]);
    printf("\n");
  }

  printf("\nvertex and edge labels\n\n");
  igraph_gspan(&graphs, &vertex_colors, &edge_colors,
               &igraph_db_mib_support, /*min_supp=*/ 1, /*max_edges=*/ 5,
               &result_graphs, &result_vertex_colors, &result_edge_colors, &result_supps);
  for (i = 0; i < igraph_vector_ptr_size(&result_graphs); i++) {
    printf("supp=%ld\n", (long int) VECTOR(result_supps)[i]);
    igraph_i_print((igraph_t *) VECTOR(result_graphs)[i],
                   (igraph_vector_int_t *) VECTOR(result_vertex_colors)[i],
                   (igraph_vector_int_t *) VECTOR(result_edge_colors)[i]);
    printf("\n");
  }

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

