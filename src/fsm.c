/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   IGraph library AcGM implementation.
   Copyright (C) 2015  Erik Scharwaechter <erik.scharwaechter@rwth-aachen.de>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <stdlib.h>
#include "igraph_fsm.h"
#include "igraph_spmatrix.h"
#include "igraph_interface.h"

/* graph1 is the larger graph, graph2 is the smaller graph */
int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_isocompat_t *node_compat_fn,
			   igraph_isocompat_t *edge_compat_fn,
			   igraph_integer_t *support) {
  igraph_bool_t iso;
  if (igraph_subisomorphic_vf2(graph1, graph2, vertex_color1, vertex_color2,
		edge_color1, edge_color2, &iso, NULL, NULL, node_compat_fn,
		edge_compat_fn, NULL)) {
    return 1;
  }
  if (iso) {
    *support = 1;
  } else {
    *support = 0;
  }
  return 0;
}


igraph_bool_t igraph_i_mib_isohandler(const igraph_vector_t *map12,
				      const igraph_vector_t *map21, void *arg) {
  igraph_spmatrix_t *target_hits = (igraph_spmatrix_t *) arg;
  long int vcount2 = igraph_vector_size(map21);
  long int i;
  for (i = 0; i < vcount2; i++) {
    igraph_spmatrix_set(target_hits, i, VECTOR(*map21)[i], 1);
  }
  return 1;
}


/* graph1 is the larger graph, graph2 is the smaller graph */
int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_isocompat_t *node_compat_fn,
		       igraph_isocompat_t *edge_compat_fn,
		       igraph_integer_t *support) {
  igraph_vector_t map21, target_counts;
  igraph_spmatrix_t target_hits;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);

  /* naive implementation: iterates over all embeddings and collects target nodes */
  igraph_vector_init(&map21, 0);
  igraph_spmatrix_init(&target_hits, vcount2, vcount1);
  if (igraph_subisomorphic_function_vf2(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
		edge_color2, NULL, &map21, (igraph_isohandler_t *) igraph_i_mib_isohandler,
		node_compat_fn, edge_compat_fn, (void *) &target_hits)) {
    igraph_spmatrix_destroy(&target_hits);
    igraph_vector_destroy(&map21);
    return 1;
  }

  igraph_vector_init(&target_counts, vcount2);
  igraph_spmatrix_rowsums(&target_hits, &target_counts);
  *support = igraph_vector_min(&target_counts);

  igraph_vector_destroy(&target_counts);
  igraph_spmatrix_destroy(&target_hits);
  igraph_vector_destroy(&map21);
  return 0;
}

int igraph_acgm(const igraph_vector_ptr_t *graphdb, igraph_support_measure_t *supp_fn,
		igraph_real_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values) {
  return 0;
}

