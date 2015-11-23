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
#include "igraph_matrix.h"
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
		edge_color1, edge_color2, 1, &iso, NULL, NULL, node_compat_fn,
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
  igraph_matrix_t *target_hits = (igraph_matrix_t *) arg;
  long int vcount2 = igraph_vector_size(map21);
  long int i;
  for (i = 0; i < vcount2; i++) {
    igraph_matrix_set(target_hits, i, VECTOR(*map21)[i], 1.);
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
  igraph_matrix_t target_hits;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);

  /* naive implementation: iterates over all embeddings and collects target nodes */
  igraph_vector_init(&map21, 0);
  igraph_matrix_init(&target_hits, vcount2, vcount1);
  igraph_matrix_null(&target_hits);
  if (igraph_subisomorphic_function_vf2(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
		edge_color2, 1, NULL, &map21, (igraph_isohandler_t *) igraph_i_mib_isohandler,
		node_compat_fn, edge_compat_fn, (void *) &target_hits)) {
    igraph_matrix_destroy(&target_hits);
    igraph_vector_destroy(&map21);
    return 1;
  }

  igraph_vector_init(&target_counts, vcount2);
  igraph_matrix_rowsum(&target_hits, &target_counts);
  *support = igraph_vector_min(&target_counts);

  igraph_vector_destroy(&target_counts);
  igraph_matrix_destroy(&target_hits);
  igraph_vector_destroy(&map21);
  return 0;
}


// graph1 is the larger graph, graph2 is the smaller graph
// map21 may contain fixed assignments, only entries set to -1 are filled
// no range check and no check for duplicate fixed assignments is done
int igraph_i_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_bool_t induced,
			   igraph_vector_t *map21,
			   igraph_bool_t *iso) {
  long int vcount1 = igraph_vcount(graph1);
  long int vcount2 = igraph_vcount(graph2);
  long int fixed_count = 0;
  long int i, j, k, stack_pos;
  int end = 0, success = 1;
  long int pattern_node, other_pattern_node, target_node, other_target_node;
  igraph_vector_t indeg1, indeg2, outdeg1, outdeg2, inneighs1, inneighs2, outneighs1, outneighs2;
  igraph_vector_t node_ordering;
  igraph_vector_t partial_solution_stack;
  igraph_integer_t eid1, eid2;

  *iso = 0;

  // create a static ordering of the pattern nodes: fixed assignments come first,
  // followed by all other nodes in arbitrary order (here: reversed for simplicity)
  IGRAPH_CHECK(igraph_vector_init(&node_ordering, vcount2));
  if (map21 == NULL) {
    // no fixed nodes specified
    for (i = 0; i < vcount2; i++) {
      VECTOR(node_ordering)[i] = vcount2-1-i;
    }
    fixed_count = 0;
  } else {
    // possibly fixed nodes specified (but can also be -1 everywhere)
    i = 0; // index in node_ordering
    j = vcount2-1; // reverse index in node_ordering
    for (k = 0; k < vcount2; k++) {
      if ((int)VECTOR(*map21)[k] == -1) {
      	VECTOR(node_ordering)[j] = k;
      	j--;
      } else {
      	VECTOR(node_ordering)[i] = k;
      	i++;
      }
    }
    fixed_count = i;
  }

  // STEP 1: check the fixed assignments for consistency and add to partial solution

  IGRAPH_CHECK(igraph_vector_init(&partial_solution_stack, vcount2)); // filled with zeros
  stack_pos = 0;

  IGRAPH_CHECK(igraph_vector_init(&indeg1, 0));
  IGRAPH_CHECK(igraph_vector_init(&indeg2, 0));
  IGRAPH_CHECK(igraph_vector_init(&outdeg1, 0));
  IGRAPH_CHECK(igraph_vector_init(&outdeg2, 0));

  IGRAPH_CHECK(igraph_vector_init(&inneighs1, 0));
  IGRAPH_CHECK(igraph_vector_init(&inneighs2, 0));
  IGRAPH_CHECK(igraph_vector_init(&outneighs1, 0));
  IGRAPH_CHECK(igraph_vector_init(&outneighs2, 0));

  // iterate over all fixed nodes
  for (i = 0; !end && i < fixed_count; i++) {
    pattern_node = VECTOR(node_ordering)[i];
    target_node = VECTOR(*map21)[pattern_node];

    // check colors
    if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			    != VECTOR(*vertex_color2)[pattern_node])) {
      end = 1;
      break;
    }

    // check degrees
    IGRAPH_CHECK(igraph_degree(graph1, &indeg1, igraph_vss_1(target_node),
	    IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &indeg2, igraph_vss_1(pattern_node),
	    IGRAPH_IN, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph1, &outdeg1, igraph_vss_1(target_node),
	    IGRAPH_OUT, IGRAPH_LOOPS));
    IGRAPH_CHECK(igraph_degree(graph2, &outdeg2, igraph_vss_1(pattern_node),
	    IGRAPH_OUT, IGRAPH_LOOPS));
    if (!end && ((VECTOR(indeg1)[0] < VECTOR(indeg2)[0])
	  || (VECTOR(outdeg1)[0] < VECTOR(outdeg2)[0]))) {
      end = 1;
      break;
    }

    // check edges to other fixed assignments
    IGRAPH_CHECK(igraph_neighbors(graph1, &inneighs1, (igraph_integer_t) target_node, IGRAPH_IN));
    IGRAPH_CHECK(igraph_neighbors(graph2, &inneighs2, (igraph_integer_t) pattern_node, IGRAPH_IN));
    IGRAPH_CHECK(igraph_neighbors(graph1, &outneighs1, (igraph_integer_t) target_node, IGRAPH_OUT));
    IGRAPH_CHECK(igraph_neighbors(graph2, &outneighs2, (igraph_integer_t) pattern_node, IGRAPH_OUT));
    for (j = 0; !end && j < i; j++) {
      other_pattern_node = VECTOR(node_ordering)[j];
      other_target_node = VECTOR(*map21)[pattern_node];
      if (igraph_vector_binsearch2(&inneighs2, other_pattern_node)) {
	// there is an edge pattern_node <- other_pattern_node
	if (igraph_vector_binsearch2(&inneighs1, other_target_node)) {
	  // there is an edge target_node <- other_target_node
	  IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) other_target_node,
			  (igraph_integer_t) target_node, 1, 1));
	  IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) other_pattern_node,
			  (igraph_integer_t) pattern_node, 1, 1));
	  if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
	      VECTOR(*edge_color2)[(long int)eid2]) {
	    end = 1;
	  }
	} else {
	  // edge between target nodes is missing
	  end = 1;
	}
      }
      if (!end && igraph_vector_binsearch2(&outneighs2, other_pattern_node)) {
	// there is an edge pattern_node -> other_pattern_node
	if (igraph_vector_binsearch2(&outneighs1, other_target_node)) {
	  // there is an edge target_node -> other_target_node
	  IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) target_node,
			  (igraph_integer_t) other_target_node, 1, 1));
	  IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) pattern_node,
			  (igraph_integer_t) other_pattern_node, 1, 1));
	  if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
	      VECTOR(*edge_color2)[(long int)eid2]) {
	    end = 1;
	  }
	} else {
	  // edge between target nodes is missing
	  end = 1;
	}
      }
      if (induced) {
	if (!end && igraph_vector_binsearch2(&inneighs1, other_target_node)) {
	  // there is an edge target_node <- other_target_node
	  if (igraph_vector_binsearch2(&inneighs2, other_pattern_node)) {
	    // there is an edge pattern_node <- other_pattern_node
	    IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) other_target_node,
			    (igraph_integer_t) target_node, 1, 1));
	    IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) other_pattern_node,
			    (igraph_integer_t) pattern_node, 1, 1));
	    if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		VECTOR(*edge_color2)[(long int)eid2]) {
	      end = 1;
	    }
	  } else {
	    // edge between pattern nodes is missing
	    end = 1;
	  }
	}
	if (!end && igraph_vector_binsearch2(&outneighs1, other_target_node)) {
	  // there is an edge target_node -> other_target_node
	  if (igraph_vector_binsearch2(&outneighs2, other_pattern_node)) {
	    // there is an edge pattern_node -> other_pattern_node
	    IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) target_node,
			    (igraph_integer_t) other_target_node, 1, 1));
	    IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) pattern_node,
			    (igraph_integer_t) other_pattern_node, 1, 1));
	    if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		VECTOR(*edge_color2)[(long int)eid2]) {
	      end = 1;
	    }
	  } else {
	    // edge between target nodes is missing
	    end = 1;
	  }
	}
      } // if induced
    } // for j (other fixed nodes)
    VECTOR(partial_solution_stack)[stack_pos] = target_node;
    stack_pos++;
  } // for i (fixed nodes)

  // STEP 2: fill the other assignments with DFS
  if (!end) {
    while (stack_pos < vcount2 && stack_pos >= fixed_count) {
      success = 1;
      pattern_node = VECTOR(node_ordering)[stack_pos];
      target_node = VECTOR(partial_solution_stack)[stack_pos];

      // check colors
      if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			      != VECTOR(*vertex_color2)[pattern_node])) {
	success = 0;
      }

      // check whether target node has been matched before
      for (i = 0; success && i < stack_pos; i++) {
	other_target_node = VECTOR(partial_solution_stack)[i];
	if (other_target_node == target_node) {
	  success = 0;
	}
      }

      // check degrees
      if (success) {
	IGRAPH_CHECK(igraph_degree(graph1, &indeg1, igraph_vss_1(target_node),
	      IGRAPH_IN, IGRAPH_LOOPS));
	IGRAPH_CHECK(igraph_degree(graph2, &indeg2, igraph_vss_1(pattern_node),
	    IGRAPH_IN, IGRAPH_LOOPS));
	IGRAPH_CHECK(igraph_degree(graph1, &outdeg1, igraph_vss_1(target_node),
	    IGRAPH_OUT, IGRAPH_LOOPS));
	IGRAPH_CHECK(igraph_degree(graph2, &outdeg2, igraph_vss_1(pattern_node),
	    IGRAPH_OUT, IGRAPH_LOOPS));
	if ((VECTOR(indeg1)[0] < VECTOR(indeg2)[0])
	      || (VECTOR(outdeg1)[0] < VECTOR(outdeg2)[0])) {
	  success = 0;
	}
      }

      // check edges to already matched nodes
      IGRAPH_CHECK(igraph_neighbors(graph1, &inneighs1, (igraph_integer_t) target_node, IGRAPH_IN));
      IGRAPH_CHECK(igraph_neighbors(graph2, &inneighs2, (igraph_integer_t) pattern_node, IGRAPH_IN));
      IGRAPH_CHECK(igraph_neighbors(graph1, &outneighs1, (igraph_integer_t) target_node, IGRAPH_OUT));
      IGRAPH_CHECK(igraph_neighbors(graph2, &outneighs2, (igraph_integer_t) pattern_node, IGRAPH_OUT));
      for (i = 0; success && i < stack_pos; i++) {
	other_pattern_node = VECTOR(node_ordering)[i];
	other_target_node = VECTOR(partial_solution_stack)[i];
	if (igraph_vector_binsearch2(&inneighs2, other_pattern_node)) {
	  // there is an edge pattern_node <- other_pattern_node
	  if (igraph_vector_binsearch2(&inneighs1, other_target_node)) {
	    // there is an edge target_node <- other_target_node
	    IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) other_target_node,
			    (igraph_integer_t) target_node, 1, 1));
	    IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) other_pattern_node,
			    (igraph_integer_t) pattern_node, 1, 1));
	    if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		VECTOR(*edge_color2)[(long int)eid2]) {
	      success = 0;
	    }
	  } else {
	    // edge between target nodes is missing
	    success = 0;
	  }
	}
	if (success && igraph_vector_binsearch2(&outneighs2, other_pattern_node)) {
	  // there is an edge pattern_node -> other_pattern_node
	  if (igraph_vector_binsearch2(&outneighs1, other_target_node)) {
	    // there is an edge target_node -> other_target_node
	    IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) target_node,
			    (igraph_integer_t) other_target_node, 1, 1));
	    IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) pattern_node,
			    (igraph_integer_t) other_pattern_node, 1, 1));
	    if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		VECTOR(*edge_color2)[(long int)eid2]) {
	      success = 0;
	    }
	  } else {
	    // edge between target nodes is missing
	    success = 0;
	  }
	}
	if (induced) {
	  if (success && igraph_vector_binsearch2(&inneighs1, other_target_node)) {
	    // there is an edge target_node <- other_target_node
	    if (igraph_vector_binsearch2(&inneighs2, other_pattern_node)) {
	      // there is an edge pattern_node <- other_pattern_node
	      IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) other_target_node,
			      (igraph_integer_t) target_node, 1, 1));
	      IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) other_pattern_node,
			      (igraph_integer_t) pattern_node, 1, 1));
	      if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		  VECTOR(*edge_color2)[(long int)eid2]) {
		success = 0;
	      }
	    } else {
	      // edge between pattern nodes is missing
	      success = 0;
	    }
	  }
	  if (success && igraph_vector_binsearch2(&outneighs1, other_target_node)) {
	    // there is an edge target_node -> other_target_node
	    if (igraph_vector_binsearch2(&outneighs2, other_pattern_node)) {
	      // there is an edge pattern_node -> other_pattern_node
	      IGRAPH_CHECK(igraph_get_eid(graph1, &eid1, (igraph_integer_t) target_node,
			      (igraph_integer_t) other_target_node, 1, 1));
	      IGRAPH_CHECK(igraph_get_eid(graph2, &eid2, (igraph_integer_t) pattern_node,
			      (igraph_integer_t) other_pattern_node, 1, 1));
	      if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		  VECTOR(*edge_color2)[(long int)eid2]) {
		success = 0;
	      }
	    } else {
	      // edge between target nodes is missing
	      success = 0;
	    }
	  }
	} // if induced
      } // for j (other fixed nodes)

      if (success) {
	if (stack_pos == vcount2-1) {
	  /* successful matching, finish processing */
	  break;
	}
	// match the next pattern node according to the node ordering
	stack_pos++;
	VECTOR(partial_solution_stack)[stack_pos] = 0;
      } else {
	// perform backtracking if all candidates have been tried
	while (VECTOR(partial_solution_stack)[stack_pos] == vcount1-1) {
	  VECTOR(partial_solution_stack)[stack_pos] = 0;
	  stack_pos--;
	}
	if (stack_pos >= fixed_count) {
	  VECTOR(partial_solution_stack)[stack_pos]++; // try next parent
	}
      }
    }
  }

  if (!end && success && stack_pos == vcount2-1) {
    *iso = 1;
    if (map21 != NULL) {
      for (i = 0; i < vcount2; i++) {
	VECTOR(*map21)[(long int) VECTOR(node_ordering)[i]] = VECTOR(partial_solution_stack)[i];
      }
    }
  }

  igraph_vector_destroy(&node_ordering);
  igraph_vector_destroy(&partial_solution_stack);
  igraph_vector_destroy(&indeg1);
  igraph_vector_destroy(&indeg2);
  igraph_vector_destroy(&outdeg1);
  igraph_vector_destroy(&outdeg2);
  igraph_vector_destroy(&inneighs1);
  igraph_vector_destroy(&inneighs2);
  igraph_vector_destroy(&outneighs1);
  igraph_vector_destroy(&outneighs2);

  return 0;
}


/* graph1 is the larger graph, graph2 is the smaller graph */
int igraph_mib_support_fast(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_isocompat_t *node_compat_fn,
		       igraph_isocompat_t *edge_compat_fn,
		       igraph_integer_t *support) {
  igraph_vector_t target_counts, map21;
  igraph_bool_t iso;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);
  long int i, j;

  IGRAPH_CHECK(igraph_vector_init(&target_counts, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&map21, vcount2));
  for (i = 0; i < vcount2; i++) {
    for (j = 0; j < vcount1; j++) {
      igraph_vector_fill(&map21, -1.);
      VECTOR(map21)[i] = j; // force assignment: pattern node i -> target node j
      iso = 0;
      if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
	      edge_color2, 1, &map21, &iso)) {
        igraph_vector_destroy(&target_counts);
        igraph_vector_destroy(&map21);
        return 1;
      }
      if (iso) {
	VECTOR(target_counts)[i] = VECTOR(target_counts)[i] + 1;
      }
    }
  }

  *support = igraph_vector_min(&target_counts);
  igraph_vector_destroy(&target_counts);
  igraph_vector_destroy(&map21);
  return 0;
}

int igraph_acgm(const igraph_vector_ptr_t *graphdb, igraph_support_measure_t *supp_fn,
		igraph_real_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values) {
  return 0;
}

