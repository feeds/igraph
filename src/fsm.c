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
#include "igraph_stack.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_components.h"


// graph1 is the larger graph, graph2 is the smaller graph
int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_isocompat_t *node_compat_fn,
			   igraph_isocompat_t *edge_compat_fn,
			   igraph_bool_t induced,
			   igraph_integer_t *support) {
  igraph_bool_t iso;
  if (igraph_subisomorphic_vf2(graph1, graph2, vertex_color1, vertex_color2,
		edge_color1, edge_color2, induced, &iso, NULL, NULL, node_compat_fn,
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
/* naive implementation: iterates over all embeddings and collects target nodes */
int igraph_mib_support_slow(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_isocompat_t *node_compat_fn,
		       igraph_isocompat_t *edge_compat_fn,
		       igraph_bool_t induced,
		       igraph_integer_t *support) {
  igraph_vector_t map21, target_counts;
  igraph_matrix_t target_hits;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);

  igraph_vector_init(&map21, 0);
  igraph_matrix_init(&target_hits, vcount2, vcount1);
  igraph_matrix_null(&target_hits);
  if (igraph_subisomorphic_function_vf2(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
		edge_color2, induced, NULL, &map21, (igraph_isohandler_t *) igraph_i_mib_isohandler,
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
// Can handle a single fixed assignment (pattern node, target node) passed as a length-2 vector
// NOTE: Only works for connected pattern graphs!
//
// Algorithm:
//    1) build a DFS ordering of the pattern nodes (intuition: when matching the next node
//       starting from a partial solution, we only have to consider the neighbors of all
//       previously matched target nodes
//    2) match all pattern nodes in the order specified by the DFS ordering, using the neighbors
//       of their DFS parents as candidates, while maintaining the subgraph isomorphism
//       properties for every partial solution (matching node labels, matching degrees, matching
//       edges, in that order)
//
int igraph_i_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_bool_t induced,
			   igraph_vector_t *fixed,
			   igraph_bool_t *iso) {
  long int vcount1 = igraph_vcount(graph1);
  long int vcount2 = igraph_vcount(graph2);
  long int i, j, fixed_count, partial_solution_pos, pred;
  long int pattern_node, other_pattern_node, target_node, other_target_node;
  long int indeg1, indeg2, outdeg1, outdeg2;
  igraph_integer_t eid1, eid2;
  int end, success;
  igraph_vector_t node_ordering;
  igraph_vector_t pred_idx;
  igraph_vector_t visited;
  igraph_vector_t state_target_idx;
  igraph_vector_t state_target_node;
  igraph_stack_t dfs_node_stack;
  igraph_stack_t dfs_pred_stack;
  igraph_bool_t conn, directed;

  *iso = 0;
  end = 0;
  success = 1;
  directed = igraph_is_directed(graph1);

  // TODO: currently works only for connected patterns
  IGRAPH_CHECK(igraph_is_connected(graph2, &conn, IGRAPH_WEAK));
  if (!conn) {
    *iso = 0;
    return 0;
  }

  // create a static ordering of the pattern nodes by DFS
  // if a fixed assignment is given, use this node as root, otherwise take the one with index 0
  IGRAPH_CHECK(igraph_stack_init(&dfs_node_stack, vcount2*vcount2));
  IGRAPH_CHECK(igraph_stack_init(&dfs_pred_stack, vcount2*vcount2));
  IGRAPH_CHECK(igraph_vector_init(&node_ordering, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&pred_idx, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&visited, vcount2));
  if (fixed == NULL) {
    IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, 0));
    fixed_count = 0;
  } else {
    IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, VECTOR(*fixed)[0]));
    fixed_count = 1;
  }
  IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, -1)); // first node has no predecessor
  i = 0;
  while (!igraph_stack_empty(&dfs_node_stack)) {
    pattern_node = (long int) igraph_stack_pop(&dfs_node_stack);
    pred = (long int) igraph_stack_pop(&dfs_pred_stack);
    if (VECTOR(visited)[pattern_node] == 0) {
      // insert current node into ordering and set predecessor index
      VECTOR(node_ordering)[i] = pattern_node;
      VECTOR(pred_idx)[i] = pred;

      // add neighbors to stack
      for (j = 0; j < DEGREE(*graph2, pattern_node); j++) {
	IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, NEIGHBOR(*graph2, pattern_node, j)));
	IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, i));
      }

      VECTOR(visited)[pattern_node] = 1;
      i++;
    }
  }

  // initialize the state representation
  IGRAPH_CHECK(igraph_vector_init(&state_target_idx, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&state_target_node, vcount2));

  // STEP 1: check the fixed assignment for consistency and add to partial solution

  if (fixed != NULL) {
    pattern_node = VECTOR(*fixed)[0];
    target_node = VECTOR(*fixed)[1];

    // check color
    if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			    != VECTOR(*vertex_color2)[pattern_node])) {
      end = 1;
    }

    // check degree
    if (!end) {
      indeg1 = IN_DEGREE(*graph1, target_node);
      indeg2 = IN_DEGREE(*graph2, pattern_node);
      outdeg1 = OUT_DEGREE(*graph1, target_node);
      outdeg2 = OUT_DEGREE(*graph2, pattern_node);
      if (directed) {
	if ((indeg1 < indeg2) || (outdeg1 < outdeg2)) {
	  end = 1;
	}
      } else {
	if (indeg1+outdeg1 < indeg2+outdeg2) {
	  end = 1;
	}
      }
    }

    // initialize first node with fixed assignment
    VECTOR(state_target_idx)[0] = target_node;
    VECTOR(state_target_node)[0] = target_node;
    //printf("fixed %ld -> %ld, match? %d\n", pattern_node, target_node, !end);
  }

  // STEP 2: fill the other assignments with DFS

  if (!end && vcount2 > fixed_count) {
    // initialize first free assignment
    partial_solution_pos = fixed_count;
    VECTOR(state_target_idx)[partial_solution_pos] = 0;
    while (1) {
      success = 1;
      pattern_node = VECTOR(node_ordering)[partial_solution_pos];
      if (partial_solution_pos == 0) {
	// target index is actual target node
	VECTOR(state_target_node)[partial_solution_pos] = VECTOR(state_target_idx)[partial_solution_pos];
      } else {
	// target index is index in parent's neighborhood
	VECTOR(state_target_node)[partial_solution_pos] = NEIGHBOR(*graph1,
	    (long int) VECTOR(state_target_node)[(long int) VECTOR(pred_idx)[partial_solution_pos]],
	    (long int) VECTOR(state_target_idx)[partial_solution_pos]);
      }
      target_node = VECTOR(state_target_node)[partial_solution_pos];
      //igraph_vector_print(&state_target_node);

      // check colors
      if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			      != VECTOR(*vertex_color2)[pattern_node])) {
	success = 0;
      }

      // check whether target node has been matched before
      for (i = 0; success && i < partial_solution_pos; i++) {
	other_target_node = (long int) VECTOR(state_target_node)[i];
	if (other_target_node == target_node) {
	  success = 0;
	}
      }

      // check degrees
      if (success) {
	indeg1 = IN_DEGREE(*graph1, target_node);
	indeg2 = IN_DEGREE(*graph2, pattern_node);
	outdeg1 = OUT_DEGREE(*graph1, target_node);
	outdeg2 = OUT_DEGREE(*graph2, pattern_node);
	if (directed) {
	  if ((indeg1 < indeg2) || (outdeg1 < outdeg2)) {
	    success = 0;
	  }
	} else {
	  if (indeg1+outdeg1 < indeg2+outdeg2) {
	    success = 0;
	  }
	}
      }

      // check edges to already matched nodes
      for (i = 0; success && i < partial_solution_pos; i++) {
	other_pattern_node = VECTOR(node_ordering)[i];
	other_target_node = VECTOR(state_target_node)[i];

	igraph_get_eid(graph1, &eid1, other_target_node, target_node, 1, 0);
	igraph_get_eid(graph2, &eid2, other_pattern_node, pattern_node, 1, 0);
	if (eid2 > -1) {
	  if (eid1 == -1) {
	    success = 0;
	  } else {
	    if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		  VECTOR(*edge_color2)[(long int)eid2]) {
	      success = 0;
	    }
	  }
	} else {
	  if (induced && eid1 > -1) {
	    success = 0;
	  }
	}

	if (directed) {
	  igraph_get_eid(graph1, &eid1, target_node, other_target_node, 1, 0);
	  igraph_get_eid(graph2, &eid2, pattern_node, other_pattern_node, 1, 0);
	  if (eid2 > -1) {
	    if (eid1 == -1) {
	      success = 0;
	    } else {
	      if (edge_color1 && VECTOR(*edge_color1)[(long int)eid1] !=
		    VECTOR(*edge_color2)[(long int)eid2]) {
		success = 0;
	      }
	    }
	  } else {
	    if (induced && eid1 > -1) {
	      success = 0;
	    }
	  }
	}
      } // for i (other fixed nodes)

      if (success) {
	// partial solution is consistent
	if (partial_solution_pos == vcount2-1) {
	  // partial solution is a full solution, successful finish
	  break;
	}
	// initialize next position
	partial_solution_pos++;
	VECTOR(state_target_idx)[partial_solution_pos] = 0;
      } else {
	// partial solution has failed

	while ((partial_solution_pos > 0)
		&& (VECTOR(state_target_idx)[partial_solution_pos]
		    == DEGREE(*graph1, (long int) VECTOR(state_target_node)[(long int)
				VECTOR(pred_idx)[partial_solution_pos]])-1)) {
	  // all nodes from parent's neighborhood have been tried, perform backtracking
	  partial_solution_pos--;
	}

	if (partial_solution_pos == 0) {
	  if (fixed_count == 0) {
	    // special case: first position, no fixed nodes -> all target nodes are candidates
	    if (VECTOR(state_target_idx)[0] == vcount1-1) {
	      // all target nodes have been tried, no candidates left
	      break;
	    } else {
	      // try next target node
	      VECTOR(state_target_idx)[0] += 1;
	    }
	  } else if (fixed_count == 1) {
	    // special case: first position, fixed node -> cannot change assignment
	    break;
	  }
	} else {
	  // there are node candidates left in the parent's neighborhood, proceed to next node
	  VECTOR(state_target_idx)[partial_solution_pos] += 1;
	}
      }
    } // DFS
  } // if (!end)

  if (!end && success) {
    //printf("   success!\n");
    *iso = 1;
  }

  igraph_stack_destroy(&dfs_node_stack);
  igraph_stack_destroy(&dfs_pred_stack);
  igraph_vector_destroy(&visited);
  igraph_vector_destroy(&node_ordering);
  igraph_vector_destroy(&pred_idx);
  igraph_vector_destroy(&state_target_idx);
  igraph_vector_destroy(&state_target_node);

  return 0;
}


// graph1 is the larger graph, graph2 is the smaller graph.
// node_compat_fn and edge_compat_fn are unused.
int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_isocompat_t *node_compat_fn,
		       igraph_isocompat_t *edge_compat_fn,
		       igraph_bool_t induced,
		       igraph_integer_t *support) {
  igraph_vector_t target_counts, fixed;
  igraph_bool_t iso;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);
  long int i, j;

  // TODO: consider automorphisms of the pattern graph. If i and j are isomorphic,
  // and we found all matchings for i, we can reuse them for j.

  IGRAPH_CHECK(igraph_vector_init(&target_counts, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&fixed, 2));
  for (i = 0; i < vcount2; i++) {
    VECTOR(fixed)[0] = i; // force assignment: pattern node i
    for (j = 0; j < vcount1; j++) {
      VECTOR(fixed)[1] = j; // force assignment: target node j
      iso = 0;
      if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
	      edge_color2, induced, &fixed, &iso)) {
        igraph_vector_destroy(&target_counts);
        igraph_vector_destroy(&fixed);
        return 1;
      }
      if (iso) {
	VECTOR(target_counts)[i] = VECTOR(target_counts)[i] + 1;
      }
    }
  }

  //igraph_vector_print(&target_counts);
  *support = igraph_vector_min(&target_counts);
  igraph_vector_destroy(&target_counts);
  igraph_vector_destroy(&fixed);
  return 0;
}

int igraph_acgm(const igraph_vector_ptr_t *graphdb, igraph_support_measure_t *supp_fn,
		igraph_real_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values) {
  return 0;
}

