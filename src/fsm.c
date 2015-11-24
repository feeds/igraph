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
/* naive implementation: iterates over all embeddings and collects target nodes */
int igraph_mib_support_slow(const igraph_t *graph1,
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
// Can handle a single fixed assignment (pattern node, target node) passed as a length-2 vector
// NOTE: Only works for connected pattern graphs!
// Algorithm:
//    1) build a DFS ordering of the pattern nodes (intuition: when matching the next node
//       starting from a partial solution, we only have to consider the neighbors of all
//       previously matched target nodes
//    2) match all pattern nodes in the order specified by the DFS ordering, using the neighbors
//       of already matched target nodes as candidates, while maintaining the subgraph isomorphism
//       properties for every partial solution (matching node labels, matching degrees, matching
//       edges, in that order)
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
  long int i, j, fixed_count, partial_solution_pos;
  int end = 0, success = 1;
  long int pattern_node, other_pattern_node, target_node, other_target_node;
  igraph_vector_t indeg1, indeg2, outdeg1, outdeg2;
  igraph_vector_t inneighs1, inneighs2, outneighs1, outneighs2, neighs;
  igraph_vector_t node_ordering;
  igraph_vector_t visited;
  igraph_vector_t state_target_idx;
  igraph_vector_t state_nbrhood_idx;
  igraph_vector_ptr_t state_nbrhood_ptr;
  igraph_stack_t dfs_stack;
  igraph_integer_t eid1, eid2;

  *iso = 0;

  // create a static ordering of the pattern nodes by DFS
  // if a fixed assignment is given, use this node as root, otherwise take the one with index 0
  IGRAPH_CHECK(igraph_stack_init(&dfs_stack, vcount2*vcount2)); // allocate enough memory (V^2)
  IGRAPH_CHECK(igraph_vector_init(&node_ordering, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&visited, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&neighs, 0));
  if (fixed == NULL) {
    IGRAPH_CHECK(igraph_stack_push(&dfs_stack, 0));
    fixed_count = 0;
  } else {
    IGRAPH_CHECK(igraph_stack_push(&dfs_stack, VECTOR(*fixed)[0]));
    fixed_count = 1;
  }
  i = 0;
  while (!igraph_stack_empty(&dfs_stack)) {
    pattern_node = (long int) igraph_stack_pop(&dfs_stack);
    if (VECTOR(visited)[pattern_node] == 0) {
      // insert current node into ordering
      VECTOR(node_ordering)[i] = pattern_node;
      VECTOR(visited)[pattern_node] = 1;
      i++;

      // add neighbors to stack
      IGRAPH_CHECK(igraph_neighbors(graph2, &neighs, (igraph_integer_t)pattern_node, IGRAPH_ALL));
      for (j = 0; j < igraph_vector_size(&neighs); j++) {
	IGRAPH_CHECK(igraph_stack_push(&dfs_stack, VECTOR(neighs)[j]));
      }
    }
  }

  IGRAPH_CHECK(igraph_vector_init(&state_target_idx, vcount2)); // filled with zeros
  IGRAPH_CHECK(igraph_vector_init(&state_nbrhood_idx, vcount2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&state_nbrhood_ptr, vcount2));
  for (i = 1; i < vcount2; i++) { // we don't need the pointer at pos 0
    VECTOR(state_nbrhood_ptr)[i] = igraph_Calloc(1, igraph_vector_t);
    IGRAPH_CHECK(igraph_vector_init((igraph_vector_t *) VECTOR(state_nbrhood_ptr)[i], 0));
  }

  IGRAPH_CHECK(igraph_vector_init(&indeg1, 0));
  IGRAPH_CHECK(igraph_vector_init(&indeg2, 0));
  IGRAPH_CHECK(igraph_vector_init(&outdeg1, 0));
  IGRAPH_CHECK(igraph_vector_init(&outdeg2, 0));

  IGRAPH_CHECK(igraph_vector_init(&inneighs1, 0));
  IGRAPH_CHECK(igraph_vector_init(&inneighs2, 0));
  IGRAPH_CHECK(igraph_vector_init(&outneighs1, 0));
  IGRAPH_CHECK(igraph_vector_init(&outneighs2, 0));

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
      IGRAPH_CHECK(igraph_degree(graph1, &indeg1, igraph_vss_1(target_node),
	      IGRAPH_IN, IGRAPH_LOOPS));
      IGRAPH_CHECK(igraph_degree(graph2, &indeg2, igraph_vss_1(pattern_node),
	      IGRAPH_IN, IGRAPH_LOOPS));
      IGRAPH_CHECK(igraph_degree(graph1, &outdeg1, igraph_vss_1(target_node),
	      IGRAPH_OUT, IGRAPH_LOOPS));
      IGRAPH_CHECK(igraph_degree(graph2, &outdeg2, igraph_vss_1(pattern_node),
	      IGRAPH_OUT, IGRAPH_LOOPS));
      if ((VECTOR(indeg1)[0] < VECTOR(indeg2)[0]) || (VECTOR(outdeg1)[0] < VECTOR(outdeg2)[0])) {
	end = 1;
      }
    }
    // initialize first node with fixed assignment
    VECTOR(state_target_idx)[0] = target_node;
    VECTOR(state_nbrhood_idx)[0] = 0;
    VECTOR(state_nbrhood_ptr)[0] = NULL; // unused, always use full set of nodes at pos 0
    if (vcount2 > 1) {
      igraph_neighbors(graph1, (igraph_vector_t *)VECTOR(state_nbrhood_ptr)[1],
		       (igraph_integer_t) target_node, IGRAPH_ALL);
    }
  }

  // STEP 2: fill the other assignments with DFS

  if (!end && vcount2 > fixed_count) {
    // initialize first free assignment
    partial_solution_pos = fixed_count;
    VECTOR(state_target_idx)[partial_solution_pos] = 0;
    VECTOR(state_nbrhood_idx)[partial_solution_pos] = 1;
    while (1) {
      success = 1;
      pattern_node = VECTOR(node_ordering)[partial_solution_pos];
      if (partial_solution_pos == 0) {
	// target index is actual target node
	target_node = VECTOR(state_target_idx)[partial_solution_pos];
      } else {
	// target index is index in some parent's neighborhood
	// TODO: collapse again
	igraph_vector_t *nbrhood = (igraph_vector_t *)VECTOR(state_nbrhood_ptr)[
				  (int)VECTOR(state_nbrhood_idx)[partial_solution_pos]];
	target_node = VECTOR(*nbrhood)[(int)VECTOR(state_target_idx)[partial_solution_pos]];
      }

      // check colors
      if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			      != VECTOR(*vertex_color2)[pattern_node])) {
	success = 0;
      }

      // check whether target node has been matched before
      if (partial_solution_pos > 0 && target_node == (int)VECTOR(state_target_idx)[0]) {
	success = 0;
      }
      for (i = 1; success && i < partial_solution_pos; i++) {
	other_target_node = VECTOR(*(igraph_vector_t *)VECTOR(state_nbrhood_ptr)[
				    (int)VECTOR(state_nbrhood_idx)[i]]
			    )[(int)VECTOR(state_target_idx)[i]];
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
      for (i = 0; success && i < partial_solution_pos; i++) {
	other_pattern_node = VECTOR(node_ordering)[i];
	if (i == 0) {
	  // target index is actual target node
	  other_target_node = VECTOR(state_target_idx)[i];
	} else {
	  // target index is index in some parent's neighborhood
	  other_target_node = VECTOR(*(igraph_vector_t *)VECTOR(state_nbrhood_ptr)[
				    (int)VECTOR(state_nbrhood_idx)[i]]
			    )[(int)VECTOR(state_target_idx)[i]];
	}
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
	// partial solution is consistent
	if (partial_solution_pos == vcount2-1) {
	  // partial solution is a full solution
	  break;
	}
	// match the next pattern node according to the node ordering
	partial_solution_pos++;
	VECTOR(state_target_idx)[partial_solution_pos] = 0;
	VECTOR(state_nbrhood_idx)[partial_solution_pos] = 1;
	igraph_neighbors(graph1, (igraph_vector_t *)VECTOR(state_nbrhood_ptr)[partial_solution_pos],
			  (igraph_integer_t) target_node, IGRAPH_ALL);
      } else {
	// partial solution has failed

	while ((partial_solution_pos > 0)
		&& (VECTOR(state_target_idx)[partial_solution_pos]
		    == igraph_vector_size((igraph_vector_t *)VECTOR(state_nbrhood_ptr)[
					(int)VECTOR(state_nbrhood_idx)[partial_solution_pos]])-1)
		&& (VECTOR(state_nbrhood_idx)[partial_solution_pos] == partial_solution_pos)) {
	  // all nodes from all neighborhoods have been tried, perform backtracking
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
	  if (VECTOR(state_target_idx)[partial_solution_pos]
	        == igraph_vector_size((igraph_vector_t *)VECTOR(state_nbrhood_ptr)[
	      			  (int)VECTOR(state_nbrhood_idx)[partial_solution_pos]])-1) {
	    // all nodes from current neighborhood have been tried, proceed to next neighborhood
	    VECTOR(state_nbrhood_idx)[partial_solution_pos] += 1;
	    VECTOR(state_target_idx)[partial_solution_pos] = 0;
	  } else {
	    // there are node candidates left in the current neighborhood, proceed to next node
	    VECTOR(state_target_idx)[partial_solution_pos] += 1;
	  }
	}
      }
    } // DFS
  } // if (!end)

  if (!end && success) {
    *iso = 1;
    //if (map21 != NULL) {
    //  for (i = 0; i < vcount2; i++) {
    //    VECTOR(*map21)[(long int) VECTOR(node_ordering)[i]] = VECTOR(partial_solution_stack)[i];
    //  }
    //}
  }

  igraph_stack_destroy(&dfs_stack);
  igraph_vector_destroy(&visited);
  igraph_vector_destroy(&neighs);
  igraph_vector_destroy(&node_ordering);
  igraph_vector_destroy(&state_target_idx);
  igraph_vector_destroy(&state_nbrhood_idx);
  for (i = 1; i < vcount2; i++) {
    igraph_vector_destroy((igraph_vector_t *) VECTOR(state_nbrhood_ptr)[i]);
    igraph_free(VECTOR(state_nbrhood_ptr)[i]);
  }
  igraph_vector_ptr_destroy(&state_nbrhood_ptr);
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
int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_isocompat_t *node_compat_fn,
		       igraph_isocompat_t *edge_compat_fn,
		       igraph_integer_t *support) {
  igraph_vector_t target_counts, fixed;
  igraph_bool_t iso;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);
  long int i, j;

  IGRAPH_CHECK(igraph_vector_init(&target_counts, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&fixed, 2));
  for (i = 0; i < vcount2; i++) {
    VECTOR(fixed)[0] = i; // force assignment: pattern node i
    for (j = 0; j < vcount1; j++) {
      VECTOR(fixed)[1] = j; // force assignment: target node j
      iso = 0;
      if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
	      edge_color2, 1, &fixed, &iso)) {
        igraph_vector_destroy(&target_counts);
        igraph_vector_destroy(&fixed);
        return 1;
      }
      if (iso) {
	VECTOR(target_counts)[i] = VECTOR(target_counts)[i] + 1;
      }
    }
  }

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

