/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   IGraph library: frequent subgraph mining algorithms
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
#include "igraph_constructors.h"


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
		       igraph_bool_t induced,
		       igraph_integer_t *support,
		       igraph_integer_t min_supp) {
  igraph_vector_t map21, target_counts;
  igraph_matrix_t target_hits;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);

  igraph_vector_init(&map21, 0);
  igraph_matrix_init(&target_hits, vcount2, vcount1);
  igraph_matrix_null(&target_hits);
  if (igraph_subisomorphic_function_vf2(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
		edge_color2, induced, NULL, &map21, (igraph_isohandler_t *) igraph_i_mib_isohandler,
		NULL, NULL, (void *) &target_hits)) {
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
int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_bool_t induced,
		       igraph_integer_t *support,
		       igraph_integer_t min_supp) {
  igraph_vector_t target_counts, fixed;
  igraph_bool_t iso;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);
  long int i, j, automorphic_node;

  IGRAPH_CHECK(igraph_vector_init(&fixed, 2));

  // find all automorphic pattern nodes
  igraph_matrix_t automorphic_nodes;
  IGRAPH_CHECK(igraph_matrix_init(&automorphic_nodes, vcount2, vcount2));
  for (i = 0; i < vcount2; i++) {
    VECTOR(fixed)[0] = i; // force assignment: pattern node i
    for (j = 0; j < i; j++) {
      VECTOR(fixed)[1] = j; // force assignment: pattern node j
      iso = 0;
      if (igraph_i_subisomorphic(graph2, graph2, vertex_color2, vertex_color2, edge_color2,
	      edge_color2, induced, &fixed, &iso)) {
	igraph_vector_destroy(&fixed);
	igraph_matrix_destroy(&automorphic_nodes);
	return 1;
      }
      if (iso) {
	MATRIX(automorphic_nodes, i, j) = 1;
      }
    }
  }

  // test all possible pairs (pattern node i, target node j) for isomorphism
  IGRAPH_CHECK(igraph_vector_init(&target_counts, vcount2));
  for (i = 0; i < vcount2; i++) {
    VECTOR(fixed)[0] = i; // force assignment: pattern node i

    // check if this node is isomorphic to a previously checked node
    automorphic_node = -1;
    for (j = 0; j < i; j++) {
      if (MATRIX(automorphic_nodes, i, j) == 1) {
	automorphic_node = j;
	break;
      }
    }

    if (automorphic_node >= 0) {
      // we found an automorphic node, reuse its result
      VECTOR(target_counts)[i] = VECTOR(target_counts)[j];
    } else {
      // no automorphic node found, test all possible target assignments
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

      // early termination: the support can only be smaller than or equal to VECTOR(target_counts)[i].
      // if that value is already smaller than min_supp, we don't need to continue.
      if (min_supp >= 0 && VECTOR(target_counts)[i] < min_supp) {
	*support = 0;
	igraph_vector_destroy(&target_counts);
	igraph_vector_destroy(&fixed);
	igraph_matrix_destroy(&automorphic_nodes);
	return 0;
      }
    }
  }

  //igraph_vector_print(&target_counts);
  *support = igraph_vector_min(&target_counts);
  igraph_vector_destroy(&target_counts);
  igraph_vector_destroy(&fixed);
  igraph_matrix_destroy(&automorphic_nodes);
  return 0;
}

// graph1 is the larger graph, graph2 is the smaller graph
int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_bool_t induced,
			   igraph_integer_t *support,
			   igraph_integer_t min_supp) {
  igraph_bool_t iso;
  if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2,
		edge_color1, edge_color2, induced, NULL, &iso)) {
    return 1;
  }
  if (iso) {
    *support = 1;
  } else {
    *support = 0;
  }
  return 0;
}


int igraph_db_mib_support(const igraph_vector_ptr_t *graphs,
			  const igraph_vector_ptr_t *vertex_colors,
			  const igraph_vector_ptr_t *edge_colors,
			  const igraph_t *pattern,
			  const igraph_vector_int_t *pattern_vcolors,
			  const igraph_vector_int_t *pattern_ecolors,
			  igraph_bool_t induced,
			  igraph_integer_t *support) {
  long int i;
  igraph_integer_t gsupp;
  *support = 0;
  if (vertex_colors != NULL) {
    if (edge_colors != NULL) {
      for (i = 0; i < igraph_vector_ptr_size(graphs); i++) {
	igraph_mib_support((igraph_t *) VECTOR(*graphs)[i], pattern,
			   (igraph_vector_int_t *) VECTOR(*vertex_colors)[i], pattern_vcolors,
			   (igraph_vector_int_t *) VECTOR(*edge_colors)[i], pattern_ecolors,
			   /*induced=*/ 0, &gsupp, 1);
	*support += gsupp;
      }
    } else {
      for (i = 0; i < igraph_vector_ptr_size(graphs); i++) {
	igraph_mib_support((igraph_t *) VECTOR(*graphs)[i], pattern,
			   (igraph_vector_int_t *) VECTOR(*vertex_colors)[i], pattern_vcolors,
			   NULL, NULL, /*induced=*/ 0, &gsupp, 1);
	*support += gsupp;
      }
    }
  } else {
    if (edge_colors != NULL) {
      for (i = 0; i < igraph_vector_ptr_size(graphs); i++) {
	igraph_mib_support((igraph_t *) VECTOR(*graphs)[i], pattern, NULL, NULL,
			   (igraph_vector_int_t *) VECTOR(*edge_colors)[i], pattern_ecolors,
			   /*induced=*/ 0, &gsupp, 1);
	*support += gsupp;
      }
    } else {
      for (i = 0; i < igraph_vector_ptr_size(graphs); i++) {
	igraph_mib_support((igraph_t *) VECTOR(*graphs)[i], pattern, NULL, NULL,
			   NULL, NULL, /*induced=*/ 0, &gsupp, 1);
	*support += gsupp;
      }
    }
  }
  return 0;
}


int igraph_db_shallow_support(const igraph_vector_ptr_t *graphs,
			      const igraph_vector_ptr_t *vertex_colors,
			      const igraph_vector_ptr_t *edge_colors,
			      const igraph_t *pattern,
			      const igraph_vector_int_t *pattern_vcolors,
			      const igraph_vector_int_t *pattern_ecolors,
			      igraph_bool_t induced,
			      igraph_integer_t *support) {
  return 0;
}


int igraph_acgm(const igraph_vector_ptr_t *graphs, const igraph_vector_ptr_t *vertex_colors,
		const igraph_vector_ptr_t *edge_colors, igraph_db_support_measure_t *supp_fn,
		igraph_integer_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values) {
  return 0;
}


// assert: #graphs == #vertex_colors == #edge_colors
// assert: 0 <= node and edge colors <= MAX_COLOR
// assert: min_supp > 0
// assert: undirected graphs
int igraph_gspan(const igraph_vector_ptr_t *graphs, const igraph_vector_ptr_t *vertex_colors,
		 const igraph_vector_ptr_t *edge_colors, igraph_db_support_measure_t *db_supp_measure,
		 igraph_integer_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		 igraph_vector_t *support_values) {

  long int graph_count = igraph_vector_ptr_size(graphs);
  long int i, j, k;
  igraph_vector_int_t vcolor_freq, ecolor_freq; // frequencies of all colors
  igraph_vector_int_t freq_ecolors, freq_vcolors; // lists of all frequent colors
  igraph_vector_int_t *vcolor, *ecolor;
  igraph_integer_t supp;
  igraph_t *g;

  // FIND FREQUENT VERTEX AND EDGE COLORS

  // count all color occurrences
  igraph_vector_int_init(&vcolor_freq, MAX_COLOR+1);
  igraph_vector_int_init(&ecolor_freq, MAX_COLOR+1);
  if (vertex_colors != NULL || edge_colors != NULL) {
    for (i = 0; i < graph_count; i++) {
      g = (igraph_t *) VECTOR(*graphs)[i];
      if (vertex_colors != NULL) {
	vcolor = (igraph_vector_int_t *) VECTOR(*vertex_colors)[i];
	for (j = 0; j < igraph_vcount(g); j++) {
	  VECTOR(vcolor_freq)[VECTOR(*vcolor)[j]] += 1;
	}
      }
      if (edge_colors != NULL) {
	ecolor = (igraph_vector_int_t *) VECTOR(*edge_colors)[i];
	for (j = 0; j < igraph_ecount(g); j++) {
	  VECTOR(ecolor_freq)[VECTOR(*ecolor)[j]] += 1;
	}
      }
    }
  }

  // keep only frequent colors
  igraph_vector_int_init(&freq_vcolors, MAX_COLOR+1);
  igraph_vector_int_init(&freq_ecolors, MAX_COLOR+1);
  igraph_vector_int_fill(&freq_vcolors, -1);
  igraph_vector_int_fill(&freq_ecolors, -1);
  if (vertex_colors != NULL) {
    j = 0;
    for (i = 0; i <= MAX_COLOR; i++) {
      if (VECTOR(vcolor_freq)[i] >= min_supp) {
	VECTOR(freq_vcolors)[j] = i;
	j += 1;
      }
    }
  }
  if (edge_colors != NULL) {
    j = 0;
    for (i = 0; i <= MAX_COLOR; i++) {
      if (VECTOR(ecolor_freq)[i] >= min_supp) {
	VECTOR(freq_ecolors)[j] = i;
	j += 1;
      }
    }
  }

  // FIND FREQUENT 1-EDGE GRAPHS

  igraph_t pattern_graph;
  igraph_vector_int_t pattern_vcolor, pattern_ecolor;
  igraph_full(&pattern_graph, 2, /*directed=*/ 0, /*self-loops=*/ 0);
  igraph_vector_int_init(&pattern_vcolor, 2);
  igraph_vector_int_init(&pattern_ecolor, 1);

  if (vertex_colors != NULL) {
    for (i = 0; VECTOR(freq_vcolors)[i] != -1; i++) {
      for (j = 0; j <= i && VECTOR(freq_vcolors)[j] != -1; j++) {
	VECTOR(pattern_vcolor)[0] = VECTOR(freq_vcolors)[i];
	VECTOR(pattern_vcolor)[1] = VECTOR(freq_vcolors)[j];
	if (edge_colors != NULL) {
	  for (k = 0; VECTOR(freq_ecolors)[k] != -1; k++) {
	    // VC[i] -- EC[k] -- VC[j]
	    VECTOR(pattern_ecolor)[0] = VECTOR(freq_ecolors)[k];
	    db_supp_measure(graphs, vertex_colors, edge_colors, &pattern_graph,
			    &pattern_vcolor, &pattern_ecolor, /*induced=*/ 0, &supp);
	    if (supp >= min_supp) {
	      printf("supp=%2ld   ", (long int) supp);
	      igraph_i_print(&pattern_graph, &pattern_vcolor, &pattern_ecolor);
	    }
	  }
	} else {
	  // VC[i] -- VC[j]
	  db_supp_measure(graphs, vertex_colors, NULL, &pattern_graph,
			  &pattern_vcolor, NULL, /*induced=*/ 0, &supp);
	  if (supp >= min_supp) {
	    printf("supp=%2ld   ", (long int) supp);
	    igraph_i_print(&pattern_graph, &pattern_vcolor, NULL);
	  }
	}
      }
    }
  } else {
    if (edge_colors != NULL) {
      for (k = 0; VECTOR(freq_ecolors)[k] != -1; k++) {
	// O -- EC[k] -- O
	VECTOR(pattern_ecolor)[0] = VECTOR(freq_ecolors)[k];
	db_supp_measure(graphs, NULL, edge_colors, &pattern_graph,
			NULL, &pattern_ecolor, /*induced=*/ 0, &supp);
	if (supp >= min_supp) {
	  printf("supp=%2ld   ", (long int) supp);
	  igraph_i_print(&pattern_graph, NULL, &pattern_ecolor);
	}
      }
    } else {
      // O -- O
      db_supp_measure(graphs, NULL, NULL, &pattern_graph, NULL, NULL, /*induced=*/ 0, &supp);
      if (supp >= min_supp) {
	printf("supp=%2ld   ", (long int) supp);
	igraph_i_print(&pattern_graph, NULL, NULL);
      }
    }
  }

  // TODO: recursively expand all frequent 1-edge graphs by pattern growth

  igraph_vector_int_destroy(&vcolor_freq);
  igraph_vector_int_destroy(&ecolor_freq);
  igraph_vector_int_destroy(&freq_vcolors);
  igraph_vector_int_destroy(&freq_ecolors);
  igraph_vector_int_destroy(&pattern_vcolor);
  igraph_vector_int_destroy(&pattern_ecolor);
  igraph_destroy(&pattern_graph);

  return 0;
}
