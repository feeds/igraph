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
#include <string.h> // memcpy
#include "igraph_fsm.h"
#include "igraph_matrix.h"
#include "igraph_stack.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_components.h"
#include "igraph_constructors.h"


// ------------- HELPER FUNCTIONS -------------

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
    return 1;
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


// ------------- SUPPORT MEASURES -------------


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

      // early termination: the support can only be smaller than or equal to
      // VECTOR(target_counts)[i].
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
		edge_color1, edge_color2, induced, /*fixed=*/ NULL, &iso)) {
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


// ------------- GSPAN -------------

// Linked list (internal)

typedef struct igraph_llist_item_t {
  void *data;
  struct igraph_llist_item_t *next;
  struct igraph_llist_item_t *prev;
} igraph_llist_item_t;

typedef struct igraph_llist_t {
  igraph_llist_item_t *first;
  igraph_llist_item_t *last;
} igraph_llist_t;

int igraph_i_llist_init(igraph_llist_t *llist);
void igraph_i_llist_destroy(igraph_llist_t *llist);
int igraph_i_llist_push_back(igraph_llist_t *llist, void *data);

int igraph_i_llist_init(igraph_llist_t *llist) {
  llist->first = NULL;
  llist->last = NULL;
  return 0;
}

void igraph_i_llist_destroy(igraph_llist_t *llist) {
  igraph_llist_item_t *cur, *prev;
  cur = llist->first;
  while (cur != NULL) {
    prev = cur;
    cur = prev->next;
    igraph_free(prev);
  }
}

int igraph_i_llist_push_back(igraph_llist_t *llist, void *data) {
  igraph_llist_item_t *new = igraph_Calloc(1, igraph_llist_item_t);
  if (new == NULL) {
    return 1;
  }
  new->data = data;
  new->next = NULL;
  new->prev = llist->last;

  if (llist->last == NULL) {
    llist->last = new;
    llist->first = new;
  } else {
    llist->last->next = new;
    llist->last = new;
  }
  return 0;
}

int igraph_i_llist_push_front(igraph_llist_t *llist, void *data) {
  igraph_llist_item_t *new = igraph_Calloc(1, igraph_llist_item_t);
  if (new == NULL) {
    return 1;
  }
  new->data = data;
  new->next = llist->first;
  new->prev = NULL;

  if (llist->first == NULL) {
    llist->first = new;
    llist->last = new;
  } else {
    llist->first->prev = new;
    llist->first = new;
  }
  return 0;
}


// DFS code related (internal)

typedef struct igraph_dfscode_edge_t {
  long int i; // source node
  long int j; // target node
  long int l_i; // source label
  long int l_ij; // edge label
  long int l_j; // target label
} igraph_dfscode_edge_t;

typedef struct igraph_dfscode_t {
  igraph_dfscode_edge_t *stor_begin;
  long int last_edge;
  long int max_edges;
} igraph_dfscode_t;

int igraph_i_dfscode_init(igraph_dfscode_t *dfscode, long int max_edges);
int igraph_i_dfscode_init_copy(igraph_dfscode_t *dfscode_to, const igraph_dfscode_t *dfscode_from);
void igraph_i_dfscode_destroy(igraph_dfscode_t *dfscode);
void igraph_i_dfscode_print(const igraph_dfscode_t *dfscode);
long int igraph_i_dfscode_size(const igraph_dfscode_t *dfscode);
int igraph_i_dfscode_push_back(igraph_dfscode_t *dfscode, const igraph_dfscode_edge_t *edge);
igraph_dfscode_edge_t igraph_i_dfscode_pop_back(igraph_dfscode_t *dfscode);
int igraph_i_dfscode_edge_compare(const igraph_dfscode_edge_t *a, const igraph_dfscode_edge_t *b);
int igraph_i_dfscode_compare(const igraph_dfscode_t *a, const igraph_dfscode_t *b);
int igraph_i_dfscode_to_graph(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors);
igraph_bool_t igraph_i_dfscode_is_canonical(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors);
igraph_bool_t igraph_i_dfscode_is_canonical_rec(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		long int dfscode_pos, igraph_vector_int_t *ordering, long int ordering_pos,
		igraph_vector_int_t *visited_nodes);
int igraph_i_dfscode_extend(const igraph_vector_ptr_t *graphs,
		const igraph_vector_ptr_t *vertex_colors, const igraph_vector_ptr_t *edge_colors,
		igraph_db_support_measure_t *db_supp_measure,
		igraph_integer_t min_supp, igraph_integer_t max_edges,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors,
		igraph_dfscode_t *seed_dfscode, igraph_llist_t *result_graph_list,
		igraph_llist_t *result_vcolor_list, igraph_llist_t *result_ecolor_list,
		igraph_llist_t *result_supp_list);

int igraph_i_dfscode_init(igraph_dfscode_t *dfscode, long int max_edges) {
  dfscode->stor_begin = igraph_Calloc(max_edges, igraph_dfscode_edge_t);
  if (dfscode->stor_begin == NULL) {
    return 1;
  }
  dfscode->last_edge = -1;
  dfscode->max_edges = max_edges;
  return 0;
}

int igraph_i_dfscode_init_copy(igraph_dfscode_t *dfscode_to, const igraph_dfscode_t *dfscode_from) {
  IGRAPH_CHECK(igraph_i_dfscode_init(dfscode_to, dfscode_from->max_edges));
  memcpy(dfscode_to->stor_begin, dfscode_from->stor_begin,
	      sizeof(igraph_dfscode_edge_t)*dfscode_from->last_edge+1);
  return 0;
}

void igraph_i_dfscode_destroy(igraph_dfscode_t *dfscode) {
  igraph_free(dfscode->stor_begin);
}

void igraph_i_dfscode_print(const igraph_dfscode_t *dfscode) {
  long int i;
  for (i = 0; i < dfscode->last_edge+1; i++) {
    printf("(%ld, %ld, %ld, %ld, %ld) ", VECTOR(*dfscode)[i].i, VECTOR(*dfscode)[i].j,
	VECTOR(*dfscode)[i].l_i, VECTOR(*dfscode)[i].l_ij, VECTOR(*dfscode)[i].l_j);
  }
  printf("\n");
}

long int igraph_i_dfscode_size(const igraph_dfscode_t *dfscode) {
  return dfscode->last_edge+1;
}

int igraph_i_dfscode_push_back(igraph_dfscode_t *dfscode, const igraph_dfscode_edge_t *edge) {
  if (dfscode->last_edge == dfscode->max_edges+1) {
    return 1;
  }
  dfscode->last_edge += 1;
  VECTOR(*dfscode)[dfscode->last_edge] = *edge;
  return 0;
}

// user has to make sure that the code is non-empty
igraph_dfscode_edge_t igraph_i_dfscode_pop_back(igraph_dfscode_t *dfscode) {
  dfscode->last_edge -= 1;
  return VECTOR(*dfscode)[dfscode->last_edge+1];
}

// definition from CloseGraph paper (Yan & Han 2003)
int igraph_i_dfscode_edge_compare(const igraph_dfscode_edge_t *a, const igraph_dfscode_edge_t *b) {
  // first priority: DFS edge ordering (i,j)
  if ((a->i < a->j) && (b->i < b->j)) {
    // a and b are forward edges
    if ((a->j < b->j) || ((a->i > b->i) && (a->j == b->j)))
      return -1; // a < b
  }
  if ((a->i > a->j) && (b->i > b->j)) {
    // a and b are backward edges
    if ((a->i < b->i) || ((a->i == b->i) && (a->j < b->j))) {
      return -1; // a < b
    }
  }
  if ((a->i > a->j) && (b->i < b->j)) {
    // a is backward, b is forward edge
    if (a->i < b->j) {
      return -1; // a < b
    }
  }
  if ((a->i < a->j) && (b->i > b->j)) {
    // a is forward, b is backward edge
    if (a->j <= b->i) {
      return -1; // a < b
    }
  }

  if (a->i == b->i && a->j == b->j) {
    // second priority: label of node i
    if (a->l_i < b->l_i) {
      return -1; // a < b
    }
    if (a->l_i == b->l_j) {
      // third priority: label of edge (i,j)
      if (a->l_ij < b->l_ij) {
	return -1; // a < b
      }
      if (a->l_ij == b->l_ij) {
	// fourth priority: label of node j
	if (a->l_j < b->l_j) {
	  return -1; // a < b
	}
	if (a->l_j == b->l_j) {
	  // all entries are equal
	  return 0; // a == b
	}
      }
    }
  }

  return 1; // a > b
}

// DFS lexicographic order
int igraph_i_dfscode_compare(const igraph_dfscode_t *a, const igraph_dfscode_t *b) {
  long int min_len = (a->last_edge < b->last_edge) ? a->last_edge+1 : b->last_edge+1;
  long int i;
  int cmp;
  for (i = 0; i < min_len; i++) {
    cmp = igraph_i_dfscode_edge_compare(&VECTOR(*a)[i], &VECTOR(*b)[i]);
    if (cmp == -1) {
      return -1; // a < b
    }
    if (cmp == 0) {
      continue;
    }
    return 1; // a > b
  }

  // the first min_len edges are equal, compare lengths
  if (a->last_edge < b->last_edge) {
    return -1; // a < b
  }
  if (a->last_edge == b->last_edge) {
    return 0; // a == b
  }
  return 1; // a > b
}

int igraph_i_dfscode_to_graph(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors) {
  igraph_vector_t edges;
  long int i;
  long int rightmost_vertex = ((VECTOR(*dfscode)[dfscode->last_edge].i
				< VECTOR(*dfscode)[dfscode->last_edge].j)
			      ? VECTOR(*dfscode)[dfscode->last_edge].j
			      : VECTOR(*dfscode)[dfscode->last_edge].i);

  IGRAPH_CHECK(igraph_empty(graph, rightmost_vertex+1, IGRAPH_UNDIRECTED));
  IGRAPH_CHECK(igraph_vector_int_init(vertex_colors, rightmost_vertex+1));
  IGRAPH_CHECK(igraph_vector_int_init(edge_colors, igraph_i_dfscode_size(dfscode)));
  IGRAPH_CHECK(igraph_vector_init(&edges, 2*igraph_i_dfscode_size(dfscode)));

  for (i = 0; i < igraph_i_dfscode_size(dfscode); i++) {
    VECTOR(edges)[2*i] = VECTOR(*dfscode)[i].i;
    VECTOR(edges)[2*i+1] = VECTOR(*dfscode)[i].j;
    VECTOR(*edge_colors)[i] = VECTOR(*dfscode)[i].l_ij;
    VECTOR(*vertex_colors)[VECTOR(*dfscode)[i].i] = VECTOR(*dfscode)[i].l_i;
    VECTOR(*vertex_colors)[VECTOR(*dfscode)[i].j] = VECTOR(*dfscode)[i].l_j;
  }

  igraph_add_edges(graph, &edges, 0);
  igraph_vector_destroy(&edges);

  return 0;
}

// Mentioned briefly in [Yan & Han 2002], DFS version of [Borgelt 2006]:
// Perform DFS using all possible nodes as the root node to check whether a smaller
// code than dfscode can be constructed. As soon as we have found a smaller prefix, we
// can terminate.
// Original paper says: almost the same as enumerating all automorphisms, but taking
// advantage of the prefix property for pruning.
// TODO: only works for undirected graphs!!!
igraph_bool_t igraph_i_dfscode_is_canonical(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors) {
  long int root;
  igraph_vector_int_t visited_nodes, ordering;
  igraph_bool_t is_smaller;

  if (igraph_i_dfscode_size(dfscode) == 0) {
    // always minimal
    return 1;
  }
  if (igraph_i_dfscode_size(dfscode) == 1) {
    // only minimal if the first edge connects nodes 0 and 1
    if ((VECTOR(*dfscode)[0].i != 0) || (VECTOR(*dfscode)[0].j != 1)) {
      return 0;
    }
  }

  IGRAPH_CHECK(igraph_vector_int_init(&visited_nodes, igraph_vcount(graph)));
  IGRAPH_CHECK(igraph_vector_int_init(&ordering, igraph_vcount(graph)));

  for (root = 0; root < igraph_vcount(graph); root++) {
    if (vertex_colors) {
      if (VECTOR(*vertex_colors)[root] < VECTOR(*dfscode)[0].l_i) {
	// new root has a smaller label than dfscode root. dfscode not minimal.
	igraph_vector_int_destroy(&ordering);
	igraph_vector_int_destroy(&visited_nodes);
	return 0;
      }
      if (VECTOR(*vertex_colors)[root] > VECTOR(*dfscode)[0].l_i) {
	// new root has a larger label than dfscode root. it cannot give a smaller DFS code.
	continue;
      }
    }

    // recursively build all DFS orderings starting from this root node
    VECTOR(ordering)[0] = root;
    igraph_vector_int_fill(&visited_nodes, -1);
    VECTOR(visited_nodes)[root] = 0;
    is_smaller = igraph_i_dfscode_is_canonical_rec(dfscode, graph, vertex_colors, edge_colors,
						   0, &ordering, 0, &visited_nodes);
    if (!is_smaller) {
      // dfscode is not smaller than or equal to all codes grown by DFS from this root node,
      // hence it is not minimal
      igraph_vector_int_destroy(&ordering);
      igraph_vector_int_destroy(&visited_nodes);
      return 0;
    }
    VECTOR(visited_nodes)[root] = -1;
  }

  // dfscode is smaller than or equal to all codes grown from any root node,
  // hence it is minimal
  igraph_vector_int_destroy(&ordering);
  igraph_vector_int_destroy(&visited_nodes);
  return 1;
}


// Check if the provided dfscode is smaller than or equal to all codes that can be grown
// from the provided partial DFS ordering. The code is assumed to be equal up to dfscode_pos.
// Returns TRUE (1) if the provided dfscode is smaller than or equal to all codes, FALSE (0)
// if a code is found that is smaller than dfscode.
igraph_bool_t igraph_i_dfscode_is_canonical_rec(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		long int dfscode_pos, igraph_vector_int_t *ordering, long int ordering_pos,
		igraph_vector_int_t *visited_nodes) {
  long int i, j, k, neigh, backward_neigh, backward_edge_count, ext_node, ext_node_pos;
  int cmp;
  igraph_dfscode_edge_t new_edge;
  igraph_integer_t eid;
  igraph_bool_t is_less_eq;

  if (ordering_pos == igraph_vector_int_size(ordering)-1) {
    // DFS ordering is complete (no extension possible) and we have
    // not found a smaller code
    return 1;
  }

  // find the next extendible node: walk up the DFS tree and look for unvisited neighbors
  ext_node_pos = -1;
  for (i = ordering_pos; i >= 0; i--) {
    for (j = 0; j < DEGREE(*graph, VECTOR(*ordering)[i]); j++) {
      if (VECTOR(*visited_nodes)[NEIGHBOR(*graph, VECTOR(*ordering)[i], j)] == -1) {
	ext_node_pos = i;
	i = -1;
	break;
      }
    }
  }

  // iterate over all neighbors of the extendible node
  // TODO: we could sort neighbors by the edge label (smallest first) here to speed up processing
  ext_node = VECTOR(*ordering)[ext_node_pos];
  for (i = 0; i < DEGREE(*graph, ext_node); i++) {
    neigh = NEIGHBOR(*graph, ext_node, i);
    if (VECTOR(*visited_nodes)[neigh] >= 0) {
      // skip visited neighbors
      continue;
    }

    // check if visiting the current neighbor would result in an edge code smaller
    // than the one from dfscode

    new_edge = (igraph_dfscode_edge_t) {.i = ext_node_pos, .j = ordering_pos+1};
    if (vertex_colors) {
      new_edge.l_i = VECTOR(*vertex_colors)[ext_node];
      new_edge.l_j = VECTOR(*vertex_colors)[neigh];
    }
    if (edge_colors) {
      igraph_get_eid(graph, &eid, ext_node, neigh, 1, 0);
      new_edge.l_ij = VECTOR(*edge_colors)[eid];
    }

    cmp = igraph_i_dfscode_edge_compare(&new_edge, &VECTOR(*dfscode)[dfscode_pos]);
    if (cmp == -1) {
      // new edge is smaller, dfscode not minimal
      //printf("a\n");
      return 0;
    } else if (cmp == 1) {
      // new edge is larger, cannot result in smaller code at later point in time,
      // try next neigh
      // TODO: if we sort the neighbors as noted above, we can return 1 here
      continue;
    }

    // new edge is equal to existing one, insert into ordering
    VECTOR(*ordering)[ordering_pos+1] = neigh;
    VECTOR(*visited_nodes)[neigh] = ordering_pos+1;

    // check all backward edges originating from neigh
    // TODO: should we check the right-most path only?
    backward_edge_count = 0;
    for (j = 0; j < ext_node_pos; j++) {
      backward_neigh = VECTOR(*ordering)[j];

      // check if backward_neigh is actually a neighbor of neigh
      for (k = 0; k < DEGREE(*graph, neigh); k++) {
	if (NEIGHBOR(*graph, neigh, k) == backward_neigh) {
	  // it is, create the backward edge
	  new_edge = (igraph_dfscode_edge_t) {.i = ordering_pos+1, .j = j};
	  if (vertex_colors) {
	    new_edge.l_i = VECTOR(*vertex_colors)[neigh];
	    new_edge.l_j = VECTOR(*vertex_colors)[backward_neigh];
	  }
	  if (edge_colors) {
	    igraph_get_eid(graph, &eid, neigh, backward_neigh, /*directed=*/ 1, /*error=*/ 0);
	    new_edge.l_ij = VECTOR(*edge_colors)[eid];
	  }

	  // check if the backward edge is smaller than the next dfscode position
	  cmp = igraph_i_dfscode_edge_compare(&new_edge,
		    &VECTOR(*dfscode)[dfscode_pos+backward_edge_count+1]);
	  if (cmp == -1) {
	    // new edge is smaller, dfscode not minimal
	    //igraph_vector_int_print(ordering);
	    //printf("b -- bwEdge (%ld,%ld) < (%ld,%ld)\n", neigh, backward_neigh,
	    //           VECTOR(*dfscode)[dfscode_pos+backward_edge_count+1].i,
	    //           VECTOR(*dfscode)[dfscode_pos+backward_edge_count+1].j);
	    return 0;
	  } else if (cmp == 1) {
	    // new edge is larger, cannot result in smaller code at later point in time,
	    // try next choice for neigh (outmost loop)
	    j = ordering_pos+1; // breaks loop over j and continues loop over i
	    break;
	  }

	  backward_edge_count += 1;
	  break;
	}
      }
    } // for j (backward edges)

    if (j == ordering_pos+1) {
      // adding a backward edge resulted in a larger code than dfscode, try next neigh
      // TODO: if we sort the neighbors as noted above, we can return 1 here
      continue;
    }

    // recursively continue search for DFS ordering
    is_less_eq = igraph_i_dfscode_is_canonical_rec(dfscode, graph, vertex_colors, edge_colors,
						   dfscode_pos+backward_edge_count+1,
						   ordering, ordering_pos+1, visited_nodes);
    if (!is_less_eq) {
      // we can grow a DFS code from the current node ordering that is smaller than dfscode,
      // hence dfscode is not minimal
      //printf("c\n");
      return 0;
    }
    VECTOR(*visited_nodes)[neigh] = -1;

    // dfscode is less than or equal to all DFS codes grown from the current node ordering
    // try the next assignment for position ordering_pos

  } // for i (neigh)

  // dfscode is less than or equal to all DFS codes grown from the current node ordering
  // for all assignments for position ordering_pos, perform backtracking
  return 1;
}




// TODO: Implement speed-up for labelled graphs:
// Check that a new edge is larger than all existing edges incident to the target
// node (for backward extension) or the source node (for forward extension from any
// node other than right-most vertex). Otherwise, the resulting DFS code can not be
// minimal. See Yan & Han (2002), Section 5.1 (2), and Borgelt (2006), Section 4.
int igraph_i_dfscode_extend(const igraph_vector_ptr_t *graphs,
		const igraph_vector_ptr_t *vertex_colors, const igraph_vector_ptr_t *edge_colors,
		igraph_db_support_measure_t *db_supp_measure,
		igraph_integer_t min_supp, igraph_integer_t max_edges,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors,
		igraph_dfscode_t *seed_dfscode, igraph_llist_t *result_graph_list,
		igraph_llist_t *result_vcolor_list, igraph_llist_t *result_ecolor_list,
		igraph_llist_t *result_supp_list) {
  long int i, j, cur_color, rightmost_vertex_color;
  long int cur_vertex, prev_vertex, rightmost_vertex, rightmost_vertex_pred;
  igraph_stack_int_t rightmost_path, rightmost_path_colors;
  igraph_dfscode_edge_t new_edge;
  igraph_t *seed_graph;
  igraph_vector_int_t *seed_vcolors;
  igraph_vector_int_t *seed_ecolors;
  igraph_integer_t *seed_supp;

  // create graph from DFS code
  seed_graph = igraph_Calloc(1, igraph_t);
  seed_vcolors = igraph_Calloc(1, igraph_vector_int_t);
  seed_ecolors = igraph_Calloc(1, igraph_vector_int_t);
  seed_supp = igraph_Calloc(1, igraph_integer_t);
  IGRAPH_CHECK(igraph_i_dfscode_to_graph(seed_dfscode, seed_graph, seed_vcolors, seed_ecolors));

  igraph_i_dfscode_print(seed_dfscode);
  if (!igraph_i_dfscode_is_canonical(seed_dfscode, seed_graph, seed_vcolors, seed_ecolors)) {
    // seed not in canonical form, prune
    printf("   not canonical\n");
    igraph_destroy(seed_graph);
    igraph_vector_int_destroy(seed_vcolors);
    igraph_vector_int_destroy(seed_ecolors);
    igraph_free(seed_graph);
    igraph_free(seed_vcolors);
    igraph_free(seed_ecolors);
    igraph_free(seed_supp);
    return 0;
  }

  // compute seed support
  db_supp_measure(graphs, vertex_colors, edge_colors, seed_graph,
		  seed_vcolors, seed_ecolors, /*induced=*/ 0, seed_supp);
  if (*seed_supp < min_supp) {
    // infrequent seed, free memory and prune
    igraph_destroy(seed_graph);
    igraph_vector_int_destroy(seed_vcolors);
    igraph_vector_int_destroy(seed_ecolors);
    igraph_free(seed_graph);
    igraph_free(seed_vcolors);
    igraph_free(seed_ecolors);
    igraph_free(seed_supp);
    return 0;
  } else {
    // frequent seed, add to result
    igraph_i_llist_push_back(result_graph_list, seed_graph);
    igraph_i_llist_push_back(result_vcolor_list, seed_vcolors);
    igraph_i_llist_push_back(result_ecolor_list, seed_ecolors);
    igraph_i_llist_push_back(result_supp_list, seed_supp);
  }

  if (igraph_i_dfscode_size(seed_dfscode) == max_edges) {
    // pattern growth limit reached, prune.
    // we do not free the memory, because the graph was added to the result list
    return 0;
  }

  // determine right-most path and last backward extension of the right-most vertex (if any)
  rightmost_vertex = ((VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i
			< VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j)
		      ? VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j
		      : VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i);
  rightmost_vertex_color = ((VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i
			      < VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j)
			    ? VECTOR(*seed_dfscode)[seed_dfscode->last_edge].l_j
			    : VECTOR(*seed_dfscode)[seed_dfscode->last_edge].l_i);
  igraph_stack_int_init(&rightmost_path, seed_dfscode->max_edges+1);
  igraph_stack_int_init(&rightmost_path_colors, seed_dfscode->max_edges+1);
  igraph_stack_int_push(&rightmost_path, rightmost_vertex);
  igraph_stack_int_push(&rightmost_path_colors, rightmost_vertex_color);
  rightmost_vertex_pred = 0;
  prev_vertex = rightmost_vertex;
  for (i = igraph_i_dfscode_size(seed_dfscode)-1; i >= 0; i--) {
    if ((VECTOR(*seed_dfscode)[i].i < VECTOR(*seed_dfscode)[i].j)
	  && (VECTOR(*seed_dfscode)[i].j == prev_vertex)) {
      igraph_stack_int_push(&rightmost_path, VECTOR(*seed_dfscode)[i].i);
      igraph_stack_int_push(&rightmost_path_colors, VECTOR(*seed_dfscode)[i].l_i);

      if (prev_vertex == rightmost_vertex)
	rightmost_vertex_pred = VECTOR(*seed_dfscode)[i].i;
      prev_vertex = VECTOR(*seed_dfscode)[i].i;
    }
  }

  // perform extensions
  while (!igraph_stack_int_empty(&rightmost_path)) {
    cur_vertex = igraph_stack_int_pop(&rightmost_path);
    cur_color = igraph_stack_int_pop(&rightmost_path_colors);

    // forward extension (to new vertex)
    new_edge = (igraph_dfscode_edge_t) {.i = cur_vertex, .j=rightmost_vertex+1,
					.l_i = cur_color, .l_ij = 0, .l_j = 0};
    if (vertex_colors != NULL) {
      for (i = 0; VECTOR(*freq_vcolors)[i] != -1; i++) {
	new_edge.l_j = VECTOR(*freq_vcolors)[i];
	if (edge_colors != NULL) {
	  for (j = 0; VECTOR(*freq_ecolors)[j] != -1; j++) {
	    new_edge.l_ij = VECTOR(*freq_ecolors)[j];
	    igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	    igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
				    min_supp, max_edges, freq_vcolors, freq_ecolors,
				    seed_dfscode, result_graph_list, result_vcolor_list,
				    result_ecolor_list, result_supp_list);
	    igraph_i_dfscode_pop_back(seed_dfscode);
	  }
	} else {
	  igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	  igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
				  min_supp, max_edges, freq_vcolors, freq_ecolors,
				  seed_dfscode, result_graph_list, result_vcolor_list,
				  result_ecolor_list, result_supp_list);
	  igraph_i_dfscode_pop_back(seed_dfscode);
	}
      }
    } else {
      if (edge_colors != NULL) {
	for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
	  new_edge.l_ij = VECTOR(*freq_ecolors)[i];
	  igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	  igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
				  min_supp, max_edges, freq_vcolors, freq_ecolors,
				  seed_dfscode, result_graph_list, result_vcolor_list,
				  result_ecolor_list, result_supp_list);
	  igraph_i_dfscode_pop_back(seed_dfscode);
	}
      } else {
	igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
				min_supp, max_edges, freq_vcolors, freq_ecolors,
				seed_dfscode, result_graph_list, result_vcolor_list,
				result_ecolor_list, result_supp_list);
	igraph_i_dfscode_pop_back(seed_dfscode);
      }
    }

    // backward extension (from right-most vertex to current node on right-most path)
    if (cur_vertex == rightmost_vertex) {
      // no self-loops
      continue;
    }
    if (cur_vertex == rightmost_vertex_pred) {
      // this edge already exists as a forward edge
      continue;
    }
    if ((VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i
	  > VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j)
	&& (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i == rightmost_vertex)
	&& (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j >= cur_vertex)) {
      // last edge was a backward edge starting from rightmost_vertex,
      // and this backward edge ended AFTER (or at) cur_vertex. a backward extension to
      // cur_vertex would result in a non-minimal DFS code (or a duplicate edge).
      continue;
    }
    new_edge = (igraph_dfscode_edge_t) {.i = rightmost_vertex, .j=cur_vertex,
					.l_i = rightmost_vertex_color, .l_ij = 0,
					.l_j = cur_color};
    if (edge_colors != NULL) {
      for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
	new_edge.l_ij = VECTOR(*freq_ecolors)[i];
	igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
				min_supp, max_edges, freq_vcolors, freq_ecolors,
				seed_dfscode, result_graph_list, result_vcolor_list,
				result_ecolor_list, result_supp_list);
	igraph_i_dfscode_pop_back(seed_dfscode);
      }
    } else {
      igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
      igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
			      min_supp, max_edges, freq_vcolors, freq_ecolors,
			      seed_dfscode, result_graph_list, result_vcolor_list,
			      result_ecolor_list, result_supp_list);
      igraph_i_dfscode_pop_back(seed_dfscode);
    }
  }
  igraph_stack_int_destroy(&rightmost_path);

  return 0;
}

// public interface

// assert: #graphs == #vertex_colors == #edge_colors
// assert: 0 <= node and edge colors <= MAX_COLOR
// assert: min_supp > 0
// assert: undirected graphs
int igraph_gspan(const igraph_vector_ptr_t *graphs, const igraph_vector_ptr_t *vertex_colors,
		 const igraph_vector_ptr_t *edge_colors, igraph_db_support_measure_t *db_supp_measure,
		 igraph_integer_t min_supp, igraph_integer_t max_edges,
		 igraph_vector_ptr_t *frequent_subgraphs,
		 igraph_vector_ptr_t *frequent_subgraph_vcolors,
		 igraph_vector_ptr_t *frequent_subgraph_ecolors,
		 igraph_vector_int_t *frequent_subgraph_supps) {

  long int graph_count = igraph_vector_ptr_size(graphs);
  long int i, j, k;
  igraph_vector_int_t vcolor_freq, ecolor_freq; // frequencies of all colors
  igraph_vector_int_t freq_ecolors, freq_vcolors; // lists of all frequent colors
  igraph_vector_int_t *vcolor, *ecolor;
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

  // BUILD ALL 1-EDGE GRAPHS AS SEEDS

  igraph_dfscode_t *pattern_dfscode;
  igraph_dfscode_edge_t pattern_dfscode_edge;
  igraph_llist_t initial_patterns;
  igraph_llist_t result_graph_list, result_vcolor_list, result_ecolor_list, result_supp_list;

  IGRAPH_CHECK(igraph_i_llist_init(&initial_patterns));
  IGRAPH_CHECK(igraph_i_llist_init(&result_graph_list));
  IGRAPH_CHECK(igraph_i_llist_init(&result_vcolor_list));
  IGRAPH_CHECK(igraph_i_llist_init(&result_ecolor_list));
  IGRAPH_CHECK(igraph_i_llist_init(&result_supp_list));

  if (vertex_colors != NULL) {
    for (i = 0; VECTOR(freq_vcolors)[i] != -1; i++) {
      for (j = 0; j <= i && VECTOR(freq_vcolors)[j] != -1; j++) {
	if (edge_colors != NULL) {
	  for (k = 0; VECTOR(freq_ecolors)[k] != -1; k++) {
	    // VC[i] -- EC[k] -- VC[j]
	    pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	    pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			      .l_i = VECTOR(freq_vcolors)[i],
			      .l_ij = VECTOR(freq_ecolors)[k],
			      .l_j = VECTOR(freq_vcolors)[j]};
	    igraph_i_dfscode_init(pattern_dfscode, max_edges);
	    igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
	    igraph_i_llist_push_back(&initial_patterns, pattern_dfscode);
	  }
	} else {
	  // VC[i] -- VC[j]
	  pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	  pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			    .l_i = VECTOR(freq_vcolors)[i],
			    .l_ij = 0,
			    .l_j = VECTOR(freq_vcolors)[j]};
	  igraph_i_dfscode_init(pattern_dfscode, max_edges);
	  igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
	  igraph_i_llist_push_back(&initial_patterns, pattern_dfscode);
	}
      }
    }
  } else {
    if (edge_colors != NULL) {
      for (k = 0; VECTOR(freq_ecolors)[k] != -1; k++) {
	// O -- EC[k] -- O
	pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			  .l_i = 0,
			  .l_ij = VECTOR(freq_ecolors)[k],
			  .l_j = 0};
	igraph_i_dfscode_init(pattern_dfscode, max_edges);
	igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
	igraph_i_llist_push_back(&initial_patterns, pattern_dfscode);
      }
    } else {
      // O -- O
      pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
      pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			.l_i = 0,
			.l_ij = 0,
			.l_j = 0};
      igraph_i_dfscode_init(pattern_dfscode, max_edges);
      igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
      igraph_i_llist_push_back(&initial_patterns, pattern_dfscode);
    }
  }

  // RECURSIVELY EXPAND ALL FREQUENT 1-EDGE GRAPHS BY PATTERN GROWTH

  igraph_llist_item_t *item;
  for (item = initial_patterns.first; item != NULL; item = item->next) {
    pattern_dfscode = (igraph_dfscode_t *) item->data;
    IGRAPH_CHECK(igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, db_supp_measure,
		    min_supp, max_edges, &freq_vcolors, &freq_ecolors,
		    pattern_dfscode, &result_graph_list,
		    &result_vcolor_list, &result_ecolor_list, &result_supp_list));
  }

  // PREPARE RESULT SET

  // count patterns in result list
  long int pattern_count = 0;
  for (item = result_graph_list.first; item != NULL; item = item->next)
    pattern_count += 1;

  // store result in the provided containers (if provided...)
  // the user has to free the allocated memory
  if (frequent_subgraphs != NULL) {
    igraph_vector_ptr_resize(frequent_subgraphs, pattern_count);
    for (item = result_graph_list.first, i = 0; item != NULL; item = item->next, i++)
      VECTOR(*frequent_subgraphs)[i] = item->data;
  }
  if (frequent_subgraph_vcolors != NULL) {
    igraph_vector_ptr_resize(frequent_subgraph_vcolors, pattern_count);
    for (item = result_vcolor_list.first, i = 0; item != NULL; item = item->next, i++)
      VECTOR(*frequent_subgraph_vcolors)[i] = item->data;
  }
  if (frequent_subgraph_ecolors != NULL) {
    igraph_vector_ptr_resize(frequent_subgraph_ecolors, pattern_count);
    for (item = result_ecolor_list.first, i = 0; item != NULL; item = item->next, i++)
      VECTOR(*frequent_subgraph_ecolors)[i] = item->data;
  }
  if (frequent_subgraph_supps != NULL) {
    igraph_vector_int_resize(frequent_subgraph_supps, pattern_count);
    for (item = result_supp_list.first, i = 0; item != NULL; item = item->next, i++)
      VECTOR(*frequent_subgraph_supps)[i] = *(long int *)item->data;
  }

  // CLEAN UP

  igraph_vector_int_destroy(&vcolor_freq);
  igraph_vector_int_destroy(&ecolor_freq);
  igraph_vector_int_destroy(&freq_vcolors);
  igraph_vector_int_destroy(&freq_ecolors);

  if (frequent_subgraphs == NULL) {
    // the user doesn't want the result, so we have to free the allocated memory
    for (item = result_graph_list.first; item != NULL; item = item->next) {
      igraph_destroy((igraph_t *) item->data);
    }
  }
  if (frequent_subgraph_vcolors == NULL) {
    // same here
    for (item = result_vcolor_list.first; item != NULL; item = item->next) {
      igraph_vector_int_destroy((igraph_vector_int_t *) item->data);
    }
  }
  if (frequent_subgraph_ecolors == NULL) {
    // same here
    for (item = result_ecolor_list.first; item != NULL; item = item->next) {
      igraph_vector_int_destroy((igraph_vector_int_t *) item->data);
    }
  }

  // free memory for temporary structures
  for (item = initial_patterns.first; item != NULL; item = item->next) {
    igraph_i_dfscode_destroy((igraph_dfscode_t *) item->data);
  }
  for (item = result_supp_list.first; item != NULL; item = item->next) {
    igraph_free(item->data);
  }

  igraph_i_llist_destroy(&initial_patterns);
  igraph_i_llist_destroy(&result_supp_list);
  igraph_i_llist_destroy(&result_graph_list);
  igraph_i_llist_destroy(&result_vcolor_list);
  igraph_i_llist_destroy(&result_ecolor_list);

  return 0;
}
