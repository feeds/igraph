/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   Frequent subgraph mining algorithms
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

/*
 * TECHNICAL TODO:
 * - Better use vertex/edge selectors and iterators. Problem: igraph_vit_create() with
 *   igraph_vs_adj() as selector uses the inefficient igraph_neighbors() function, and
 *   with igraph_vs_vector() performs an O(N) check for invalid vertex ids in the vector.
 * - Replace manual linked list iterations by iteration functions
 */

#include <stdlib.h>
#include <string.h> // memcpy
#include <limits.h> // INT_MAX
#include <omp.h>
#include "igraph_fsm.h"
#include "igraph_matrix.h"
#include "igraph_stack.h"
#include "igraph_list.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_structural.h"


// ------------- RUNTIME STATISTICS -------------

static long int igraph_fsm_stats_subiso_success_count = 0;
static long int igraph_fsm_stats_subiso_fail_count = 0;
static long int igraph_fsm_stats_subiso_failed_edge_existence_count = 0;
static long int igraph_fsm_stats_subiso_failed_edge_color_count = 0;
static long int igraph_fsm_stats_subiso_failed_edge_timestamp_count = 0;
static long int igraph_fsm_stats_subiso_failed_node_duplicate_count = 0;
static long int igraph_fsm_stats_subiso_failed_node_degree_count = 0;
static long int igraph_fsm_stats_subiso_failed_node_color_count = 0;
static long int igraph_fsm_stats_aggregated_support_count = 0;
static long int igraph_fsm_stats_mibsupport_count = 0;
static long int igraph_fsm_stats_mibsupport_subiso_success_count = 0;
static long int igraph_fsm_stats_mibsupport_subiso_fail_count = 0;
static long int igraph_fsm_stats_shallowsuppport_count = 0;
static long int igraph_fsm_stats_egobasedsuppport_count = 0;
static long int igraph_fsm_stats_noncanonical_count = 0;
static long int igraph_fsm_stats_infrequent_count = 0;
static long int igraph_fsm_stats_frequent_count = 0;
static long int igraph_fsm_stats_lfrminer_invalid_count = 0;


// ------------- HELPER FUNCTIONS -------------

void igraph_print_stats(const igraph_t *g) {
  igraph_bool_t simple, multi;
  igraph_is_simple(g, &simple);
  igraph_has_multiple(g, &multi);
  printf("vcount %ld ecount %ld directed %d simple %d multi %d c %.1f\n",
      (long int) igraph_vcount(g), (long int) igraph_ecount(g), (int)igraph_is_directed(g),
      simple, multi, 2.*igraph_ecount(g)/igraph_vcount(g));
}

void igraph_print(const igraph_t *g, const igraph_vector_int_t *vcolors,
		    const igraph_vector_int_t *ecolors, const igraph_vector_int_t *etimes) {
  long int i;
  for (i = 0; i < igraph_ecount(g); i++) {
    printf("%ld(%d) --%d[%d]--%s %ld(%d)\n",
		(long int) VECTOR(g->from)[i],
		((vcolors != NULL) ? VECTOR(*vcolors)[(long int) VECTOR(g->from)[i]] : 0),
		((ecolors != NULL) ? VECTOR(*ecolors)[i] : 0),
		((etimes != NULL) ? VECTOR(*etimes)[i] : 0),
		(igraph_is_directed(g) ? ">": ""),
		(long int) VECTOR(g->to)[i],
		((vcolors != NULL) ? VECTOR(*vcolors)[(long int) VECTOR(g->to)[i]] : 0));
  }
}


int igraph_write_colored_graph(igraph_t *g, igraph_vector_int_t *vcolors,
      igraph_vector_int_t *ecolors, igraph_vector_int_t *etimes, FILE *f) {
  long int i;
  igraph_integer_t from, to;
  for (i = 0; i < igraph_vcount(g); i++) {
    if (vcolors != NULL)
      fprintf(f, "v %ld %d\n", i, VECTOR(*vcolors)[i]);
    else
      fprintf(f, "v %ld\n", i);
  }
  for (i = 0; i < igraph_ecount(g); i++) {
    IGRAPH_CHECK(igraph_edge(g, i, &from, &to));
    fprintf(f, "e %ld %ld", (long int) from, (long int) to);
    if (ecolors != NULL) {
      fprintf(f, " %d", VECTOR(*ecolors)[i]);
      if (etimes != NULL) {
	fprintf(f, " %d", VECTOR(*etimes)[i]);
      }
    } else {
      if (etimes != NULL) {
	fprintf(f, " %d", VECTOR(*etimes)[i]);
      }
    }
    fprintf(f, "\n");
  }
  return 0;
}


int igraph_write_colored_graph_gz(igraph_t *g, igraph_vector_int_t *vcolors,
      igraph_vector_int_t *ecolors, igraph_vector_int_t *etimes, gzFile f) {
  long int i;
  igraph_integer_t from, to;
  for (i = 0; i < igraph_vcount(g); i++) {
    if (vcolors != NULL)
      gzprintf(f, "v %ld %d\n", i, VECTOR(*vcolors)[i]);
    else
      gzprintf(f, "v %ld\n", i);
  }
  for (i = 0; i < igraph_ecount(g); i++) {
    IGRAPH_CHECK(igraph_edge(g, i, &from, &to));
    gzprintf(f, "e %ld %ld", (long int) from, (long int) to);
    if (ecolors != NULL) {
      gzprintf(f, " %d", VECTOR(*ecolors)[i]);
      if (etimes != NULL) {
	gzprintf(f, " %d", VECTOR(*etimes)[i]);
      }
    } else {
      if (etimes != NULL) {
	gzprintf(f, " %d", VECTOR(*etimes)[i]);
      }
    }
    gzprintf(f, "\n");
  }
  return 0;
}


// graph1 is the larger graph, graph2 is the smaller graph
// Can handle a single fixed assignment (pattern node, target node) passed as a length-2 vector
// 
// Unconnected components can be allowed, but should not be 
// used if it can be avoided. In that setting, a fixed assignment is prohibited.
//
// DFS version of the RI algorithm by Bonnici et al.: "A subgraph isomorphism algorithm
// and its application to biochemical data", in: BMC Bioinformatics, vol. 14(7), 2013.
// We assume that DFS performs better than the original GreatestConstraintFirst search
// strategy, because it guarantees (for connected graphs) that every pattern node has a
// well-defined parent. In the original paper, parents can be undefined, which means that
// all target nodes are candidates.
//
// Algorithm:
//    1) build a DFS ordering of the pattern nodes (intuition: when matching the next node
//       starting from a partial solution, we only have to consider the neighbors of the
//       current node's parent, and avoid an open search over all possible target assignments)
//    2) match all pattern nodes in the order specified by the DFS ordering, using the neighbors
//       of their DFS parents as candidates, while maintaining the subgraph isomorphism
//       properties for every partial solution (matching node labels, matching degrees, matching
//       edges, in that order)
//
// The implementation makes heavy use of the internals of the igraph_t data type.
int igraph_i_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamp1,
			   const igraph_vector_int_t *edge_timestamp2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *fixed,
         igraph_bool_t allow_unconnected,
			   igraph_bool_t *iso) {
  long int vcount1 = igraph_vcount(graph1);
  long int vcount2 = igraph_vcount(graph2);
  long int i, j, k, fixed_count, partial_solution_pos, pred;
  long int pattern_node, other_pattern_node, target_node, other_target_node;
  long int indeg1, indeg2, outdeg1, outdeg2;
  igraph_integer_t eid1, eid2;
  int end, success;

  igraph_vector_t membership;
  igraph_vector_t cluster_size;
  igraph_integer_t number_of_clusters; 
  igraph_vector_int_t cluster_max_deg;
  igraph_vector_int_t cluster_start_node;

  igraph_vector_t node_ordering;
  igraph_vector_t pred_idx;
  igraph_vector_t visited;
  igraph_vector_t state_target_idx;
  igraph_vector_t state_target_node;
  igraph_stack_t dfs_node_stack;
  igraph_stack_t dfs_pred_stack;
  igraph_bool_t directed;

  // initialize different variants
  long int germ_delta = 0;
  long int lfrminer_se_timestamp = 0;
  if (((variant == IGRAPH_GSPAN_LFRMINER) || (variant == IGRAPH_GSPAN_GERM))
	&& (edge_timestamp1 == NULL)) {
    IGRAPH_ERROR("GERM and LFR-Miner need edge timestamps!", IGRAPH_EINVAL);
  }

  *iso = 0;
  end = 0;
  success = 1;
  directed = igraph_is_directed(graph1);

  // STEP -1: Identify weakly connected clusters
  IGRAPH_CHECK(igraph_vector_init(&membership, 0));
  IGRAPH_CHECK(igraph_vector_init(&cluster_size, 0));

  // generate (weakly connected) clusters
  igraph_clusters(graph2, &membership, &cluster_size, &number_of_clusters,
                  IGRAPH_WEAK);
 
  // if you only care about connected graphs and have several clusters. Abort.
  if(!allow_unconnected && number_of_clusters > 1){
    *iso = 0;
    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&cluster_size);
    return 1;
  }

  // STEP 0: create a static ordering of the pattern nodes by DFS
  // if a fixed assignment is given, use this node as root, otherwise take the one with
  // the largest degree (heuristic for better pruning from RI algorithm)
  IGRAPH_CHECK(igraph_stack_init(&dfs_node_stack, vcount2*vcount2));
  IGRAPH_CHECK(igraph_stack_init(&dfs_pred_stack, vcount2*vcount2));
  IGRAPH_CHECK(igraph_vector_init(&node_ordering, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&pred_idx, vcount2));
  IGRAPH_CHECK(igraph_vector_init(&visited, vcount2));

  IGRAPH_CHECK(igraph_vector_int_init(&cluster_max_deg,number_of_clusters));
  IGRAPH_CHECK(igraph_vector_int_init(&cluster_start_node,number_of_clusters));
 
  igraph_vector_int_fill(&cluster_max_deg,-1);


  if(allow_unconnected && fixed != NULL){ 
    // ignore fixed s/e-vertexes in multiple-cluster setting (might be good idea to allow this)
    printf("Fixed start and end edges are not allowed in the multiple-cluster setting.\n");
    return 1;
  }

  if(fixed != NULL){
    IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, fixed[0]));
    number_of_clusters = 1;
    fixed_count = 1;
  } else {
    // find max-degree nodes for all clusters
    int cluster_id;
    for (i = 0; i < vcount2; i++) {
      cluster_id = VECTOR(membership)[i];

      if (DEGREE(*graph2, i) > VECTOR(cluster_max_deg)[cluster_id]) { 
        VECTOR(cluster_max_deg)[cluster_id] = DEGREE(*graph2, i);
        VECTOR(cluster_start_node)[cluster_id] = i;
      } 
    }

    // TODO (!!) performance: order by cluster size?
    IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, VECTOR(cluster_start_node)[0])); // push vertex of first cluster first
    fixed_count = 0;
  }

  // assert: one node is already pushed to node_stack
  IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, -1)); // first node has no predecessor

  i = 0; // start DFS
  for(k = 0; k < number_of_clusters; k++){ // start a DFS for each cluster
    
    if(k != 0){ // fetch a new starting node
      IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, VECTOR(cluster_start_node)[k]));
      IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, -1)); // node from new cluster doesnt have a predecessor
    }     

    // DFS for one cluster
    while (!igraph_stack_empty(&dfs_node_stack)) {
      pattern_node = (long int) igraph_stack_pop(&dfs_node_stack);
      pred = (long int) igraph_stack_pop(&dfs_pred_stack);
      if (VECTOR(visited)[pattern_node] == 1)
        continue;

      // insert current node into ordering and set predecessor index
      VECTOR(node_ordering)[i] = pattern_node;
      VECTOR(pred_idx)[i] = pred;

      // add neighbors to stack
      // TODO: in the order specified by RI's scoring functions
      for (j = 0; j < DEGREE(*graph2, pattern_node); j++) {
        if (VECTOR(visited)[NEIGHBOR(*graph2, pattern_node, j)] == 1) {
          continue;
        }
        if ((variant == IGRAPH_GSPAN_LFRMINER)
  	    && VECTOR(*vertex_color2)[NEIGHBOR(*graph2, pattern_node, j)] == 1) {
  	// only the e node has label 1, and it is treated specially (see below)
  	continue;
        }
        IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, NEIGHBOR(*graph2, pattern_node, j)));
        IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, i));
      }

      if ((variant == IGRAPH_GSPAN_LFRMINER) && (i == 0)) {
        // we just added the s node at the first position of the DFS node ordering,
        // and added its neighbors (except for e node) to the stack.
        // The next node in the ordering must be the e node (which is a neighbor), so we
        // put in on top of the stack.
        for (j = 0; j < DEGREE(*graph2, pattern_node); j++) {
  	if (VECTOR(*vertex_color2)[NEIGHBOR(*graph2, pattern_node, j)] == 1) {
  	  IGRAPH_CHECK(igraph_stack_push(&dfs_node_stack, NEIGHBOR(*graph2, pattern_node, j)));
  	  IGRAPH_CHECK(igraph_stack_push(&dfs_pred_stack, i));
  	  break;
  	}
        }
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
    pattern_node = fixed[0];
    target_node = fixed[1];

    // check color
    if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			    != VECTOR(*vertex_color2)[pattern_node])) {
      end = 1;
      igraph_fsm_stats_subiso_failed_node_color_count++;
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
	  igraph_fsm_stats_subiso_failed_node_degree_count++;
	}
      } else {
	if (indeg1+outdeg1 < indeg2+outdeg2) {
	  end = 1;
	  igraph_fsm_stats_subiso_failed_node_degree_count++;
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
      if (partial_solution_pos == 1) {
	// this will be the first matching edge, initialize data for GERM and LFR-Miner
	germ_delta = -1;
	lfrminer_se_timestamp = -1;
      }

      success = 1;

      pattern_node = VECTOR(node_ordering)[partial_solution_pos]; // get the node from working dfs pos
      if (partial_solution_pos == 0 
          || (long int) VECTOR(pred_idx)[partial_solution_pos] == -1) { // fb change here?
	// target index is actual target node
	VECTOR(state_target_node)[partial_solution_pos] = VECTOR(state_target_idx)[partial_solution_pos];
      } else {
	// target index is index in parent's neighborhood
	VECTOR(state_target_node)[partial_solution_pos] = NEIGHBOR(*graph1,
	    (long int) VECTOR(state_target_node)[(long int) VECTOR(pred_idx)[partial_solution_pos]],
	    (long int) VECTOR(state_target_idx)[partial_solution_pos]);
      }
      target_node = VECTOR(state_target_node)[partial_solution_pos];

      //////////////////////////////////////////////
      //  CHECK CONSTRAINTS
      /////////////////////////////////////////////

      // check colors
      if (vertex_color1 && (VECTOR(*vertex_color1)[target_node]
			      != VECTOR(*vertex_color2)[pattern_node])) {
	success = 0;
	igraph_fsm_stats_subiso_failed_node_color_count++;
      }

      // check whether target node has been matched before
      for (i = 0; success && i < partial_solution_pos; i++) {
	other_target_node = (long int) VECTOR(state_target_node)[i];
	if (other_target_node == target_node) {
	  success = 0;
	  igraph_fsm_stats_subiso_failed_node_duplicate_count++;
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
	    igraph_fsm_stats_subiso_failed_node_degree_count++;
	  }
	} else {
	  if (indeg1+outdeg1 < indeg2+outdeg2) {
	    success = 0;
	    igraph_fsm_stats_subiso_failed_node_degree_count++;
	  }
	}
      }

      // check edges to already matched nodes
      for (i = 0; success && i < partial_solution_pos; i++) {
	other_pattern_node = VECTOR(node_ordering)[i];
	other_target_node = VECTOR(state_target_node)[i];

	igraph_get_eid(graph2, &eid2, other_pattern_node, pattern_node, /*directed*/1, /*err*/0);
	if (eid2 > -1) {
	  igraph_get_eid(graph1, &eid1, other_target_node, target_node, /*directed*/1, /*err*/0);
	  if (eid1 == -1) {
	    success = 0;
	    igraph_fsm_stats_subiso_failed_edge_existence_count++;
	  } else {
	    // check edge timestamps (for some variants only)
	    if (variant == IGRAPH_GSPAN_GERM) {
	      // edge timestamps have to match with a fixed time gap that is
	      // determined by the first matched edge
	      if ((i == 0) && (partial_solution_pos == 1)) {
		// this is the first edge we are matching, store the time gap
		germ_delta = (VECTOR(*edge_timestamp1)[(long int)eid1]
				  - VECTOR(*edge_timestamp2)[(long int)eid2]);
		if (germ_delta < 0) {
		  success = 0;
		  igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
		}
	      }
	      if (success && (VECTOR(*edge_timestamp1)[(long int)eid1] !=
		    VECTOR(*edge_timestamp2)[(long int)eid2] + germ_delta)) {
		success = 0;
		igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
	      }
	    } else if (variant == IGRAPH_GSPAN_LFRMINER) {
	      // edge timestamps have to be smaller than the timestamp of the (s, e) edge
	      if ((i == 0) && (partial_solution_pos == 1)) {
		// this is the (s,e) edge, because we forced node e to be at pos 1 in ordering!
		// store edge timestamp
		lfrminer_se_timestamp = VECTOR(*edge_timestamp1)[(long int)eid1];
	      } else if (lfrminer_se_timestamp < 0) {
		// no (s,e) edge matched before? should never happen
		printf("strange\n");
		success = 0;
	      } else if (VECTOR(*edge_timestamp1)[(long int)eid1] >= lfrminer_se_timestamp) {
		success = 0;
		igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
	      }
	    }

	    // check edge color
	    if (success && edge_color1 && (VECTOR(*edge_color1)[(long int)eid1] !=
	          VECTOR(*edge_color2)[(long int)eid2])) {
	      success = 0;
	      igraph_fsm_stats_subiso_failed_edge_color_count++;
	    }
	  }
	} else {
	  if ((variant == IGRAPH_GSPAN_LFRMINER) && (i == 0) && (partial_solution_pos == 1)) {
	    // no (s,e) edge present, fail
	    success = 0;
	    igraph_fsm_stats_subiso_failed_edge_existence_count++;
	  }
	  if (success && induced) {
	    igraph_get_eid(graph1, &eid1, other_target_node, target_node, /*directed*/1, /*err*/0);
	    if (eid1 > -1) {
	      success = 0;
	      igraph_fsm_stats_subiso_failed_edge_existence_count++;
	    }
	  }
	}

	if (success && directed) { // check the other edge direction, too
	  igraph_get_eid(graph2, &eid2, pattern_node, other_pattern_node, /*directed*/1, /*err*/0);
	  if (eid2 > -1) {
	    igraph_get_eid(graph1, &eid1, target_node, other_target_node, /*directed*/1, /*err*/0);
	    if (eid1 == -1) {
	      success = 0;
	      igraph_fsm_stats_subiso_failed_edge_existence_count++;
	    } else {
	      // check edge timestamps (for some variants only)
	      if (variant == IGRAPH_GSPAN_GERM) {
		// edge timestamps have to match with a fixed time gap that is
		// determined by the first matched edge
		if ((i == 0) && (partial_solution_pos == 1) && (germ_delta < 0)) {
		  // this is the first edge we are matching, store the time gap
		  germ_delta = (VECTOR(*edge_timestamp1)[(long int)eid1]
				    - VECTOR(*edge_timestamp2)[(long int)eid2]);
		  if (germ_delta < 0) {
		    success = 0;
		    igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
		  }
		}
		if (success && (VECTOR(*edge_timestamp1)[(long int)eid1] !=
		      VECTOR(*edge_timestamp2)[(long int)eid2] + germ_delta)) {
		  success = 0;
		  igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
		}
	      } else if (variant == IGRAPH_GSPAN_LFRMINER) {
		// edge timestamps have to be smaller than the timestamp of the (s, e) edge
		if (lfrminer_se_timestamp < 0) {
		  // no (s,e) matched before. we should never end up here.
		  printf("strange\n");
		  success = 0;
		} else if (VECTOR(*edge_timestamp1)[(long int)eid1] >= lfrminer_se_timestamp) {
		  success = 0;
		  igraph_fsm_stats_subiso_failed_edge_timestamp_count++;
		}
	      }

	      // check edge color
	      if (success && edge_color1 && (VECTOR(*edge_color1)[(long int)eid1] !=
	            VECTOR(*edge_color2)[(long int)eid2])) {
		success = 0;
		igraph_fsm_stats_subiso_failed_edge_color_count++;
	      }
	    }
	  } else if (induced) {
	    igraph_get_eid(graph1, &eid1, target_node, other_target_node, /*directed*/1, /*err*/0);
	    if (eid1 > -1) {
	      success = 0;
	      igraph_fsm_stats_subiso_failed_edge_existence_count++;
	    }
	  }
	}
      } // for i (other matched nodes)

      /////////////////////////////////////
      // NEXT NODE SELECTION
      ////////////////////////////////////


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


  // TODO: perform backtracking IF
  // necessary condition: partial_solution_pos > 0
  // if :normal node: all child nodes must have been tried
  // if :first (cluster) node: all nodes from big graph have been tried.
	while ((partial_solution_pos > 0) // and if NOT -1 ! -> consider other nodes as candidates, like with 1
		&& ((( VECTOR(pred_idx)[partial_solution_pos] > -1 ) 
        && (VECTOR(state_target_idx)[partial_solution_pos] == DEGREE(*graph1, (long int) VECTOR(state_target_node)[
                                                    (long int) VECTOR(pred_idx)[partial_solution_pos]])-1))
        || // or
        ((VECTOR(pred_idx)[partial_solution_pos] == -1) 
          && VECTOR(state_target_idx)[partial_solution_pos] == vcount1-1))
        ){
	  // all nodes from parent's neighborhood have been tried (or all nodes have been tried for cluster starter)
    // perform backtracking
	  partial_solution_pos--;
	}

  // TODO: here we chack only if he first node is without a predecessor -> aply that learning for other ones as well
  // ? how can we make sure we dont visit a previously visited edge
  // ? We work with NO induced graphs! for this setting (also without induced?)
	if (partial_solution_pos == 0 
    ||(long int) VECTOR(pred_idx)[partial_solution_pos] == -1 // fb changed
    ) {  
	  if (fixed_count == 0) {
	    // special case: first position, no fixed nodes -> all target nodes are candidates
	    if (VECTOR(state_target_idx)[partial_solution_pos] == vcount1-1) {
	      
        // backtrack to next unfully elaborated cluster
        while((partial_solution_pos > 0) && 
          (((long int) VECTOR(pred_idx)[partial_solution_pos] != -1)
               || (VECTOR(state_target_idx)[partial_solution_pos] == vcount1-1))){
          partial_solution_pos--;
        }

        if(partial_solution_pos == 0 && VECTOR(state_target_idx)[partial_solution_pos] == vcount1-1){
          // all target nodes have been tried, no candidates left
          break; 
        } else {
          // try new start node for cluster
          VECTOR(state_target_idx)[partial_solution_pos] += 1;
        }
        
	      
	    } else {
	      // try next target node
	      VECTOR(state_target_idx)[partial_solution_pos] += 1;
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
    *iso = 1;
    igraph_fsm_stats_subiso_success_count++;
  } else {
    *iso = 0;
    igraph_fsm_stats_subiso_fail_count++;
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

int igraph_subisomorphic_evomine(const igraph_t *graph1, const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamp1,
			   const igraph_vector_int_t *edge_timestamp2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *fixed,
         igraph_bool_t allow_unconnected,
			   igraph_bool_t *iso) {
  return igraph_i_subisomorphic(graph1, graph2,
			   vertex_color1,
			   vertex_color2,
			   edge_color1,
			   edge_color2,
			   edge_timestamp1,
			   edge_timestamp2,
			   induced,
			   variant,
			   variant_data,
			   fixed,
         allow_unconnected,
			   iso);
}

int igraph_isomorphic_evomine(const igraph_t *graph1, const igraph_t *graph2,
         const igraph_vector_int_t *vertex_color1,
         const igraph_vector_int_t *vertex_color2,
         const igraph_vector_int_t *edge_color1,
         const igraph_vector_int_t *edge_color2,
         /* igraph_bool_t induced, */
         /*igraph_gspan_variant_t variant,*/ // TODO erik: only normal, right? 
         /* void *variant_data, */
         /* long int *fixed, */
         igraph_bool_t *iso){

  *iso = 0;
  long int vcount1 = igraph_vcount(graph1);
  long int vcount2 = igraph_vcount(graph2);

  if(vcount2 != vcount1){ // no isomorphism
    return 0;
  }

  long int ecount1 = igraph_ecount(graph1);
  long int ecount2 = igraph_ecount(graph2);

  if(ecount2 != ecount1){ // no isomorphism
    return 0;
  }

  // the graphs counts are identical 
  // if one is the subgraph of the other, it will be a full isomorphism
  return igraph_i_subisomorphic(graph1, graph2,
                                       vertex_color1,
                                       vertex_color2,
                                       edge_color1,
                                       edge_color2,
                                       NULL,
                                       NULL,
                                       1,
                                       0,
                                       NULL,
                                       NULL,
                                       1,
                                       iso);

}



// ------------- SUPPORT MEASURES -------------


// graph1 is the larger graph, graph2 is the smaller graph.
int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       const igraph_vector_int_t *edge_timestamp1,
		       const igraph_vector_int_t *edge_timestamp2,
		       igraph_bool_t induced,
		       igraph_gspan_variant_t variant,
		       void *variant_data,
		       long int *support,
		       igraph_integer_t min_supp) {
  long int fixed[2];
  igraph_bool_t iso;
  int abort_status;
  long int vcount1 = igraph_vcount(graph1), vcount2 = igraph_vcount(graph2);
  long int i, j, automorphic_node;
  long int cur_supp, cur_count;

  // find all automorphic pattern nodes
  igraph_matrix_t automorphic_nodes;
  IGRAPH_CHECK(igraph_matrix_init(&automorphic_nodes, vcount2, vcount2));
  for (i = 0; i < vcount2; i++) {
    fixed[0] = i; // force assignment: pattern node i
    for (j = 0; j < i; j++) {
      fixed[1] = j; // force assignment: pattern node j
      iso = 0;
      if (igraph_i_subisomorphic(graph2, graph2, vertex_color2, vertex_color2, edge_color2,
	      edge_color2, edge_timestamp2, edge_timestamp2, induced, variant /*for GERM*/,
	      /*variant_data=*/ NULL, fixed, 0, &iso)) {
	igraph_matrix_destroy(&automorphic_nodes);
	return 1;
      }
      if (iso) {
	MATRIX(automorphic_nodes, i, j) = 1;
      }
    }
  }

  // test all possible pairs (pattern node i, target node j) for isomorphism
  cur_supp = -1;
  abort_status = 0; // 0=run, 1=early termination, 2=fail
  #pragma omp parallel for private(j, automorphic_node, cur_count, iso, fixed)
  for (i = 0; i < vcount2; i++) {
    if (abort_status > 0)
      continue; // for OpenMP

    fixed[0] = i; // force assignment: pattern node i

    // check if this node is isomorphic to a previously checked node
    automorphic_node = -1;
    for (j = 0; j < i; j++) {
      if (MATRIX(automorphic_nodes, i, j) == 1) {
	automorphic_node = j;
	break;
      }
    }
    if (automorphic_node >= 0) {
      continue;
    }

    // no automorphic node found, test all possible target assignments
    cur_count = 0;
    for (j = 0; j < vcount1; j++) {
      fixed[1] = j; // force assignment: target node j
      iso = 0;
      if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2, edge_color1,
	      edge_color2, edge_timestamp1, edge_timestamp2, induced, variant, variant_data,
	      fixed,0, &iso)) {
	#pragma omp critical (status)
	abort_status = 2;
	#pragma omp flush(abort_status)
      }
      if (iso) {
	igraph_fsm_stats_mibsupport_subiso_success_count++;
        cur_count++;

        // early termination: if we have already found at least min_supp
        // assignments, we can continue with the next pattern node
        // NOTE: the reported support values will always be equal to min_supp
        //if (min_supp >= 0 && cur_count >= min_supp) {
        //  break;
        //}

        // early pruning: if we have already found more target node assignments for
        // the current pattern node than for some previous pattern node, the current
        // pattern node is irrelevant for the support and we can stop trying more
        // candidates
        if (cur_supp >= 0 && cur_count >= cur_supp) {
          //printf("prune\n");
          break;
        }
      } else {
	igraph_fsm_stats_mibsupport_subiso_fail_count++;
      }
    }

    // keep track of current minimum number of target node assignments
    #pragma omp critical (update_support)
    {
      if (cur_supp < 0 || cur_count < cur_supp) {
	cur_supp = cur_count;
      }

      // early termination: the support can only be smaller than or equal to
      // VECTOR(target_counts)[i].
      // if that value is already smaller than min_supp, we don't need to continue.
      if (min_supp >= 0 && cur_supp < min_supp) {
	//printf("terminate\n");
	#pragma omp critical (status)
	abort_status = 1; // early termination
	#pragma omp flush(abort_status)
      }
    }
  }

  igraph_matrix_destroy(&automorphic_nodes);

  if (abort_status == 1) {
    // early termination: a pattern node was discovered that can be matched to less
    // than min_supp target nodes, so the pattern cannot be frequent. however, there
    // might be nodes with even less target nodes, which we did not check.
    *support = 0;
    return 0;
  }

  if (abort_status == 2) {
    // an error occurred during subisomorphism check
    return 1;
  }

  *support = cur_supp;
  igraph_fsm_stats_mibsupport_count++;

  return 0;
}


// graph1 is the larger graph, graph2 is the smaller graph
//
// NOTE: Only use this in conjunction with IGRAPH_GSPAN_LFRMINER!
// The code assumes graph2 contains a single node with color 0 (start node),
// a single node with color 1 (end node), and all other nodes should have color 2.
// graph1 should have no node labels, and edge labels that encode timestamps.
int igraph_egobased_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamp1,
			   const igraph_vector_int_t *edge_timestamp2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *support,
			   igraph_integer_t min_supp) {
  long int i;
  long int fixed[2], start_node;
  igraph_bool_t iso, abort_status;

  // determine start node in graph2 (node label 0)
  start_node = -1;
  for (i = 0; i < igraph_vcount(graph2); i++) {
    if (VECTOR(*vertex_color2)[i] == 0) {
      start_node = i;
      break;
    }
  }
  if (start_node < 0) {
    IGRAPH_ERROR("no start node in pattern (set vcolor = 2 for exactly one node))", IGRAPH_EINVAL);
  }

  // check for all possible target nodes whether they can be used as
  // a start node for the pattern in graph2
  *support = 0;
  abort_status = 0;
  #pragma omp parallel for private(fixed, iso)
  for (i = 0; i < igraph_vcount(graph1); i++) {
    if (abort_status)
      continue;

    fixed[0] = start_node;
    fixed[1] = i;
    if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2,
		  edge_color1, edge_color2, edge_timestamp1, edge_timestamp2, induced,
		  variant, variant_data, fixed,0, &iso)) {
      abort_status = 1;
      #pragma omp flush(abort_status)
    }
    if (iso) {
      #pragma omp atomic
      (*support)++;
    }

    // early pruning
    //if ((*support+igraph_vcount(graph1)-i-1) < min_supp) {
      // support cannot become larger than min_supp anymore
    //  return 0;
    //}
  }

  if (abort_status) {
    // an error occurred during subisomorphism check
    return 1;
  }

  igraph_fsm_stats_egobasedsuppport_count++;
  return 0;
}


// graph1 is the larger graph, graph2 is the smaller graph
int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamp1,
			   const igraph_vector_int_t *edge_timestamp2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *support,
			   igraph_integer_t min_supp) {
  igraph_bool_t iso;
  if (igraph_i_subisomorphic(graph1, graph2, vertex_color1, vertex_color2,
		edge_color1, edge_color2, edge_timestamp1, edge_timestamp2, induced,
		variant, variant_data, /*fixed=*/ NULL, 0,&iso)) {
    return 1;
  }
  if (iso) {
    *support = 1;
  } else {
    *support = 0;
  }
  igraph_fsm_stats_shallowsuppport_count++;
  return 0;
}


int igraph_aggregated_db_support(const igraph_vector_ptr_t *graphs,
			  const igraph_vector_ptr_t *vertex_colors,
			  const igraph_vector_ptr_t *edge_colors,
			  const igraph_vector_ptr_t *edge_times,
			  const igraph_t *pattern,
			  const igraph_vector_int_t *pattern_vcolors,
			  const igraph_vector_int_t *pattern_ecolors,
			  const igraph_vector_int_t *pattern_etimes,
			  igraph_bool_t induced,
			  igraph_gspan_variant_t variant,
			  void *variant_data,
			  igraph_support_measure_t single_graph_support,
			  long int *support,
			  igraph_integer_t min_supp) {
  long int i;
  long int gsupp;
  igraph_vector_int_t *vcolors = NULL, *ecolors = NULL, *etimes = NULL;
  *support = 0;

  for (i = 0; i < igraph_vector_ptr_size(graphs); i++) {
    if (vertex_colors != NULL)
      vcolors = (igraph_vector_int_t *) VECTOR(*vertex_colors)[i];
    if (edge_colors != NULL)
      ecolors = (igraph_vector_int_t *) VECTOR(*edge_colors)[i];
    if (edge_times != NULL)
      etimes = (igraph_vector_int_t *) VECTOR(*edge_times)[i];

    single_graph_support((igraph_t *) VECTOR(*graphs)[i], pattern,
			  vcolors, pattern_vcolors,
			  ecolors, pattern_ecolors,
			  etimes, pattern_etimes,
			  induced, variant, variant_data, &gsupp, /*min_supp=*/ 0);
    *support += gsupp;
  }

  igraph_fsm_stats_aggregated_support_count++;
  return 0;
}


// ------------- GSPAN -------------

// DFS code related (internal)

typedef struct igraph_dfscode_edge_t {
  long int i; // source node
  long int j; // target node
  long int d; // direction, where 0 means e=(i,j) and 1 means e=(j,i)
  long int l_i; // source label
  long int l_ij; // edge label
  long int t_ij; // edge timestamp
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
igraph_bool_t igraph_i_dfscode_contains_edge(igraph_dfscode_t *dfscode, long int v1, long int v2);
int igraph_i_dfscode_to_graph(const igraph_dfscode_t *dfscode, igraph_bool_t directed,
		igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		igraph_vector_int_t *edge_timestamps);
igraph_bool_t igraph_i_dfscode_is_canonical(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		igraph_vector_int_t *edge_timestamps);
igraph_bool_t igraph_i_dfscode_is_canonical_rec(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		igraph_vector_int_t *edge_timestamps,
		long int dfscode_pos, igraph_vector_int_t *ordering, long int ordering_pos,
		igraph_vector_int_t *visited_nodes);
int igraph_i_dfscode_extend(const igraph_vector_ptr_t *graphs,
		const igraph_vector_ptr_t *vertex_colors, const igraph_vector_ptr_t *edge_colors,
		const igraph_vector_ptr_t *edge_times,
		igraph_support_measure_t *single_graph_support,
		igraph_integer_t min_supp, igraph_integer_t max_edges,
		igraph_gspan_variant_t variant, void *variant_data,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors,
		igraph_dfscode_t *seed_dfscode, igraph_llist_ptr_t *result_graph_list,
		igraph_llist_ptr_t *result_vcolor_list,
		igraph_llist_ptr_t *result_ecolor_list, igraph_llist_ptr_t *result_etimes_list,
		igraph_llist_long_t *result_supp_list);
int igraph_i_vector_ptr_to_vector_int_minmax(const igraph_vector_ptr_t *vec_ptr, long int *min,
		long int *max);
int igraph_i_frequent_colors(const igraph_vector_ptr_t *vertex_colors,
		const igraph_vector_ptr_t *edge_colors, long int min_supp,
		long int max_vcolor, long int max_ecolor,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors);
int igraph_i_build_seeds_default(igraph_bool_t has_vcolors, igraph_bool_t has_ecolors,
				 const igraph_vector_int_t *freq_vcolors,
				 const igraph_vector_int_t *freq_ecolors,
				 igraph_integer_t max_edges,
				 igraph_gspan_variant_t variant,
				 void *variant_data,
				 igraph_llist_ptr_t *initial_patterns);
int igraph_i_build_seeds_lfrminer(igraph_bool_t has_ecolors,
				 const igraph_vector_int_t *freq_ecolors,
				 igraph_integer_t max_edges,
				 igraph_llist_ptr_t *initial_patterns);
igraph_bool_t igraph_i_lfrminer_valid(igraph_t *g);

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
    printf("(%ld,%ld,%s,%ld,%ld,%ld,%ld) ", VECTOR(*dfscode)[i].i, VECTOR(*dfscode)[i].j,
	((VECTOR(*dfscode)[i].d == 0) ? ">" : "<"),
	VECTOR(*dfscode)[i].l_i,
	VECTOR(*dfscode)[i].l_ij, VECTOR(*dfscode)[i].t_ij,
	VECTOR(*dfscode)[i].l_j);
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

// definition from CloseGraph paper (Yan & Han 2003),
// extended by edge directions and timestamps
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
    // check direction of edge
    if (a->d < b->d) {
      return -1; // a < b
    }
    if (a->d == b->d) {
      // check label of node i
      if (a->l_i < b->l_i) {
	return -1; // a < b
      }
      if (a->l_i == b->l_j) {
	// check label of edge (i,j)
	if (a->l_ij < b->l_ij) {
	  return -1; // a < b
	}
	if (a->l_ij == b->l_ij) {
	  // check timestamp of edge (i,j)
	  if (a->t_ij < b->t_ij) {
	    return -1; // a < b
	  }
	  if (a->t_ij == b->t_ij) {
	    // check label of node j
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


igraph_bool_t igraph_i_dfscode_contains_edge(igraph_dfscode_t *dfscode, long int v1, long int v2) {
  long int i;
  for (i = 0; i <= dfscode->last_edge; i++) {
    if (((VECTOR(*dfscode)[i].i == v1) && (VECTOR(*dfscode)[i].j == v2))
	|| ((VECTOR(*dfscode)[i].i == v2) && (VECTOR(*dfscode)[i].j == v1))) {
      return 1;
    }
  }
  return 0;
}


int igraph_i_dfscode_to_graph(const igraph_dfscode_t *dfscode, igraph_bool_t directed,
		igraph_t *graph, igraph_vector_int_t *vertex_colors,
		igraph_vector_int_t *edge_colors, igraph_vector_int_t *edge_timestamps) {
  igraph_vector_t edges;
  long int i;
  long int rightmost_vertex = ((VECTOR(*dfscode)[dfscode->last_edge].i
				< VECTOR(*dfscode)[dfscode->last_edge].j)
			      ? VECTOR(*dfscode)[dfscode->last_edge].j
			      : VECTOR(*dfscode)[dfscode->last_edge].i);

  IGRAPH_CHECK(igraph_empty(graph, rightmost_vertex+1,
		    (directed ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED)));
  IGRAPH_CHECK(igraph_vector_int_init(vertex_colors, rightmost_vertex+1));
  IGRAPH_CHECK(igraph_vector_int_init(edge_colors, igraph_i_dfscode_size(dfscode)));
  IGRAPH_CHECK(igraph_vector_int_init(edge_timestamps, igraph_i_dfscode_size(dfscode)));
  IGRAPH_CHECK(igraph_vector_init(&edges, 2*igraph_i_dfscode_size(dfscode)));

  for (i = 0; i < igraph_i_dfscode_size(dfscode); i++) {
    if (VECTOR(*dfscode)[i].d == 0) {
      VECTOR(edges)[2*i] = VECTOR(*dfscode)[i].i;
      VECTOR(edges)[2*i+1] = VECTOR(*dfscode)[i].j;
    } else {
      VECTOR(edges)[2*i] = VECTOR(*dfscode)[i].j;
      VECTOR(edges)[2*i+1] = VECTOR(*dfscode)[i].i;
    }
    VECTOR(*edge_colors)[i] = VECTOR(*dfscode)[i].l_ij;
    VECTOR(*edge_timestamps)[i] = VECTOR(*dfscode)[i].t_ij;
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
igraph_bool_t igraph_i_dfscode_is_canonical(const igraph_dfscode_t *dfscode, igraph_t *graph,
		igraph_vector_int_t *vertex_colors, igraph_vector_int_t *edge_colors,
		igraph_vector_int_t *edge_timestamps) {
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
						   edge_timestamps,
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
		igraph_vector_int_t *edge_timestamps,
		long int dfscode_pos, igraph_vector_int_t *ordering, long int ordering_pos,
		igraph_vector_int_t *visited_nodes) {
  long int i, j, k, neigh, backward_neigh, backward_edge_count, ext_node, ext_node_pos;
  int cmp;
  igraph_dfscode_edge_t new_edge;
  igraph_integer_t eid;
  igraph_bool_t is_less_eq, is_outneigh, directed;

  directed = igraph_is_directed(graph);

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
    if (directed) {
      is_outneigh = (i < OUT_DEGREE(*graph, ext_node));
      new_edge.d = !is_outneigh; // 0 if e=(ext_node,neigh), 1 if e=(neigh,ext_node)
    }
    if (vertex_colors) {
      new_edge.l_i = VECTOR(*vertex_colors)[ext_node];
      new_edge.l_j = VECTOR(*vertex_colors)[neigh];
    }
    igraph_get_eid(graph, &eid, ext_node, neigh, 1, 0);
    if (edge_colors) {
      new_edge.l_ij = VECTOR(*edge_colors)[eid];
    }
    if (edge_timestamps) {
      new_edge.t_ij = VECTOR(*edge_timestamps)[eid];
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
	  if (directed) {
	    is_outneigh = (k < OUT_DEGREE(*graph, neigh));
	    new_edge.d = !is_outneigh; // 0 if e=(neigh,backward_neigh), 1 otherwise
	  }
	  if (vertex_colors) {
	    new_edge.l_i = VECTOR(*vertex_colors)[neigh];
	    new_edge.l_j = VECTOR(*vertex_colors)[backward_neigh];
	  }
	  igraph_get_eid(graph, &eid, neigh, backward_neigh, /*directed=*/ 1, /*error=*/ 0);
	  if (edge_colors) {
	    new_edge.l_ij = VECTOR(*edge_colors)[eid];
	  }
	  if (edge_timestamps) {
	    new_edge.t_ij = VECTOR(*edge_timestamps)[eid];
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
						   edge_timestamps,
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


// g is valid if the start node (0) and end node (1) are connected to all other nodes
igraph_bool_t igraph_i_lfrminer_valid(igraph_t *g) {
  igraph_vector_int_t hits;
  long int i;

  if (igraph_is_directed(g)) {
    igraph_vector_int_init(&hits, igraph_vcount(g));

    // check start node
    VECTOR(hits)[0] = 1;
    for (i = 0; i < DEGREE(*g, 0); i++) {
      VECTOR(hits)[NEIGHBOR(*g, 0, i)] = 1;
    }
    if (igraph_vector_int_min(&hits) == 0)
      return 0;

    // check end node
    igraph_vector_int_fill(&hits, 0);
    VECTOR(hits)[1] = 1;
    for (i = 0; i < DEGREE(*g, 1); i++) {
      VECTOR(hits)[NEIGHBOR(*g, 1, i)] = 1;
    }
    if (igraph_vector_int_min(&hits) == 0)
      return 0;
    else
      return 1;
  } else { // undirected
    return ((DEGREE(*g, 0) == igraph_vcount(g)-1)
	      && (DEGREE(*g, 1) == igraph_vcount(g)-1));
  }
}


// TODO: Implement speed-up for labelled graphs:
// Check that a new edge is larger than all existing edges incident to the target
// node (for backward extension) or the source node (for forward extension from any
// node other than right-most vertex). Otherwise, the resulting DFS code can not be
// minimal. See Yan & Han (2002), Section 5.1 (2), and Borgelt (2006), Section 4.
int igraph_i_dfscode_extend(const igraph_vector_ptr_t *graphs,
		const igraph_vector_ptr_t *vertex_colors,
		const igraph_vector_ptr_t *edge_colors,
		const igraph_vector_ptr_t *edge_times,
		igraph_support_measure_t *single_graph_support,
		igraph_integer_t min_supp, igraph_integer_t max_edges,
		igraph_gspan_variant_t variant,
		void *variant_data,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors,
		igraph_dfscode_t *seed_dfscode, igraph_llist_ptr_t *result_graph_list,
		igraph_llist_ptr_t *result_vcolor_list,
		igraph_llist_ptr_t *result_ecolor_list,
		igraph_llist_ptr_t *result_etimes_list,
		igraph_llist_long_t *result_supp_list) {
  long int i, j, d, t, cur_color, rightmost_vertex_color;
  long int cur_vertex, prev_vertex, rightmost_vertex;
  igraph_stack_int_t rightmost_path, rightmost_path_colors;
  igraph_dfscode_edge_t new_edge, rightmost_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1};
  igraph_t *seed_graph;
  igraph_vector_int_t *seed_vcolors;
  igraph_vector_int_t *seed_ecolors;
  igraph_vector_int_t *seed_etimes;
  long int seed_supp;
  igraph_bool_t directed;

  if ((graphs == NULL) || (igraph_vector_ptr_size(graphs) == 0)) {
    IGRAPH_ERROR("no graph database specified", IGRAPH_EINVAL);
  }
  directed = igraph_is_directed((igraph_t *) VECTOR(*graphs)[0]);

  // create graph from DFS code
  seed_graph = igraph_Calloc(1, igraph_t);
  seed_vcolors = igraph_Calloc(1, igraph_vector_int_t);
  seed_ecolors = igraph_Calloc(1, igraph_vector_int_t);
  seed_etimes = igraph_Calloc(1, igraph_vector_int_t);
  IGRAPH_CHECK(igraph_i_dfscode_to_graph(seed_dfscode, directed,
			    seed_graph, seed_vcolors, seed_ecolors, seed_etimes));

  // check if DFS code is canonical
  igraph_i_dfscode_print(seed_dfscode);
  if (!igraph_i_dfscode_is_canonical(seed_dfscode, seed_graph, seed_vcolors,
	    seed_ecolors, seed_etimes)) {
    printf("   not canonical\n");
    igraph_destroy(seed_graph);
    igraph_vector_int_destroy(seed_vcolors);
    igraph_vector_int_destroy(seed_ecolors);
    igraph_vector_int_destroy(seed_etimes);
    igraph_free(seed_graph);
    igraph_free(seed_vcolors);
    igraph_free(seed_ecolors);
    igraph_free(seed_etimes);
    igraph_fsm_stats_noncanonical_count++;
    return 0;
  }

  // compute seed support
  igraph_aggregated_db_support(graphs, vertex_colors, edge_colors, edge_times, seed_graph,
		  seed_vcolors, seed_ecolors, seed_etimes, /*induced=*/ 0, variant, variant_data,
		  single_graph_support,
		  &seed_supp, min_supp);
  printf("   supp: %ld\n", seed_supp);
  if (seed_supp < min_supp) {
    // infrequent seed, free memory and prune
    igraph_destroy(seed_graph);
    igraph_vector_int_destroy(seed_vcolors);
    igraph_vector_int_destroy(seed_ecolors);
    igraph_vector_int_destroy(seed_etimes);
    igraph_free(seed_graph);
    igraph_free(seed_vcolors);
    igraph_free(seed_ecolors);
    igraph_free(seed_etimes);
    igraph_fsm_stats_infrequent_count++;
    return 0;
  } else if ((variant == IGRAPH_GSPAN_LFRMINER) && !igraph_i_lfrminer_valid(seed_graph)) {
    // LFR-Miner mode: pattern is frequent but invalid, since not all nodes are connected
    // to start node 0 and end node 1. It can, however, be extended to a valid pattern
    // (otherwise it would not have been generated due to our extension modification).
    // Hence, free memory, don't add it to the result, but continue processing.
    igraph_destroy(seed_graph);
    igraph_vector_int_destroy(seed_vcolors);
    igraph_vector_int_destroy(seed_ecolors);
    igraph_vector_int_destroy(seed_etimes);
    igraph_free(seed_graph);
    igraph_free(seed_vcolors);
    igraph_free(seed_ecolors);
    igraph_free(seed_etimes);
    igraph_fsm_stats_lfrminer_invalid_count++;
    igraph_fsm_stats_frequent_count++;
  } else {
    // frequent seed, add to result
    igraph_llist_ptr_push_back(result_graph_list, seed_graph);
    igraph_llist_ptr_push_back(result_vcolor_list, seed_vcolors);
    igraph_llist_ptr_push_back(result_ecolor_list, seed_ecolors);
    igraph_llist_ptr_push_back(result_etimes_list, seed_etimes);
    igraph_llist_long_push_back(result_supp_list, seed_supp);
    igraph_fsm_stats_frequent_count++;
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
  prev_vertex = rightmost_vertex;
  for (i = igraph_i_dfscode_size(seed_dfscode)-1; i >= 0; i--) {
    if ((VECTOR(*seed_dfscode)[i].i < VECTOR(*seed_dfscode)[i].j)
	  && (VECTOR(*seed_dfscode)[i].j == prev_vertex)) {
      igraph_stack_int_push(&rightmost_path, VECTOR(*seed_dfscode)[i].i);
      igraph_stack_int_push(&rightmost_path_colors, VECTOR(*seed_dfscode)[i].l_i);

      if (prev_vertex == rightmost_vertex) {
	rightmost_edge = VECTOR(*seed_dfscode)[i];
      }
      prev_vertex = VECTOR(*seed_dfscode)[i].i;
    }
  }

  // iterate over all vertices from the right-most path and perform the extensions
  while (!igraph_stack_int_empty(&rightmost_path)) {
    cur_vertex = igraph_stack_int_pop(&rightmost_path);
    cur_color = igraph_stack_int_pop(&rightmost_path_colors);

    // FORWARD EXTENSIONS
    // from right-most path to new vertex

    new_edge = (igraph_dfscode_edge_t) {.i = cur_vertex, .j = rightmost_vertex+1, .d = 0,
					.l_i = cur_color, .l_ij = 0, .l_j = 0};
    switch(variant) {
      case IGRAPH_GSPAN_LFRMINER:
	// edge timestamps are not used in the patterns, only in the target graph!

	// enforce additional connectivity constraints to avoid invalid patterns
	// TODO: augment variant_data with STRICT flag to turn on/off connectivity constraints
	if ((cur_vertex == 0) && (rightmost_vertex >= 1)) {
	  break;
	}
	if (cur_vertex >= 1) {
	  for (i = 2; i <= rightmost_vertex; i++) {
	    if (!(igraph_i_dfscode_contains_edge(seed_dfscode, i, 0)
		    && igraph_i_dfscode_contains_edge(seed_dfscode, i, 1))) {
	      break;
	    }
	  }
	  if (i <= rightmost_vertex) {
	    break;
	  }
	}

	// extend with edge(s) to new standard node (label 2)
	new_edge.l_j = 2;
	for (d = 0; d <= (long int) directed; d++) {
	  // extend with all possible edge directions
	  new_edge.d = d;
	  if (edge_colors != NULL) {
	    // extend with all frequent edge colors
	    for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
	      new_edge.l_ij = VECTOR(*freq_ecolors)[i];
	      igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	      igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				single_graph_support,
				min_supp, max_edges, variant, variant_data,
				freq_vcolors, freq_ecolors,
				seed_dfscode, result_graph_list, result_vcolor_list,
				result_ecolor_list, result_etimes_list, result_supp_list);
	      igraph_i_dfscode_pop_back(seed_dfscode);
	    }
	  } else {
	    // no edge color
	    igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	    igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				    single_graph_support,
				    min_supp, max_edges, variant, variant_data,
				    freq_vcolors, freq_ecolors,
				    seed_dfscode, result_graph_list, result_vcolor_list,
				    result_ecolor_list, result_etimes_list, result_supp_list);
	    igraph_i_dfscode_pop_back(seed_dfscode);
	  }
	}
	break;

      case IGRAPH_GSPAN_GERM:
	for (d = 0; d <= (long int) directed; d++) {
	  // all possible edge directions
	  new_edge.d = d;
	  for (t = 0; t <= ((igraph_germ_data_t *)variant_data)->max_rel_timestamp; t++) {
	    // all possible relative edge timestamps
	    new_edge.t_ij = t;
	    if (vertex_colors != NULL) {
	      for (i = 0; VECTOR(*freq_vcolors)[i] != -1; i++) {
		// extend to all possible source vertex colors
		new_edge.l_j = VECTOR(*freq_vcolors)[i];
		if (edge_colors != NULL) {
		  // extend with all frequent edge colors
		  for (j = 0; VECTOR(*freq_ecolors)[j] != -1; j++) {
		    new_edge.l_ij = VECTOR(*freq_ecolors)[j];
		    igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		    igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
					    single_graph_support,
					    min_supp, max_edges, variant, variant_data,
					    freq_vcolors, freq_ecolors,
					    seed_dfscode, result_graph_list, result_vcolor_list,
					    result_ecolor_list, result_etimes_list, result_supp_list);
		    igraph_i_dfscode_pop_back(seed_dfscode);
		  }
		} else {
		  // no edge color
		  igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		  igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
					  single_graph_support,
					  min_supp, max_edges, variant, variant_data,
					  freq_vcolors, freq_ecolors,
					  seed_dfscode, result_graph_list, result_vcolor_list,
					  result_ecolor_list, result_etimes_list, result_supp_list);
		  igraph_i_dfscode_pop_back(seed_dfscode);
		}
	      }
	    } else { // no vertex colors
	      if (edge_colors != NULL) {
		// extend with all frequent edge colors
		for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
		  new_edge.l_ij = VECTOR(*freq_ecolors)[i];
		  igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		  igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				    single_graph_support,
				    min_supp, max_edges, variant, variant_data,
				    freq_vcolors, freq_ecolors,
				    seed_dfscode, result_graph_list, result_vcolor_list,
				    result_ecolor_list, result_etimes_list, result_supp_list);
		  igraph_i_dfscode_pop_back(seed_dfscode);
		}
	      } else {
		// no edge color
		igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
					single_graph_support,
					min_supp, max_edges, variant, variant_data,
					freq_vcolors, freq_ecolors,
					seed_dfscode, result_graph_list, result_vcolor_list,
					result_ecolor_list, result_etimes_list, result_supp_list);
		igraph_i_dfscode_pop_back(seed_dfscode);
	      }
	    } // if (vertex_colors)
	  } // edge timestamps
	} // edge directions
	break;

      case IGRAPH_GSPAN_EVOMINE:
      case IGRAPH_GSPAN_DEFAULT:
      default:
	for (d = 0; d <= (long int) directed; d++) {
	  new_edge.d = d;
	  if (vertex_colors != NULL) {
	    for (i = 0; VECTOR(*freq_vcolors)[i] != -1; i++) {
	      new_edge.l_j = VECTOR(*freq_vcolors)[i]; // extend to all possible vertex colors
	      if (edge_colors != NULL) {
		// extend with all frequent edge colors
		for (j = 0; VECTOR(*freq_ecolors)[j] != -1; j++) {
		  new_edge.l_ij = VECTOR(*freq_ecolors)[j];
		  igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		  igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
					  single_graph_support,
					  min_supp, max_edges, variant, variant_data,
					  freq_vcolors, freq_ecolors,
					  seed_dfscode, result_graph_list, result_vcolor_list,
					  result_ecolor_list, result_etimes_list, result_supp_list);
		  igraph_i_dfscode_pop_back(seed_dfscode);
		}
	      } else {
		// no edge color
		igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
					single_graph_support,
					min_supp, max_edges, variant, variant_data,
					freq_vcolors, freq_ecolors,
					seed_dfscode, result_graph_list, result_vcolor_list,
					result_ecolor_list, result_etimes_list, result_supp_list);
		igraph_i_dfscode_pop_back(seed_dfscode);
	      }
	    }
	  } else { // no vertex colors
	    if (edge_colors != NULL) {
	      // extend with all frequent edge colors
	      for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
		new_edge.l_ij = VECTOR(*freq_ecolors)[i];
		igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				  single_graph_support,
				  min_supp, max_edges, variant, variant_data,
				  freq_vcolors, freq_ecolors,
				  seed_dfscode, result_graph_list, result_vcolor_list,
				  result_ecolor_list, result_etimes_list, result_supp_list);
		igraph_i_dfscode_pop_back(seed_dfscode);
	      }
	    } else {
	      // no edge color
	      igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	      igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				      single_graph_support,
				      min_supp, max_edges, variant, variant_data,
				      freq_vcolors, freq_ecolors,
				      seed_dfscode, result_graph_list, result_vcolor_list,
				      result_ecolor_list, result_etimes_list, result_supp_list);
	      igraph_i_dfscode_pop_back(seed_dfscode);
	    }
	  } // if (vertex_colors)
	} // directions
	break;
    }

    // BACKWARD EXTENSION
    // from right-most vertex to current node on right-most path

    // prune invalid cases
    if (cur_vertex == rightmost_vertex) {
      // no self-loops
      continue;
    }
    if ((cur_vertex == rightmost_edge.i) && !directed) {
      // this edge already exists as a forward edge
      continue;
    }
    if ((VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i
	  > VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j)
	&& (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i == rightmost_vertex)) {
      // last edge was a backward edge starting from rightmost_vertex
      if (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j > cur_vertex) {
	// the last edge ended AFTER cur_vertex. a backward extension to
	// cur_vertex would result in a non-minimal DFS code.
	continue;
      }
      if (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j == cur_vertex && !directed) {
	// the last edge ended at cur_vertex. a backward extension to
	// cur_vertex would result in a duplicate edge.
	// TODO: can this happen?
	continue;
      }
    }

    new_edge = (igraph_dfscode_edge_t) {.i = rightmost_vertex, .j = cur_vertex, .d = 0,
					.l_i = rightmost_vertex_color, .l_ij = 0, .t_ij = 0,
					.l_j = cur_color};

    // extend in all possible directions
    for (d = 0; d <= (long int) directed; d++) {
      new_edge.d = d;

      // more pruning for invalid directed cases
      if (directed && (cur_vertex == rightmost_edge.i) && (d == 1-rightmost_edge.d)) {
	// this edge already exists as a forward edge with the same direction
	// NOTE: forward edge direction d corresponds to backward edge direction 1-d
	continue;
      }
      if (directed && (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i
	    > VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j)
	  && (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].i == rightmost_vertex)
	  && (VECTOR(*seed_dfscode)[seed_dfscode->last_edge].j == cur_vertex)
	  && (d == VECTOR(*seed_dfscode)[seed_dfscode->last_edge].d)) {
	// this edge already exists as a backward edge with the same direction
	continue;
      }

      // actual backward extension
      switch (variant) {
	case IGRAPH_GSPAN_GERM:
	  // extend with all possible relative edge timestamps
	  for (t = 0; t <= ((igraph_germ_data_t *)variant_data)->max_rel_timestamp; t++) {
	    new_edge.t_ij = t;
	    if (edge_colors != NULL) {
	      // extend with all possible edge colors
	      for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
		new_edge.l_ij = VECTOR(*freq_ecolors)[i];
		igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
		igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
			    single_graph_support,
			    min_supp, max_edges, variant, variant_data,
			    freq_vcolors, freq_ecolors,
			    seed_dfscode, result_graph_list, result_vcolor_list,
			    result_ecolor_list, result_etimes_list, result_supp_list);
		igraph_i_dfscode_pop_back(seed_dfscode);
	      }
	    } else {
	      // extend with unlabelled edge
	      igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	      igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
				single_graph_support,
				min_supp, max_edges, variant, variant_data,
				freq_vcolors, freq_ecolors,
				seed_dfscode, result_graph_list, result_vcolor_list,
				result_ecolor_list, result_etimes_list, result_supp_list);
	      igraph_i_dfscode_pop_back(seed_dfscode);
	    }
	  }
	  break;

	case IGRAPH_GSPAN_LFRMINER:
	  // enforce additional connectivity constraint to avoid invalid patterns
	  // TODO: augment variant_data with STRICT flag to turn on/off connectivity constraints
	  if (cur_vertex >= 2) {
	    if (!(igraph_i_dfscode_contains_edge(seed_dfscode, rightmost_vertex, 0)
	          && igraph_i_dfscode_contains_edge(seed_dfscode, rightmost_vertex, 1))) {
	      break;
	    }
	  }
	  // NO BREAK HERE! we perform the standard backward extension from below!
	case IGRAPH_GSPAN_EVOMINE:
	case IGRAPH_GSPAN_DEFAULT:
	default:
	  if (edge_colors != NULL) {
	    // extend with all possible edge colors
	    for (i = 0; VECTOR(*freq_ecolors)[i] != -1; i++) {
	      new_edge.l_ij = VECTOR(*freq_ecolors)[i];
	      igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	      igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
			  single_graph_support,
			  min_supp, max_edges, variant, variant_data,
			  freq_vcolors, freq_ecolors,
			  seed_dfscode, result_graph_list, result_vcolor_list,
			  result_ecolor_list, result_etimes_list, result_supp_list);
	      igraph_i_dfscode_pop_back(seed_dfscode);
	    }
	  } else {
	    // extend with unlabelled edge
	    igraph_i_dfscode_push_back(seed_dfscode, &new_edge);
	    igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
			      single_graph_support,
			      min_supp, max_edges, variant, variant_data,
			      freq_vcolors, freq_ecolors,
			      seed_dfscode, result_graph_list, result_vcolor_list,
			      result_ecolor_list, result_etimes_list, result_supp_list);
	    igraph_i_dfscode_pop_back(seed_dfscode);
	  }
	  break;
      } // switch
    } // directions
  } // loop over right-most path
  igraph_stack_int_destroy(&rightmost_path);
  igraph_stack_int_destroy(&rightmost_path_colors);

  return 0;
}


// LFR-Miner is currently only implemented for unlabelled nodes!
// Edge labels supported. Node labels mark s and e.
// All other nodes created during extension will have the color 2.
int igraph_i_build_seeds_lfrminer(igraph_bool_t has_ecolors,
				 const igraph_vector_int_t *freq_ecolors,
				 igraph_integer_t max_edges,
				 igraph_llist_ptr_t *initial_patterns) {
  igraph_dfscode_t *pattern_dfscode;
  igraph_dfscode_edge_t pattern_dfscode_edge;
  long int k;

  if (has_ecolors) {
    for (k = 0; VECTOR(*freq_ecolors)[k] != -1; k++) {
      // s --EC[k]--> e
      pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
      pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			.l_i = 0,
			.d = 0,
			.l_ij = VECTOR(*freq_ecolors)[k],
			.l_j = 1};
      igraph_i_dfscode_init(pattern_dfscode, max_edges);
      igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
      igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode);
    }
  } else {
    // s --> e
    pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
    pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
		      .l_i = 0,
		      .d = 0,
		      .l_ij = 0,
		      .l_j = 1};
    igraph_i_dfscode_init(pattern_dfscode, max_edges);
    igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge);
    igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode);
  }

  return 0;
}


// NOTE: When run in GERM mode, initial edges always have a timestamp of 0,
//       because the timestamp field in igraph_dfscode_edge_t is initialized with 0
int igraph_i_build_seeds_default(igraph_bool_t has_vcolors, igraph_bool_t has_ecolors,
				 const igraph_vector_int_t *freq_vcolors,
				 const igraph_vector_int_t *freq_ecolors,
				 igraph_integer_t max_edges,
				 igraph_gspan_variant_t variant,
				 void *variant_data,
				 igraph_llist_ptr_t *initial_patterns) {
  igraph_dfscode_t *pattern_dfscode;
  igraph_dfscode_edge_t pattern_dfscode_edge;
  long int i, j, k;

  if (has_vcolors) {
    for (i = 0; VECTOR(*freq_vcolors)[i] != -1; i++) {
      for (j = i; VECTOR(*freq_vcolors)[j] != -1; j++) {
	if (has_ecolors) {
	  for (k = 0; VECTOR(*freq_ecolors)[k] != -1; k++) {
	    // VC[i] -- EC[k] -- VC[j]
	    pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	    pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			.l_i = VECTOR(*freq_vcolors)[i],
			.d = 0,
			.l_ij = VECTOR(*freq_ecolors)[k],
			.l_j = VECTOR(*freq_vcolors)[j]};
	    IGRAPH_CHECK(igraph_i_dfscode_init(pattern_dfscode, max_edges));
	    IGRAPH_CHECK(igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge));
	    IGRAPH_CHECK(igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode));
	  }
	} else {
	  // VC[i] -- VC[j]
	  pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	  pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			    .l_i = VECTOR(*freq_vcolors)[i],
			    .d = 0,
			    .l_ij = 0,
			    .l_j = VECTOR(*freq_vcolors)[j]};
	  IGRAPH_CHECK(igraph_i_dfscode_init(pattern_dfscode, max_edges));
	  IGRAPH_CHECK(igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge));
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode));
	}
      }
    }
  } else { // no vcolors
    if (has_ecolors) {
      for (k = 0; VECTOR(*freq_ecolors)[k] != -1; k++) {
	// v -- EC[k] -- v
	pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
	pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			  .l_i = 0,
			  .d = 0,
			  .l_ij = VECTOR(*freq_ecolors)[k],
			  .l_j = 0};
	IGRAPH_CHECK(igraph_i_dfscode_init(pattern_dfscode, max_edges));
	IGRAPH_CHECK(igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge));
	IGRAPH_CHECK(igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode));
      }
    } else {
      // v -- v
      pattern_dfscode = igraph_Calloc(1, igraph_dfscode_t);
      pattern_dfscode_edge = (igraph_dfscode_edge_t) {.i = 0, .j = 1,
			.l_i = 0,
			.d = 0,
			.l_ij = 0,
			.l_j = 0};
      IGRAPH_CHECK(igraph_i_dfscode_init(pattern_dfscode, max_edges));
      IGRAPH_CHECK(igraph_i_dfscode_push_back(pattern_dfscode, &pattern_dfscode_edge));
      IGRAPH_CHECK(igraph_llist_ptr_push_back(initial_patterns, pattern_dfscode));
    }
  }

  return 0;
}


int igraph_i_vector_ptr_to_vector_int_minmax(const igraph_vector_ptr_t *vec_ptr,
      long int *min, long int *max) {
  long int i, n;
  int cur_max, cur_min;
  *min = INT_MAX;
  *max = -1;
  n = igraph_vector_ptr_size(vec_ptr);

  for (i = 0; i < n; i++) {
    IGRAPH_CHECK(igraph_vector_int_minmax((igraph_vector_int_t *) VECTOR(*vec_ptr)[i],
		    &cur_min, &cur_max));
    if (cur_max > *max) {
      *max = cur_max;
    }
    if (cur_min < *min) {
      *min = cur_min;
    }
  }
  return 0;
}


int igraph_i_frequent_colors(const igraph_vector_ptr_t *vertex_colors,
		const igraph_vector_ptr_t *edge_colors, long int min_supp,
		long int max_vcolor, long int max_ecolor,
		igraph_vector_int_t *freq_vcolors, igraph_vector_int_t *freq_ecolors) {
  long int i, j, graph_count;
  igraph_vector_long_t vcolor_freq, ecolor_freq;
  igraph_vector_int_t *colors;

  if ((vertex_colors == NULL) && (edge_colors == NULL)) {
    IGRAPH_CHECK(igraph_vector_int_resize(freq_vcolors, 1));
    IGRAPH_CHECK(igraph_vector_int_resize(freq_ecolors, 1));
    igraph_vector_int_fill(freq_vcolors, -1);
    igraph_vector_int_fill(freq_ecolors, -1);
    return 0;
  } else if (vertex_colors != NULL) {
    graph_count = igraph_vector_ptr_size(vertex_colors);
  } else {
    graph_count = igraph_vector_ptr_size(edge_colors);
  }

  if (vertex_colors != NULL) {
    IGRAPH_CHECK(igraph_vector_long_init(&vcolor_freq, max_vcolor+1));
    IGRAPH_CHECK(igraph_vector_int_resize(freq_vcolors, max_vcolor+2));
    // the last element is just for guaranteed loop termination below
  } else {
    IGRAPH_CHECK(igraph_vector_int_resize(freq_vcolors, 1));
  }
  if (edge_colors != NULL) {
    IGRAPH_CHECK(igraph_vector_long_init(&ecolor_freq, max_ecolor+1));
    IGRAPH_CHECK(igraph_vector_int_resize(freq_ecolors, max_ecolor+2));
  } else {
    IGRAPH_CHECK(igraph_vector_int_resize(freq_ecolors, 1));
  }

  // count all color occurrences
  for (i = 0; i < graph_count; i++) {
    if (vertex_colors != NULL) {
      colors = (igraph_vector_int_t *) VECTOR(*vertex_colors)[i];
      for (j = 0; j < igraph_vector_int_size(colors); j++) {
        VECTOR(vcolor_freq)[VECTOR(*colors)[j]] += 1;
      }
    }
    if (edge_colors != NULL) {
      colors = (igraph_vector_int_t *) VECTOR(*edge_colors)[i];
      for (j = 0; j < igraph_vector_int_size(colors); j++) {
        VECTOR(ecolor_freq)[VECTOR(*colors)[j]] += 1;
      }
    }
  }

  // keep only frequent colors
  igraph_vector_int_fill(freq_vcolors, -1);
  igraph_vector_int_fill(freq_ecolors, -1);

  if (vertex_colors != NULL) {
    printf("node color frequencies: ");
    igraph_vector_long_print(&vcolor_freq);
    j = 0;
    for (i = 0; i <= max_vcolor; i++) {
      if (VECTOR(vcolor_freq)[i] >= min_supp) {
	VECTOR(*freq_vcolors)[j] = i;
	j += 1;
      }
    }
  }
  if (edge_colors != NULL) {
    printf("edge color frequencies: ");
    igraph_vector_long_print(&ecolor_freq);
    j = 0;
    for (i = 0; i <= max_ecolor; i++) {
      if (VECTOR(ecolor_freq)[i] >= min_supp) {
	VECTOR(*freq_ecolors)[j] = i;
	j += 1;
      }
    }
  }

  if (vertex_colors != NULL)
    igraph_vector_long_destroy(&vcolor_freq);
  if (edge_colors != NULL)
    igraph_vector_long_destroy(&ecolor_freq);

  return 0;
}



// public interface

// assert: #graphs == #vertex_colors == #edge_colors
// assert: 0 <= node and edge colors
// assert: min_supp > 0
// assert: undirected graphs
int igraph_gspan(const igraph_vector_ptr_t *graphs, const igraph_vector_ptr_t *vertex_colors,
		 const igraph_vector_ptr_t *edge_colors, const igraph_vector_ptr_t *edge_times,
		 igraph_support_measure_t *single_graph_support,
		 igraph_integer_t min_supp, igraph_integer_t max_edges,
		 igraph_gspan_variant_t variant,
		 igraph_vector_ptr_t *frequent_subgraphs,
		 igraph_vector_ptr_t *frequent_subgraph_vcolors,
		 igraph_vector_ptr_t *frequent_subgraph_ecolors,
		 igraph_vector_ptr_t *frequent_subgraph_etimes,
		 igraph_vector_long_t *frequent_subgraph_supps) {
  long int max_vcolor, max_ecolor, min_vcolor, min_ecolor, min_etime, max_etime;
  void *variant_data;
  igraph_vector_int_t freq_ecolors, freq_vcolors;
  igraph_llist_ptr_t initial_patterns;
  igraph_llist_item_ptr_t *item_ptr;
  igraph_llist_ptr_t result_graph_list, result_vcolor_list, result_ecolor_list, result_etimes_list;
  igraph_llist_long_t result_supp_list;

  // FIND FREQUENT VERTEX AND EDGE COLORS

  // determine minimum/maximum colors and timestamps
  if (vertex_colors != NULL)
    IGRAPH_CHECK(igraph_i_vector_ptr_to_vector_int_minmax(vertex_colors, &min_vcolor, &max_vcolor));
  if (edge_colors != NULL)
    IGRAPH_CHECK(igraph_i_vector_ptr_to_vector_int_minmax(edge_colors, &min_ecolor, &max_ecolor));
  if (edge_times != NULL)
    IGRAPH_CHECK(igraph_i_vector_ptr_to_vector_int_minmax(edge_times, &min_etime, &max_etime));

  // determine frequent vertex and edge colors
  IGRAPH_CHECK(igraph_vector_int_init(&freq_vcolors, 0));
  IGRAPH_CHECK(igraph_vector_int_init(&freq_ecolors, 0));
  IGRAPH_CHECK(igraph_i_frequent_colors(vertex_colors, edge_colors, min_supp,
			max_vcolor, max_ecolor,
			&freq_vcolors, &freq_ecolors));
  printf("GSPAN frequent ecolors: ");
  igraph_vector_int_print(&freq_ecolors);
  printf("GSPAN frequent vcolors: ");
  igraph_vector_int_print(&freq_vcolors);

  if (variant == IGRAPH_GSPAN_GERM) {
    variant_data = igraph_Calloc(1, igraph_germ_data_t);
    ((igraph_germ_data_t *)variant_data)->max_rel_timestamp = max_etime-min_etime;
  } else if (variant == IGRAPH_GSPAN_EVOMINE) {
    variant_data = igraph_Calloc(1, igraph_evomine_data_t);
    ((igraph_evomine_data_t *)variant_data)->max_vcolor = max_vcolor;
    ((igraph_evomine_data_t *)variant_data)->max_ecolor = max_ecolor;
  } else {
    variant_data = NULL;
  }

  // clean termination for dry run
  if (max_edges == 0) {
    return 0;
  }

  // BUILD ALL 1-EDGE GRAPHS AS SEEDS

  IGRAPH_CHECK(igraph_llist_ptr_init(&initial_patterns));

  switch (variant) {
    case IGRAPH_GSPAN_LFRMINER:
      IGRAPH_CHECK(igraph_i_build_seeds_lfrminer((edge_colors != NULL), &freq_ecolors,
				    max_edges, &initial_patterns));
      break;
    case IGRAPH_GSPAN_EVOMINE:
    case IGRAPH_GSPAN_GERM:
    case IGRAPH_GSPAN_DEFAULT:
    default:
      IGRAPH_CHECK(igraph_i_build_seeds_default((vertex_colors != NULL), (edge_colors != NULL),
		      &freq_vcolors, &freq_ecolors, max_edges, variant, variant_data,
		      &initial_patterns));
      break;
  }

  // RECURSIVELY EXPAND ALL FREQUENT 1-EDGE GRAPHS BY PATTERN GROWTH

  IGRAPH_CHECK(igraph_llist_ptr_init(&result_graph_list));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_vcolor_list));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_ecolor_list));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_etimes_list));
  IGRAPH_CHECK(igraph_llist_long_init(&result_supp_list));

  for (item_ptr = initial_patterns.first; item_ptr != NULL; item_ptr = item_ptr->next) {
    IGRAPH_CHECK(igraph_i_dfscode_extend(graphs, vertex_colors, edge_colors, edge_times,
		    single_graph_support,
		    min_supp, max_edges, variant, variant_data, &freq_vcolors, &freq_ecolors,
		    (igraph_dfscode_t *) item_ptr->data, &result_graph_list,
		    &result_vcolor_list, &result_ecolor_list, &result_etimes_list,
		    &result_supp_list));
  }

  // PREPARE RESULT SET

  // store result in the provided containers (if provided...)
  // the user has to free the allocated memory
  if (frequent_subgraphs != NULL)
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_graph_list, frequent_subgraphs, 0));
  if (frequent_subgraph_vcolors != NULL)
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_vcolor_list, frequent_subgraph_vcolors, 0));
  if (frequent_subgraph_ecolors != NULL)
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_ecolor_list, frequent_subgraph_ecolors, 0));
  if (frequent_subgraph_etimes != NULL)
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_etimes_list, frequent_subgraph_etimes, 0));
  if (frequent_subgraph_supps != NULL)
    IGRAPH_CHECK(igraph_llist_long_to_vector(&result_supp_list, frequent_subgraph_supps, 0));

  // CLEAN UP

  if (variant_data != NULL)
    igraph_free(variant_data);

  igraph_vector_int_destroy(&freq_vcolors);
  igraph_vector_int_destroy(&freq_ecolors);

  if (frequent_subgraphs == NULL) {
    // the user doesn't want the result, so we have to free the allocated memory
    for (item_ptr = result_graph_list.first; item_ptr != NULL; item_ptr = item_ptr->next) {
      igraph_destroy((igraph_t *) item_ptr->data);
    }
  }
  if (frequent_subgraph_vcolors == NULL) {
    // same here
    for (item_ptr = result_vcolor_list.first; item_ptr != NULL; item_ptr = item_ptr->next) {
      igraph_vector_int_destroy((igraph_vector_int_t *) item_ptr->data);
    }
  }
  if (frequent_subgraph_ecolors == NULL) {
    // same here
    for (item_ptr = result_ecolor_list.first; item_ptr != NULL; item_ptr = item_ptr->next) {
      igraph_vector_int_destroy((igraph_vector_int_t *) item_ptr->data);
    }
  }
  if (frequent_subgraph_etimes == NULL) {
    // same here
    for (item_ptr = result_etimes_list.first; item_ptr != NULL; item_ptr = item_ptr->next) {
      igraph_vector_int_destroy((igraph_vector_int_t *) item_ptr->data);
    }
  }

  // free memory for temporary structures
  for (item_ptr = initial_patterns.first; item_ptr != NULL; item_ptr = item_ptr->next) {
    igraph_i_dfscode_destroy((igraph_dfscode_t *) item_ptr->data);
  }

  igraph_llist_ptr_destroy(&initial_patterns);
  igraph_llist_long_destroy(&result_supp_list);
  igraph_llist_ptr_destroy(&result_graph_list);
  igraph_llist_ptr_destroy(&result_vcolor_list);
  igraph_llist_ptr_destroy(&result_ecolor_list);
  igraph_llist_ptr_destroy(&result_etimes_list);

  // STATISTICS

  printf("\nGSPAN STATISTICS\n"
	  "total subisomorphism checks: %ld\n"
	  "   successful: %ld\n"
	  "   failed: %ld\n"
	  "   failed edge existence: %ld\n"
	  "   failed edge color: %ld\n"
	  "   failed edge timestamp: %ld\n"
	  "   failed node degree: %ld\n"
	  "   failed node duplicate: %ld\n"
	  "   failed node color: %ld\n"
	  "DB support computations: %ld\n"
	  "   MIB support computations: %ld\n"
	  "   MIB subisomorphism checks: %ld\n"
	  "      MIB successful subisomorphism checks: %ld\n"
	  "      MIB failed subisomorphism checks: %ld\n"
	  "   Shallow support computations: %ld\n"
	  "   Ego-based support computations: %ld\n"
	  "Infrequent pattern candidates: %ld\n"
	  "Frequent patterns: %ld\n"
	  "   Invalid intermediate LF patterns: %ld\n"
	  "Non-canonical pattern candidates: %ld\n",
	  igraph_fsm_stats_subiso_success_count+igraph_fsm_stats_subiso_fail_count,
	  igraph_fsm_stats_subiso_success_count,
	  igraph_fsm_stats_subiso_fail_count,
	  igraph_fsm_stats_subiso_failed_edge_existence_count,
	  igraph_fsm_stats_subiso_failed_edge_color_count,
	  igraph_fsm_stats_subiso_failed_edge_timestamp_count,
	  igraph_fsm_stats_subiso_failed_node_degree_count,
	  igraph_fsm_stats_subiso_failed_node_duplicate_count,
	  igraph_fsm_stats_subiso_failed_node_color_count,
	  igraph_fsm_stats_aggregated_support_count,
	  igraph_fsm_stats_mibsupport_count,
	  (igraph_fsm_stats_mibsupport_subiso_success_count
		+ igraph_fsm_stats_mibsupport_subiso_fail_count),
	  igraph_fsm_stats_mibsupport_subiso_success_count,
	  igraph_fsm_stats_mibsupport_subiso_fail_count,
	  igraph_fsm_stats_shallowsuppport_count,
	  igraph_fsm_stats_egobasedsuppport_count,
	  igraph_fsm_stats_infrequent_count,
	  igraph_fsm_stats_frequent_count,
	  igraph_fsm_stats_lfrminer_invalid_count,
	  igraph_fsm_stats_noncanonical_count);

  return 0;
}
