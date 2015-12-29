/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   Algorithms for dynamic graphs
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
#include "igraph_dynamic.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_list.h"
#include "igraph_memory.h"
#include "igraph_structural.h"


int igraph_i_compute_joint_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_t *changed_nodes, igraph_vector_t *neighborhood);

int igraph_i_compute_union_graph_projection(igraph_t *graph1,
	  igraph_vector_int_t *vcolors1, igraph_vector_int_t *ecolors1,
	  igraph_t *graph2,
	  igraph_vector_int_t *vcolors2, igraph_vector_int_t *ecolors2,
	  igraph_vs_t node_selector,
	  igraph_t *union_graph,
	  igraph_vector_int_t *union_graph_vcolors,
	  igraph_vector_int_t *union_graph_ecolors,
	  igraph_integer_t max_vcolor, igraph_integer_t max_ecolor);


// input file should be in the format:
//    v $vid [...]
//    e $vid1 $vid2 $creation [$deletion [...]]
// where -1 for $deletion means that the edge is never deleted and [...] is ignored. Edges
// must be sorted by $creation. Node ids start with 0. Example:
//    v 0
//    v 1
//    ...
//    e 0 3 0 -1
//    e 1 41 0 3
//    e 1 8 1 -1
//    ...
int igraph_read_dynamic_velist(igraph_vector_ptr_t *graphs, FILE *instream) {
  char buf[32];
  long int v1, v2, max_vid, timestamp, last_timestamp;
  char ve;
  igraph_t *graph, *graph_copy;
  igraph_llist_int_t add_edges;
  igraph_llist_ptr_t graph_list;
  igraph_vector_t edges;

  IGRAPH_CHECK(igraph_llist_int_init(&add_edges));
  IGRAPH_CHECK(igraph_llist_ptr_init(&graph_list));
  IGRAPH_CHECK(igraph_vector_init(&edges, 0));

  // read vertices
  max_vid = -1;
  last_timestamp = -1;
  while (fgets(buf, 32, instream)) {
    if (sscanf(buf, "%c %ld %ld %ld", &ve, &v1, &v2, &timestamp) < 1) {
      // ignore lines that cannot be parsed
      continue;
    }
    if (ve == 'v') {
      if (v1 > max_vid) {
	max_vid = v1;
      }
    } else if (ve == 'e') {
      // we reached the first edge
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, v1));
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, v2));
      last_timestamp = timestamp;
      break;
    } else {
      // ignore lines not starting with v or e
    }
  }

  // read more edges
  while (fgets(buf, 32, instream)) {
    if (sscanf(buf, "%c %ld %ld %ld", &ve, &v1, &v2, &timestamp) < 1) {
      // ignore lines that cannot be parsed
      continue;

    }
    if (ve != 'e') {
      // ignore lines that do not define an edge
      continue;
    }

    if (timestamp > last_timestamp) {
      // we have processed all edges from this timestamp, construct graph

      // initialize new graph from previous one (or empty one, if this is the first timestamp)
      graph = igraph_Calloc(1, igraph_t);
      if (igraph_llist_ptr_size(&graph_list) == 0) {
	IGRAPH_CHECK(igraph_empty(graph, max_vid+1, /*directed=*/ 0));
      } else {
	IGRAPH_CHECK(igraph_copy(graph, (igraph_t *) igraph_llist_ptr_back(&graph_list)));
      }

      // if current and the last timestamp are not consecutive, add copies
      // of the current graph to the list to fill the gap
      while (last_timestamp < timestamp-1) {
	graph_copy = igraph_Calloc(1, igraph_t);
	IGRAPH_CHECK(igraph_copy(graph_copy, graph));
	IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph_copy));
	last_timestamp++;
      }

      // add edges to the current graph
      IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges));
      IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

      // append to graph list
      IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph));

      // destroy old edge list and create new empty ones
      igraph_llist_int_destroy(&add_edges);
      IGRAPH_CHECK(igraph_llist_int_init(&add_edges));

      last_timestamp = timestamp;
    } else if (timestamp < last_timestamp) {
      IGRAPH_ERROR("edges not sorted by timestamp", IGRAPH_PARSEERROR);
      return 1;
    }

    IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, v1));
    IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, v2));
  }

  // process last timestamp
  graph = igraph_Calloc(1, igraph_t);
  if (igraph_llist_ptr_size(&graph_list) == 0) {
    IGRAPH_CHECK(igraph_empty(graph, max_vid+1, /*directed=*/ 0));
  } else {
    IGRAPH_CHECK(igraph_copy(graph, (igraph_t *) igraph_llist_ptr_back(&graph_list)));
  }
  while (last_timestamp < timestamp-1) {
    graph_copy = igraph_Calloc(1, igraph_t);
    IGRAPH_CHECK(igraph_copy(graph_copy, graph));
    IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph_copy));
    last_timestamp++;
  }
  IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));
  IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph));

  // convert list into vector for output
  IGRAPH_CHECK(igraph_llist_ptr_to_vector(&graph_list, graphs));

  // TODO: process deleted edges

  igraph_llist_int_destroy(&add_edges);
  igraph_llist_ptr_destroy(&graph_list);
  igraph_vector_destroy(&edges);

  return 0;
}


// assert: all nodes from seed_nodes appear in graph1 and graph2
// assert: undirected graphs
//
// runtime: O(N_seeds*k_max + N_seeds*k_max*log(N_seeds*k_max) + (N_seeds*k_max+N_seeds))
int igraph_i_compute_joint_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_t *seed_nodes, igraph_vector_t *neighborhood) {
  igraph_llist_t node_list;
  long int i, j, k, changed_count, neigh, neigh_size;
  igraph_vector_t neighborhood_dups;

  changed_count = igraph_vector_size(seed_nodes);

  igraph_llist_init(&node_list);
  igraph_vector_init(&neighborhood_dups, 0);

  for (i = 0; i < changed_count; i++) {
    // collect neighbors in graph1
    for (j = 0; j < DEGREE(*graph1, i); j++) {
      neigh = NEIGHBOR(*graph1, i, j);
      if (!igraph_vector_binsearch2(seed_nodes, neigh)) { // filter to avoid duplicates
	igraph_llist_push_back(&node_list, neigh);
      }
    }
    // collect neighbors in graph2
    for (j = 0; j < DEGREE(*graph2, i); j++) {
      neigh = NEIGHBOR(*graph2, i, j);
      if (!igraph_vector_binsearch2(seed_nodes, neigh)) {
	igraph_llist_push_back(&node_list, neigh);
      }
    }
  }
  igraph_llist_to_vector(&node_list, &neighborhood_dups); // contains duplicates!
  neigh_size = igraph_vector_size(&neighborhood_dups);

  // sort neighborhood and merge it with the seed nodes, while removing duplicates

  igraph_vector_resize(neighborhood, changed_count+neigh_size);
  igraph_vector_sort(&neighborhood_dups);
  i = 0; // index in seed_nodes
  j = 0; // index in neighborhood_dups
  k = 0; // index in neighborhood (merge result)

  while (i < changed_count || j < neigh_size) {
    // if we reached the end of one of the two vectors, fill it with the other
    if (i == changed_count) {
      VECTOR(*neighborhood)[k] = VECTOR(neighborhood_dups)[j];
      while ((VECTOR(neighborhood_dups)[j] == VECTOR(*neighborhood)[k]) && (j < neigh_size))
	j++;
      k++;
      continue;
    }
    if (j == neigh_size) {
      VECTOR(*neighborhood)[k] = VECTOR(*seed_nodes)[i];
      k++; i++;
      continue;
    }

    // both vectors are not exhausted yet
    if (VECTOR(*seed_nodes)[i] < VECTOR(neighborhood_dups)[j]) {
      VECTOR(*neighborhood)[k] = VECTOR(*seed_nodes)[i];
      k++; i++;
    } else if (VECTOR(*seed_nodes)[i] == VECTOR(neighborhood_dups)[j]) {
      VECTOR(*neighborhood)[k] = VECTOR(*seed_nodes)[i];
      while ((VECTOR(neighborhood_dups)[j] == VECTOR(*neighborhood)[k]) && (j < neigh_size))
	j++;
      k++; i++;
    } else if (VECTOR(*seed_nodes)[i] > VECTOR(neighborhood_dups)[j]) {
      VECTOR(*neighborhood)[k] = VECTOR(neighborhood_dups)[j];
      while ((VECTOR(neighborhood_dups)[j] == VECTOR(*neighborhood)[k]) && (j < neigh_size))
	j++;
      k++;
    }
  }
  igraph_vector_resize(neighborhood, k);
  igraph_llist_destroy(&node_list);

  return 0;
}


// NOTE: we assume that the (undirected) edges in graph1 and graph2 are stored in exactly
// the same way, i.e. if the edge (i,j) is stored in the out-neighbors list of i in graph1,
// it should also be stored in the out-neighbors list of i in graph2 (and not in the
// in-neighbors of j). If that does not hold, we have two edges connecting i and j in
// the union graph.
//
// assert: ecolors and vcolors must be from the range [1, max_ecolor] or [1, max_vcolor], resp.!
int igraph_i_compute_union_graph_projection(igraph_t *graph1,
	  igraph_vector_int_t *vcolors1, igraph_vector_int_t *ecolors1,
	  igraph_t *graph2,
	  igraph_vector_int_t *vcolors2, igraph_vector_int_t *ecolors2,
	  igraph_vs_t node_selector,
	  igraph_t *union_graph,
	  igraph_vector_int_t *union_graph_vcolors,
	  igraph_vector_int_t *union_graph_ecolors,
	  igraph_integer_t max_vcolor,
	  igraph_integer_t max_ecolor) {
  long int i, j, k, eid1, eid2, orig_node;
  igraph_integer_t vcount;
  igraph_llist_int_t edge_list, ecolors_list;
  igraph_vector_int_t bw_index;
  igraph_vector_t edges;
  igraph_vit_t vit;

  igraph_bool_t has_vcolors = ((vcolors1 != NULL) && (vcolors2 != NULL)
				  && (union_graph_vcolors != NULL));
  igraph_bool_t has_ecolors = ((ecolors1 != NULL) && (ecolors2 != NULL)
				  && (union_graph_ecolors != NULL));

  // TODO: make sure that we only get vs_1, vs_vector, or vs_all
  IGRAPH_CHECK(igraph_vs_size(graph1, &node_selector, &vcount)); // same result for graph2

  IGRAPH_CHECK(igraph_empty(union_graph, vcount, /*directed=*/ 0));
  if (has_vcolors) {
    IGRAPH_CHECK(igraph_vector_int_resize(union_graph_vcolors, vcount));
  }

  // create backward index for vertex id lookup
  // TODO: only need this for vs_vector
  igraph_vit_create(graph1, node_selector, &vit);
  igraph_vector_int_init(&bw_index, igraph_vcount(graph1));
  igraph_vector_int_fill(&bw_index, -1);
  for (i=0; (i<vcount) && (!IGRAPH_VIT_END(vit)); i++, IGRAPH_VIT_NEXT(vit)) {
    VECTOR(bw_index)[(long int) IGRAPH_VIT_GET(vit)] = i;
  }

  IGRAPH_CHECK(igraph_llist_int_init(&edge_list));
  IGRAPH_CHECK(igraph_llist_int_init(&ecolors_list));
  for (i=0, IGRAPH_VIT_RESET(vit); (i<vcount) && (!IGRAPH_VIT_END(vit)); i++, IGRAPH_VIT_NEXT(vit)) {
    orig_node = (long int) IGRAPH_VIT_GET(vit);

    // set new vertex color to the base-(max_vcolor+1) number (vcolor1,vcolor2)
    if (has_vcolors) {
      VECTOR(*union_graph_vcolors)[i] = ((max_vcolor+1)*VECTOR(*vcolors1)[orig_node]
		+ VECTOR(*vcolors2)[orig_node]);
    }

    // fill edge list with edges from both input graphs
    j = 0; // neighbor index in graph1
    k = 0; // neighbor index in graph2
    while ((j < OUT_DEGREE(*graph1, orig_node)) || (k < OUT_DEGREE(*graph2, orig_node))) {
      // if we reached the end of one of the two vectors, fill it with the other
      if (j == OUT_DEGREE(*graph1, orig_node)) {
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)] >= 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)]);
	  if (has_ecolors) {
	    eid2 = OUT_NEIGH_TO_EID(*graph2, orig_node, k);
	    igraph_llist_int_push_back(&ecolors_list, VECTOR(*ecolors2)[eid2]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b01);
	  }
	}
	k++;
	continue;
      }
      if (k == OUT_DEGREE(*graph2, orig_node)) {
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] >= 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]);
	  if (has_ecolors) {
	    eid1 = OUT_NEIGH_TO_EID(*graph1, orig_node, j);
	    igraph_llist_int_push_back(&ecolors_list, (max_ecolor+1)*VECTOR(*ecolors1)[eid1]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b10);
	  }
	}
	j++;
	continue;
      }

      // both vectors are not exhausted yet
      if (OUT_NEIGHBOR(*graph1, orig_node, j) < OUT_NEIGHBOR(*graph2, orig_node, k)) {
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] >= 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]);
	  if (has_ecolors) {
	    eid1 = OUT_NEIGH_TO_EID(*graph1, orig_node, j);
	    igraph_llist_int_push_back(&ecolors_list, (max_ecolor+1)*VECTOR(*ecolors1)[eid1]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b10);
	  }
	}
	j++;
      } else if (OUT_NEIGHBOR(*graph1, orig_node, j) == OUT_NEIGHBOR(*graph2, orig_node, k)) {
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] >= 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]);
	  if (has_ecolors) {
	    eid1 = OUT_NEIGH_TO_EID(*graph1, orig_node, j);
	    eid2 = OUT_NEIGH_TO_EID(*graph2, orig_node, k);
	    igraph_llist_int_push_back(&ecolors_list, (max_ecolor+1)*VECTOR(*ecolors1)[eid1]
			    + VECTOR(*ecolors2)[eid2]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b11);
	  }
	}
	j++; k++;
      } else { // (OUT_NEIGHBOR(*graph1, orig_node, j) > OUT_NEIGHBOR(*graph2, orig_node, k))
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)] >= 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)]);
	  if (has_ecolors) {
	    eid2 = OUT_NEIGH_TO_EID(*graph2, orig_node, k);
	    igraph_llist_int_push_back(&ecolors_list, VECTOR(*ecolors2)[eid2]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b01);
	  }
	}
	k++;
      }
    }
  }

  // add edges to graph and set edge colors
  igraph_vector_init(&edges, 0);
  igraph_llist_int_to_vector_real(&edge_list, &edges);
  IGRAPH_CHECK(igraph_add_edges(union_graph, &edges, 0));
  igraph_llist_int_to_vector(&ecolors_list, union_graph_ecolors);

  //printf("edges\n");
  //igraph_vector_print(&edges);
  //printf("ecolors\n");
  //igraph_vector_int_print(union_graph_ecolors);
  //printf("ug nodes %d edges %d\n", igraph_vcount(union_graph), igraph_ecount(union_graph));

  // clean up
  igraph_llist_int_destroy(&edge_list);
  igraph_llist_int_destroy(&ecolors_list);
  igraph_vector_int_destroy(&bw_index);
  igraph_vector_destroy(&edges);

  return 0;
}


// assert: all graphs have the same number of nodes, and node IDs correspond with each other
// assert: undirected graph
int igraph_compute_dynamic_union_graph_projection(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_projection_t proj_type,
      igraph_vector_ptr_t *result, igraph_vector_ptr_t *result_vcolors,
      igraph_vector_ptr_t *result_ecolors,
      igraph_integer_t max_vcolor, igraph_integer_t max_ecolor) {
  long int T = igraph_vector_ptr_size(graphs);
  long int N = igraph_vcount((igraph_t *) VECTOR(*graphs)[0]);
  long int t, i, j, eid1, eid2;
  igraph_llist_int_t node_list;
  igraph_llist_ptr_t result_list, result_vcolors_list, result_ecolors_list;
  igraph_vector_t changed_nodes, neighborhood;
  igraph_t *union_graph;
  igraph_vector_int_t *union_graph_vcolors, *union_graph_ecolors;
  igraph_vs_t node_selector;

  if (proj_type == IGRAPH_PROJECTION_NEIGHBORS) {
    igraph_vector_init(&changed_nodes, 0);
    igraph_vector_init(&neighborhood, 0);
  }
  igraph_llist_ptr_init(&result_list);
  igraph_llist_ptr_init(&result_vcolors_list);
  igraph_llist_ptr_init(&result_ecolors_list);

  // iterate over all timestamps
  for (t = 0; t < T-1; t++) {
    switch (proj_type) {
    case IGRAPH_PROJECTION_FULL:
      IGRAPH_CHECK(igraph_vs_all(&node_selector));
      break;
    case IGRAPH_PROJECTION_EVENT:
      /* not implemented */
      return 1;
    case IGRAPH_PROJECTION_NEIGHBORS:
      // determine the nodes that have changed between the current and the next timestep
      igraph_llist_int_init(&node_list);
      for (i = 0; i < N; i++) {
	// check if node color has changed
	if ((vcolors != NULL) && (VECTOR(*(igraph_vector_int_t *) VECTOR(*vcolors)[t])[i]
	      != VECTOR(*(igraph_vector_int_t *) VECTOR(*vcolors)[t+1])[i])) {
	  igraph_llist_int_push_back(&node_list, i);
	  continue;
	}

	// check if node degree has changed
	if (DEGREE(*(igraph_t *) VECTOR(*graphs)[t], i)
	      != DEGREE(*(igraph_t *) VECTOR(*graphs)[t+1], i)) {
	  igraph_llist_int_push_back(&node_list, i);
	  continue;
	}

	// check if incident edges have changed
	for (j = 0; j < DEGREE(*(igraph_t *) VECTOR(*graphs)[t], i); j++) {
	  if (NEIGHBOR(*(igraph_t *) VECTOR(*graphs)[t], i, j)
		!= NEIGHBOR(*(igraph_t *) VECTOR(*graphs)[t+1], i, j)) {
	    igraph_llist_int_push_back(&node_list, i);
	    break;
	  }

	  // check edge label
	  if (ecolors != NULL) {
	    eid1 = NEIGH_TO_EID(*(igraph_t *) VECTOR(*graphs)[t], i, j);
	    eid2 = NEIGH_TO_EID(*(igraph_t *) VECTOR(*graphs)[t+1], i, j);
	    if (VECTOR(*(igraph_vector_int_t *) VECTOR(*ecolors)[t])[eid1]
		  != VECTOR(*(igraph_vector_int_t *) VECTOR(*ecolors)[t+1])[eid2]) {
	      igraph_llist_int_push_back(&node_list, i);
	      break;
	    }
	  }
	}
      }
      igraph_llist_int_to_vector_real(&node_list, &changed_nodes); // sorted by construction
      igraph_llist_int_destroy(&node_list);

      // compute joint 1-hop neighborhood at timesteps t and t+1
      igraph_i_compute_joint_neighborhood((igraph_t *) VECTOR(*graphs)[t],
	    (igraph_t *) VECTOR(*graphs)[t+1], &changed_nodes, &neighborhood);

      IGRAPH_CHECK(igraph_vs_vector(&node_selector, &neighborhood));
      break;
    } // switch (proj_type)

    // compute the union graph projection
    union_graph = igraph_Calloc(1, igraph_t);
    union_graph_vcolors = igraph_Calloc(1, igraph_vector_int_t);
    union_graph_ecolors = igraph_Calloc(1, igraph_vector_int_t);
    igraph_vector_int_init(union_graph_vcolors, 0);
    igraph_vector_int_init(union_graph_ecolors, 0);
    if (vcolors != NULL) {
      if (ecolors != NULL) {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	     node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	     max_vcolor, max_ecolor);
      } else {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t], NULL,
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1], NULL,
	     node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	     max_vcolor, max_ecolor);
      }
    } else {
      if (ecolors != NULL) {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     NULL, (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	     node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	     max_vcolor, max_ecolor);
      } else {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t], NULL, NULL,
	     (igraph_t *) VECTOR(*graphs)[t+1], NULL, NULL,
	     node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	     max_vcolor, max_ecolor);
      }
    }

    //printf("ug %ld: vcount %d ecount %d\n", t, igraph_vcount(union_graph), igraph_ecount(union_graph));
    igraph_llist_ptr_push_back(&result_list, union_graph);
    igraph_llist_ptr_push_back(&result_vcolors_list, union_graph_vcolors);
    igraph_llist_ptr_push_back(&result_ecolors_list, union_graph_ecolors);
  } // for t = 1...T

  // return result
  // TODO: free memory if ptrs == NULL
  if (result != NULL) {
    igraph_llist_ptr_to_vector(&result_list, result);
  }
  if (result_vcolors != NULL) {
    igraph_llist_ptr_to_vector(&result_vcolors_list, result_vcolors);
  }
  if (result_ecolors != NULL) {
    igraph_llist_ptr_to_vector(&result_ecolors_list, result_ecolors);
  }

  // clean up
  if (proj_type == IGRAPH_PROJECTION_NEIGHBORS) {
    igraph_vector_destroy(&changed_nodes);
    igraph_vector_destroy(&neighborhood);
  }
  igraph_llist_ptr_destroy(&result_list);
  igraph_llist_ptr_destroy(&result_vcolors_list);
  igraph_llist_ptr_destroy(&result_ecolors_list);

  return 0;
}

