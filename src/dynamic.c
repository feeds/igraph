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


int igraph_i_compute_dynamic_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_int_t *changed_nodes, igraph_vector_int_t *neighborhood);

int igraph_i_compute_union_graph_projection(igraph_t *graph1,
	  igraph_vector_int_t *vcolors1, igraph_vector_int_t *ecolors1,
	  igraph_t *graph2,
	  igraph_vector_int_t *vcolors2, igraph_vector_int_t *ecolors2,
	  igraph_vector_int_t *neighborhood,
	  igraph_t *union_graph,
	  igraph_vector_int_t *union_graph_vcolors,
	  igraph_vector_int_t *union_graph_ecolors);


int igraph_read_dynamic_velist(igraph_vector_ptr_t *graphs, FILE *instream) {
  char buf[32];
  long int v1, v2, max_vid, timestamp, last_timestamp;
  char ve;
  igraph_t *graph;
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

    if (timestamp != last_timestamp) {
      // we have processed all edges from this timestamp, construct graph

      // initialize new graph from previous one (or empty one, if this is the first timestamp)
      graph = igraph_Calloc(1, igraph_t);
      if (igraph_llist_ptr_size(&graph_list) == 0) {
	IGRAPH_CHECK(igraph_empty(graph, max_vid+1, /*directed=*/ 0));
      } else {
	IGRAPH_CHECK(igraph_copy(graph, (igraph_t *) igraph_llist_ptr_back(&graph_list)));
      }

      // add edges
      IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges));
      IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

      // append to graph list
      IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph));

      // destroy old edge list and create new empty ones
      igraph_llist_int_destroy(&add_edges);
      IGRAPH_CHECK(igraph_llist_int_init(&add_edges));

      last_timestamp = timestamp;
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
  IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));
  IGRAPH_CHECK(igraph_llist_ptr_push_back(&graph_list, graph));

  // convert list into vector for output
  IGRAPH_CHECK(igraph_llist_ptr_to_vector(&graph_list, graphs));

  igraph_llist_int_destroy(&add_edges);
  igraph_llist_ptr_destroy(&graph_list);
  igraph_vector_destroy(&edges);

  return 0;
}


// assert: all nodes from seed_nodes appear in graph1 and graph2
// assert: undirected graphs
//
// runtime: O(N_seeds*k_max + N_seeds*k_max*log(N_seeds*k_max) + (N_seeds*k_max+N_seeds))
int igraph_i_compute_dynamic_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_int_t *seed_nodes, igraph_vector_int_t *neighborhood) {
  igraph_llist_int_t node_list;
  long int i, j, k, changed_count, neigh, neigh_size;
  igraph_vector_int_t neighborhood_dups;

  changed_count = igraph_vector_int_size(seed_nodes);

  igraph_llist_int_init(&node_list);
  igraph_vector_int_init(&neighborhood_dups, 0);

  for (i = 0; i < changed_count; i++) {
    // collect neighbors in graph1
    for (j = 0; j < DEGREE(*graph1, i); j++) {
      neigh = NEIGHBOR(*graph1, i, j);
      if (!igraph_vector_int_binsearch2(seed_nodes, neigh)) { // filter to avoid duplicates
	igraph_llist_int_push_back(&node_list, neigh);
      }
    }
    // collect neighbors in graph2
    for (j = 0; j < DEGREE(*graph2, i); j++) {
      neigh = NEIGHBOR(*graph2, i, j);
      if (!igraph_vector_int_binsearch2(seed_nodes, neigh)) {
	igraph_llist_int_push_back(&node_list, neigh);
      }
    }
  }
  igraph_llist_int_to_vector(&node_list, &neighborhood_dups); // contains duplicates!
  neigh_size = igraph_vector_int_size(&neighborhood_dups);

  // sort neighborhood and merge it with the seed nodes, while removing duplicates

  igraph_vector_int_resize(neighborhood, changed_count+neigh_size);
  igraph_vector_int_sort(&neighborhood_dups);
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
  igraph_vector_int_resize(neighborhood, k);
  igraph_llist_int_destroy(&node_list);

  return 0;
}


int igraph_i_compute_union_graph_projection(igraph_t *graph1,
	  igraph_vector_int_t *vcolors1, igraph_vector_int_t *ecolors1,
	  igraph_t *graph2,
	  igraph_vector_int_t *vcolors2, igraph_vector_int_t *ecolors2,
	  igraph_vector_int_t *neighborhood,
	  igraph_t *union_graph,
	  igraph_vector_int_t *union_graph_vcolors,
	  igraph_vector_int_t *union_graph_ecolors) {
  // TODO
  return 0;
}


// assert: all graphs have the same number of nodes, and node IDs correspond with each other
// assert: undirected graph
int igraph_compute_dynamic_neighborhoods(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_vector_ptr_t *result, igraph_vector_ptr_t *result_vcolors,
      igraph_vector_ptr_t *result_ecolors) {
  long int T = igraph_vector_ptr_size(graphs);
  long int N = igraph_vcount((igraph_t *) VECTOR(*graphs)[0]);
  long int t, i, j, changed_count;
  igraph_llist_int_t node_list;
  igraph_llist_ptr_t result_list, result_vcolors_list, result_ecolors_list;
  igraph_vector_int_t changed_nodes, neighborhood;
  igraph_t *union_graph;
  igraph_vector_int_t *union_graph_vcolors, *union_graph_ecolors;

  igraph_vector_int_init(&changed_nodes, 0);
  igraph_vector_int_init(&neighborhood, 0);
  igraph_llist_ptr_init(&result_list);
  igraph_llist_ptr_init(&result_vcolors_list);
  igraph_llist_ptr_init(&result_ecolors_list);

  // iterate over all timestamps
  for (t = 0; t < T-1; t++) {
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
	  continue;
	}

	// TODO: check edge labels, geteid needed?
      }
    }
    changed_count = igraph_llist_int_size(&node_list);
    igraph_llist_int_to_vector(&node_list, &changed_nodes); // sorted by construction
    igraph_llist_int_destroy(&node_list);

    // compute 1-hop neighborhood at timesteps t and t+1
    igraph_i_compute_dynamic_neighborhood((igraph_t *) VECTOR(*graphs)[t],
	  (igraph_t *) VECTOR(*graphs)[t+1], &changed_nodes, &neighborhood);

    // compute the union graph projection
    union_graph = igraph_Calloc(1, igraph_t);
    union_graph_vcolors = igraph_Calloc(1, igraph_vector_int_t);
    union_graph_ecolors = igraph_Calloc(1, igraph_vector_int_t);
    if (vcolors != NULL) {
      if (ecolors != NULL) {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	     &neighborhood, union_graph, union_graph_vcolors, union_graph_ecolors);
      } else {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t], NULL,
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1], NULL,
	     &neighborhood, union_graph, union_graph_vcolors, union_graph_ecolors);
      }
    } else {
      if (ecolors != NULL) {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	     NULL, (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	     (igraph_t *) VECTOR(*graphs)[t+1],
	     (igraph_vector_int_t *) VECTOR(*vcolors)[t+1],
	     (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	     &neighborhood, union_graph, union_graph_vcolors, union_graph_ecolors);
      } else {
	igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t], NULL, NULL,
	     (igraph_t *) VECTOR(*graphs)[t+1], NULL, NULL,
	     &neighborhood, union_graph, union_graph_vcolors, union_graph_ecolors);
      }
    }

    igraph_llist_ptr_push_back(&result_list, union_graph);
    igraph_llist_ptr_push_back(&result_vcolors_list, union_graph_vcolors);
    igraph_llist_ptr_push_back(&result_ecolors_list, union_graph_ecolors);

    printf("timestep %ld: %ld changed, %ld neighs\n", t, changed_count,
		igraph_vector_int_size(&neighborhood));
  } // for t = 1...T

  // return result
  igraph_llist_ptr_to_vector(&result_list, result);
  igraph_llist_ptr_to_vector(&result_vcolors_list, result_vcolors);
  igraph_llist_ptr_to_vector(&result_ecolors_list, result_ecolors);

  // clean up
  igraph_vector_int_destroy(&changed_nodes);
  igraph_vector_int_destroy(&neighborhood);
  igraph_llist_ptr_destroy(&result_list);
  igraph_llist_ptr_destroy(&result_vcolors_list);
  igraph_llist_ptr_destroy(&result_ecolors_list);

  return 0;
}
