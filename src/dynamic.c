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
#include "igraph_dynamic.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_list.h"
#include "igraph_memory.h"


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

