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

//#include <string.h> // memcpy
//#include "igraph_fsm.h"
//#include "igraph_matrix.h"
//#include "igraph_stack.h"
//#include "igraph_components.h"
//#include "igraph_constructors.h"

// TODO: right now, it reads all edges into a single graph
int igraph_read_dynamic_velist(igraph_vector_ptr_t *graphs, FILE *instream) {
  char buf[32];
  long int i1, i2, max_vid;
  char ve;
  igraph_t *graph;
  igraph_llist_int_t add_edges;
  igraph_llist_int_t del_edges;
  igraph_vector_t edges;

  IGRAPH_CHECK(igraph_llist_int_init(&add_edges));
  IGRAPH_CHECK(igraph_llist_int_init(&del_edges));

  // read vertices
  max_vid = -1;
  while (fgets(buf, 32, instream)) {
    sscanf(buf, "%c %ld %ld", &ve, &i1, &i2);
    if (ve == 'v') {
      if (i1 > max_vid) {
	max_vid = i1;
      }
    } else if (ve == 'e') {
      // we reached the first edge
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i1));
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i2));
    } else {
      // ignore lines not starting with v or e
    }
  }

  // read more edges
  while (fgets(buf, 32, instream)) {
    sscanf(buf, "%c %ld %ld", &ve, &i1, &i2);
    if (ve == 'e') {
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i1));
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i2));
    } else {
      // ignore lines not starting with v or e
    }
  }

  // initialize new graph
  graph = igraph_Calloc(1, igraph_t);
  IGRAPH_CHECK(igraph_empty(graph, max_vid+1, 0));

  // add edges
  IGRAPH_CHECK(igraph_vector_init(&edges, 0));
  IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));
  igraph_vector_destroy(&edges);

  // append to graph list
  IGRAPH_CHECK(igraph_vector_ptr_resize(graphs, 1));
  VECTOR(*graphs)[0] = graph;

  igraph_llist_int_destroy(&add_edges);
  igraph_llist_int_destroy(&del_edges);
  return 0;
}

