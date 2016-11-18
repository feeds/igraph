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
#include <limits.h>
#include <time.h>
#include "igraph_fsm.h"
#include "igraph_dynamic.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_list.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_games.h"


int igraph_i_compute_joint_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_t *changed_nodes, igraph_vector_t *neighborhood);

int igraph_i_compute_dynamic_node_selectors_full(long int T, igraph_vector_ptr_t *node_selectors);
int igraph_i_compute_dynamic_node_selectors_neighbors(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_vector_ptr_t *node_selectors);
int igraph_i_compute_dynamic_node_selectors_event(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_vector_ptr_t *node_selectors);

int igraph_i_compute_union_graph_projection(igraph_t *graph1,
	  igraph_vector_int_t *vcolors1, igraph_vector_int_t *ecolors1,
	  igraph_t *graph2,
	  igraph_vector_int_t *vcolors2, igraph_vector_int_t *ecolors2,
	  igraph_vs_t node_selector,
	  igraph_t *union_graph,
	  igraph_vector_int_t *union_graph_vcolors,
	  igraph_vector_int_t *union_graph_ecolors,
	  igraph_integer_t max_vcolor, igraph_integer_t max_ecolor);


// assert: input graph must have edge timestamps!
// assert: edge colors > 0
// assert: maximum edge color appears in first timestep
//
// input file should be in the format:
//    v $vid [...]
//    e $vid1 $vid2 [$color] $creation [$deletion [...]]
// where -1 for $deletion means that the edge is never deleted and [...] is ignored. Edges
// must be sorted by $creation. Node ids start with 0. Example (unlabelled edges):
//    v 0
//    v 1
//    ...
//    e 0 3 0 -1
//    e 1 41 0 3
//    e 1 8 1 -1
//    ...
int igraph_read_and_project_dynamic_velist(FILE *instream, igraph_bool_t directed,
      igraph_bool_t has_vcolors, igraph_bool_t has_ecolors, igraph_bool_t has_etimesdel,
      igraph_projection_t proj_type, long int timestep_limit,
      igraph_vector_ptr_t *result_graphs, igraph_vector_ptr_t *result_vcolors,
      igraph_vector_ptr_t *result_ecolors, gzFile fgz) {
  char buf[32];
  long int i, i1, i2, i3, i4, i5, max_vid, timestamp, last_timestamp, n_fields, processed_graphs;
  long int max_ecolor = 0;
  long int tid = 0;
  igraph_integer_t eid;
  char ve;
  igraph_t *graph1=NULL, *graph2=NULL;
  igraph_llist_int_t add_edges, add_ecolors;
  igraph_llist_ptr_t del_edges; // holds a list of integer lists l, with l->first = del_time,
                                // followed by all edges (v1, v2) to delete at del_time
  igraph_llist_ptr_t result_list, result_vcolors_list, result_ecolors_list;
  igraph_vector_t edges, deletions;
  igraph_vector_int_t ecolors;
  igraph_vector_ptr_t tmp_ug_graphs, tmp_ug_vcolors, tmp_ug_ecolors;
  igraph_vector_ptr_t tmp_db, tmp_db_vcolors, tmp_db_ecolors;
  igraph_llist_item_ptr_t *item_ptr;
  igraph_llist_item_int_t *item_int;

  if (has_etimesdel && has_ecolors) {
    IGRAPH_ERROR("velist reader with edge colors and deletions not implemented",
	IGRAPH_UNIMPLEMENTED);
  }

  IGRAPH_CHECK(igraph_llist_int_init(&add_edges));
  IGRAPH_CHECK(igraph_llist_int_init(&add_ecolors));
  IGRAPH_CHECK(igraph_llist_ptr_init(&del_edges));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_list));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_vcolors_list));
  IGRAPH_CHECK(igraph_llist_ptr_init(&result_ecolors_list));
  IGRAPH_CHECK(igraph_vector_init(&edges, 0));
  IGRAPH_CHECK(igraph_vector_init(&deletions, 0));
  IGRAPH_CHECK(igraph_vector_int_init(&ecolors, 0));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db, 2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db_vcolors, 2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db_ecolors, 2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_ug_graphs, 0));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_ug_vcolors, 0));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_ug_ecolors, 0));

  // read vertices
  max_vid = -1;
  last_timestamp = -1;
  processed_graphs = 0;
  while (fgets(buf, 32, instream)) {
    i1 = -1; i2 = -1; i3 = -1; i4 = -1; i5 = -1;
    n_fields = sscanf(buf, "%c %ld %ld %ld %ld %ld", &ve, &i1, &i2, &i3, &i4, &i5);
    if (n_fields < 2) {
      // ignore lines that cannot be parsed
      continue;
    }
    if (ve == 'v') {
      if (i1 > max_vid) {
	max_vid = i1;
      }
      if (has_vcolors) {
	// TODO: store i2 in vcolor list (static!)
      }
    } else if (ve == 'e') {
      // we reached the first edge
      if (n_fields < 3) {
	IGRAPH_ERROR("edge definition incomplete, node ids missing", IGRAPH_PARSEERROR);
      }
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i1));
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i2));
      if (has_ecolors) {
	IGRAPH_CHECK(igraph_llist_int_push_back(&add_ecolors, i3));
	if (i3 > max_ecolor)
	  max_ecolor = i3;
	last_timestamp = i4;
	if (has_etimesdel && (i5 > last_timestamp)) {
	  // init edge list for new deletion timestamp and mark edge for deletion at time i5
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&del_edges,igraph_Calloc(1, igraph_llist_int_t)));
	  IGRAPH_CHECK(igraph_llist_int_init((igraph_llist_int_t *) del_edges.first->data));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i5));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i1));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i2));
	}
      } else {
	last_timestamp = i3;
	if (has_etimesdel && (i4 > last_timestamp)) {
	  // init edge list for new deletion timestamp and mark edge (eid 0) for deletion at time i4
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&del_edges,igraph_Calloc(1, igraph_llist_int_t)));
	  IGRAPH_CHECK(igraph_llist_int_init((igraph_llist_int_t *) del_edges.first->data));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i4));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i1));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.first->data,i2));
	}
      }
      break;
    } else {
      // ignore lines not starting with v or e
    }
  }

  // read more edges
  timestamp = -1;
  while (fgets(buf, 32, instream)) {
    i1 = -1; i2 = -1; i3 = -1; i4 = -1; i5 = -1;
    n_fields = sscanf(buf, "%c %ld %ld %ld %ld %ld", &ve, &i1, &i2, &i3, &i4, &i5);
    if (n_fields < 2) {
      // ignore lines that cannot be parsed
      continue;

    }
    if (ve != 'e') {
      // ignore lines that do not define an edge
      continue;
    }

    timestamp = (has_ecolors ? i4 : i3);
    if (timestamp > last_timestamp) {
      // we have processed all edges from this timestamp, construct graph

      // initialize new graph from previous one (or empty one, if this is the first timestamp)
      if (graph1 == NULL) {
	graph1 = igraph_Calloc(1, igraph_t);
	IGRAPH_CHECK(igraph_empty(graph1, max_vid+1, directed));
	if (proj_type == IGRAPH_PROJECTION_NONE) {
	  if (fgz != NULL) {
	    // output to file
	  } else {
	    // add to result list
	    IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph1));
	  }
	}
      } else if (graph2 == NULL) {
	graph2 = igraph_Calloc(1, igraph_t);
	IGRAPH_CHECK(igraph_copy(graph2, graph1));
	if (proj_type == IGRAPH_PROJECTION_NONE) {
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph2));
	}
      } else {
	if (proj_type == IGRAPH_PROJECTION_NONE) {
	  graph1 = graph2;
	  graph2 = igraph_Calloc(1, igraph_t);
	  IGRAPH_CHECK(igraph_copy(graph2, graph1));
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph2));
	} else {
	  igraph_destroy(graph1);
	  IGRAPH_CHECK(igraph_copy(graph1, graph2));
	}
      }

      // prepare edge insertions
      IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges, 0));
      IGRAPH_CHECK(igraph_llist_int_to_vector(&add_ecolors, &ecolors, 0));

      // prepare edge deletions
      IGRAPH_CHECK(igraph_vector_resize(&deletions, 0));
      if (has_etimesdel) {
	// check if there are remaining edge deletions to perform at this timestamp
	// NOTE: if there are timestamps with only edge deletions, these deletions are
	//       peformed at the next timestamp that has insertions
	for (item_ptr = del_edges.first; item_ptr != NULL; item_ptr = item_ptr->next) {
	  if (((igraph_llist_int_t *) item_ptr->data)->first->data <= last_timestamp) {
	    if (item_ptr != NULL) {
	      igraph_vector_resize(&deletions, igraph_vector_size(&deletions)
		    + (igraph_llist_int_size((igraph_llist_int_t *) item_ptr->data)-1)/2);
	      // retrieve eids
	      for (item_int = ((igraph_llist_int_t *)item_ptr->data)->first->next, i = 0;
		 	  item_int != NULL;
		 	  item_int = item_int->next->next, i++) {
		IGRAPH_CHECK(igraph_get_eid(((graph2 != NULL) ? graph2 : graph1),
		      &eid, item_int->data, item_int->next->data, 1, 1));
		VECTOR(deletions)[igraph_vector_size(&deletions)
		    - (igraph_llist_int_size((igraph_llist_int_t *) item_ptr->data)-1)/2+i] = eid;
	      }
	      // destroy deletion list (and add dummy)
	      igraph_llist_int_destroy((igraph_llist_int_t *)item_ptr->data);
	      igraph_llist_int_init((igraph_llist_int_t *)item_ptr->data);
	      igraph_llist_int_push_back((igraph_llist_int_t *)item_ptr->data, INT_MAX);
	    }
	  }
	}
      }

      // perform insertions/deletions
      // TODO: ecolor deletions
      if (graph2 == NULL) {
	IGRAPH_CHECK(igraph_delete_edges(graph1, igraph_ess_vector(&deletions)));
	IGRAPH_CHECK(igraph_add_edges(graph1, &edges, 0));
	printf("graph %ld (+%ld, -%ld): ", last_timestamp, igraph_vector_size(&edges)/2,
		  igraph_vector_size(&deletions));
	igraph_print_stats(graph1);
      } else {
	IGRAPH_CHECK(igraph_delete_edges(graph2, igraph_ess_vector(&deletions)));
	IGRAPH_CHECK(igraph_add_edges(graph2, &edges, 0));
	printf("graph %ld (+%ld, -%ld): ", last_timestamp, igraph_vector_size(&edges)/2,
		  igraph_vector_size(&deletions));
	igraph_print_stats(graph2);
      }

      // compute the projection
      if ((proj_type != IGRAPH_PROJECTION_NONE) && (graph2 != NULL)) {
	VECTOR(tmp_db)[0] = graph1;
	VECTOR(tmp_db)[1] = graph2;
	if (has_ecolors) {
	  VECTOR(tmp_db_ecolors)[0] = &ecolors;
	  VECTOR(tmp_db_ecolors)[1] = &ecolors; // they are static
	}
	if (fgz == NULL) {
	  IGRAPH_CHECK(igraph_compute_dynamic_union_graph_projection(&tmp_db,
		  /*vcolors=*/ NULL, (has_ecolors ? &tmp_db_ecolors : NULL), proj_type,
		  &tmp_ug_graphs, /*tmp_ug_vcolors=*/ NULL, &tmp_ug_ecolors,
		  /*max_vcolor=*/ 0, (has_ecolors ? max_ecolor : 0)));
	  for (i = 0; i < igraph_vector_ptr_size(&tmp_ug_graphs); i++) {
	    IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, VECTOR(tmp_ug_graphs)[i]));
	    //IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_vcolors_list,VECTOR(tmp_ug_vcolors)[i]));
	    IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_ecolors_list, VECTOR(tmp_ug_ecolors)[i]));
	  }
	} else {
	  IGRAPH_CHECK(igraph_write_dynamic_union_graph_projection(&tmp_db,
		  /*vcolors=*/ NULL, (has_ecolors ? &tmp_db_ecolors : NULL), proj_type,
		  /*max_vcolor=*/ 0, (has_ecolors ? max_ecolor : 0), fgz, &tid));
	}
      }

      // destroy old edge list and create new empty ones
      igraph_llist_int_destroy(&add_edges);
      IGRAPH_CHECK(igraph_llist_int_init(&add_edges));

      if (last_timestamp < timestamp-1) {
	// current and last timestamp are not consecutive
	// NOTE: this is only a problem if we have timesteps that contain only edge deletions
	//       and no insertions
	printf("+++ gap in edge timestamp, ignoring\n");
      }
      last_timestamp = timestamp;
      processed_graphs++;
      if ((timestep_limit >= 0) && (processed_graphs >= timestep_limit)) {
	break;
      }
    } else if (timestamp < last_timestamp) {
      IGRAPH_ERROR("edges not sorted by timestamp", IGRAPH_PARSEERROR);
      return 1;
    }

    IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i1));
    IGRAPH_CHECK(igraph_llist_int_push_back(&add_edges, i2));

    if (has_ecolors) {
      IGRAPH_CHECK(igraph_llist_int_push_back(&add_ecolors, i3));
      if (i3 > max_ecolor)
	max_ecolor = i3;
      if (has_etimesdel && (i5 > last_timestamp)) {
	// check if the deletion timestamp i5 already exists in list
	for (item_ptr = del_edges.first; item_ptr != NULL; item_ptr = item_ptr->next) {
	  if (((igraph_llist_int_t *) item_ptr->data)->first->data == i5) {
	    break;
	  }
	}
	if (item_ptr != NULL) {
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) item_ptr->data,i1));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) item_ptr->data,i2));
	} else {
	  // init edge list for new deletion timestamp i5
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&del_edges,igraph_Calloc(1, igraph_llist_int_t)));
	  IGRAPH_CHECK(igraph_llist_int_init((igraph_llist_int_t *) del_edges.last->data));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i5));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i1));
	  IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i2));
	}
      }
    } else if (has_etimesdel && (i4 > last_timestamp)) {
      // check if the deletion timestamp i4 already exists in list
      for (item_ptr = del_edges.first; item_ptr != NULL; item_ptr = item_ptr->next) {
	if (((igraph_llist_int_t *) item_ptr->data)->first->data == i4) {
	  break;
	}
      }
      if (item_ptr != NULL) {
	IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) item_ptr->data,i1));
	IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) item_ptr->data,i2));
      } else {
	// init edge list for new deletion timestamp i4
	IGRAPH_CHECK(igraph_llist_ptr_push_back(&del_edges,igraph_Calloc(1, igraph_llist_int_t)));
	IGRAPH_CHECK(igraph_llist_int_init((igraph_llist_int_t *) del_edges.last->data));
	IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i4));
	IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i1));
	IGRAPH_CHECK(igraph_llist_int_push_back((igraph_llist_int_t *) del_edges.last->data,i2));
      }
    }
  }

  // process last timestamp
  if ((timestep_limit < 0) || (processed_graphs < timestep_limit)) {
    if (graph1 == NULL) {
      graph1 = igraph_Calloc(1, igraph_t);
      IGRAPH_CHECK(igraph_empty(graph1, max_vid+1, directed));
      if (proj_type == IGRAPH_PROJECTION_NONE) {
	IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph1));
      } else {
	IGRAPH_ERROR("one timestep not enough to compute any projections!", IGRAPH_EINVAL);
      }
    } else if (graph2 == NULL) {
      graph2 = igraph_Calloc(1, igraph_t);
      IGRAPH_CHECK(igraph_copy(graph2, graph1));
      if (proj_type == IGRAPH_PROJECTION_NONE) {
	IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph2));
      }
    } else {
      if (proj_type == IGRAPH_PROJECTION_NONE) {
	graph1 = graph2;
	graph2 = igraph_Calloc(1, igraph_t);
	IGRAPH_CHECK(igraph_copy(graph2, graph1));
	IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, graph2));
      } else {
	igraph_destroy(graph1);
	IGRAPH_CHECK(igraph_copy(graph1, graph2));
      }
    }
    IGRAPH_CHECK(igraph_llist_int_to_vector_real(&add_edges, &edges, 0));
    IGRAPH_CHECK(igraph_llist_int_to_vector(&add_ecolors, &ecolors, 0));
    IGRAPH_CHECK(igraph_vector_resize(&deletions, 0));
    if (has_etimesdel) {
      for (item_ptr = del_edges.first; item_ptr != NULL; item_ptr = item_ptr->next) {
	if (((igraph_llist_int_t *) item_ptr->data)->first->data <= last_timestamp) {
	  if (item_ptr != NULL) {
	    igraph_vector_resize(&deletions, igraph_vector_size(&deletions)
		  + (igraph_llist_int_size((igraph_llist_int_t *) item_ptr->data)-1)/2);
	    // retrieve eids
	    for (item_int = ((igraph_llist_int_t *)item_ptr->data)->first->next, i = 0;
			item_int != NULL;
			item_int = item_int->next->next, i++) {
	      IGRAPH_CHECK(igraph_get_eid(((graph2 != NULL) ? graph2 : graph1),
		    &eid, item_int->data, item_int->next->data, 1, 1));
	      VECTOR(deletions)[igraph_vector_size(&deletions)
		  - (igraph_llist_int_size((igraph_llist_int_t *) item_ptr->data)-1)/2+i] = eid;
	    }
	    // destroy deletion list (and add dummy)
	    igraph_llist_int_destroy((igraph_llist_int_t *)item_ptr->data);
	    igraph_llist_int_init((igraph_llist_int_t *)item_ptr->data);
	    igraph_llist_int_push_back((igraph_llist_int_t *)item_ptr->data, INT_MAX);
	  }
	}
      }
    }
    IGRAPH_CHECK(igraph_delete_edges(graph2, igraph_ess_vector(&deletions)));
    IGRAPH_CHECK(igraph_add_edges(graph2, &edges, 0));
    // TODO: ecolor deletions
    printf("graph %ld (+%ld, -%ld): ", timestamp, igraph_vector_size(&edges)/2,
	      igraph_vector_size(&deletions));
    igraph_print_stats(graph2);
    if (proj_type != IGRAPH_PROJECTION_NONE) {
      VECTOR(tmp_db)[0] = graph1;
      VECTOR(tmp_db)[1] = graph2;
      if (has_ecolors) {
	VECTOR(tmp_db_ecolors)[0] = &ecolors;
	VECTOR(tmp_db_ecolors)[1] = &ecolors; // they are static
      }
      if (fgz == NULL) {
	IGRAPH_CHECK(igraph_compute_dynamic_union_graph_projection(&tmp_db,
		/*vcolors=*/ NULL, (has_ecolors ? &tmp_db_ecolors : NULL), proj_type,
		&tmp_ug_graphs, /*tmp_ug_vcolors=*/ NULL, &tmp_ug_ecolors,
		/*max_vcolor=*/ 0, (has_ecolors ? max_ecolor : 0)));
	for (i = 0; i < igraph_vector_ptr_size(&tmp_ug_graphs); i++) {
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_list, VECTOR(tmp_ug_graphs)[i]));
	  //IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_vcolors_list,VECTOR(tmp_ug_vcolors)[i]));
	  IGRAPH_CHECK(igraph_llist_ptr_push_back(&result_ecolors_list, VECTOR(tmp_ug_ecolors)[i]));
	}
      } else {
	IGRAPH_CHECK(igraph_write_dynamic_union_graph_projection(&tmp_db,
		/*vcolors=*/ NULL, (has_ecolors ? &tmp_db_ecolors : NULL), proj_type,
		/*max_vcolor=*/ 0, (has_ecolors ? max_ecolor : 0), fgz, &tid));
      }
    }
  }

  // convert lists into vector for output
  if (fgz == NULL) {
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_list, result_graphs, 0));
    //IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_vcolors_list, result_vcolors, 0));
    IGRAPH_CHECK(igraph_llist_ptr_to_vector(&result_ecolors_list, result_ecolors, 0));
  }

  printf("+++ ecolors %ld min %d max %d\n", igraph_vector_int_size(&ecolors),
      igraph_vector_int_min(&ecolors), igraph_vector_int_max(&ecolors));

  igraph_llist_int_destroy(&add_edges);
  igraph_llist_int_destroy(&add_ecolors);
  igraph_llist_ptr_destroy(&del_edges); // TODO: free all lists
  igraph_llist_ptr_destroy(&result_list);
  igraph_llist_ptr_destroy(&result_vcolors_list);
  igraph_llist_ptr_destroy(&result_ecolors_list);
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&deletions);
  igraph_vector_ptr_destroy(&tmp_db);
  igraph_vector_ptr_destroy(&tmp_db_vcolors);
  igraph_vector_ptr_destroy(&tmp_db_ecolors);
  igraph_vector_ptr_destroy(&tmp_ug_graphs);
  igraph_vector_ptr_destroy(&tmp_ug_vcolors);
  igraph_vector_ptr_destroy(&tmp_ug_ecolors);

  return 0;
}


// assert: seed_nodes is sorted ascendingly
// assert: all nodes from seed_nodes appear in graph1 and graph2
// assert: undirected graphs (why? seems to work with directed graphs)
//
// naive algorithm, differs from thesis!
// runtime: O(N_seeds*k_max + N_seeds*k_max*log(N_seeds*k_max) + (N_seeds*k_max+N_seeds))
int igraph_i_compute_joint_neighborhood(igraph_t *graph1, igraph_t *graph2,
      igraph_vector_t *seed_nodes, igraph_vector_t *neighborhood) {
  igraph_llist_t node_list;
  long int i, j, k, changed_count, neigh, neigh_size, node;
  igraph_vector_t neighborhood_dups;

  changed_count = igraph_vector_size(seed_nodes);

  igraph_llist_init(&node_list);
  igraph_vector_init(&neighborhood_dups, 0);

  for (i = 0; i < changed_count; i++) {
    node = VECTOR(*seed_nodes)[i];
    // collect neighbors in graph1
    for (j = 0; j < DEGREE(*graph1, node); j++) {
      neigh = NEIGHBOR(*graph1, node, j);
      if (!igraph_vector_binsearch2(seed_nodes, neigh)) { // filter to avoid duplicates
	igraph_llist_push_back(&node_list, neigh);
      }
    }
    // collect neighbors in graph2
    for (j = 0; j < DEGREE(*graph2, node); j++) {
      neigh = NEIGHBOR(*graph2, node, j);
      if (!igraph_vector_binsearch2(seed_nodes, neigh)) {
	igraph_llist_push_back(&node_list, neigh);
      }
    }
  }
  igraph_llist_to_vector(&node_list, &neighborhood_dups, 0); // contains duplicates!
  igraph_vector_sort(&neighborhood_dups);
  neigh_size = igraph_vector_size(&neighborhood_dups);

  // sort neighborhood and merge it with the seed nodes, while removing duplicates

  igraph_vector_resize(neighborhood, changed_count+neigh_size); // max. possible size
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
// assert: node_selector (if created from vector) must be sorted ascendingly
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
  igraph_bool_t directed = igraph_is_directed(graph1);

  // make sure that we only get supported vertex selectors
  if (!((igraph_vs_type(&node_selector) == IGRAPH_VS_1)
      || (igraph_vs_type(&node_selector) == IGRAPH_VS_ALL)
      || (igraph_vs_type(&node_selector) == IGRAPH_VS_SEQ)
      || (igraph_vs_type(&node_selector) == IGRAPH_VS_VECTOR)
      || (igraph_vs_type(&node_selector) == IGRAPH_VS_VECTORPTR))) {
    IGRAPH_ERROR("invalid node selector", IGRAPH_EINVAL);
  }
  IGRAPH_CHECK(igraph_vs_size(graph1, &node_selector, &vcount)); // same result for graph2

  IGRAPH_CHECK(igraph_empty(union_graph, vcount, directed));
  if (has_vcolors) {
    IGRAPH_CHECK(igraph_vector_int_resize(union_graph_vcolors, vcount));
  }

  // create backward index for vertex id lookup
  // TODO: we don't need this for VS_1 and VS_ALL
  igraph_vit_create(graph1, node_selector, &vit);
  igraph_vector_int_init(&bw_index, igraph_vcount(graph1));
  for (i=0; (i<vcount) && (!IGRAPH_VIT_END(vit)); i++, IGRAPH_VIT_NEXT(vit)) {
    // NOTE: we add 1 to the vertex id i to distinguish uninitialized elements from
    // bw_index (value 0) from initialized elements. we have to subtract 1 again when using
    // the backward index
    VECTOR(bw_index)[(long int) IGRAPH_VIT_GET(vit)] = i+1;
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
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)] > 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,
		VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)]-1);
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
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] > 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,
		VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]-1);
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
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] > 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,
		VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]-1);
	  if (has_ecolors) {
	    eid1 = OUT_NEIGH_TO_EID(*graph1, orig_node, j);
	    igraph_llist_int_push_back(&ecolors_list, (max_ecolor+1)*VECTOR(*ecolors1)[eid1]);
	  } else {
	    igraph_llist_int_push_back(&ecolors_list, 0b10);
	  }
	}
	j++;
      } else if (OUT_NEIGHBOR(*graph1, orig_node, j) == OUT_NEIGHBOR(*graph2, orig_node, k)) {
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)] > 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,
		VECTOR(bw_index)[OUT_NEIGHBOR(*graph1,orig_node,j)]-1);
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
	if (VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)] > 0) {
	  igraph_llist_int_push_back(&edge_list, i);
	  igraph_llist_int_push_back(&edge_list,
		VECTOR(bw_index)[OUT_NEIGHBOR(*graph2,orig_node,k)]-1);
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
  igraph_llist_int_to_vector_real(&edge_list, &edges, 0);
  IGRAPH_CHECK(igraph_add_edges(union_graph, &edges, 0));
  igraph_llist_int_to_vector(&ecolors_list, union_graph_ecolors, 0);

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


// node_selectors points to a length T-1 igraph_vector_ptr object that contains for each
// timestep an igraph_llist_ptr_t object with a list of node selectors (igraph_vs_t *)
int igraph_i_compute_dynamic_node_selectors_full(long int T, igraph_vector_ptr_t *node_selectors) {
  long int t;
  igraph_vs_t *node_selector;
  IGRAPH_CHECK(igraph_vector_ptr_resize(node_selectors, T-1));
  for (t = 0; t < T-1; t++) {
    VECTOR(*node_selectors)[t] = igraph_Calloc(1, igraph_llist_ptr_t);
    IGRAPH_CHECK(igraph_llist_ptr_init((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t]));
    node_selector = igraph_Calloc(1, igraph_vs_t);
    IGRAPH_CHECK(igraph_vs_all(node_selector));
    IGRAPH_CHECK(igraph_llist_ptr_push_back((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t],
		  node_selector));
  }
  return 0;
}


int igraph_i_compute_dynamic_node_selectors_neighbors(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_vector_ptr_t *node_selectors) {
  igraph_llist_int_t node_list;
  igraph_vector_t changed_nodes, neighborhood;
  igraph_vs_t *node_selector;
  long int t, i, j, eid1, eid2;
  long int N = igraph_vcount((igraph_t *) VECTOR(*graphs)[0]);
  long int T = igraph_vector_ptr_size(graphs);

  IGRAPH_CHECK(igraph_vector_init(&changed_nodes, 0));
  IGRAPH_CHECK(igraph_vector_init(&neighborhood, 0));

  // iterate over all timestamps and collect all changed nodes
  IGRAPH_CHECK(igraph_vector_ptr_resize(node_selectors, T-1));
  for (t = 0; t < T-1; t++) {
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
    igraph_llist_int_to_vector_real(&node_list, &changed_nodes, 0); // sorted by construction
    igraph_llist_int_destroy(&node_list);

    // compute joint 1-hop neighborhood at timesteps t and t+1
    igraph_i_compute_joint_neighborhood((igraph_t *) VECTOR(*graphs)[t],
          (igraph_t *) VECTOR(*graphs)[t+1], &changed_nodes, &neighborhood);

    // create node selector list for this timestep
    VECTOR(*node_selectors)[t] = igraph_Calloc(1, igraph_llist_ptr_t);
    IGRAPH_CHECK(igraph_llist_ptr_init((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t]));
    node_selector = igraph_Calloc(1, igraph_vs_t);
    IGRAPH_CHECK(igraph_vs_vector_copy(node_selector, &neighborhood));
    igraph_llist_ptr_push_back((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t], node_selector);
  }

  igraph_vector_destroy(&changed_nodes);
  igraph_vector_destroy(&neighborhood);

  return 0;
}

int igraph_i_compute_dynamic_node_selectors_event(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_vector_ptr_t *node_selectors) {
  igraph_vector_t node_event, edge_event, neighborhood;
  igraph_vs_t *node_selector;
  long int t, i, j, k, nei1, nei2, eid1, eid2;
  long int N = igraph_vcount((igraph_t *) VECTOR(*graphs)[0]);
  long int T = igraph_vector_ptr_size(graphs);
  igraph_bool_t new_event;

  IGRAPH_CHECK(igraph_vector_init(&node_event, 1));
  IGRAPH_CHECK(igraph_vector_init(&edge_event, 2));
  IGRAPH_CHECK(igraph_vector_init(&neighborhood, 0));

  // iterate over all timestamps and collect node and edge events
  IGRAPH_CHECK(igraph_vector_ptr_resize(node_selectors, T-1));
  for (t = 0; t < T-1; t++) {
    VECTOR(*node_selectors)[t] = igraph_Calloc(1, igraph_llist_ptr_t);
    IGRAPH_CHECK(igraph_llist_ptr_init((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t]));

    for (i = 0; i < N; i++) {
      // check for node event
      if ((vcolors != NULL) && (VECTOR(*(igraph_vector_int_t *) VECTOR(*vcolors)[t])[i]
            != VECTOR(*(igraph_vector_int_t *) VECTOR(*vcolors)[t+1])[i])) {
	VECTOR(node_event)[0] = i;
	igraph_i_compute_joint_neighborhood((igraph_t *) VECTOR(*graphs)[t],
		  (igraph_t *) VECTOR(*graphs)[t+1], &node_event, &neighborhood);
	node_selector = igraph_Calloc(1, igraph_vs_t);
	IGRAPH_CHECK(igraph_vs_vector_copy(node_selector, &neighborhood));
	igraph_llist_ptr_push_back((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t], node_selector);
      }

      // check for edge events, i.e. check which incident edges have changed
      j = 0; // neighbor index at timestep t
      k = 0; // neighbor index at timestep t+1
      while ((j < OUT_DEGREE(*(igraph_t *) VECTOR(*graphs)[t], i)) ||
		  (k < OUT_DEGREE(*(igraph_t *) VECTOR(*graphs)[t+1], i))) {
	nei1 = OUT_NEIGHBOR(*(igraph_t *) VECTOR(*graphs)[t], i, j);
	nei2 = OUT_NEIGHBOR(*(igraph_t *) VECTOR(*graphs)[t+1], i, k);
	new_event = 0;

	if (j == OUT_DEGREE(*(igraph_t *) VECTOR(*graphs)[t], i)) {
	  // neighbors at timestep t exhausted, all remaining neighbors at time t+1
	  // make edge insertion events
	  VECTOR(edge_event)[0] = i;
	  VECTOR(edge_event)[1] = nei2;
	  new_event = 1;
	  k++;
	} else if (k == OUT_DEGREE(*(igraph_t *) VECTOR(*graphs)[t+1], i)) {
	  // neighbors at timestep t+1 exhausted, all remaining neighbors at time t
	  // make edge deletion events
	  VECTOR(edge_event)[0] = i;
	  VECTOR(edge_event)[1] = nei1;
	  new_event = 1;
	  j++;
	} else {
	  // neighbors remaining at both timesteps
          if (nei1 < nei2) {
	    // edge event: link to nei1 has been removed
	    VECTOR(edge_event)[0] = i;
	    VECTOR(edge_event)[1] = nei1;
	    new_event = 1;
	    j++;
          } else if (nei1 == nei2) {
	    // neighbors are equal, check the edge label
            if (ecolors != NULL) {
              eid1 = NEIGH_TO_EID(*(igraph_t *) VECTOR(*graphs)[t], i, j);
              eid2 = NEIGH_TO_EID(*(igraph_t *) VECTOR(*graphs)[t+1], i, k);
              if (VECTOR(*(igraph_vector_int_t *) VECTOR(*ecolors)[t])[eid1]
		    != VECTOR(*(igraph_vector_int_t *) VECTOR(*ecolors)[t+1])[eid2]) {
		// edge event: link to nei1/nei2 has been relabelled
		VECTOR(edge_event)[0] = i;
		VECTOR(edge_event)[1] = nei1;
		new_event = 1;
              }
            }
	    j++;
	    k++;
	  } else if (nei1 > nei2) {
	    // edge event: link to nei2 has been added
	    VECTOR(edge_event)[0] = i;
	    VECTOR(edge_event)[1] = nei2;
	    new_event = 1;
	    k++;
	  }
	}

	// create node selector for edge event
	if (new_event) {
	  igraph_vector_sort(&edge_event);
	  igraph_i_compute_joint_neighborhood((igraph_t *) VECTOR(*graphs)[t],
		    (igraph_t *) VECTOR(*graphs)[t+1], &edge_event, &neighborhood);
	  node_selector = igraph_Calloc(1, igraph_vs_t);
	  //igraph_vector_print(&edge_event);
	  //igraph_vector_print(&neighborhood);
	  IGRAPH_CHECK(igraph_vs_vector_copy(node_selector, &neighborhood));
	  igraph_llist_ptr_push_back((igraph_llist_ptr_t *) VECTOR(*node_selectors)[t],
				      node_selector);
	  new_event = 0;
	}
      } // loop over all neighbors of i at timesteps t and t+1
    } // for all nodes i
  } // for all timesteps t

  igraph_vector_destroy(&node_event);
  igraph_vector_destroy(&edge_event);
  igraph_vector_destroy(&neighborhood);

  return 0;
}

// assert: all graphs have the same number of nodes, and node IDs correspond with each other
int igraph_write_dynamic_union_graph_projection(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_projection_t proj_type,
      igraph_integer_t max_vcolor, igraph_integer_t max_ecolor,
      gzFile fgz, long int *tid) {
  long int T = igraph_vector_ptr_size(graphs);
  long int t;
  igraph_t union_graph;
  igraph_vector_int_t union_graph_vcolors, union_graph_ecolors;
  igraph_vs_t *node_selector;
  igraph_vector_ptr_t node_selectors;
  igraph_llist_ptr_t *node_selector_list;
  igraph_llist_item_ptr_t *item_ptr;

  igraph_vector_ptr_init(&node_selectors, 0);

  // depending on the projection type, compute the node selectors for all timesteps
  switch (proj_type) {
    case IGRAPH_PROJECTION_EVENT:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_event(graphs,
					      vcolors, ecolors, &node_selectors));
      break;
    case IGRAPH_PROJECTION_NEIGHBORS:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_neighbors(graphs,
					      vcolors, ecolors, &node_selectors));
      break;
    case IGRAPH_PROJECTION_FULL:
    default:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_full(T, &node_selectors));
      break;
  }

  // compute the union graph projection for all node selectors at each timestep
  for (t = 0; t < T-1; t++) {
    node_selector_list = (igraph_llist_ptr_t *) VECTOR(node_selectors)[t];
    for (item_ptr = node_selector_list->first; item_ptr != NULL; item_ptr = item_ptr->next) {
      node_selector = (igraph_vs_t *) item_ptr->data;
      igraph_vector_int_init(&union_graph_vcolors, 0);
      igraph_vector_int_init(&union_graph_ecolors, 0);
      if (vcolors != NULL) {
	if (ecolors != NULL) {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t],
	       (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	       (igraph_t *) VECTOR(*graphs)[t+1],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t+1],
	       (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	       *node_selector, &union_graph, &union_graph_vcolors, &union_graph_ecolors,
	       max_vcolor, max_ecolor);
	} else {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t], NULL,
	       (igraph_t *) VECTOR(*graphs)[t+1],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t+1], NULL,
	       *node_selector, &union_graph, &union_graph_vcolors, &union_graph_ecolors,
	       max_vcolor, max_ecolor);
	}
      } else {
	if (ecolors != NULL) {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	       NULL, (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	       (igraph_t *) VECTOR(*graphs)[t+1],
	       NULL,
	       (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	       *node_selector, &union_graph, &union_graph_vcolors, &union_graph_ecolors,
	       max_vcolor, max_ecolor);
	} else {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t], NULL, NULL,
	       (igraph_t *) VECTOR(*graphs)[t+1], NULL, NULL,
	       *node_selector, &union_graph, &union_graph_vcolors, &union_graph_ecolors,
	       max_vcolor, max_ecolor);
	}
      }

      //printf("ug %ld: vcount %d ecount %d\n", t, igraph_vcount(union_graph), igraph_ecount(union_graph));
      gzprintf(fgz, "t # %ld\n", *tid);
      igraph_write_colored_graph_gz(&union_graph, (vcolors != NULL) ? &union_graph_vcolors : NULL,
		&union_graph_ecolors, /*etimes*/ NULL, fgz);
      (*tid)++;

      igraph_vector_int_destroy(&union_graph_vcolors);
      igraph_vector_int_destroy(&union_graph_ecolors);
      igraph_destroy(&union_graph);

      igraph_vs_destroy(node_selector);
      igraph_free(node_selector);
    } // for all node selectors at timestep t

    igraph_llist_ptr_destroy(node_selector_list);
    igraph_free(node_selector_list);
  } // for t = 1...T

  igraph_vector_ptr_destroy(&node_selectors);

  return 0;
}

// assert: all graphs have the same number of nodes, and node IDs correspond with each other
int igraph_compute_dynamic_union_graph_projection(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_projection_t proj_type,
      igraph_vector_ptr_t *result, igraph_vector_ptr_t *result_vcolors,
      igraph_vector_ptr_t *result_ecolors,
      igraph_integer_t max_vcolor, igraph_integer_t max_ecolor) {
  long int T = igraph_vector_ptr_size(graphs);
  long int t;
  igraph_llist_ptr_t result_list, result_vcolors_list, result_ecolors_list;
  igraph_t *union_graph;
  igraph_vector_int_t *union_graph_vcolors, *union_graph_ecolors;
  igraph_vs_t *node_selector;
  igraph_vector_ptr_t node_selectors;
  igraph_llist_ptr_t *node_selector_list;
  igraph_llist_item_ptr_t *item_ptr;

  igraph_llist_ptr_init(&result_list);
  igraph_llist_ptr_init(&result_vcolors_list);
  igraph_llist_ptr_init(&result_ecolors_list);
  igraph_vector_ptr_init(&node_selectors, 0);

  // depending on the projection type, compute the node selectors for all timesteps
  switch (proj_type) {
    case IGRAPH_PROJECTION_EVENT:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_event(graphs,
					      vcolors, ecolors, &node_selectors));
      break;
    case IGRAPH_PROJECTION_NEIGHBORS:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_neighbors(graphs,
					      vcolors, ecolors, &node_selectors));
      break;
    case IGRAPH_PROJECTION_FULL:
    default:
      IGRAPH_CHECK(igraph_i_compute_dynamic_node_selectors_full(T, &node_selectors));
      break;
  }

  // compute the union graph projection for all node selectors at each timestep
  for (t = 0; t < T-1; t++) {
    node_selector_list = (igraph_llist_ptr_t *) VECTOR(node_selectors)[t];
    for (item_ptr = node_selector_list->first; item_ptr != NULL; item_ptr = item_ptr->next) {
      node_selector = (igraph_vs_t *) item_ptr->data;
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
	       *node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	       max_vcolor, max_ecolor);
	} else {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t], NULL,
	       (igraph_t *) VECTOR(*graphs)[t+1],
	       (igraph_vector_int_t *) VECTOR(*vcolors)[t+1], NULL,
	       *node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	       max_vcolor, max_ecolor);
	}
      } else {
	if (ecolors != NULL) {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t],
	       NULL, (igraph_vector_int_t *) VECTOR(*ecolors)[t],
	       (igraph_t *) VECTOR(*graphs)[t+1],
	       NULL,
	       (igraph_vector_int_t *) VECTOR(*ecolors)[t+1],
	       *node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	       max_vcolor, max_ecolor);
	} else {
	  igraph_i_compute_union_graph_projection((igraph_t *) VECTOR(*graphs)[t], NULL, NULL,
	       (igraph_t *) VECTOR(*graphs)[t+1], NULL, NULL,
	       *node_selector, union_graph, union_graph_vcolors, union_graph_ecolors,
	       max_vcolor, max_ecolor);
	}
      }

      //printf("ug %ld: vcount %d ecount %d\n", t, igraph_vcount(union_graph), igraph_ecount(union_graph));
      igraph_llist_ptr_push_back(&result_list, union_graph);
      igraph_llist_ptr_push_back(&result_vcolors_list, union_graph_vcolors);
      igraph_llist_ptr_push_back(&result_ecolors_list, union_graph_ecolors);

      igraph_vs_destroy(node_selector);
      igraph_free(node_selector);
    } // for all node selectors at timestep t

    igraph_llist_ptr_destroy(node_selector_list);
    igraph_free(node_selector_list);
  } // for t = 1...T

  // return result
  // TODO: free memory if ptrs == NULL
  if (result != NULL) {
    igraph_llist_ptr_to_vector(&result_list, result, 0);
  }
  if (result_vcolors != NULL) {
    igraph_llist_ptr_to_vector(&result_vcolors_list, result_vcolors, 0);
  }
  if (result_ecolors != NULL) {
    igraph_llist_ptr_to_vector(&result_ecolors_list, result_ecolors, 0);
  }

  // clean up
  igraph_vector_ptr_destroy(&node_selectors);
  igraph_llist_ptr_destroy(&result_list);
  igraph_llist_ptr_destroy(&result_vcolors_list);
  igraph_llist_ptr_destroy(&result_ecolors_list);

  return 0;
}


int igraph_read_transactions_velist(FILE *instream, igraph_bool_t directed,
	igraph_bool_t has_vcolors,
	igraph_bool_t has_ecolors,
	igraph_bool_t has_etimes,
	igraph_vector_ptr_t *graphs,
	igraph_vector_ptr_t *vcolors,
	igraph_vector_ptr_t *ecolors,
	igraph_vector_ptr_t *etimes,
	igraph_vector_long_t *supps) {
  char buf[32];
  int timestamp;
  long int i1, i2, i3, i4, n_fields, max_vid = 0;
  long int supp = -1;
  igraph_llist_ptr_t result_graphs_list;
  igraph_llist_ptr_t result_vcolors_list;
  igraph_llist_ptr_t result_ecolors_list;
  igraph_llist_ptr_t result_etimes_list;
  igraph_llist_long_t result_supps_list;
  igraph_t *graph;
  igraph_llist_t edge_list;
  igraph_llist_int_t vcolor_list;
  igraph_llist_int_t ecolor_list;
  igraph_llist_int_t etimes_list;
  igraph_vector_t edges;
  igraph_vector_int_t *vcolor;
  igraph_vector_int_t *ecolor;
  igraph_vector_int_t *etime;

  igraph_llist_ptr_init(&result_graphs_list);
  igraph_llist_ptr_init(&result_vcolors_list);
  igraph_llist_ptr_init(&result_ecolors_list);
  igraph_llist_ptr_init(&result_etimes_list);
  igraph_llist_long_init(&result_supps_list);
  igraph_vector_init(&edges, 0);

  // read and parse first tid
  if (!fgets(buf, 32, instream)) {
    IGRAPH_ERROR("could not read from file", IGRAPH_PARSEERROR);
  }
  if (sscanf(buf, "t # %d %ld", &timestamp, &supp) < 1) {
    IGRAPH_ERROR("invalid file format, missing tid", IGRAPH_PARSEERROR);
  }

  // read first vertex
  if (!fgets(buf, 32, instream)) {
    IGRAPH_ERROR("could not read from file", IGRAPH_PARSEERROR);
  }

  do {
    igraph_llist_init(&edge_list);
    igraph_llist_int_init(&vcolor_list);
    igraph_llist_int_init(&ecolor_list);
    igraph_llist_int_init(&etimes_list);

    // parse vertices
    do {
      i1 = -1; i2 = -1;
      n_fields = sscanf(buf, "v %ld %ld", &i1, &i2);
      if (n_fields < 1) {
	break; // we reached the first edge
      }
      if (i1 > max_vid) {
	max_vid = i1;
      }
      if (has_vcolors) {
	igraph_llist_int_push_back(&vcolor_list, i2);
      }
    } while (fgets(buf, 32, instream));

    // parse edges
    do {
      i1 = -1; i2 = -1; i3 = -1;
      n_fields = sscanf(buf, "e %ld %ld %ld %ld", &i1, &i2, &i3, &i4);
      if (n_fields < 2) {
	break; // we reached the next tid
      }
      igraph_llist_push_back(&edge_list, i1);
      igraph_llist_push_back(&edge_list, i2);
      if (has_ecolors) {
	igraph_llist_int_push_back(&ecolor_list, i3);
	if (has_etimes) {
	  igraph_llist_int_push_back(&etimes_list, i4);
	}
      } else {
	if (has_etimes) {
	  igraph_llist_int_push_back(&etimes_list, i3);
	}
      }
    } while (fgets(buf, 32, instream));

    // allocate graph
    graph = igraph_Calloc(1, igraph_t);
    igraph_empty(graph, max_vid+1, directed);
    igraph_llist_to_vector(&edge_list, &edges, 0);
    igraph_add_edges(graph, &edges, 0);
    igraph_llist_ptr_push_back(&result_graphs_list, graph);
    if (has_vcolors) {
      vcolor = igraph_Calloc(1, igraph_vector_int_t);
      igraph_vector_int_init(vcolor, 0);
      igraph_llist_int_to_vector(&vcolor_list, vcolor, 0);
      igraph_llist_ptr_push_back(&result_vcolors_list, vcolor);
    }
    if (has_ecolors) {
      ecolor = igraph_Calloc(1, igraph_vector_int_t);
      igraph_vector_int_init(ecolor, 0);
      igraph_llist_int_to_vector(&ecolor_list, ecolor, 0);
      igraph_llist_ptr_push_back(&result_ecolors_list, ecolor);
    }
    if (has_etimes) {
      etime = igraph_Calloc(1, igraph_vector_int_t);
      igraph_vector_int_init(etime, 0);
      igraph_llist_int_to_vector(&etimes_list, etime, 0);
      igraph_llist_ptr_push_back(&result_etimes_list, etime);
    }
    if (supps != NULL) {
      igraph_llist_long_push_back(&result_supps_list, supp);
    }
    igraph_llist_destroy(&edge_list);
    igraph_llist_int_destroy(&vcolor_list);
    igraph_llist_int_destroy(&ecolor_list);
    igraph_llist_int_destroy(&etimes_list);

    // parse tid
    if (sscanf(buf, "t # %d %ld", &timestamp, &supp) < 1) {
      break;
    }
    max_vid = -1;
  } while (fgets(buf, 32, instream));

  igraph_llist_ptr_to_vector(&result_graphs_list, graphs, 0);
  if (has_vcolors) {
    igraph_llist_ptr_to_vector(&result_vcolors_list, vcolors, 0);
  }
  if (has_ecolors) {
    igraph_llist_ptr_to_vector(&result_ecolors_list, ecolors, 0);
  }
  if (has_etimes) {
    igraph_llist_ptr_to_vector(&result_etimes_list, etimes, 0);
  }
  if (supps != NULL) {
    igraph_llist_long_to_vector(&result_supps_list, supps, 0);
  }

  igraph_llist_ptr_destroy(&result_graphs_list);
  igraph_llist_ptr_destroy(&result_vcolors_list);
  igraph_llist_ptr_destroy(&result_ecolors_list);
  igraph_llist_ptr_destroy(&result_etimes_list);
  igraph_llist_long_destroy(&result_supps_list);

  return 0;
}

int igraph_read_and_project_transactions_velist(FILE *instream, igraph_bool_t directed,
	igraph_bool_t has_vcolors, igraph_bool_t has_ecolors,
	igraph_integer_t max_vcolor, igraph_integer_t max_ecolor,
	igraph_projection_t proj_type, igraph_integer_t timestep_limit,
	gzFile fgz) {
  char buf[32];
  int timestamp;
  long int i1, i2, i3, n_fields, max_vid = 0;
  long int tid = 0;

  igraph_t *graph1 = NULL;
  igraph_vector_int_t *vcolor1 = NULL;
  igraph_vector_int_t *ecolor1 = NULL;
  igraph_t *graph2 = NULL;
  igraph_vector_int_t *vcolor2 = NULL;
  igraph_vector_int_t *ecolor2 = NULL;

  igraph_llist_t edge_list;
  igraph_llist_int_t vcolor_list;
  igraph_llist_int_t ecolor_list;

  igraph_t *graph = NULL;
  igraph_vector_t edges;
  igraph_vector_int_t *vcolor = NULL;
  igraph_vector_int_t *ecolor = NULL;

  igraph_vector_ptr_t tmp_db, tmp_db_vcolors, tmp_db_ecolors;

  IGRAPH_CHECK(igraph_vector_init(&edges, 0));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db, 2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db_vcolors, 2));
  IGRAPH_CHECK(igraph_vector_ptr_init(&tmp_db_ecolors, 2));

  // read and parse first timestamp
  if (!fgets(buf, 32, instream)) {
    IGRAPH_ERROR("could not read from file", IGRAPH_PARSEERROR);
  }
  if (sscanf(buf, "t # %d", &timestamp) < 1) {
    IGRAPH_ERROR("invalid file format, missing tid", IGRAPH_PARSEERROR);
  }

  // read first vertex
  if (!fgets(buf, 32, instream)) {
    IGRAPH_ERROR("could not read from file", IGRAPH_PARSEERROR);
  }

  do {
    igraph_llist_init(&edge_list);
    igraph_llist_int_init(&vcolor_list);
    igraph_llist_int_init(&ecolor_list);

    // parse vertices
    do {
      i1 = -1; i2 = -1;
      n_fields = sscanf(buf, "v %ld %ld", &i1, &i2);
      if (n_fields < 1) {
	break; // we reached the first edge
      }
      if (i1 > max_vid) {
	max_vid = i1;
      }
      if (has_vcolors) {
	igraph_llist_int_push_back(&vcolor_list, i2);
      }
    } while (fgets(buf, 32, instream));

    // parse edges
    do {
      i1 = -1; i2 = -1; i3 = -1;
      n_fields = sscanf(buf, "e %ld %ld %ld", &i1, &i2, &i3);
      if (n_fields < 2) {
	break; // we reached the next tid
      }
      igraph_llist_push_back(&edge_list, i1);
      igraph_llist_push_back(&edge_list, i2);
      if (has_ecolors) {
	igraph_llist_int_push_back(&ecolor_list, i3);
      }
    } while (fgets(buf, 32, instream));

    // allocate graph
    graph = igraph_Calloc(1, igraph_t);
    igraph_empty(graph, max_vid+1, directed);
    igraph_llist_to_vector(&edge_list, &edges, 0);
    igraph_add_edges(graph, &edges, 0);
    if (has_vcolors) {
      vcolor = igraph_Calloc(1, igraph_vector_int_t);
      igraph_vector_int_init(vcolor, 0);
      igraph_llist_int_to_vector(&vcolor_list, vcolor, 0);
    }
    if (has_ecolors) {
      ecolor = igraph_Calloc(1, igraph_vector_int_t);
      igraph_vector_int_init(ecolor, 0);
      igraph_llist_int_to_vector(&ecolor_list, ecolor, 0);
    }

    // add to length-2 ring buffer
    if (graph1 == NULL) {
      graph1 = graph;
      vcolor1 = vcolor;
      ecolor1 = ecolor;
    } else {
      if (graph2 == NULL) {
	graph2 = graph;
	vcolor2 = vcolor;
	ecolor2 = ecolor;
      } else {
	igraph_destroy(graph1);
	igraph_free(graph1);
	if (has_vcolors) {
	  igraph_vector_int_destroy(vcolor1);
	  igraph_free(vcolor1);
	}
	if (has_ecolors) {
	  igraph_vector_int_destroy(ecolor1);
	  igraph_free(ecolor1);
	}
	graph1 = graph2;
	vcolor1 = vcolor2;
	ecolor1 = ecolor2;
	graph2 = graph;
	vcolor2 = vcolor;
	ecolor2 = ecolor;
      }
    }

    // destroy helpers
    igraph_llist_destroy(&edge_list);
    igraph_llist_int_destroy(&vcolor_list);
    igraph_llist_int_destroy(&ecolor_list);

    // compute projection
    if ((proj_type != IGRAPH_PROJECTION_NONE) && (graph2 != NULL)) {
      VECTOR(tmp_db)[0] = graph1;
      VECTOR(tmp_db)[1] = graph2;
      if (has_vcolors) {
	VECTOR(tmp_db_vcolors)[0] = vcolor1;
	VECTOR(tmp_db_vcolors)[1] = vcolor2;
      }
      if (has_ecolors) {
	VECTOR(tmp_db_ecolors)[0] = ecolor1;
	VECTOR(tmp_db_ecolors)[1] = ecolor2;
      }
      IGRAPH_CHECK(igraph_write_dynamic_union_graph_projection(&tmp_db,
	      (has_vcolors ? &tmp_db_vcolors : NULL),
	      (has_ecolors ? &tmp_db_ecolors : NULL),
	      proj_type,
	      (has_vcolors ? max_vcolor : 0),
	      (has_ecolors ? max_ecolor : 0),
	      fgz, &tid));
    }

    // parse timestamp
    if (sscanf(buf, "t # %d", &timestamp) < 1) {
      break;
    }
    if ((timestep_limit >= 0) && (timestamp >= timestep_limit)) {
      break;
    }
    max_vid = -1;
  } while (fgets(buf, 32, instream));

  if (graph1 != NULL) {
    igraph_destroy(graph1);
    igraph_free(graph1);
    if (has_vcolors) {
      igraph_vector_int_destroy(vcolor1);
      igraph_free(vcolor1);
    }
    if (has_ecolors) {
      igraph_vector_int_destroy(ecolor1);
      igraph_free(ecolor1);
    }
  }
  if (graph2 != NULL) {
    igraph_destroy(graph2);
    igraph_free(graph2);
    if (has_vcolors) {
      igraph_vector_int_destroy(vcolor2);
      igraph_free(vcolor2);
    }
    if (has_ecolors) {
      igraph_vector_int_destroy(ecolor2);
      igraph_free(ecolor2);
    }
  }

  igraph_vector_ptr_destroy(&tmp_db);
  igraph_vector_ptr_destroy(&tmp_db_vcolors);
  igraph_vector_ptr_destroy(&tmp_db_ecolors);

  return 0;
}


int igraph_write_avm(long int N, long int T, int avg_degree,
	    double opinion_prior, double rewiring_p, int initial_graph_generator,
	    FILE *outstream) {
  igraph_t graph;
  igraph_vector_int_t opinions;
  long int i, t, eid, v_resolv, v_other, v_new;
  igraph_integer_t v1, v2;

  srand(time(NULL));

  // initialize graph
  if (initial_graph_generator == 0) {
    IGRAPH_CHECK(igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNM, N, (int)(avg_degree*N/2.),
		/*directed*/ 0, /*loops*/ 0));
  } else {
    IGRAPH_CHECK(igraph_barabasi_game(&graph, N, /*power=*/1.0,
		/*outedges per node=*/(avg_degree/2), /*outseq=*/ NULL,
		/*outpref=*/ 1, /*A=*/ 1.0, /*directed=*/ 0,
		/*algo=*/IGRAPH_BARABASI_PSUMTREE, /*startfrom=*/ NULL));
  }

  // initialize opinions
  IGRAPH_CHECK(igraph_vector_int_init(&opinions, N));
  for (i = 0; i < N; i++) {
    VECTOR(opinions)[i] = 1+((rand()/(double)RAND_MAX)<opinion_prior);
  }

  // write first graph
  fprintf(outstream, "t # 0\n");
  igraph_write_colored_graph(&graph, &opinions, NULL, NULL, outstream);

  // evolve network
  for (t = 0; t < T; t++) {
    // randomly choose an edge that connects nodes with different opinions
    do {
      eid = rand() % igraph_ecount(&graph);
      IGRAPH_CHECK(igraph_edge(&graph, eid, &v1, &v2));
      printf("+++ sampling edge %ld\n", eid);
    } while (VECTOR(opinions)[v1] == VECTOR(opinions)[v2]); // TODO: possibly infinite

    // randomly choose one of the two nodes as the resolver
    if (rand() % 2) {
      v_resolv = v1;
      v_other = v2;
    } else {
      v_resolv = v2;
      v_other = v1;
    }

    // randomly choose an action to take (rewiring or adoption)
    if ((rand()/(double)RAND_MAX)<rewiring_p) {
      // rewiring
      do {
	v_new = rand() % N;
      } while ((v_new == v_resolv) || (VECTOR(opinions)[v_new] != VECTOR(opinions)[v_resolv]));//inf
      printf("time %ld: %ld rewires from %ld to %ld\n", t+1, v_resolv,
		  v_other, v_new);
      IGRAPH_CHECK(igraph_delete_edges(&graph, igraph_ess_1(eid)));
      IGRAPH_CHECK(igraph_add_edge(&graph, v_resolv, v_new));
    } else {
      // adoption
      printf("time %ld: %ld adopts opinion %d from %ld\n", t+1, v_resolv,
		  VECTOR(opinions)[v_other], v_other);
      VECTOR(opinions)[v_resolv] = VECTOR(opinions)[v_other];
    }

    // write graph
    fprintf(outstream, "t # %ld\n", t+1);
    igraph_write_colored_graph(&graph, &opinions, NULL, NULL, outstream);
  }

  return 0;
}

// seperate an evolutionary graph pattern g into it's two timesteps a -> b
void igraph_seperate_graph_pattern(
    const igraph_t *g, 
    const igraph_vector_int_t *vcolors, const igraph_vector_int_t *ecolors,
    int max_vcolor, int max_ecolor, 
    igraph_t *a_g, igraph_t *b_g
    ,igraph_vector_int_t *a_vc, igraph_vector_int_t *b_vc
    ,igraph_vector_int_t *a_ecolors, igraph_vector_int_t *b_ecolors
    ) {
 
  long int i;
  long int N = igraph_vcount(g);
  igraph_integer_t from, to;

  if (ecolors == NULL) {
    printf("ERROR: need edge colors, forgot -c -m MAX_ECOLOR flags?\n");
    return;
  }

  // vertices
  igraph_vector_int_init(a_vc,N);
  igraph_vector_int_init(b_vc,N);

  // if existent, add node-colors
  // otherwise, label all nodes '1'
  if (vcolors != NULL){
    for (i = 0; i < igraph_vcount(g); i++) {
        VECTOR(*a_vc)[i] = VECTOR(*vcolors)[i]/(max_vcolor+1);
        VECTOR(*b_vc)[i] = VECTOR(*vcolors)[i]%(max_vcolor+1);
    }
  } else {
    igraph_real_t one = 1;
    igraph_vector_int_fill(a_vc,one);
    igraph_vector_int_fill(b_vc,one); 
  }
  
  // graphs
  igraph_empty(a_g, N, igraph_is_directed(g));
  igraph_empty(b_g, N, igraph_is_directed(g));

  // edges
  igraph_vector_t a_edges, b_edges;
  igraph_llist_t a_edge_list, b_edge_list;

  igraph_llist_int_t a_ecolors_list;
  igraph_llist_int_t b_ecolors_list;

  igraph_vector_init(&a_edges, 0);
  igraph_vector_init(&b_edges, 0);

  igraph_llist_init(&a_edge_list);
  igraph_llist_init(&b_edge_list);

  igraph_llist_int_init(&a_ecolors_list);
  igraph_llist_int_init(&b_ecolors_list);


  for (i = 0; i < igraph_ecount(g); i++) {
    igraph_edge(g, i, &from, &to);

    if (VECTOR(*ecolors)[i] / (max_ecolor+1) > 0) {
       int a_color = VECTOR(*ecolors)[i]/(max_ecolor+1);
       igraph_llist_int_push_back(&a_ecolors_list, a_color);
       igraph_llist_push_back(&a_edge_list, from);
       igraph_llist_push_back(&a_edge_list, to);
    }


    if (VECTOR(*ecolors)[i] % (max_ecolor+1) > 0) {
      int b_color = VECTOR(*ecolors)[i]%(max_ecolor+1);
      igraph_llist_int_push_back(&b_ecolors_list,b_color);
      igraph_llist_push_back(&b_edge_list, from);
      igraph_llist_push_back(&b_edge_list, to);
    }
  }

  igraph_llist_to_vector(&a_edge_list, &a_edges, 0);
  igraph_llist_to_vector(&b_edge_list, &b_edges, 0);

  // "return"
  igraph_add_edges(a_g, &a_edges, 0);  
  igraph_add_edges(b_g, &b_edges, 0);  

  igraph_vector_int_init(a_ecolors,0);
  igraph_vector_int_init(b_ecolors,0);
  igraph_llist_int_to_vector(&a_ecolors_list,a_ecolors, 0);
  igraph_llist_int_to_vector(&b_ecolors_list,b_ecolors, 0);

}
