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

#ifndef IGRAPH_DYNAMIC_H
#define IGRAPH_DYNAMIC_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include <stdio.h>
#include <zlib.h>
#include "igraph_vector_ptr.h"
#include "igraph_fsm.h"

__BEGIN_DECLS

typedef enum igraph_projection_t {
  IGRAPH_PROJECTION_NONE,
  IGRAPH_PROJECTION_FULL,
  IGRAPH_PROJECTION_NEIGHBORS,
  IGRAPH_PROJECTION_EVENT
} igraph_projection_t;

typedef enum igraph_event_type_t {
  IGRAPH_EVENT_NODE_RELABEL,
  IGRAPH_EVENT_EDGE_INSERTION,
  IGRAPH_EVENT_EDGE_DELETION,
  IGRAPH_EVENT_EDGE_RELABEL
} igraph_event_type_t;

typedef struct igraph_event_t {
  igraph_event_type_t type;
  long int v1, v2, label;
} igraph_event_t;


int igraph_write_avm(long int N, long int T, int avg_degree,
	    double opinion_prior, double rewiring_p, int initial_graph_generator,
	    FILE *outstream);

int igraph_read_transactions_velist(FILE *instream, igraph_bool_t directed,
	igraph_bool_t has_vcolors,
	igraph_bool_t has_ecolors,
	igraph_bool_t has_etimes,
	igraph_vector_ptr_t *graphs,
	igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
	igraph_vector_ptr_t *etimes,
	igraph_vector_long_t *supps);
int igraph_read_and_project_transactions_velist(FILE *instream, igraph_bool_t directed,
	igraph_bool_t has_vcolors, igraph_bool_t has_ecolors,
	igraph_integer_t max_vcolor, igraph_integer_t max_ecolor,
	igraph_projection_t proj_type, igraph_integer_t timestep_limit,
	gzFile fgz);

int igraph_read_and_project_dynamic_velist(FILE *instream, igraph_bool_t directed,
      igraph_bool_t has_vcolors, igraph_bool_t has_ecolors, igraph_bool_t has_etimesdel,
      igraph_projection_t proj_type, long int timestep_limit,
      igraph_vector_ptr_t *graphs, igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      gzFile fgz);

int igraph_compute_dynamic_union_graph_projection(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_projection_t proj_type,
      igraph_vector_ptr_t *result, igraph_vector_ptr_t *result_vcolors,
      igraph_vector_ptr_t *result_ecolors,
      igraph_integer_t max_vcolor, igraph_integer_t max_ecolor);

int igraph_write_dynamic_union_graph_projection(igraph_vector_ptr_t *graphs,
      igraph_vector_ptr_t *vcolors, igraph_vector_ptr_t *ecolors,
      igraph_projection_t proj_type,
      igraph_integer_t max_vcolor, igraph_integer_t max_ecolor,
      gzFile fgz, long int *tid);

void igraph_evomine_stream(igraph_t *initial_network, igraph_event_t *stream);

void igraph_seperate_graph_pattern(
    const igraph_t *g, 
    const igraph_vector_int_t *vcolors, const igraph_vector_int_t *ecolors,
    int max_vcolor, int max_ecolor, 
    igraph_t *a_g, igraph_t *b_g
    ,igraph_vector_int_t *a_vc, igraph_vector_int_t *b_vc
    ,igraph_vector_int_t *a_ecolors, igraph_vector_int_t *b_ecolors
    );

__END_DECLS

#endif
