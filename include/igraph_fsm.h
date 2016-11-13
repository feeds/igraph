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

#ifndef IGRAPH_FSM_H
#define IGRAPH_FSM_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_vector_ptr.h"
#include "igraph_topology.h"

#include <zlib.h>

__BEGIN_DECLS

typedef enum igraph_gspan_variant_t {
  IGRAPH_GSPAN_DEFAULT, // works on databases of graphs with or without node/edge labels
  IGRAPH_GSPAN_EVOMINE, // assumes the database contains union graphs (with label strings)
  IGRAPH_GSPAN_GERM, // edge labels = timestamps, vertex labels supported
  IGRAPH_GSPAN_LFRMINER // edge labels = timestamps, no vertex labels supported
} igraph_gspan_variant_t;

typedef struct igraph_lfrminer_data_t {
  void *unused;
} igraph_lfrminer_data_t;

typedef struct igraph_germ_data_t {
  long int max_rel_timestamp;
} igraph_germ_data_t;

typedef struct igraph_evomine_data_t {
  long int max_vcolor;
  long int max_ecolor;
} igraph_evomine_data_t;

// HELPERS

void igraph_print(const igraph_t *g, const igraph_vector_int_t *vcolors,
		    const igraph_vector_int_t *ecolors, const igraph_vector_int_t *etimes);

void igraph_print_stats(const igraph_t *g);

int igraph_write_colored_graph(igraph_t *g, igraph_vector_int_t *vcolors,
      igraph_vector_int_t *ecolors, igraph_vector_int_t *etimes, FILE *f);

int igraph_write_colored_graph_gz(igraph_t *g, igraph_vector_int_t *vcolors,
      igraph_vector_int_t *ecolors, igraph_vector_int_t *etimes, gzFile f);

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
			   igraph_bool_t *iso);

int igraph_isomorphic_evomine(const igraph_t *graph1, const igraph_t *graph2,
         const igraph_vector_int_t *vertex_color1,
         const igraph_vector_int_t *vertex_color2,
         const igraph_vector_int_t *edge_color1,
         const igraph_vector_int_t *edge_color2,
         igraph_bool_t *iso);

// O(1) access to in- and out-neighbors and degrees

#ifndef OUT_NEIGHBOR
# define OUT_NEIGHBOR(g, v, i) ((long int) VECTOR((g).to)[(long int) VECTOR((g).oi)[(long int) VECTOR((g).os)[(v)] + i]])
#endif

#ifndef IN_NEIGHBOR
# define IN_NEIGHBOR(g, v, i) ((long int) VECTOR((g).from)[(long int) VECTOR((g).ii)[(long int) VECTOR((g).is)[(v)] + i]])
#endif	

#ifndef NEIGHBOR
# define NEIGHBOR(g, v, i) (((i) < VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)]) ? OUT_NEIGHBOR((g),(v),(i)) : IN_NEIGHBOR((g),(v),((i)-(long int)(VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)]))))
#endif

#ifndef IN_DEGREE
# define IN_DEGREE(g, v) ((long int) (VECTOR((g).is)[(v)+1]-VECTOR((g).is)[(v)]))
#endif

#ifndef OUT_DEGREE
# define OUT_DEGREE(g, v) ((long int) (VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)]))
#endif

#ifndef DEGREE
# define DEGREE(g, v) (IN_DEGREE((g), (v)) + OUT_DEGREE((g), (v)))
#endif

#ifndef IN_NEIGH_TO_EID
# define IN_NEIGH_TO_EID(g, v, i) ((long int) VECTOR((g).ii)[(long int) VECTOR((g).is)[(v)] + i])
#endif

#ifndef OUT_NEIGH_TO_EID
# define OUT_NEIGH_TO_EID(g, v, i) ((long int) VECTOR((g).oi)[(long int) VECTOR((g).os)[(v)] + i])
#endif

#ifndef NEIGH_TO_EID
# define NEIGH_TO_EID(g, v, i) (((i) < VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)]) ? OUT_NEIGH_TO_EID((g),(v),(i)) : IN_NEIGH_TO_EID((g),(v),((i)-(long int)(VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)]))))
#endif


/* -------------------------------------------------- */
/* Frequent subgraph mining                           */
/* -------------------------------------------------- */

// Single graph support measures

typedef int igraph_support_measure_t(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamps1,
			   const igraph_vector_int_t *edge_timestamps2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *support,
			   igraph_integer_t min_supp);

int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamps1,
			   const igraph_vector_int_t *edge_timestamps2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *support,
			   igraph_integer_t min_supp);

int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       const igraph_vector_int_t *edge_timestamps1,
		       const igraph_vector_int_t *edge_timestamps2,
		       igraph_bool_t induced,
		       igraph_gspan_variant_t variant,
		       void *variant_data,
		       long int *support,
		       igraph_integer_t min_supp);

int igraph_egobased_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   const igraph_vector_int_t *edge_timestamps1,
			   const igraph_vector_int_t *edge_timestamps2,
			   igraph_bool_t induced,
			   igraph_gspan_variant_t variant,
			   void *variant_data,
			   long int *support,
			   igraph_integer_t min_supp);

// Graph database support measures

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
			  igraph_integer_t min_supp);

// Mining algorithms

int igraph_gspan(const igraph_vector_ptr_t *graphs, const igraph_vector_ptr_t *vertex_colors,
		const igraph_vector_ptr_t *edge_colors, const igraph_vector_ptr_t *edge_times,
		igraph_support_measure_t *single_graph_support,
		igraph_integer_t min_supp, igraph_integer_t max_edges,
		igraph_gspan_variant_t variant,
		igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_ptr_t *frequent_subgraph_vcolors,
		igraph_vector_ptr_t *frequent_subgraph_ecolors,
		igraph_vector_ptr_t *frequent_subgraph_etimes,
		igraph_vector_long_t *frequent_subgraph_supps);

__END_DECLS

#endif
