/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   IGraph library AcGM implementation.
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
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

__BEGIN_DECLS

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
# define IN_DEGREE(g, v) ((long int) VECTOR((g).is)[(v)+1]-VECTOR((g).is)[(v)])
#endif

#ifndef OUT_DEGREE
# define OUT_DEGREE(g, v) ((long int) VECTOR((g).os)[(v)+1]-VECTOR((g).os)[(v)])
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


/* -------------------------------------------------- */
/* Frequent subgraph mining                           */
/* -------------------------------------------------- */

typedef int igraph_support_measure_t(const igraph_t *graph1,
				     const igraph_t *graph2,
				     const igraph_vector_int_t *vertex_color1,
				     const igraph_vector_int_t *vertex_color2,
				     const igraph_vector_int_t *edge_color1,
				     const igraph_vector_int_t *edge_color2,
				     igraph_bool_t induced,
				     igraph_integer_t *support,
				     igraph_integer_t min_supp);

int igraph_shallow_support(const igraph_t *graph1,
			   const igraph_t *graph2,
			   const igraph_vector_int_t *vertex_color1,
			   const igraph_vector_int_t *vertex_color2,
			   const igraph_vector_int_t *edge_color1,
			   const igraph_vector_int_t *edge_color2,
			   igraph_bool_t induced,
			   igraph_integer_t *support,
			   igraph_integer_t min_supp);

int igraph_mib_support(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_bool_t induced,
		       igraph_integer_t *support,
		       igraph_integer_t min_supp);

int igraph_mib_support_slow(const igraph_t *graph1,
		       const igraph_t *graph2,
		       const igraph_vector_int_t *vertex_color1,
		       const igraph_vector_int_t *vertex_color2,
		       const igraph_vector_int_t *edge_color1,
		       const igraph_vector_int_t *edge_color2,
		       igraph_bool_t induced,
		       igraph_integer_t *support,
		       igraph_integer_t min_supp);

int igraph_acgm(const igraph_vector_ptr_t *graphdb, igraph_support_measure_t *supp_fn,
		igraph_integer_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values);

int igraph_gspan(const igraph_vector_ptr_t *graphdb, igraph_support_measure_t *supp_fn,
		igraph_integer_t min_supp, igraph_vector_ptr_t *frequent_subgraphs,
		igraph_vector_t *support_values);

__END_DECLS

#endif
