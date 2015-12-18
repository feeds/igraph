/* vim:set ts=8 sw=2 sts=2 noet:  */
/* 
   Linked list data type
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

#ifndef IGRAPH_LIST_H
#define IGRAPH_LIST_H

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
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Plain linked list                                  */
/* -------------------------------------------------- */

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_LONG
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_INT
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_CHAR
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_PTR
#include "igraph_pmt.h"
#include "igraph_list_pmt.h"
#include "igraph_pmt_off.h"
#undef BASE_PTR

int igraph_llist_int_to_vector_real(igraph_llist_int_t *llist, igraph_vector_t *vector);

__END_DECLS

#endif
