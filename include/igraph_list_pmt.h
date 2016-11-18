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

/**
 * Linked list data type.
 * \ingroup internal
 */

typedef struct TYPE(igraph_llist_item) {
  BASE data;
  struct TYPE(igraph_llist_item) *next;
  struct TYPE(igraph_llist_item) *prev;
} TYPE(igraph_llist_item);

typedef struct TYPE(igraph_llist) {
  TYPE(igraph_llist_item) *first;
  TYPE(igraph_llist_item) *last;
  long int len;
} TYPE(igraph_llist);

int FUNCTION(igraph_llist,init)(TYPE(igraph_llist) *llist);
void FUNCTION(igraph_llist,destroy)(TYPE(igraph_llist) *llist);
void FUNCTION(igraph_llist,print)(TYPE(igraph_llist) *llist);
int FUNCTION(igraph_llist,push_back)(TYPE(igraph_llist) *llist, BASE data);
long int FUNCTION(igraph_llist,size)(TYPE(igraph_llist) *llist);
BASE FUNCTION(igraph_llist,back)(TYPE(igraph_llist) *llist);
BASE FUNCTION(igraph_llist,first)(TYPE(igraph_llist) *llist);
int FUNCTION(igraph_llist,to_vector)(TYPE(igraph_llist) *llist, TYPE(igraph_vector) *vector,
      int skip);
int FUNCTION(igraph_llist, insert_at)(TYPE(igraph_llist) * llist, BASE data, long int position);
BASE FUNCTION(igraph_llist, item_at)(TYPE(igraph_llist) * llist, long int position);

