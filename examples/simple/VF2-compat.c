/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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

#include <igraph.h>
#include <stdlib.h>

/* ----------------------------------------------------------- */

/* Vertices/edges with the same parity match */
igraph_bool_t compat_parity(const igraph_t *graph1, 
			    const igraph_t *graph2,
			    const igraph_integer_t g1_num,
			    const igraph_integer_t g2_num,
			    void *arg) {
  return (g1_num % 2) == (g2_num %2);
}

/* Nothing vertex/edge 0 in graph1 */
igraph_bool_t compat_not0(const igraph_t *graph1, 
			  const igraph_t *graph2,
			  const igraph_integer_t g1_num,
			  const igraph_integer_t g2_num,
			  void *arg) {
  return g1_num != 0;
}

int match_rings() {

  igraph_t r1, r2;
  igraph_bool_t iso;
  igraph_integer_t count;
  igraph_ring(&r1, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
  igraph_ring(&r2, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);

  igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
			&iso, /*map12=*/ 0, /*map21=*/ 0, 
			/*node_compat_fn=*/ 0, /*edge_compat_fn=*/ 0, 
			/*arg=*/ 0);
  if (!iso) { exit(1); }

  igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
			&iso, /*map12=*/ 0, /*map21=*/ 0, 
			compat_parity, /*edge_compat_fn=*/ 0, /*arg=*/ 0);
  if (!iso) { exit(2); }  

  igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
  			&iso, /*map12=*/ 0, /*map21=*/ 0,
  			compat_not0, /*edge_compat_fn=*/ 0, /*arg=*/ 0);
  if (iso) { exit(3); }

  /* ------- */

  igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
			&iso, /*map12=*/ 0, /*map21=*/ 0, 
			/*node_compat_fn=*/ 0, compat_parity, /*arg=*/ 0);
  if (!iso) { exit(4); }  

  igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
			&iso, /*map12=*/ 0, /*map21=*/ 0, 
			/*node_compat_fn=*/ 0, compat_not0, /*arg=*/ 0);
  if (iso) { exit(5); }  

  /* ------- */

  igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
				&count, /*node_compat_fn=*/ 0, 
				/*edge_compat_fn=*/ 0, /*arg=*/ 0);
  
  if (count != 20) { exit(6); }

  igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
				&count, compat_parity, /*edge_compat_fn=*/ 0,
				/*arg=*/ 0);

  if (count != 10) { exit(7); }

  igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
				&count, compat_not0, /*edge_compat_fn=*/ 0,
				/*arg=*/ 0);

  if (count != 0) { exit(8); }  

  /* ------- */ 

  igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
				&count, /*node_compat_fn=*/ 0, compat_parity,
				/*arg=*/ 0);

  if (count != 10) { exit(9); }

  igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ 0, 0, 0, 0,
				&count, /*node_compat_fn=*/ 0, compat_not0,
				/*arg=*/ 0);

  if (count != 0) { exit(10); }    

  igraph_destroy(&r1);
  igraph_destroy(&r2);
  return 0;
}

int match_rings_open_closed() {
  igraph_t ro, rc;
  igraph_bool_t iso;
  igraph_integer_t count;
  igraph_ring(&ro, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 0);
  igraph_ring(&rc, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);  

  igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
			   &iso, /*map12=*/ 0, /*map21=*/ 0, 
			   /*node_compat_fn=*/ 0, /*edge_compat_fn=*/ 0, 
			   /*arg=*/ 0);
  if (!iso) { exit(31); }

  igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
			   &iso, /*map12=*/ 0, /*map21=*/ 0, 
			   compat_parity, /*edge_compat_fn=*/ 0, 
			   /*arg=*/ 0);
  if (!iso) { exit(32); }

  igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
			   &iso, /*map12=*/ 0, /*map21=*/ 0, 
			   compat_not0, /*edge_compat_fn=*/ 0, 
			   /*arg=*/ 0);
  if (iso) { exit(33); }
  
  /* ------- */ 

  igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
			   &iso, /*map12=*/ 0, /*map21=*/ 0, 
			   /*node_compat_fn=*/ 0, compat_parity,
			   /*arg=*/ 0);
  if (!iso) { exit(34); }

  igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
			   &iso, /*map12=*/ 0, /*map21=*/ 0, 
			   /*node_compat_fn=*/ 0, compat_not0,
			   /*arg=*/ 0);
  if (!iso) { exit(35); }  

  /* ------- */

  igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
  				   &count, /*node_compat_fn=*/ 0,
  				   /*edge_compat_fn=*/ 0, /*arg=*/ 0);
  
  if (count != 20) { exit(36); }

  igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
  				   &count, compat_parity,
  				   /*edge_compat_fn=*/ 0, /*arg=*/ 0);

  if (count != 10) { exit(37); }

  igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
  				   &count, compat_not0, /*edge_compat_fn=*/ 0,
  				   /*arg=*/ 0);

  if (count != 0) { exit(38); }

  /* ------- */ 

  igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
				   &count, /*node_compat_fn=*/ 0, 
				   compat_parity, /*arg=*/ 0);
  
  if (count != 10) { exit(39); }

  igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ 0, 0, 0, 0, /* induced= */ 0,
				   &count, /*node_compat_fn=*/ 0, compat_not0,
				   /*arg=*/ 0);

  if (count != 2) { exit(40); }

  igraph_destroy(&ro);
  igraph_destroy(&rc);
  return 0;  
}

int matching_functions() {
  igraph_t g1, g2, g3;
  igraph_vector_t map21, edges1, edges2, edges3;
  igraph_bool_t iso;

  igraph_real_t my_map21[3] = {0, 1, 2};
  igraph_real_t my_edges1[8] = {0,1, 0,2, 0,3, 1,2};
  igraph_real_t my_edges2[4] = {0,1, 0,2};
  igraph_real_t my_edges3[6] = {0,1, 0,2, 1,2};

  igraph_vector_init_copy(&map21, my_map21, 3);
  igraph_vector_init_copy(&edges1, my_edges1, 8);
  igraph_vector_init_copy(&edges2, my_edges2, 4);
  igraph_vector_init_copy(&edges3, my_edges3, 6);

  igraph_create(&g1, &edges1, 4, 0);
  igraph_create(&g2, &edges2, 3, 0);
  igraph_create(&g3, &edges3, 3, 0);

  igraph_is_subisomorphism(0, &map21, &g1, &g2, 0, 0, 0, 0, 0, 0, 0, 0, &iso);
  if (iso) { exit(41); }

  igraph_is_subisomorphism(0, &map21, &g1, &g3, 0, 0, 0, 0, 0, 0, 0, 0, &iso);
  if (!iso) { exit(42); }

  return 0;
}

/* ----------------------------------------------------------- */

int main() {
  matching_functions();
  match_rings();
  match_rings_open_closed();
  return 0;
}
