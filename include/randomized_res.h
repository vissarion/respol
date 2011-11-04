// Copyright 2011 National and Kapodistrian University of Athens, Greece.
//
// This file is part of respol.
//
// Respol is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// Respol is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with respol, see
// <http://www.gnu.org/licenses/>.

#include <CGAL/constructions_d.h>

// generate a random vector in the cone defined by the vectors
// Pdets[i]-inner_p, for every i in idx
Vector_d generate_random_vector(HD& Pdets,
																std::vector<size_t> idx,
																PPoint_d inner_p){
		CGAL::Random rng;
		std::cout << PD << " " << (Pdets[0]).end()-(Pdets[0]).begin() << 
		           " " << Pdets[0] << std::endl;
		
		Vector_d r(PD,CGAL::NULL_VECTOR);
		for (std::vector<size_t>::const_iterator idx_it=idx.begin();
		     idx_it!=idx.end(); ++idx_it){
		  PPoint_d p(PD,Pdets[*idx_it].begin(),Pdets[*idx_it].end());
		  Vector_d v=p-inner_p;
		  v*=Field(rng.get_int(0,100));
		  r+=v;
		  std::cout << v << "-->" << r << std::endl;
		}
  	return r;
}

// 
int run_rand_vec(HD& Pdets){
	PPoint_d a(PD,Pdets[0].begin(),Pdets[0].end());
	PPoint_d b(PD,Pdets[1].begin(),Pdets[1].end());
	
	PPoint_d inner_p = CGAL::midpoint(a,b);
  size_t idx_v[] = {0,2,3,5};
  std::vector<size_t> idx (idx_v, idx_v + sizeof(idx_v) / sizeof(size_t) );
  generate_random_vector(Pdets,idx,inner_p);
  
}
