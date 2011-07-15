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

#include <fstream>
#include <iostream>
#include <list>
#include <cassert>
#include <ctime>

// these have no meaning;  cayley_trick function will change them
// according to the input
int D;      	 	// this is the dimension of the supports
int CD;	   	// this is the Cayley space + 1 for lifting
int PD;				//this is the dimension of the projection
//int vec[]={0,1,2,10,11,};
//int PD = 6;//sizeof(vec) / sizeof(int);				

//#include <../include/cgal_chd.h>
//#include <../include/cgal_chd_hornus.h>
//#include <../include/cgal_chd_hornus_with_cgal_det.h>
#include <../include/cgal_chd_hornus_cellinfo.h>
//#include <../include/cgal_chd_hornus_cellinfo_with_cgal_det.h>

//////////////////////////////////////////////////////////////////
// main

int main(const int argc, const char** argv) {
	double tstart1, tstop1, tstart2, tstop2, tstartall, tstopall;
	
	// start clocking
	tstartall = (double)clock()/(double)CLOCKS_PER_SEC;
	
	// mi will be the cardinalities of the support sets
	// n is sum{mi}
 	// construct an empty pointset
	std::vector<int> mi;
 	int n;
 	std::vector<std::vector<Field> > pointset;

	// initialize all the above
 	// read input (pointset, mi, n), apply cayley trick 
 	cayley_trick(pointset, mi, n);
	
	// this is the dimension of the resultant (and secondary) polytope 
 	int RD = n - 2*D -1;
 		
	// compute the big matrix
	// you don't have to homogenize!
  HD dets(pointset.begin(),pointset.end());
  
	// this is the big matrix for dimension PD it is empty at 
	// the beginning and we add the points when they are computed
	HD Pdets;

	// define the projection
	// attention! proj SHOULD BE SORTED
	vector<int> proj = proj_first_coord(PD,n,mi);
	//vector<int> proj = proj_more_coord(PD,n,mi);
  //std::cout << proj << std::endl;
	//vector<int> proj = full_proj(PD,n,mi);
	//vector<int> proj (vec, vec + sizeof(vec) / sizeof(int) );
	
	// the data structure to hold the res polytope
	int numof_triangs=0, numof_init_Res_vertices;
	Triangulation Res(PD);
	
	//compute the res polytope
	pair<int,int> num_of_triangs = compute_res_faster(pointset,n,mi,RD,proj,dets,Pdets,Res);
	
	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
	
	// print the result i.e. the proj of the Resultant polytope 
	print_res_vertices(Res);
	
	//print_res_facets_number(Res);
	
	// print some statistics
	//print_statistics(num_of_triangs.first, 
	//                 num_of_triangs.second,
	//                 Res.number_of_vertices(), 
	//                tstopall-tstartall,
	//                 volume(Res,Pdets));
	
	print_statistics_small(CD-1, 
												 PD,
												 pointset.size(),
												 num_of_triangs.first+num_of_triangs.second,
												 Res.number_of_vertices(), 
												 tstopall-tstartall,
	//											 volume(Res,Pdets));
	         0);
	//Pdets.print_matrix(cout);
	
	return 0;
}
