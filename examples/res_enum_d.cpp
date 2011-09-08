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

#include <cassert>
#include <ctime>

// these values have no meaning;  cayley_trick function will change them
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
//#include <../include/cgal_chd_hornus_cellinfo_placing.h>
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
	std::vector<int> mi,proj;
	int n;
 	std::vector<std::vector<Field> > pointset;

	// initialize all the above
 	// read input (pointset, mi, n), apply cayley trick
 	cayley_trick(pointset, mi, proj, n);

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
	//std::vector<int> proj = proj_first_coord(PD,n,mi);
	//std::vector<int> proj = proj_more_coord(PD,n,mi);
  //std::cout << proj << std::endl;
	//std::vector<int> proj = full_proj(PD,n,mi);
	//std::vector<int> proj (vec, vec + sizeof(vec) / sizeof(int) );

	// the data structure to hold the res polytope
	int numof_triangs=0, numof_init_Res_vertices;
	Triangulation Res(PD);

	//compute the res polytope
	std::pair<int,int> num_of_triangs =
    compute_res_faster(pointset,n,mi,RD,proj,dets,Pdets,Res);

	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;

	// print the result i.e. the proj of the Resultant polytope
	print_res_vertices(Res,std::cout);
 
	//print_res_facets_number(Res);

  // print some statistics
  #ifdef PRINT_INFO
  print_statistics (num_of_triangs.first, 
                   num_of_triangs.second,
                   Res.number_of_vertices(), 
                   compute_extreme_res_vertices_maple(Res), 
                   tstopall-tstartall, // overall time
                   dets.get_determinant_time()+
                   Pdets.get_determinant_time(), // determinant time
                   volume(Res,Pdets));
  #else
  print_statistics_small(CD-1, 
                         PD,
                         pointset.size(),
                         num_of_triangs.first+num_of_triangs.second,
                         Res.number_of_vertices(),
                         compute_extreme_res_vertices_maple(Res), 
                         tstopall-tstartall, // overall time
                         dets.get_determinant_time()+
                         Pdets.get_determinant_time(), // determinant time
                         volume(Res,Pdets));
  #endif
  //Pdets.print_matrix(std::cout);

  return 0;
}
// vim: ts=2:expandtab
