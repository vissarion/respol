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
#include <cgal_chd_hornus_cellinfo.h>
#include <randomized_res.h>
#include <simplices_augmentation.h>
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
 	// read input (pointset, mi, n), apply cayley trick, define projection
 	read_pointset(pointset, mi, proj, n);
	int initial_pointset_size = pointset.size();
	
	// remove spesialized redundant (non-extreme) points
	#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY
	compute_extreme_points(pointset,mi,proj);
	#endif
  
  // compute the cayley points set
	cayley_trick(pointset, mi);
	
	// this is the dimension of the resultant (and secondary) polytope
 	int RD = n - 2*D -1;
  
	// compute the big matrix
	// you don't have to homogenize!
  HD dets(pointset.begin(),pointset.end());
  //dets.print_matrix(std::cout);
  
	// this is the big matrix for dimension PD it is empty at
	// the beginning and we add the points when they are computed
	HD Pdets;

		// the data structure to hold the res polytope
	int numof_triangs=0, numof_init_Res_vertices;
  Triangulation Res(PD);
  
  //////////////////////////////////////////////////////////////////////
	//COMPUTE THE RES POLYTOPE
	
	std::pair<int,int> num_of_triangs =
    compute_res(pointset,n,mi,RD,proj,dets,Pdets,Res);
  
  //std::pair<int,int> num_of_triangs =
  //InnerQwithsimplices(pointset,n,mi,RD,proj,dets,Pdets,Res);
  
	//std::pair<int,int> num_of_triangs =
  //RandomizedInnerQ(pointset,n,mi,RD,proj,dets,Pdets,Res);
  
  //std::pair<int,int> num_of_triangs =
  //compute_res_rand_uniform(pointset,n,mi,RD,proj,dets,Pdets,Res,100);
  
  //////////////////////////////////////////////////////////////////////
  
	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
  
  double recompute_time = compute_Res_offline(Pdets,Res);
  
	// print the result i.e. the proj of the Resultant polytope
 	//#ifdef PRINT_INFO	
	//print_res_vertices(Res,std::cout);
  //#endif
	//print_res_facets_number(Res);
  generate_polymake_scripts(Res);
  
  // print some statistics
/*  #ifdef PRINT_INFO
  print_statistics (Res.current_dimension(),
                   num_of_triangs.first, 
                   num_of_triangs.second,
                   Res.number_of_vertices(), 
#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY                   
                   count_extreme_vertices(Res),
#else
                   -1,
#endif                   
                   tstopall-tstartall, // overall time
                   dets.get_determinant_time()+
                   Pdets.get_determinant_time(), // determinant time
                   volume(Res,Pdets));
  #else
*/
  print_statistics_small(CD-1, 
                         PD,
                         Res.current_dimension(),
                         initial_pointset_size,
                         pointset.size(),
                         num_of_triangs.first+num_of_triangs.second,
                         Res.number_of_vertices(),
#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY                   
                         count_extreme_vertices(Res),
#else
                         -1,
#endif
                         tstopall-tstartall, // overall time
                         conv_time, // Res convex hull time
                         recompute_time, // Res convex hull offline time
                         dets.get_determinant_time()+
                         Pdets.get_determinant_time(), // determinant time
                         volume(Res,Pdets));
  //#endif
  
  //std::cout << "convex hull time = " << conv_time << std::endl;
  Pdets.print_matrix(std::cout);
  
  #ifdef PRINT_INFO
    //recompute_Res(Res);
  #endif
  
  #ifdef USE_LRSLIB
  ch();
  #endif
  
  return 0;
}
// vim: ts=2:expandtab
