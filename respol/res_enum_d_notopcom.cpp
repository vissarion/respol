// Copyright 2011-2012 National and Kapodistrian University of Athens,
// Greece.
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

const int D = 3;      	 	// this is the dimension of the supports
const int CD = 2*D+1;	   	// this is the Cayley space + 1 for lifting
const int PD = D+1;				// this is the dimension of the projection

#include <../include/cgal_chd_notopcom.h>


//////////////////////////////////////////////////////////////////
// main

int main(const int argc, const char** argv) {
	double tstart1, tstop1, tstart2, tstop2, tstartall, tstopall;
	
	// start clocking
	tstartall = (double)clock()/(double)CLOCKS_PER_SEC;
	
	// mi will be the cards of the support sets and m is the total card
  std::vector<int> mi;
 	int m;
 	
 	// construct an empty pointset and an index
 	std::vector<std::vector<Field> > pointset;
 	map<std::vector<Field>,int> points_index;
 	
 	//read input (points, mi, m), apply cayley trick 
 	// (now you have the pointset), make the index
 	cayley_trick(pointset, points_index, mi, m);
	
	// compute the big matrix
	// BUT first homogenize the pointset, dirty way..
	vector<vector<Field> > homo_pointset;
	homogenize(pointset,homo_pointset);
	std::cout << homo_pointset << std::endl;
  HD dets(homo_pointset.begin(),homo_pointset.end());
	
	//define the projection
	vector<int> proj = proj_first_coord(PD,m,mi);
	
	// the data structure to hold the res polytope
	int numof_triangs=0, numof_init_Res_vertices;
	Convex_hull_d CH(PD);
	
	//compute the res polytope
	compute_res(pointset,points_index,m,mi,proj,dets,numof_triangs, numof_init_Res_vertices,CH);
	//compute_res_fast(pointset,points_index,m,mi,proj,dets,numof_triangs, numof_init_Res_vertices,CH);
	
	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
	
	// print the vertices of the res polytope
	for (Vertex_iterator_d vit = CH.vertices_begin(); vit != CH.vertices_end(); vit++)
		std::cout << vit->point() << " ";
	std::cout << std::endl;
	
	// print some statistics
	print_statistics(numof_triangs, numof_init_Res_vertices, CH.number_of_vertices(), tstopall-tstartall);
	
	return 0;
}
