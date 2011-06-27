#include <fstream>
#include <iostream>
#include <list>
#include <cassert>
#include <ctime>

const int D = 3;      	 	// this is the dimension of the supports
const int CD = 2*D+1;	   	// this is the Cayley space + 1 for lifting
const int PD = D+1;				// this is the dimension of the projection

//#define PRINT_INFO
//#include <../include/cgal_chd.h>
#include <../include/cgal_chd_hornus.h>

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
	// you don't have to homogenize!
  HD dets(pointset.begin(),pointset.end());
	//
	HD Pdets;
	
	// define the projection
	// attention! proj SHOULD BE SORTED
	vector<int> proj = proj_first_coord(PD,m,mi);
	
	// the data structure to hold the res polytope
	int numof_triangs=0, numof_init_Res_vertices;
	Triangulation CH(PD);
	
	//compute the res polytope
	//compute_res(pointset,points_index,m,mi,proj,dets,numof_triangs, numof_init_Res_vertices,CH);
	compute_res_faster(pointset,points_index,m,mi,proj,dets,Pdets,numof_triangs, numof_init_Res_vertices,CH);
	
	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
	
	//print_res_vertices(CH);
	
	// print some statistics
	print_statistics(numof_triangs, numof_init_Res_vertices, CH.number_of_vertices(), tstopall-tstartall);
	
	return 0;
}
