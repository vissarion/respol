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

#include "respol_config.h"
#include <iostream>
#include <fstream>
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
#include <res_enum_functions.h>
//#include <randomized_res.h>
//#include <simplices_augmentation.h>
//#include <../include/cgal_chd_hornus_cellinfo_placing.h>
//#include <../include/cgal_chd_hornus_cellinfo_with_cgal_det.h>

//////////////////////////////////////////////////////////////////
// main

int main(const int argc,const char** argv){

  // Parse command-line options.
  //if (argc==1){
//		std::cerr<<"Use -h option to see how to use respol."<<std::endl;
//		exit(-1);
//	}		
  ResPol::config c;
  c.verbose=0;
  c.read_from_file=false;
  c.output_f_vector=false;
  c.polytope_type=0; // resultant polytope
  for(int i=1;i<argc;++i){
    bool correct=false;
    if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
      std::cerr<<
        "-f, --f-vector\t\toutput the f-vector (require polymake)\n"<<
        "-h, --help\t\tshow this message\n"<<
        "-i, --input file\tset input file (read input from stdin otherwise)\n"<<
        "-r, --resultant\t\tcompute the resultant polytope (default)\n"<<
        "-d, --discriminant\t\tcompute the discriminant polytope (require tropli)\n"<<
        "-s, --secondary\t\tcompute the secondary polytope\n"<<
        "-v, --verbose n\t\tset verbosity level to 0 (default), 1 , 2, 3,"<<
        " 4(returns the vertices of computed polytope)\n";
      exit(-1);
    }
    if(!strcmp(argv[i],"-f")||!strcmp(argv[i],"--f-vector")){
      c.output_f_vector=true;
      correct=true;
    }
    if(!strcmp(argv[i],"-i")||!strcmp(argv[i],"--input")){
      c.read_from_file=true;
      c.inp.open(argv[++i],std::ifstream::in);
      correct=true;
    }
    if(!strcmp(argv[i],"-r")||!strcmp(argv[i],"--resultant")){
      c.polytope_type=0;
      correct=true;
    }
    if(!strcmp(argv[i],"-s")||!strcmp(argv[i],"--secondary")){
      c.polytope_type=1;
      correct=true;
    }
    if(!strcmp(argv[i],"-d")||!strcmp(argv[i],"--discriminant")){
      c.polytope_type=2;
      correct=true;
    }
    if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")){
      c.verbose=atoi(argv[++i]);
      correct=true;
    }
    if(correct==false){
      std::cerr<<"unknown parameter \'"<<argv[i]<<
        "\', try "<<argv[0]<<" --help"<<std::endl;
      exit(-2);
    }
  }

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
 	if(c.read_from_file)
    read_pointset(c.inp, pointset, mi, proj, n);
  else
    read_pointset(std::cin, pointset, mi, proj, n);
	int initial_pointset_size = pointset.size();
	
	// remove spesialized redundant (non-extreme) points
	#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY
	compute_extreme_points(pointset,mi,proj,c);
	#endif
  
  // compute the cayley points set
	cayley_trick(pointset, mi);
	int RD = n - 2*D -1;
	
	/*
	// this is the dimension of the resultant (or secondary or discriminant) polytope
 	int RD;
  if(c.polytope_type==0){
    // compute the cayley points set for the resultant case
          cayley_trick(pointset, mi);
          // this is the dimension of the resultant
          RD = n - 2*D -1;
  }else{
          CD = D+1;
          // this is the dimension of the secondary and discriminant (no Cayley trick)
          RD = n - D -1;
  }
  */
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
    compute_res(pointset,n,mi,RD,proj,dets,Pdets,Res,c);
  
  //RandomizedInnerQ(pointset,n,mi,RD,proj,dets,Pdets,Res);
  
  //std::pair<int,int> num_of_triangs =
  //InnerQwithsimplices(pointset,n,mi,RD,proj,dets,Pdets,Res,c);
  
	//std::pair<int,int> num_of_triangs =
  //RandomizedInnerQ(pointset,n,mi,RD,proj,dets,Pdets,Res);
  
  //std::pair<int,int> num_of_triangs =
  //compute_res_rand_uniform(pointset,n,mi,RD,proj,dets,Pdets,Res,100);
  
  //////////////////////////////////////////////////////////////////////
  
	// stop clocking
	tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
  
  //double recompute_time = compute_Res_offline(Pdets,Res);
  double recompute_time =-1;//= recompute_Res(Res);
  
	// print the result i.e. the proj of the Resultant polytope
 	//#ifdef PRINT_INFO	
	//print_res_vertices(Res,std::cout);
  //#endif
	//print_res_facets_number(Res);
  
  std::ofstream outfile;
  outfile.open("res_vertices.txt");
  print_output(Res,outfile);
  
  // print some statistics
  
  switch (c.verbose) {
  case 0:
      // for test-suite
      std::cout <<  Res.number_of_vertices() << std::endl;
      break;
      break;
  case 1:
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
                             volume(Res,Pdets),
                             Res,
                             c);
      break;
  case 2:
      pretty_print_statistics(CD-1,
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
                              volume(Res,Pdets),
                              Res,
                              c);
      break;
  default:
      break;

  }

  if (c.verbose != 0) {
      // print the vertices of the polytope
      std::cout<< "\nThe vertices of the " ;
      switch(c.polytope_type){
      case 0: std::cout << " resultant "; break;
      case 1: std::cout << " secondary "; break;
      case 2: std::cout << " discriminant "; break;
      }
      std::cout<< "polytope: " << std::endl;
      print_res_vertices(Res,std::cout);
  }
  if(c.output_f_vector){
		generate_polymake_scripts(Res);
    std::ofstream polymakefile;
    polymakefile.open("f_vector.polymake");
    print_polymake_fvector(Res,polymakefile);
    int res = system ("polymake f_vector.polymake");
    std::cout<<std::endl;
  }
  /*
  int cells, triang_facets, facets, edges, vertices;
  f_vector(Res,cells, triang_facets, facets, edges, vertices,c.verbose);
  std::cout << vertices << " " << facets << std::endl;
  */
  #ifdef USE_LRSLIB
  // TODO: call LRS functions, see test_lrs.cpp for an example
  #endif
    
  return 0;
}
// vim: ts=2:expandtab
