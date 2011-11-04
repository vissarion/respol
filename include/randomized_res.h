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

typedef std::set<size_t> cone;
typedef std::set<cone> Cones;
typedef std::set<std::vector<Field> > VPolytope;
  
VPolytope init_VRes(Triangulation& Res){
	typedef Triangulation::Vertex_iterator  VCI;
	
	VPolytope VRes;
	for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
		PPoint_d p = vit->point();
		std::vector<Field> point(p.cartesian_begin(),p.cartesian_end());
		VRes.insert(point);
	}
	return VRes;
}	  
Cones construct_init_cones(Triangulation& Res){
	typedef Triangulation::Full_cell_handle             Simplex;
  typedef std::vector<Simplex> Simplices;

	Cones C;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  Res.incident_full_cells(Res.infinite_vertex(), out);
  size_t finite_cells=0;
  for (std::vector<Simplex>::iterator cit=inf_simplices.begin();
             cit!=inf_simplices.end();
             cit++){
    cone c;
    for (int i=1; i<=Res.current_dimension(); i++){
			if (!Res.is_infinite((*cit)->vertex(i))){
			  c.insert((*cit)->vertex(i)->point().index());
			}
		}
		C.insert(c);
  }
  #ifdef PRINT_INFO
    std::cout << C << std::endl;
  #endif
  return C;
}

// generate a random vector in the cone defined by the vectors
// Pdets[i]-inner_p, for every i in idx
Vector_d generate_random_vector(HD& Pdets,
																cone idx,
																PPoint_d inner_p){
		CGAL::Random rng;
		#ifdef PRINT_INFO
		  //std::cout << PD << " " << (Pdets[0]).end()-(Pdets[0]).begin() << 
		  //         " " << Pdets[0] << std::endl;
		#endif
		Vector_d r(PD,CGAL::NULL_VECTOR);
		for (cone::const_iterator idx_it=idx.begin();
		     idx_it!=idx.end(); ++idx_it){
		  PPoint_d p(PD,Pdets[*idx_it].begin(),Pdets[*idx_it].end());
		  Vector_d v=p-inner_p;
		  v*=Field(rng.get_int(0,100));
		  r+=v;
		  #ifdef PRINT_INFO  
		    //std::cout << v << "-->" << r << std::endl;
		  #endif
		}
  	return r;
}

// compute all Res vertices left by augmenting the simplex
// to the directions of the normal vectors of the facets

int augment_Res_rand(const std::vector<std::vector<Field> >& pointset,
										 const std::vector<int>& mi,
										 int RD,
										 const std::vector<int>& proj,
										 HD& dets,
										 HD& Pdets,
										 Triangulation& Res,
										 const CTriangulation& T,
                     Cones& C){

	#ifdef PRINT_INFO
		std::cout << "\n\nRANDOMIZED AUGMENTING RESULTANT POLYTOPE" << std::endl;
  #endif
  PPoint_d a(PD,Pdets[0].begin(),Pdets[0].end());
	PPoint_d b(PD,Pdets[1].begin(),Pdets[1].end());
	PPoint_d c(PD,Pdets[0].begin(),Pdets[0].end());
	PPoint_d d(PD,Pdets[1].begin(),Pdets[1].end());
	PPoint_d e(PD,Pdets[0].begin(),Pdets[0].end());
	PPoint_d ab = CGAL::midpoint(a,b);
	PPoint_d abc = CGAL::midpoint(ab,c);
	PPoint_d abcd = CGAL::midpoint(abc,d);
	PPoint_d inner_p = CGAL::midpoint(abcd,e);
	
	
	size_t step=0;
  while (!C.empty()){
		++step;
		Cones::iterator cit = C.end();
		--cit;
		cone c = *(cit);
	  C.erase(*cit);
	  Vector_d v = generate_random_vector(Pdets,c,inner_p);
	  #ifdef PRINT_INFO
		std::cout << "\nAUGmenting step " << step << std::endl;
		std::cout << "current normal= " << v << std::endl;
		#endif
	  std::vector<Field> new_vertex =
      compute_res_vertex2(pointset,mi,RD,proj,dets,Pdets,Res,T,v);
  }
  
	return 0;
}


////////////////////////////////////////////////////////////
// the randomized algorithm
std::pair<int,int> 
     RandomizedInnerQ(const std::vector<std::vector<Field> >& pointset,
							        int m,
							        const std::vector<int>& mi,
							        int RD,
							        const std::vector<int>& proj,
							        HD& dets,
							        HD& Pdets,
							        Triangulation& Res){

	//std::cout << "cayley dim:" << CayleyTriangulation(pointset) << std::endl;

  // construct an initial triangulation of the points that will not be projected
  CTriangulation T(CD);
  StaticTriangulation(pointset,proj,T,dets);
	//std::cout << "static dim:" << T.current_dimension() << std::endl;
  
  // start by computing a simplex
  int start_triangs = initialize_Res(pointset,mi,RD,proj,dets,Pdets,Res,T);

  VPolytope VRes = init_VRes(Res);
  std::cout << VRes << std::endl;
  Cones C = construct_init_cones(Res);
  
  // augment simplex to compute the Res polytope
  int augment_triangs = augment_Res_rand(pointset,mi,RD,proj,dets,Pdets,Res,T,C);

  // number of triangulations computed
  std::pair<int,int> num_of_triangs(start_triangs,augment_triangs);

  return num_of_triangs;
}


///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///! DRAFT PAGE
///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// test function 
int run_rand_vec(HD& Pdets){
	PPoint_d a(PD,Pdets[0].begin(),Pdets[0].end());
	PPoint_d b(PD,Pdets[1].begin(),Pdets[1].end());
	
	PPoint_d inner_p = CGAL::midpoint(a,b);
  size_t idx_v[] = {0,2,3,5};
  std::set<size_t> idx (idx_v, idx_v + sizeof(idx_v) / sizeof(size_t) );
  generate_random_vector(Pdets,idx,inner_p);
}



///////////////////////////////////////////////////////////////////////
// uniformly distributed random

#ifdef RANDOM_RES
int random_compute_Res(const std::vector<std::vector<Field> >& pointset,
								 const std::vector<int>& mi,
								 int RD,
								 const std::vector<int>& proj,
								 HD& dets,
								 HD& Pdets,
								 Triangulation& Res,
								 const CTriangulation& T){

	int num_of_triangs=0;
	#ifdef PRINT_INFO
	  std::cout << "dim=" << Res.current_dimension() << std::endl;
  #endif
  // make a stack (stl vector) with normals vectors and initialize
  NV_ds normal_list_d;
  normal_list_d.random_initialize(10000);

  // compute trinagulations using normals as liftings until we  run out of  normal vectors
	std::set<std::vector<Field> > Res_vertices;
	while(!normal_list_d.empty()){
		std::cout << "normal=" << normal_list_d.back() << std::endl;
		std::vector<Field> new_vertex =
      compute_res_vertex(pointset,mi,RD,proj,dets,Pdets,Res,T,normal_list_d);
		Res_vertices.insert(new_vertex);
		num_of_triangs++;
		#ifdef PRINT_INFO
			normal_list_d.print();
			std::cout<< "current number of Res vertices: "
							 << Res_vertices.size()
							 << std::endl;
		#endif
	}
	return num_of_triangs;
}

////////////////////////////////////////////////////////////
// the uniform randomized algorithm
std::pair<int,int> compute_res_rand_uniform(
        const std::vector<std::vector<Field> >& pointset,
        int m,
        const std::vector<int>& mi,
        int RD,
        const std::vector<int>& proj,
        HD& dets,
        HD& Pdets,
        Triangulation& Res){

	//std::cout << "cayley dim:" << CayleyTriangulation(pointset) << std::endl;

  // construct an initial triangulation of the points that will not be projected
  CTriangulation T(CD);
  StaticTriangulation(pointset,proj,T,dets);
	//std::cout << "static dim:" << T.current_dimension() << std::endl;

  int start_triangs = random_compute_Res(pointset,mi,RD,proj9i,dets,Pdets,Res,T);
  std::pair<int,int> num_of_triangs(start_triangs,0);

  return num_of_triangs;
}
#endif
