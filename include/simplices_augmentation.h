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

// compute simplices that are innerQ until it compute all Q vertices

int augment_Res_simplices(const std::vector<std::vector<Field> >& pointset,
										 const std::vector<int>& mi,
										 int RD,
										 const std::vector<int>& proj,
										 HD& dets,
										 HD& Pdets,
										 Triangulation& Res,
										 const CTriangulation& T
										 ){

	#ifdef PRINT_INFO
		std::cout << "\n\nAUGMENTING RESULTANT POLYTOPE with Simplices" << std::endl;
  #endif
  typedef Triangulation::Full_cell_handle             Simplex;
  typedef std::vector<Simplex>                        Simplices;
  size_t step=0;
  
  Triangulation Res_temp(Res);
  Triangulation Res2(PD);
  
  	Simplices inf_simplices;
    std::back_insert_iterator<Simplices> out(inf_simplices);
    Res_temp.incident_full_cells(Res_temp.infinite_vertex(), out);
		std::vector<std::vector<Field> > new_vertices;
		
		for (Simplices::iterator sit=inf_simplices.begin();
			  										 sit!=inf_simplices.end();++sit){
			std::vector<PPoint_d> facet_points;
			for (int i=0; i<=Res_temp.current_dimension(); i++){
				if(!Res_temp.is_infinite((*sit)->vertex(i))){
				  facet_points.push_back((*sit)->vertex(i)->point());
				  //std::cout << (*sit)->vertex(i)->point() << " ";
				}
			}
			//std::cout << std::endl;
			PPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
			// compute a hyperplane which has in its negative side the opposite point
			PHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
			PVector_d current_vector = hp.orthogonal_direction();
			
			#ifdef PRINT_INFO
			std::cout << "\nAUGmenting step " << ++step << std::endl;
			std::cout << "number of vertices= " << Pdets.size() << std::endl;
			#endif			
			std::vector<Field> new_vertex =
	      compute_res_vertex2(pointset,mi,RD,proj,dets,Pdets,Res_temp,T,current_vector);
	    new_vertices.push_back(new_vertex);
	    int idx = Pdets.find(new_vertex);
	    PPoint_d p(PD,new_vertex.begin(),new_vertex.end());
			p.set_hash(&Pdets);
	    if (idx == -1){
			  Pdets.add_column(new_vertex);  
		    p.set_index(Pdets.size());
		    
		  } else {
		    p.set_index(idx);
		  }
		  Res2.insert(p);
	  }
}

////////////////////////////////////////////////////////////
// the simplices algorithm
std::pair<int,int> 
     InnerQwithsimplices(const std::vector<std::vector<Field> >& pointset,
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
  
  // augment simplex to compute the Res polytope
  int augment_triangs = augment_Res_simplices(pointset,mi,RD,proj,dets,Pdets,Res,T);

  // number of triangulations computed
  std::pair<int,int> num_of_triangs(start_triangs,augment_triangs);
  
  return num_of_triangs;
}

