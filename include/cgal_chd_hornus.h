#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm> 

// include d-dimensional things from CGAL
//#include <CGAL/Filtered_kernel_d.h>
//#include <CGAL/Pure_complex.h>
//#include <CGAL/Delaunay_complex.h>
//#include <CGAL/Delaunay_triangulation_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/Triangulation_vertex.h>
#include <CGAL/Triangulation_full_cell.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>

std::ostream& operator<<(std::ostream& ost, const std::vector<size_t>& V);

// for fast determinant computation
#include <../include/fast_hashed_determinant.h>
#include <../include/normal_vector_ds.h>

// for indexed points
//#include <../include/indexed_point.h>

using namespace std;

// the kernel
#include <gmpxx.h>
typedef mpq_class Field;

// Cayley space kernel for the CH of the lifted points
typedef CGAL::Cartesian_d<Field> 									CK;
//typedef Indexed_Cartesian_d<Field> 							CK;
//typedef CK::IndexedPoint_d											IndexedPoint_d;

//typedef CGAL::Triangulation_vertex<CK,size_t>			tv;
//typedef CGAL::Triangulation_full_cell<CK>					tc;
//typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<CD>,tv,tc > tds;
// not working if we'll not pas all the parameters explicitly
//e.g. CGAL::Triangulation_data_structure<CGAL::Dimension_tag<CD> > tds; gives error
//typedef CGAL::Triangulation<CK,tds> 							CTriangulation;
typedef CGAL::Triangulation<CK>			 							CTriangulation;
typedef CTriangulation::Point_d 									CPoint_d;
typedef CK::Hyperplane_d													CHyperplane_d;
typedef CK::Vector_d															CVector_d;
typedef CTriangulation::Vertex_iterator						CVertex_iterator;
typedef CTriangulation::Vertex_handle							CVertex_handle;
typedef CTriangulation::Full_cell_iterator				CSimplex_iterator;
typedef CTriangulation::Full_cell_handle					CSimplex_d;
typedef CTriangulation::Facet											CFacet;
typedef CTriangulation::Face 											CFace;
typedef CTriangulation::Facet 										CFacet;
typedef CTriangulation::Locate_type 							CLocate_type;

// projection kernel 
typedef CGAL::Cartesian_d<Field>									PK;
typedef CGAL::Triangulation<PK> 									Triangulation;
typedef PK::Vector_d															PVector_d;
typedef PK::Hyperplane_d													PHyperplane_d;
typedef Triangulation::Full_cell_handle						PSimplex_d;
typedef Triangulation::Point_d 										PPoint_d;
typedef Triangulation::Face 											PFace;
typedef Triangulation::Facet 											PFacet;
typedef Triangulation::Locate_type 								PLocate_type;
typedef Triangulation::Vertex_handle							PVertex_handle;
typedef Triangulation::Vertex_iterator						PVertex_iterator;

//misc typedefs
typedef vector<Field> 		SRvertex;
typedef vector<SRvertex> 	Resvertex;
//typedef IntegerSet 			SRvertex;
typedef set<SRvertex> 		Polytope;

// big matrix determinants typedefs
typedef FastHashedDeterminant<Field>                    HD;

typedef Normal_Vector_ds<PVector_d,Field>					NV_ds;

/////////////////////////////
//functions implementations


/////////////////////////////////////////////////////////////////
// print functions
// TODO:make them generic with template(?!)

std::ostream& operator<<(std::ostream& ost, const vector<size_t>& V) {
  for (vector<size_t>::const_iterator it=V.begin(); it!=V.end(); it++)
		ost << *it << ",";
  return ost;
}

std::ostream& operator<<(std::ostream& ost, const vector<int>& V) {
  for (vector<int>::const_iterator it=V.begin(); it!=V.end(); it++)
		ost << *it << ",";
  return ost;
}

std::ostream& operator<<(std::ostream& ost, const set<int>& V) {
  for (set<int>::const_iterator it=V.begin(); it!=V.end(); it++)
		ost << *it << ",";
  return ost;
}

std::ostream& operator<<(std::ostream& ost, const vector<set<int> >& V) {
  for (vector<set<int> >::const_iterator it=V.begin(); it!=V.end(); it++)
		ost << *it << "},{";
  return ost;
}

std::ostream& operator<<(std::ostream& ost, const SRvertex& V) {
  for (SRvertex::const_iterator it=V.begin(); it!=V.end(); it++)
		ost << *it << " ";
  return ost;
}


std::ostream& operator<<(std::ostream& ost, const Polytope& P) {
  for (Polytope::const_iterator it=P.begin(); it!=P.end(); it++)
		ost << *it << std::endl;
  return ost;
}

std::ostream& operator<<(std::ostream& ost, const Resvertex& P) {
  for (Resvertex::const_iterator it=P.begin(); it!=P.end(); it++)
		ost << *it << std::endl;
  return ost;
}

void print_vertices(vector<vector<Field> >& Poly, std::ofstream& ofs){
	ofs << "[";
	for (vector<vector<Field> >::const_iterator Polyit=Poly.begin(); Polyit!=Poly.end(); Polyit++){
		ofs << "[";
		
	  for (vector<Field>::const_iterator it=Polyit->begin(); it!=Polyit->end(); it++){
		  if (it!=Polyit->end()-1){
		    ofs << *it << ",";
		  }else{
				ofs << *it;
			}
		}
	  if (Polyit!=Poly.end()-1){
			ofs << "],";
		} else{
			ofs << "]";
		}
	}
	ofs << "]" << std::endl;
}


void print_res_vertices_with_index(Triangulation& T){
  // print the vertices of the res polytope
  std::cout << "[";
	for (PVertex_iterator vit = T.vertices_begin(); vit != T.vertices_end(); vit++)
		std::cout << vit->point() << " | " << vit->point().index() << "],[";
	std::cout << std::endl;
}

void print_res_vertices(Triangulation& T){
  // print the vertices of the res polytope
  std::cout << "[";
	for (PVertex_iterator vit = T.vertices_begin(); vit != T.vertices_end(); vit++)
		std::cout << vit->point() << "],[";
	std::cout << std::endl;
}


void print_statistics(int numoftriangs, int numofinitvertices, int numofvertices, double timeall){
	std::cout << std::endl;
	std::cout << "Number of triangs enumerated \t" << numoftriangs << std::endl;
	std::cout << "Projected Res vertices (INIT)\t" << numofinitvertices << std::endl;
	std::cout << "Projected Res vertices \t\t" << numofvertices << std::endl;
	std::cout << "Time overall   \t\t\t" << timeall << std::endl;
}

///////////////////////////////////////////////////////////
// input functions

int cayley_trick(std::vector<std::vector<Field> >& pointset, map<std::vector<Field>,int>& points_index, std::vector<int>& mi, int& m){
	int d;
	std::cin >> d;
	if (d != D){
		std::cout << "Not matching dimensions of input and compiled program!" << std::endl;
		exit(-1);
	}
	m=0;
	//std::cout << d << std::endl;
	for (int i=0; i<d+1; i++){
		int mi_temp;
		std::cin >> mi_temp;
		m += mi_temp;
		mi.push_back(mi_temp);
	}
	// compute cayley vector to augment pointset
	if (mi.size() != d+1){
		std::cout << "Input error" << std::endl;
		exit(-1);
	}
	vector<vector<Field> > cayley_vec;
	int i=0, j=-1, end=0;
	for (vector<int>::iterator vit=mi.begin() ; vit!=mi.end(); vit++){
		end += *vit;
		while (i<end){
			vector<Field> cayley_point(d,0);
			if (j>=0)	
				cayley_point[j] = 1;
			cayley_vec.push_back(cayley_point);
			i++;
		}
		j++;
	}
	string point;
  while(!std::getline(cin, point, ']').eof())
	{
		vector<Field> ipoint;
	  point.erase( std::remove( point.begin(), point.end(), '[' ), point.end() );
	  string coord;
    if (point[0]==',')
    	point.erase(0,1);
    stringstream stream(point);
    if (!point.empty()){
			while( getline(stream, coord, ',') ){
      istringstream buffer(coord);
 			Field temp;
 			buffer >> temp;
 			ipoint.push_back(temp);
	 		}
		  pointset.push_back(ipoint);
		}	  
	}
	if (m != pointset.size()){
		std::cout << "Input error" << std::endl;
		exit(-1);
	}
	// apply cayley trick and build an index
	int p_index=0;
	vector<vector<Field> >::iterator cit=cayley_vec.begin();
	for (vector<vector<Field> >::iterator vit=pointset.begin(); vit!=pointset.end(); vit++,cit++){
		vit->insert(vit->end(),cit->begin(),cit->end());
		points_index[*vit] = p_index++;
		//std::cout << *vit << "-->" << points_index[*vit] << std::endl;
	}
	// write the poinset to file for topcom test
	ofstream outfile;
  outfile.open("topcom_cayley.txt");
	print_vertices(pointset,outfile);
  return 0;
}

vector<int> proj_first_coord(const int d, int m, std::vector<int>& mi){
	//project at the first coordinate of each mi	
	vector<int> proj(d);
	proj[0]=0;
	for (int i=1; i<d; i++){
		int mm=0;
		for (int j=0; j<i; j++)
			mm+=mi[j];
		proj[i]=mm;
		//std::cout << mm << " ";
	}
	return proj;
}

vector<int> full_proj(const int d, int m, std::vector<int>& mi){
	vector<int> proj(d);
	proj[0]=0;
	for (int i=1; i<d; i++){
		proj[i]=i;
	}
	return proj;
}


// iterate through all facets determine the upper (i.e. compute determinants), 
// project and compute the r_vector (i.e. compute other determinants)
// TODO: the second determinants are minors of the first

vector<Field> project_upper_hull_r(CTriangulation& pc,
																	 map<vector<Field>, int>& points_index,
																	 HD& dets,
																	 int dcur, 
																	 int dim,
																	 std::vector<int>& mi,
																	 vector<int> proj,
																	 bool lower){
	// the result vector
	vector<Field> rho(PD,0);
	// compute the infinite simplices
	typedef std::vector<CSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  pc.incident_full_cells(pc.infinite_vertex(), out);
  
	//iterate through all infinite simplices
	for( Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
    
    //put the points of the finite facet of the infinite cell to facet_points
    std::vector<CPoint_d> facet_points;
    std::vector<size_t> det_indices;
    std::vector<Field> lifting;
    for (int i=1; i<=pc.current_dimension(); i++){
			CPoint_d p((*sit)->vertex(i)->point());
			facet_points.push_back(p);
			//std::vector<Field> ppoi;
			//for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
			//	ppoi.push_back(p[j]);
				//std::cout << *vit << ",";
			lifting.push_back(p[dim-1]);
			det_indices.push_back(p.index());
		}
		
		// compute the last coordinate of the cross product 
		// i.e. the last coordinate of an orthogonal vector to the facet
		// note that this is also the volume of the projection
		Field det = dets.homogeneous_determinant(det_indices);
		//std::cout << "det=" << det<<" ";
		bool is_upper2 = det>0;
		
		// compute a point of pc not on the facet
		CPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
    
    // distinguish between the two orthogonal vectors 
    // note that here we have to compute a determinant 
    // of one dimension higher
    if (det!=0){
	    //std::vector<Field> ppoi;
			//for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
			//	ppoi.push_back(opposite_point[j]);
			lifting.push_back(opposite_point[dim-1]);
			det_indices.push_back(opposite_point.index());
	    //std::cout << "\n"<<det_indices.size()<<" "<<lifting.size()<<std::endl;
	    Field det2 = dets.homogeneous_determinant(det_indices,lifting);
			//std::cout << " det2=" << det2;
			if (det2<0)
				is_upper2=!is_upper2;
		}	
    if (lower) is_upper2=!is_upper2;
    
    // if the facet is upper     
    if (is_upper2){
			// take the indices, sorted!
			set<int> s;
			for (std::vector<CPoint_d>::iterator vit=facet_points.begin();
			      vit!=facet_points.end(); ++vit){
				s.insert(vit->index());
		  }
			// check if the projection of the facet (i.e. a Minkowski cell)
			// is mixed (i.e. has exactly one vertex summand) 
			int d = mi.size();
		  vector<vector<int> > sets(d);
		  int current_set=0, current_m=mi[0];
		  for (set<int>::const_iterator it2=s.begin(); it2!=s.end(); it2++){
				if (*it2 >= current_m)
					current_m += mi[++current_set];
				sets[current_set].push_back(*it2);
		  } 
		  vector<int> mixed_vertices;
		  for (int i=0; i<d; i++){
		  	if (sets[i].size() == 1){
		  	  mixed_vertices.push_back(i);
		  	}
		  }
		  if (mixed_vertices.size() == 1){
		  	// find where is the mixed vertex in the vector of projection coordinates
		  	vector<int>::iterator pit = find(proj.begin(),proj.end(),sets[mixed_vertices[0]][0]);
		  	if (pit != proj.end()){
		  		rho[pit-proj.begin()] += det;
				}
		  }
		}
	}
  return rho;
}


void update_normal_list(Triangulation& ppc, NV_ds& normal_list, 
                        std::vector<PSimplex_d>::iterator sit){
	std::vector<PPoint_d> facet_points;
	for (int i=1; i<=ppc.current_dimension(); i++){
		facet_points.push_back((*sit)->vertex(i)->point());
		//std::cout << (*sit)->vertex(i)->point() << " ";
	}
	//std::cout << std::endl;
	PPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
	PHyperplane_d hp(facet_points.begin(),facet_points.end(),facet_points[0]);
	//std::cout << opposite_point << "#" << hp.has_on_positive_side(opposite_point) << "!";
	//std::cout << "dim=" << ppc.current_dimension() << " " << hp.orthogonal_direction()<< "\n";
	PVector_d vec = hp.orthogonal_direction().vector();
	if (hp.has_on_positive_side(opposite_point))
		vec = -vec;
	normal_list.put(vec);
}

// TODO: des todo pio katw
void insert_new_Rvertex(Triangulation& ppc, 
												NV_ds& normal_list, 
												SRvertex& new_vertex,
												HD& Pdets){	
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point 
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(ppc.number_of_vertices());
	new_point.set_hash(&Pdets); 
 	#ifdef PRINT_INFO  
		std::cout << "one new R-vertex found !!! "<< std::endl; 
	#endif
	PVertex_handle new_vert = ppc.insert(new_point);
	//update normal_list
	typedef std::vector<PSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  // TODO:make it more efficient
  // find only the simplices incident to the edge (new_vert,inf_vert)
	ppc.incident_full_cells(new_vert, out);
	for(Simplices::iterator sit = inf_simplices.begin(); 
                         sit != inf_simplices.end(); ++sit ){
		if (ppc.is_infinite((*sit)->vertex(0))){
			update_normal_list(ppc, normal_list, sit);
		}
	}
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for new aglorithm


// lift to zero the points that are not to be projected
// there indices are not into vector<int>& proj

vector<pair<vector<Field>,size_t> > lift_to_zero(vector<vector<Field> >& points,
													  										 vector<int>& proj){
															
  vector<pair<vector<Field>,size_t> > lifted_points;
  vector<int>::iterator lpit=proj.begin();
  for (vector<vector<Field> >::iterator pit=points.begin(); pit!=points.end(); pit++){
  	int point_index = pit-points.begin();
  	if (point_index == *lpit){
  		lpit++;
  	} else {
	  	vector<Field> lifted_point = *pit;
	  	lifted_point.push_back(0);
	  	pair<vector<Field>,size_t> lifted_point_with_index(lifted_point,point_index);
	  	lifted_points.push_back(lifted_point_with_index);
	  }
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute the Triangulation of the static points (lifted to zero)
// the points to be lifted will not be present here

void StaticTriangulation(vector<vector<Field> >& points,  
												 map<vector<Field>,int>& points_index,
												 std::vector<int>& proj, 
												 CTriangulation& CT,
												 HD& dets){
	// lift the "static" points Cayley pointset 
	// (i.e. the points that are not going to be projected)
  vector<pair<vector<Field>,size_t> > P = lift_to_zero(points,proj);

  //std::cout << "lifted points:(dim=" << dim << "\n" << P << std::endl;
  for (vector<pair<vector<Field>,size_t> >::iterator vit=P.begin();
       vit!=P.end(); vit++){
		CPoint_d p(CD,vit->first.begin(),vit->first.end());
		//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
  	p.set_index(vit->second);
  	p.set_hash(&dets);
  	CVertex_handle v = CT.insert(p);
  }
}


// lift the points that are going to be projected using nli

vector<pair<vector<Field>,size_t> > lift_to_proj(vector<vector<Field> >& points, 
																					       PVector_d nli, 
																					       std::vector<int> proj){
																			
  vector<pair<vector<Field>,size_t> > lifted_points;
  int i=0,j=0;
  for (vector<vector<Field> >::iterator vit=points.begin(); vit!=points.end(); vit++){
  	vector<Field> lifted_point = *vit;
  	if (proj[j] == i){
			//std::cout << nli << "|"<<proj<<std::endl;
  		lifted_point.push_back(nli[j]);
  		pair<vector<Field>,size_t> lifted_point_with_index(lifted_point,i);
	  	lifted_points.push_back(lifted_point_with_index);
  		j++;
  	}
  	i++;
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute lifting triangulation but not from scratch
// use a copy of CT and add new lifted points there

void LiftingTriangulationDynamic(vector<vector<Field> >& points, 
																						 PVector_d& nli, 
																						 std::vector<int>& proj, 
																						 map<vector<Field>,int>& points_index, 
																						 CTriangulation& CT,
																						 CTriangulation& LCT,
																						 HD& dets,
																						 bool lower){	
	  
  // lift the points to be projected
  vector<pair<vector<Field>,size_t> > P = lift_to_proj(points,nli,proj);

  // add the new lifted points to the copy (LCT) 
  for (vector<pair<vector<Field>,size_t> >::iterator vit=P.begin();
       vit!=P.end(); vit++){
		CPoint_d p(CD,vit->first.begin(),vit->first.end());
		//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
  	p.set_index(vit->second);
  	p.set_hash(&dets);
  	CVertex_handle v = LCT.insert(p);
  }
}

////////////////////////////////////////////////////////////
// the new (hopefully faster) algorithm
void compute_res_faster( std::vector<std::vector<Field> >& pointset, 
                 	 			 map<std::vector<Field>,int>& points_index, 
                  			 int m, 
				                 std::vector<int>& mi, 
				                 vector<int>& proj,
				                 HD& dets, 
				                 HD& Pdets,
				                 int& numof_triangs, 
				                 int& numof_init_Res_vertices, 
				                 Triangulation& T){
	
	CTriangulation CT(CD);
	StaticTriangulation(pointset,points_index,proj,CT,dets);
												 
  // make a stack (stl vector) with normals vectors and initialize
  NV_ds normal_list_d;
  int init_size = normal_list_d.size();
    
  while(!normal_list_d.empty()){
		#ifdef PRINT_INFO
		  std::cout << "dim=" << T.current_dimension() << std::endl;
			normal_list_d.print();
	  #endif
	  
	  // take the first(last more efficient with vectors??) normal and remove it from normals stack
		PVector_d nli = *(normal_list_d.begin());
		normal_list_d.erase(normal_list_d.begin());
		//PVector_d nli = normal_list_d.back();
		//normal_list_d.pop_back();
		
		// make a copy of CT
		CTriangulation LCT(CT);
		
		// outer normal vector --> upper hull projection
	  LiftingTriangulationDynamic(pointset,nli,proj,points_index,CT,LCT,dets,false); 
		
		// project LCT i.e. triangulation
	  int dcur = LCT.current_dimension();
	  vector<Field> new_vertex = project_upper_hull_r(LCT,points_index,dets,dcur,CD,mi,proj,false);
	  
	  // TODO: destroy here!
	  LCT.clear();
	  
	  #ifdef PRINT_INFO  
		  cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
		#endif
		  
	  // insert it in the complex T (if it is not already there) 
	  if (Pdets.find(new_vertex) == -1)
	  	insert_new_Rvertex(T,normal_list_d,new_vertex,Pdets);
		
		numof_triangs++;
		if (numof_triangs == init_size)
			numof_init_Res_vertices = T.number_of_vertices();
	}
	print_res_vertices(T);
}


////////////////////////////////////////////////////////////
// misc


int num_of_vertices(CTriangulation& ppc){	
	//PFace f;	
	typedef std::vector<CFace> Faces;
	Faces edges;
	std::back_insert_iterator<Faces> out(edges);
	ppc.incident_faces(ppc.infinite_vertex(), 1, out);
	// Count the number of points on the convex hull
	std::cout << "[";	
	for (Faces::iterator fit=edges.begin(); fit!=edges.end(); fit++){
		CPoint_d p = (*fit).vertex(1)->point();
		for (CPoint_d::Cartesian_const_iterator pit=p.cartesian_begin(); pit!=p.cartesian_end(); pit++){
			if (pit==p.cartesian_end()-1)
				std::cout << *pit;	
			else	
				std::cout << *pit << ",";
		}
		std::cout << "],[";
	}
	std::cout << ppc.empty() << std::endl;
	return edges.size();
}

int num_of_simplices(CTriangulation& ppc,map<vector<Field>, int>& points_index){	
	
	// Count the number of points on the convex hull
	std::cout << "{";	
	for (CSimplex_iterator fit=ppc.full_cells_begin(); fit!=ppc.full_cells_end(); fit++){
		std::vector<CPoint_d> simplex_points;
		for (int i=0; i<=ppc.current_dimension(); i++){
			CPoint_d p = fit->vertex(i)->point();
			simplex_points.push_back(p);
			std::vector<Field> ppoi;
			for (int j=0; j<ppc.current_dimension()-1; ++j) // d-1 cuts the last coordinate!!
				ppoi.push_back(p[j]);
			std::cout << points_index[ppoi] << ",";				
		}
		std::cout << "},{";
	}
	std::cout << ppc.empty() << std::endl;
	return 1;
}

////////////////////////////////////////////////////////////////////
// OLD , not used in current version
/*
void insert_new_Rvertex(Triangulation& ppc, 
												NV_ds& normal_list, 
												PPoint_d& new_point){
	PLocate_type loc_type;
  PFace f(PD);
	PFacet ft;
  ppc.locate(new_point,loc_type,f,ft); 
  if (loc_type==4 || loc_type==5){ //replace with OUTSIDE_AFFINE_HULL OUTSIDE_CONVEX_HULL
  	#ifdef PRINT_INFO  
  	  std::cout << "one new R-vertex found !!! loc_type="<< loc_type << std::endl; 
		#endif
		PVertex_handle new_vert = ppc.insert(new_point);
		//update normal_list
		typedef std::vector<PSimplex_d> Simplices;
		Simplices inf_simplices;
	  std::back_insert_iterator<Simplices> out(inf_simplices);
	  //TODO:make it more efficient
	  // find only the simplices indident to the edge (new_vert,inf_vert)
		ppc.incident_full_cells(new_vert, out);
		for(Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
			if (ppc.is_infinite((*sit)->vertex(0))){
				update_normal_list(ppc, normal_list, sit);
			}
		}
 	}else{
		;//cout << "no new R-vertex found"<< endl; 
	}
}
* 
* 
* 
void insert_new_Rvertex(Triangulation& ppc, 
												NV_ds& normal_list, 
												SRvertex& new_vertex,
												HD& Pdets){	
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point 
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(ppc.number_of_vertices());
	new_point.set_hash(&Pdets); 
	//print 
	//print_res_vertices(ppc);
	//Pdets.print_matrix(std::cout);
	// insert the point to the triangulation
	PLocate_type loc_type;
  PFace f(PD);
	PFacet ft;
  ppc.locate(new_point,loc_type,f,ft); 
  if (loc_type==4 || loc_type==5){ //replace with OUTSIDE_AFFINE_HULL OUTSIDE_CONVEX_HULL
  	#ifdef PRINT_INFO  
  	  std::cout << "one new R-vertex found !!! loc_type="<< loc_type << std::endl; 
		#endif
		PVertex_handle new_vert = ppc.insert(new_point);
		//update normal_list
		typedef std::vector<PSimplex_d> Simplices;
		Simplices inf_simplices;
	  std::back_insert_iterator<Simplices> out(inf_simplices);
	  //TODO:make it more efficient
	  // find only the simplices indident to the edge (new_vert,inf_vert)
		ppc.incident_full_cells(new_vert, out);
		for(Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
			if (ppc.is_infinite((*sit)->vertex(0))){
				update_normal_list(ppc, normal_list, sit);
			}
		}
 	}else{
		;//cout << "no new R-vertex found"<< endl; 
	}
}
* 
* 

// iterate through all facets and put all the points of upper (lower) hull facets to triang the upper hull is the default
vector<set<int> > project_lower_upper_hull(CTriangulation& pc,
                                           map<vector<Field>, int>& points_index,
                                           int dcur, 
                                           int dim, 
                                           bool lower){
	//std::cout << dim << "|" << pc.current_dimension() << std::endl;
	//SimplicialComplex triang, triang2;
	vector<set<int> > triang;
	typedef std::vector<CSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  pc.incident_full_cells(pc.infinite_vertex(), out);
	//iterate through all infinite simplices
	for( Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
    //put the points of the finite facet of the infinite cell to facet_points
    std::vector<CPoint_d> facet_points;
    for (int i=1; i<=pc.current_dimension(); i++){
			facet_points.push_back((*sit)->vertex(i)->point());
		}
		
		// compute a point of pc not on the facet
		CPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
    
    // compute a hyperplane which has in its negative side the opposite point
    CHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
        
    CVector_d vec = hp.orthogonal_direction().vector();
    //std::cout << "opp=" << hp.has_on_positive_side(vec+facet_points[0]);
    bool is_upper = vec.cartesian(CD-1) > 0;
    if (lower) is_upper=!is_upper;
    
    if (is_upper){
			set<int> s;
			for (std::vector<CPoint_d>::iterator vit=facet_points.begin();
			      vit!=facet_points.end(); ++vit){
				std::vector<Field> ppoi;
				for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
					ppoi.push_back((*vit)[j]);
				//std::cout << *vit << ",";
				s.insert(points_index[ppoi]);
		  }
			//if (s.card()==dim){
			triang.push_back(s);
			//}else{
			//  #ifdef PRINT_INFO  
			//  	std::cout << "found a low dimensional simplex. possible ERROR\n" << std::endl;
			//	#endif
			//}
		}
	}
  return triang;
}
* 
* 
* 
void update_normal_list(Triangulation& ppc, std::vector<PVector_d>& normal_list, 
                        std::vector<PSimplex_d>::iterator sit){
	std::vector<PPoint_d> facet_points;
	for (int i=1; i<=ppc.current_dimension(); i++){
		facet_points.push_back((*sit)->vertex(i)->point());
		//std::cout << (*sit)->vertex(i)->point() << " ";
	}
	//std::cout << std::endl;
	PPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
	PHyperplane_d hp(facet_points.begin(),facet_points.end(),facet_points[0]);
	//std::cout << opposite_point << "#" << hp.has_on_positive_side(opposite_point) << "!";
	//std::cout << "dim=" << ppc.current_dimension() << " " << hp.orthogonal_direction()<< "\n";
	PVector_d vec = hp.orthogonal_direction().vector();
	if (hp.has_on_positive_side(opposite_point))
		vec = -vec;
	if (std::find(normal_list.begin(),normal_list.end(),vec) 
			      == normal_list.end()){
		normal_list.push_back(vec);
	}
}
* 
* 
* 
vector<Field> compute_r_proj(const vector<set<int> >& triang, vector<vector<Field> >& points, int m, std::vector<int>& mi, int d, map<std::vector<Field>,int>& points_index,HD& dets, vector<int> proj){
	double top1, top2;
	vector<Field> r(PD,0);
	for (vector<set<int> >::const_iterator it = triang.begin(); it != triang.end();it++){
		//vector<vector<Field> > simplex;
		vector<size_t> simplex_vec;
	  vector<vector<int> > sets(d);
	  int current_set=0, current_m=mi[0];
	  for (set<int>::const_iterator it2=it->begin(); it2!=it->end(); it2++){
			if (*it2 >= current_m)
				current_m += mi[++current_set];
			sets[current_set].push_back(*it2);
			//vector<Field> point = points[*it2];
			//point.push_back(1);
	    //simplex.push_back(point);
	    simplex_vec.push_back(*it2); 
	  } 
	  vector<int> mixed_vertices;
	  for (int i=0; i<d; i++){
	  	if (sets[i].size() == 1){
	  	  mixed_vertices.push_back(i);
	  	}
	  }
	  if (mixed_vertices.size() == 1){
	  	// find where is the mixed vertex in the vector of projection coordinates
	  	vector<int>::iterator pit = find(proj.begin(),proj.end(),sets[mixed_vertices[0]][0]);
	  	if (pit != proj.end()){
	  		//std::cout << simplex_vec << std::endl;
	  	  Field volume = abs(dets.homogeneous_determinant(simplex_vec));
	  		r[pit-proj.begin()] += volume;
			}
	  }
	}
	//std::cout << r << std::endl;
	return r;
}
* 
* 
* 


// iterate through all facets and put all the points of upper (lower) hull facets to triang the upper hull is the default
vector<set<int> > project_lower_upper_hull_fast(CTriangulation& pc,
                                           map<vector<Field>, int>& points_index,
                                           HD& dets,
                                           int dcur, 
                                           int dim, 
                                           bool lower){
	//std::cout << dim << "|" << pc.current_dimension() << std::endl;
	//SimplicialComplex triang, triang2;
	vector<set<int> > triang;
	typedef std::vector<CSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  pc.incident_full_cells(pc.infinite_vertex(), out);
	//iterate through all infinite simplices
	for( Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
    //put the points of the finite facet of the infinite cell to facet_points
    std::vector<CPoint_d> facet_points;
    std::vector<size_t> det_indices;
    std::vector<Field> lifting;
    for (int i=1; i<=pc.current_dimension(); i++){
			CPoint_d p((*sit)->vertex(i)->point());
			facet_points.push_back(p);
			std::vector<Field> ppoi;
			for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
				ppoi.push_back(p[j]);
				//std::cout << *vit << ",";
			lifting.push_back(p[dim-1]);
			det_indices.push_back(points_index[ppoi]);
		}
		// compute the last coordinate of the cross product 
		// i.e. the last coordinate of an orthogonal vector to the facet
		Field det = dets.homogeneous_determinant(det_indices);
		//std::cout << "det=" << det;
		bool is_upper2 = det>0;
		
		// compute a point of pc not on the facet
		CPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
    
    // distinguish between the two orthogonal vectors 
    // note that here we have to compute a determinant 
    // of one dimension higher
    if (det!=0){
	    std::vector<Field> ppoi;
			for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
				ppoi.push_back(opposite_point[j]);
			lifting.push_back(opposite_point[dim-1]);
			det_indices.push_back(points_index[ppoi]);
	    //std::cout << "\n"<<det_indices.size()<<" "<<lifting.size()<<std::endl;
	    Field det2 = dets.homogeneous_determinant(det_indices,lifting);
			//std::cout << " det2=" << det2;
			if (det2<0)
				is_upper2=!is_upper2;
		}	
    if (lower) is_upper2=!is_upper2;
     
    if (is_upper2){
			set<int> s;
			for (std::vector<CPoint_d>::iterator vit=facet_points.begin();
			      vit!=facet_points.end(); ++vit){
				std::vector<Field> ppoi;
				for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate (the lifting)!!
					ppoi.push_back((*vit)[j]);
				//std::cout << *vit << ",";
				s.insert(points_index[ppoi]);
		  }
			triang.push_back(s);
		}
	}
  return triang;
}
*/
