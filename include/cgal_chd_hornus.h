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
//#include <CGAL/Gmpq.h>
//typedef CGAL::Gmpq Field;

// Cayley space kernel for the CH of the lifted points
typedef CGAL::Cartesian_d<Field> 									CK;
//typedef Indexed_Cartesian_d<Field> 							CK;
//typedef CK::IndexedPoint_d											IndexedPoint_d;

typedef CGAL::Triangulation_vertex<CK>						tv;
typedef CGAL::Triangulation_full_cell<CK,bool>		tc;
typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<PD>,tv,tc > tds;
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
//typedef CGAL::Triangulation<PK> 								Triangulation;
typedef CGAL::Triangulation<PK,tds> 							Triangulation;
typedef PK::Vector_d															PVector_d;
typedef PK::Hyperplane_d													PHyperplane_d;
typedef Triangulation::Full_cell_handle						PSimplex_d;
typedef Triangulation::Finite_full_cell_iterator	PSimplex_iterator;
typedef Triangulation::Point_d 										PPoint_d;
typedef Triangulation::Face 											PFace;
typedef Triangulation::Facet 											PFacet;
typedef Triangulation::Facet_iterator							PFacet_iterator;
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
// math
//Field factorial(Field n)
//{
//  return (n == Field(1) || n == Field(0)) ? Field(1) : factorial(n - 1) * n;
//}

Field factorial(Field n)
{
  return 1;
}


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

// for input topcom file construction
void print_vertices_hom(vector<vector<Field> >& Poly, std::ofstream& ofs){
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
		ofs << ",1";
	  if (Polyit!=Poly.end()-1){
			ofs << "],";
		} else{
			ofs << "]";
		}
	}
	ofs << "]" << std::endl;
}


void print_res_vertices_with_index(Triangulation& Res){
  // print the vertices of the res polytope
  std::cout << "[";
	for (PVertex_iterator vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++)
		std::cout << vit->point() << " | " << vit->point().index() << "],[";
	std::cout << std::endl;
}

void print_res_vertices(Triangulation& Res){
  // print the vertices of the res polytope
  int number_of_vertices = 0;
	std::cout << "dim=" << Res.current_dimension() << std::endl;
	for (PVertex_iterator vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
		std::cout << "[";
		for (PPoint_d::Cartesian_const_iterator cit = vit->point().cartesian_begin(); 
						cit != vit->point().cartesian_end(); cit++){
			std::cout << *cit;
			if (cit - vit->point().cartesian_begin() != vit->point().dimension()-1)
				std::cout << ",";
		}
		//std::cout << "|" << vit->point().index();
		std::cout << "]";
		if (number_of_vertices++ != Res.number_of_vertices())
			std::cout << ",";
	}
	std::cout << std::endl;
}


void print_statistics(int numoftriangs, 
											int numoftriangs2, 
											int numofvertices, 
											double timeall,
											Field volume){
	std::cout << std::endl;
	std::cout << "Num of triangs enumed (init+augment)\t" 
						<< numoftriangs+numoftriangs2 << " ("
						<< numoftriangs << "+" << numoftriangs2 
						<< ")" << std::endl;
	std::cout << "Projected Res vertices \t\t\t" << numofvertices << std::endl;
	std::cout << "Time overall   \t\t\t\t" << timeall << std::endl;
	std::cout << "Volume   \t\t\t\t" << volume 
						<< " ~ " <<  CGAL::to_double(volume) << std::endl;
}

///////////////////////////////////////////////////////////
// input functions

int cayley_trick(std::vector<std::vector<Field> >& pointset, std::vector<int>& mi, int& m){
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
	if (mi.size() != d+1){
		std::cout << "mi.size() != d+1. The number of polynomials must me one more than the dimansion!" << std::endl;
		exit(-1);
	}
	std::cout<<mi<<std::endl;
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
	//std::cout<<cayley_vec<<std::endl;
	string point;
  while(!std::getline(cin, point, ']').eof())
	{
		vector<Field> ipoint;
		//std::cout<<point<<std::endl;
		point.erase( std::remove( point.begin(), point.end(), ' ' ), point.end() );
	  point.erase( std::remove( point.begin(), point.end(), '[' ), point.end() );
	  //std::cout<<point<<std::endl;
	  string coord;
    if (point[0]==',')
    	point.erase(0,1);
    //std::cout<<point<<std::endl;
    stringstream stream(point);
    if (!point.empty()){
			while( getline(stream, coord, ',') ){
      istringstream buffer(coord);
 			Field temp;
 			buffer >> temp;
 			ipoint.push_back(temp);
	 		}
	 		//std::cout << ipoint << std::endl;
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
		//points_index[*vit] = p_index++;
		//std::cout << *vit << "-->" << points_index[*vit] << std::endl;
	}
	// write the (homogenized) poinset to file for topcom test
	ofstream outfile;
  outfile.open("topcom_cayley.txt");
	print_vertices_hom(pointset,outfile);
  return 0;
}

//////////////// projections

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

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for new algorithm


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
	  	lifted_point.push_back(Field(0));
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
												 std::vector<int>& proj, 
												 CTriangulation& T,
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
  	CVertex_handle v = T.insert(p);
  }
}


// lift the points that are going to be projected using nli

vector<pair<vector<Field>,size_t> > lift_to_proj(vector<vector<Field> >& points, 
																					       PVector_d nli, 
																					       std::vector<int> proj){
																			
  vector<pair<vector<Field>,size_t> > lifted_points;
  //int i=0;
  std::vector<int>::iterator proj_it=proj.begin();
  for (vector<vector<Field> >::iterator vit=points.begin(); 
              vit!=points.end() && proj_it!=proj.end(); vit++){
  	vector<Field> lifted_point = *vit;
  	int pindex=vit-points.begin();
  	//std::cout << vit-points.begin() << "," << proj_it-proj.begin() << std::endl;
  	if (*proj_it == pindex){
			//std::cout << nli << "|"<<proj<<std::endl;
  		lifted_point.push_back(nli[proj_it-proj.begin()]);
  		pair<vector<Field>,size_t> lifted_point_with_index(lifted_point,pindex);
	  	lifted_points.push_back(lifted_point_with_index);
	  	proj_it++;
  	}
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute lifting triangulation but not from scratch
// use a copy of T and add new lifted points there

void LiftingTriangulationDynamic(vector<vector<Field> >& points, 
																						 PVector_d& nli, 
																						 std::vector<int>& proj, 
																						 CTriangulation& T,
																						 CTriangulation& Tl,
																						 HD& dets,
																						 bool lower){	
	  
  // lift the points to be projected
  vector<pair<vector<Field>,size_t> > P = lift_to_proj(points,nli,proj);

  // add the new lifted points to the copy (Tl) 
  for (vector<pair<vector<Field>,size_t> >::iterator vit=P.begin();
       vit!=P.end(); vit++){
		CPoint_d p(CD,vit->first.begin(),vit->first.end());
  	p.set_index(vit->second);
  	p.set_hash(&dets);
  	//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
  	CVertex_handle v = Tl.insert(p);
  }
}


// iterate through all facets determine the upper (i.e. compute determinants), 
// project and compute the r_vector (i.e. compute other determinants)
// TODO: the second determinants are minors of the first

vector<Field> project_upper_hull_r(CTriangulation& pc,
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
		 // std::cout<<s<<std::endl;
		  // if there are not full dimensional cells
		  if (s.size() < CD){
				std::cout << "NOT a triangulation. Abort current computation..." << std::endl;
				vector<Field> rho_empty;
				exit(-1);
				return rho_empty;
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
	// to compute the real volume we have to divide by d!
	for (vector<Field>::iterator vit=rho.begin();vit!=rho.end();vit++){
		*vit=*vit/factorial(pc.current_dimension()-1);
	}
  return rho;
}

//////////////////////////////////////////////////////////////
// insert vertices in triangulation that holds the Resultant
// polytope

void update_normal_list(Triangulation& Res, 
												NV_ds& normal_list, 
                        std::vector<PSimplex_d> inf_simplices){
	typedef std::vector<PSimplex_d> Simplices;
	for(Simplices::iterator sit = inf_simplices.begin(); 
													sit != inf_simplices.end(); ++sit ){
		if (Res.is_infinite(*sit)){
			std::vector<PPoint_d> facet_points;
			for (int i=1; i<=Res.current_dimension(); i++){
				facet_points.push_back((*sit)->vertex(i)->point());
				//std::cout << (*sit)->vertex(i)->point() << " ";
			}
			//std::cout << std::endl;
			PPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
			// compute a hyperplane which has in its negative side the opposite point
	    PHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
	        
	    //PHyperplane_d hp(facet_points.begin(),facet_points.end(),facet_points[0]);
			//std::cout << opposite_point << "#" << hp.has_on_positive_side(opposite_point) << "!";
			//std::cout << "dim=" << Res.current_dimension() << " " << hp.orthogonal_direction()<< "\n";
			PVector_d vec = hp.orthogonal_direction().vector();
			//if (hp.has_on_positive_side(opposite_point))
			//	vec = -vec;
			//std::cout << vec << std::endl;
			normal_list.put(vec);
		}
	}
}

// TODO: des todo pio katw
void insert_new_Rvertex(Triangulation& Res, 
												NV_ds& normal_list, 
												SRvertex& new_vertex,
												HD& Pdets){	
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point 
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(Res.number_of_vertices());
	new_point.set_hash(&Pdets); 
 	#ifdef PRINT_INFO  
		std::cout << "one new R-vertex found !!! "<< std::endl; 
	#endif
	
	int prev_dim = Res.current_dimension();
	//insert it to the triangulation
	PVertex_handle new_vert = Res.insert(new_point);
	int cur_dim = Res.current_dimension();
	
	if (cur_dim != prev_dim && cur_dim>2){
		 //update normal_list
		typedef std::vector<PSimplex_d> Simplices;
		Simplices inf_simplices;
		std::back_insert_iterator<Simplices> out(inf_simplices);

		Res.incident_full_cells(Res.infinite_vertex(), out);
		std::cout<<inf_simplices.size()<<std::endl;
		update_normal_list(Res, normal_list, inf_simplices);
		
	} else {
		//update normal_list
		typedef std::vector<PSimplex_d> Simplices;
		Simplices inf_simplices;
		std::back_insert_iterator<Simplices> out(inf_simplices);
		// TODO:make it more efficient
		// find only the simplices incident to the edge (new_vert,inf_vert)
		Res.incident_full_cells(new_vert, out);
		std::cout<<inf_simplices.size()<<std::endl;
		update_normal_list(Res, normal_list, inf_simplices);
	}
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// TODO: des todo pio katw
void insert_new_Rvertex2(Triangulation& Res, 
												NV_ds& normal_list, 
												SRvertex& new_vertex,
												HD& Pdets){	
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point 
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(Res.number_of_vertices());
	new_point.set_hash(&Pdets); 
 	#ifdef PRINT_INFO  
		std::cout << "one new R-vertex found !!! "<< std::endl; 
	#endif
	
	int prev_dim = Res.current_dimension();
	//insert it to the triangulation
	PVertex_handle new_vert = Res.insert(new_point);
	int cur_dim = Res.current_dimension();
	
	normal_list.clear();
	//update normal_list
	typedef std::vector<PSimplex_d> Simplices;
	Simplices inf_simplices;
	std::back_insert_iterator<Simplices> out(inf_simplices);

	Res.incident_full_cells(Res.infinite_vertex(), out);
	std::cout<<inf_simplices.size()<<std::endl;
	update_normal_list(Res, normal_list, inf_simplices);
	
}


// compute a Res vertex by constructing a lifting triangulation 
// project and compute the volumes of some cells 

vector<Field> compute_res_vertex(std::vector<std::vector<Field> >& pointset, 
				                 std::vector<int>& mi, 
				                 int RD,
				                 vector<int>& proj,
				                 HD& dets, 
				                 HD& Pdets,
				                 Triangulation& Res,
				                 CTriangulation& T,
				                 NV_ds& normal_list_d){
													 
	// take the first(last more efficient with vectors??) normal and remove it from normals stack
	//PVector_d nli = *(normal_list_d.begin());
	//normal_list_d.erase(normal_list_d.begin());
	PVector_d nli = normal_list_d.back();
	normal_list_d.pop_back();
	//std::cout << "normal vector:" << nli << std::endl;
	
	// make a copy of T
	CTriangulation Tl(T);
	
	// outer normal vector --> upper hull projection
  LiftingTriangulationDynamic(pointset,nli,proj,T,Tl,dets,false); 
	//print_res_vertices(Tl);
	
	// project Tl i.e. triangulation
  int dcur = Tl.current_dimension();
  vector<Field> new_vertex = project_upper_hull_r(Tl,dets,dcur,CD,mi,proj,false);
  
  // TODO: destroy here!
  Tl.clear();
  
  #ifdef PRINT_INFO  
	  cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
	#endif
	  
  return new_vertex;
}


// compute Res vertices until it builts a simplex

int initialize_Res(std::vector<std::vector<Field> >& pointset, 
								 std::vector<int>& mi, 
								 int RD,
								 vector<int>& proj,
								 HD& dets, 
								 HD& Pdets,
								 Triangulation& Res,
								 CTriangulation& T){
	
	int num_of_triangs=0;								
	#ifdef PRINT_INFO
	  std::cout << "dim=" << Res.current_dimension() << std::endl;
  #endif
  // make a stack (stl vector) with normals vectors and initialize
  NV_ds normal_list_d;
  normal_list_d.simple_initialize();
  
  // compute trinagulations using normals as liftings until we compute a simplex
  // or run out of normal vectors
  int minD = (PD>RD)?RD:PD;  
  while(Res.current_dimension()<minD && !normal_list_d.empty()){
		vector<Field> new_vertex = compute_res_vertex(pointset,mi,RD,proj,dets,Pdets,Res,T,normal_list_d);
		// insert it in the complex Res (if it is not already there) 
  	if (Pdets.find(new_vertex) == -1 && new_vertex.size() != 0)
  		insert_new_Rvertex(Res,normal_list_d,new_vertex,Pdets);
		num_of_triangs++;
		#ifdef PRINT_INFO
			normal_list_d.print();
			std::cout<<"current number of Res vertices: "
							 <<Res.number_of_vertices()
							 <<" dim=" << Res.current_dimension()<<std::endl;
		#endif
	}
	return num_of_triangs;
}


// compute all Res vertices left

int augment_Res(std::vector<std::vector<Field> >& pointset, 
								 std::vector<int>& mi, 
								 int RD,
								 vector<int>& proj,
								 HD& dets, 
								 HD& Pdets,
								 Triangulation& Res,
								 CTriangulation& T){
									
	int num_of_triangs=0;
	#ifdef PRINT_INFO
		std::cout << "\n\nAUGMENTING RESULTANT POLYTOPE" << std::endl;
		//print_res_vertices(Res);
	  std::cout << "dim=" << Res.current_dimension() << std::endl;
  #endif
  // make a stack (stl vector) with normals vectors and initialize
  NV_ds normal_list_d;
  // initialize normal_list
  typedef std::vector<PSimplex_d> Simplices;
	Simplices inf_simplices;
	std::back_insert_iterator<Simplices> out(inf_simplices);

	Res.incident_full_cells(Res.infinite_vertex(), out);
	//std::cout<<inf_simplices.size()<<std::endl;
	update_normal_list(Res, normal_list_d, inf_simplices);
  
  while(!normal_list_d.empty()){
		vector<Field> new_vertex = compute_res_vertex(pointset,mi,RD,proj,dets,Pdets,Res,T,normal_list_d);
		// insert it in the complex Res (if it is not already there) 
  	if (Pdets.find(new_vertex) == -1 && new_vertex.size() != 0)
  		insert_new_Rvertex(Res,normal_list_d,new_vertex,Pdets);
		#ifdef PRINT_INFO
			normal_list_d.print();
			std::cout<<"current number of Res vertices: "<<Res.number_of_vertices()<<std::endl;
		#endif
		num_of_triangs++;
	}
	return num_of_triangs;
}

////////////////////////////////////////////////////////////
// the new (hopefully faster) algorithm
pair<int,int> compute_res_faster( std::vector<std::vector<Field> >& pointset, 
				                  			 int m, 
								                 std::vector<int>& mi, 
								                 int RD,
								                 vector<int>& proj,
								                 HD& dets, 
								                 HD& Pdets,
								                 Triangulation& Res){
	
	//std::cout << "cayley dim:" << CayleyTriangulation(pointset) << std::endl; 
	
	// construct an initial triangulation of the points that will not be projected
	CTriangulation T(CD);
	StaticTriangulation(pointset,proj,T,dets);
	std::cout << "static dim:" << T.current_dimension() << std::endl; 
												 
  // start by computing a simplex
  int start_triangs = initialize_Res(pointset,mi,RD,proj,dets,Pdets,Res,T);
  
  // augment simplex to compute the Res polytope
  int augment_triangs = augment_Res(pointset,mi,RD,proj,dets,Pdets,Res,T);
  
  // number of triangulations computed
  pair<int,int> num_of_triangs(start_triangs,augment_triangs);
  
  return num_of_triangs;
}


////////////////////////////////////////////////////////////
// misc

Field volume(Triangulation& Res, HD& Pdets){
	Field vol=0;
	for (PSimplex_iterator cit=Res.finite_full_cells_begin(); cit!=Res.finite_full_cells_end(); cit++){
		std::vector<size_t> simplex_points_indices;
		for (int i=0; i<=Res.current_dimension(); i++){
			simplex_points_indices.push_back(cit->vertex(i)->point().index());
		}
		//std::cout<<simplex_points_indices<<" --> "<<
		//abs(Pdets.homogeneous_determinant(simplex_points_indices))
		//<<std::endl;
		vol+=abs(Pdets.homogeneous_determinant(simplex_points_indices));
	}
	vol*=(1/factorial(Res.current_dimension()));
	return vol;
}


int num_of_vertices(CTriangulation& T){	
	//PFace f;	
	typedef std::vector<CFace> Faces;
	Faces edges;
	std::back_insert_iterator<Faces> out(edges);
	T.incident_faces(T.infinite_vertex(), 1, out);
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
	std::cout << T.empty() << std::endl;
	return edges.size();
}

int num_of_simplices(CTriangulation& T,map<vector<Field>, int>& points_index){	
	
	// Count the number of points on the convex hull
	std::cout << "{";	
	for (CSimplex_iterator fit=T.full_cells_begin(); fit!=T.full_cells_end(); fit++){
		std::vector<CPoint_d> simplex_points;
		for (int i=0; i<=T.current_dimension(); i++){
			CPoint_d p = fit->vertex(i)->point();
			simplex_points.push_back(p);
			//std::vector<Field> ppoi;
			//for (int j=0; j<T.current_dimension()-1; ++j) // d-1 cuts the last coordinate!!
			//	ppoi.push_back(p[j]);
			std::cout << p.index() << ",";				
		}
		std::cout << "},{";
	}
	std::cout << T.empty() << std::endl;
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
