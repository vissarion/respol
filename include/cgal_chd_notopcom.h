#include <list>
#include <vector>
#include <fstream>
#include <iostream>

//includes for CGAL
#include <CGAL/Cartesian.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Coercion_traits.h>
//linbox
#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/LA_LinBox.h>

// include 3d stuff
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Delaunay_triangulation_3.h>

// include d-dimensional things from CGAL
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Delaunay_d.h>

// for fast determinant computation
#include <../include/fast_hashed_determinant.h>

// future includes of TOPCOM
//#include "../TOPCOM-0.16.0/lib-src/Field.hh"

//#include "../../TOPCOM-0.16.2/lib-src/Field.hh"

//includes for TOPCOM hacked
//#include "../../TOPCOM-0.16.2/lib-src/SimplicialComplex.hh"
//#include "../../TOPCOM-0.16.2/lib-src/IntegerSet.hh"
//#include "../../TOPCOM-0.16.2/lib-src/PlacingTriang.hh"

using namespace std;

// typedefs for CGAL
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// kernel
#include <gmpxx.h>
typedef mpq_class Field;

typedef CGAL::Cartesian<Field>  	K;
typedef K::FT											FT;
// linbox
typedef CGAL::Linbox_rational_field<Field>    LBField;
typedef CGAL::LA_LinBox<LBField>              LA;
// d-dim types
typedef CGAL::Cartesian_d<Field>  						FCKernel;
typedef CGAL::Convex_hull_d<FCKernel> 				Convex_hull_d;
typedef CGAL::Delaunay_d<FCKernel> 						Delaunay_d;
typedef CGAL::Point_d<FCKernel>               Point_d;
typedef CGAL::Hyperplane_d<FCKernel>					Hyperplane_d;
typedef CGAL::Direction_d<FCKernel>						Direction_d;
typedef CGAL::Vector_d<FCKernel>						  Vector_d;
typedef Convex_hull_d::Vertex_handle   				Vertex_handle_d;
typedef Convex_hull_d::Vertex_iterator    		Vertex_iterator_d;
typedef Convex_hull_d::Facet_handle   				Facet_handle_d;
typedef Convex_hull_d::Facet_iterator  				Facet_iterator_d;
typedef std::vector<Direction_d>::iterator    LDIterator_3;

// typedefs
typedef vector<Field> 												SRvertex;
typedef vector<SRvertex> 											Resvertex;
typedef set<SRvertex> 												Polytope;
typedef vector<set<int> > 										Triangulation;
typedef vector<vector<Field> >								VertexSet;

// big matrix determinants typedefs
typedef FastHashedDeterminant<Field,CD>         HD;

/////////////////////////////////////////////////////////////////////
// implimentations
////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// print functions
// TODO:make them generic with template(?!)

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

void homogenize(vector<vector<Field> >& pointset, vector<vector<Field> >& homo_pointset){
  for (std::vector<std::vector<Field> >::iterator vit=pointset.begin();
	     vit!=pointset.end(); vit++){
		std::vector<Field> point = *vit;
		point.push_back(1);
		homo_pointset.push_back(point);
	}
}

/////////////////////////////////////////////////////////////
//projection functions 

SRvertex project(SRvertex vert, vector<int>& r){
	SRvertex p_vert;
	for (vector<int>::iterator it=r.begin(); it!=r.end(); it++){
		p_vert.push_back(vert[*it]);
	}
	return p_vert;
}

vector<int> proj_first_coord(const int d, int m, std::vector<int>& mi){
	//project at the first coordinate of each mi	
	vector<int> r(d);
	r[0]=0;
	for (int i=1; i<d; i++){
		int mm=0;
		for (int j=0; j<i; j++)
			mm+=mi[j];
		r[i]=mm;
		//std::cout << mm << " ";
	}
	return r;
}

vector<int> full_proj(const int d, int m, std::vector<int>& mi){
	vector<int> r(d);
	r[0]=0;
	for (int i=1; i<d; i++){
		r[i]=i;
	}
	return r;
}

////////////////////////////////////////////////////////////
// initialization functions
// i.e. construct the init lifting vectors

vector<int> int2vectord(int k, int vash, int d){
	vector<int> b;
	while (k!=0){
		b.insert(b.begin(),k%vash);
		k=k/vash;
	}
	while (b.size() != d)
		b.insert(b.begin(),0);
	return b;
}

bool is_zero(vector<Field> v){
	for (vector<Field>::iterator it=v.begin()+1; it!=v.end(); it++){
		if (*it != Field(0))
			return false;
	}	
	return true;
}

void init_normal_list(Convex_hull_d& ppc, std::vector<Vector_d>& normal_list){
	Field vec[]={-1,0,1};
  int base=sizeof(vec)/sizeof(Field);
   
  for (int i=0; i<pow(base,PD); ++i){
    vector<Field> extreme_point;
    vector<int> v=int2vectord(i,base,PD);
    //copy(v.begin(),v.end(),ostream_iterator<int>(std::cout," "));
    for (vector<int>::iterator it=v.begin(); it!=v.end(); it++){
		  extreme_point.push_back(vec[*it]);
		  //std::cout << vec[*it] << " ";
	  }
	  //std::cout << std::endl;
    if (!is_zero(extreme_point)){
      Vector_d lft(PD,extreme_point.begin(),extreme_point.end());
			normal_list.push_back(lft);
		}
	}
	//normal_list.push_back(Vector_d(0,0,0,0));
}


//////////////////////////////////////////////////////////////////
// linbox 

Field linbox_det(vector<vector<Field> >& simplex,int d){
	double lin1, lin2;
  typedef LA::Matrix  LA_Matrix;
	LA traits;
	LA_Matrix M(d,d);
	//std::cout << std::endl;
	for (int i=0; i<d; i++){
	  for (int j=0; j<d; j++){
			traits.set_matrix_entry(M,i,j,simplex[i][j]);
			//std::cout << simplex(i,j) << " ";
		}
		//std::cout << std::endl;
	}
	Field det;
	//for (int i=0; i<1000; i++)
	  //det=traits.determinant(M,CGAL::CGAL_HYBRID);
	  LinBox::det(det,M, LinBox::RingCategories::RationalTag() ,LinBox::Method::Hybrid());
	//std::cout << lin2-lin1;
	return det;
} 

//////////////////////////////////////////////////////////////////
// basic functions for the algorithm

vector<Field> compute_r_fast(const vector<set<int> >& triang, vector<vector<Field> >& points, int m, std::vector<int>& mi, int d, map<std::vector<Field>,int>& points_index,HD& dets){
	double top1, top2;
	vector<Field> r(m,0);
	for (vector<set<int> >::const_iterator it = triang.begin(); it != triang.end();it++){
		vector<vector<Field> > simplex;
		vector<size_t> simplex_vec;
	  vector<vector<int> > sets(d);
	  int current_set=0, current_m=mi[0];
	  for (set<int>::const_iterator it2=it->begin(); it2!=it->end(); it2++){
			if (*it2 >= current_m)
				current_m += mi[++current_set];
			sets[current_set].push_back(*it2);
			vector<Field> point = points[*it2];
			point.push_back(1);
	    simplex.push_back(point);
	    simplex_vec.push_back(points_index[points[*it2]]);
	  } 
	  //Field volume = abs(det(simplex));
	  //Field volume = abs(linbox_det(simplex,CD));
	  //std::cout << simplex_vec << std::endl;
	  //Field volume = abs(dets.determinant(simplex_vec));
	  
	  vector<int> mixed_vertices;
	  for (int i=0; i<d; i++){
	  	if (sets[i].size() == 1){
	  	  mixed_vertices.push_back(i);
	  	}
	  }
	  if (mixed_vertices.size() == 1){
	  	Field volume = abs(dets.determinant(simplex_vec));
	  	r[sets[mixed_vertices[0]][0]] += volume;
	  }  
	}
	std::cout << r << std::endl;
	return r;
}

vector<Field> compute_r_proj(const vector<set<int> >& triang, vector<vector<Field> >& points, int m, std::vector<int>& mi, int d, map<std::vector<Field>,int>& points_index,HD& dets, vector<int> proj){
	double top1, top2;
	vector<Field> r(PD,0);
	for (vector<set<int> >::const_iterator it = triang.begin(); it != triang.end();it++){
		vector<vector<Field> > simplex;
		vector<size_t> simplex_vec;
	  vector<vector<int> > sets(d);
	  int current_set=0, current_m=mi[0];
	  for (set<int>::const_iterator it2=it->begin(); it2!=it->end(); it2++){
			if (*it2 >= current_m)
				current_m += mi[++current_set];
			sets[current_set].push_back(*it2);
			vector<Field> point = points[*it2];
			point.push_back(1);
	    simplex.push_back(point);
	    simplex_vec.push_back(points_index[points[*it2]]);
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
	  		Field volume = abs(dets.determinant(simplex_vec));
	  		r[pit-proj.begin()] += volume;
			}
	  }
	}
	//std::cout << r << std::endl;
	return r;
}

std::vector<Direction_d> CH_normal_vectors_d(Convex_hull_d& CH){
	std::vector<Direction_d> normal_list_d;
  for (Facet_iterator_d it=CH.facets_begin(); it!=CH.facets_end(); ++it){
  	// if lower=true dir is the outer normal otherwise is the inner
  	Direction_d dir = CH.hyperplane_supporting(it).orthogonal_direction();
  	normal_list_d.push_back(dir);
  	//std::cout << dir << " ";
  }
  return normal_list_d;
}

vector<vector<Field> > lift(vector<vector<Field> >& points, Direction_d nli, std::vector<int> r){
  vector<vector<Field> > lifted_points;
  int i=0,j=0;
  for (vector<vector<Field> >::iterator vit=points.begin(); vit!=points.end(); vit++){
  	vector<Field> lv = *vit;
  	if (r[j] == i){
  		lv.push_back(nli.delta(j));
  		j++;
  	}	else { 
  		lv.push_back(0);
  	}
  	i++;
  	//std::cout << lv << "," << nli << ":" << r << "\n";
  	lifted_points.push_back(lv);
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}

// iterate through all facets and put all the points of lower hull facets to triang
vector<set<int> > project_lower_upper_hull(Convex_hull_d& CH, map<vector<Field>, 
                                     int>& points_index,
                                     int dcur, int dim, bool lower){
	vector<set<int> > triang;
	for (Facet_iterator_d it=CH.facets_begin(); it!=CH.facets_end(); ++it){
  	vector<Point_d> facet_points;
  	for (int i=0; i<dcur; ++i) //or d??, or dcur??
  		facet_points.push_back(CH.vertex_of_facet(it,i)->point());
  	// if lower=true dir is the outer normal otherwise is the inner
  	Direction_d dir = CH.hyperplane_supporting(it).orthogonal_direction();
  	dir = lower ? dir.opposite() : dir;
  	// lower or upper hull
  	bool hull_check = lower ? dir.delta(dim-1) < 0 : dir.delta(dim-1) > 0;
  	if (hull_check){ //if in lower(upper) hull
			set<int> s;
			for (vector<Point_d>::iterator vit=facet_points.begin();
			      vit!=facet_points.end(); ++vit){
				vector<Field> ppoi;
				for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate!!
					ppoi.push_back((*vit)[j]);			
				//std::cout << *vit << ",";
				s.insert(points_index[ppoi]);
			}
			//std::cout << "\n";
			if(s.size() > dim){
				//TODO:test it!!!!! 
	  		std::cout << "found a cell that is not a simplex! ERROR! ...\n";
	  	} else if(s.size() < dim){
				;//std::cout << "found a low dimensional simplex. possible ERROR\n";
	  	  //triang.insert(s);
	  	} else {
	  	  triang.push_back(s);
	  	}
  	}
  }
  return triang;
}

vector<set<int> > LiftingTriangulation(vector<vector<Field> >& points, Direction_d nli, std::vector<int> r, map<vector<Field>,int>& points_index, bool lower){
	// lift Cayley pointset
  vector<vector<Field> > P = lift(points,nli,r);
  //std::cout << "points[" << i << j << "]" << points[i][j] << std::endl;
  //std::cout << points.no() << " " << points.rank() << std::endl;
  
  // dim is the initial cayley dimension of the points -1 for 
  // removing the homogenizing ones, +1 for the lifting
  // CH in 4+1 dimensions (cayley dimension+1) here dim will be 5
  //int dim = points.rank(); 
  Convex_hull_d CH(CD);
  //std::cout << "lifted points:(dim=" << dim << "\n" << P << std::endl;
  for (int i=0; i<P.size(); i++){
		//std::cout << Point_d(CD,P[i].begin(),P[i].end()) << " point inserted" << std::endl;
  	CH.insert(Point_d(CD,P[i].begin(),P[i].end()));
  }
  
  // project CH i.e. triangulation
  int dcur = CH.current_dimension();
  vector<set<int> > t;
  t = project_lower_upper_hull(CH,points_index,dcur,CD,lower);
  #ifdef PRINT_INFO  
    std::cout << t << std::endl;
  #endif
  CH.clear(dcur);
	return t;
}

// VERY NAIVE!!!!
void insert_new_Rvertex(Convex_hull_d& CH, std::vector<Vector_d>& normal_list_d, Point_d new_point){
  int before = CH.number_of_vertices();
  CH.insert(new_point);
  int after = CH.number_of_vertices();
	if (after!=before){
		//normal_list_d.clear();
		for (Facet_iterator_d fit=CH.facets_begin(); fit!=CH.facets_end(); ++fit){
			Vector_d dir = CH.hyperplane_supporting(fit).orthogonal_direction().vector();
			if (std::find(normal_list_d.begin(),normal_list_d.end(),dir) 
			      == normal_list_d.end()){
			  if (dir.dimension() == PD) 
			  	normal_list_d.push_back(dir);
			}
	  }
	}
}

//update_CH()

////////////////////////////////////////////////////////////////////
// The algorithm!

void compute_res(std::vector<std::vector<Field> >& pointset, map<std::vector<Field>,int>& points_index, int m, std::vector<int>& mi, vector<int>& proj, HD& dets, int& numof_triangs, int& numof_init_Res_vertices, Convex_hull_d& CH){

	// make a stack (stl vector) with normals vectors
  std::vector<Vector_d> normal_list_d;
  // initialize with some normal vectors
  init_normal_list(CH, normal_list_d);
  int init_size = normal_list_d.size();
  
  while(!normal_list_d.empty()){
		#ifdef PRINT_INFO
		  std::cout << "dim=" << CH.current_dimension() << std::endl;
			std::cout << "current normal vectors:" << " ";
			typedef std::vector<Vector_d>::iterator LDIterator_d;
			for (LDIterator_d it=normal_list_d.begin(); it!=normal_list_d.end(); ++it){
		  	std::cout << *it << " ";
		  }
		  std::cout << std::endl;
	  #endif
	  
	  // take the first(last more efficient with vectors??) normal and remove it from normals stack
		Vector_d nli = *(normal_list_d.begin());
		normal_list_d.erase(normal_list_d.begin());
		
		// outer normal vector --> upper hull projection
	  vector<set<int> > t = LiftingTriangulation(pointset,nli,proj,points_index,false);
	  //std::cout << t << std::endl;
	  
	  if (t.size() != 0){//if t is indeed a triangulation!
		  // compute the new Res vertex from t
		  //std::cout << compute_r_fast(t, points, m, mi, mi.size()) << std::endl;
		  //SRvertex new_vertex = project(compute_r_fast(t, pointset, m, mi, mi.size()),proj);
		  //SRvertex new_vertex = project(compute_r_fast(t, pointset, m, mi, mi.size(),points_index,dets),proj);
		  SRvertex new_vertex = compute_r_proj(t, pointset, m, mi, mi.size(),points_index,dets,proj);
		  #ifdef PRINT_INFO  
		    cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
		  #endif
		  // insert it in the complex ppc (if it is not already there) 
		  Point_d new_point(PD,new_vertex.begin(),new_vertex.end());
		  insert_new_Rvertex(CH,normal_list_d,new_point);
	  }
		numof_triangs++;
		if (numof_triangs == init_size)
			numof_init_Res_vertices = CH.number_of_vertices();
			
	}
}

// not working :(
// I cannot copy the CH object nor delete a vertex !!!!!

/*
void compute_res_fast(std::vector<std::vector<Field> >& pointset, map<std::vector<Field>,int>& points_index, int m, std::vector<int>& mi, vector<int>& proj, HD& dets, int& numof_triangs, int& numof_init_Res_vertices, Convex_hull_d& CH){
	// the lifted static convex hull
	Convex_hull_d sCH(CD);
	
	// we have d points that will be lifted in all steps (i.e. dynamic points)
	// and n-d points that will stay in 0 (i.e.) static points
	
	// lift static points to 0 
	vector<vector<Field> > static_points;
  for (vector<vector<Field> >::iterator vit=pointset.begin(); vit!=pointset.end(); vit++){
  	vector<Field> lv = *vit;
  	lv.push_back(0);
  	static_points.push_back(lv);
  }
	
	// put static points on the convex hull
	for (vector<vector<Field> >::iterator Pit=static_points.begin(); Pit!=static_points.end(); Pit++){
		std::cout << Point_d(CD,Pit->begin(),Pit->end()) << "static point inserted" << std::endl;
  	Vertex_handle_d vertex = sCH.insert(Point_d(CD,Pit->begin(),Pit->end()));
  }
		
	// print points on sCH	
	for (Vertex_iterator_d vit = sCH.vertices_begin(); vit != sCH.vertices_end(); vit++)
		std::cout << vit->point() << "\n";
	std::cout << std::endl;
	 
	// make a stack (stl vector) with normals vectors
  std::vector<Vector_d> normal_list_d;
  // initialize with some normal vectors
  init_normal_list(CH, normal_list_d);
  int init_size = normal_list_d.size();
  
  // make a copy of static CH (i.e. dCH)
	Convex_hull_d dCH(sCH);
	
	//lift dynamic points
	Vector_d lift_test = *(normal_list_d.begin());
	vector<vector<Field> > dynamic_points;
	Vector_d::Cartesian_const_iterator lit = lift_test.cartesian_begin();
	for (vector<int>::iterator pit=proj.begin(); pit!=proj.end(); pit++){
		vector<Field> lv = pointset[*pit];
		lv.push_back(0);
		dynamic_points.push_back(lv);
	}
	
	// and insert lifted dynamic points to dCH
	for (vector<vector<Field> >::iterator Pit=dynamic_points.begin(); Pit!=dynamic_points.end(); Pit++){
		//std::cout << Point_d(CD,Pit->begin(),Pit->end()) << " point inserted" << std::endl;
  	Vertex_handle_d vertex = dCH.insert(Point_d(CD,Pit->begin(),Pit->end()));
  }
	
	// print points on CH	
	for (Vertex_iterator_d vit = dCH.vertices_begin(); vit != dCH.vertices_end(); vit++)
		std::cout << vit->point() << "\n";
	std::cout << std::endl;

  exit(0);
	
  
  while(!normal_list_d.empty()){
		#ifdef PRINT_INFO
		  std::cout << "dim=" << CH.current_dimension() << std::endl;
			std::cout << "current normal vectors:" << " ";
			typedef std::vector<Vector_d>::iterator LDIterator_d;
			for (LDIterator_d it=normal_list_d.begin(); it!=normal_list_d.end(); ++it){
		  	std::cout << *it << " ";
		  }
		  std::cout << std::endl;
	  #endif
	  
	  // take the first(last more efficient with vectors??) normal and remove it from normals stack
		Vector_d nli = *(normal_list_d.begin());
		normal_list_d.erase(normal_list_d.begin());
		
		// outer normal vector --> upper hull projection
	  vector<set<int> > t = LiftingTriangulation(pointset,nli,proj,points_index,false);
	  //std::cout << t << std::endl;
	  
	  if (t.size() != 0){//if t is indeed a triangulation!
		  // compute the new Res vertex from t
		  //std::cout << compute_r_fast(t, points, m, mi, mi.size()) << std::endl;
		  //SRvertex new_vertex = project(compute_r_fast(t, pointset, m, mi, mi.size()),proj);
		  //SRvertex new_vertex = project(compute_r_fast(t, pointset, m, mi, mi.size(),points_index,dets),proj);
		  SRvertex new_vertex = compute_r_proj(t, pointset, m, mi, mi.size(),points_index,dets,proj);
		  #ifdef PRINT_INFO  
		    cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
		  #endif
		  // insert it in the complex ppc (if it is not already there) 
		  Point_d new_point(PD,new_vertex.begin(),new_vertex.end());
		  insert_new_Rvertex(CH,normal_list_d,new_point);
	  }
		numof_triangs++;
		if (numof_triangs == init_size)
			numof_init_Res_vertices = CH.number_of_vertices();
			
	}
}
*/

//////////////////////////////////////////////////////////////////
// misc

//NOT WORKING !!!!!!! 
void insert_new_Rvertex_old(Convex_hull_d& CH, std::vector<Direction_d>& normal_list_d, Point_d new_point){
  int before = CH.number_of_vertices();
  CH.insert(new_point);
  int after = CH.number_of_vertices();
	if (after!=before){
		for (Facet_iterator_d fit=CH.facets_begin(); fit!=CH.facets_end(); ++fit){
	  	for (int i=0; i<CH.current_dimension(); i++){
				if (CH.point_of_facet(fit,i) ==  new_point){
					Direction_d dir = CH.hyperplane_supporting(fit).orthogonal_direction();
					normal_list_d.push_back(dir);
				}
			}
	  }
	}
}

Direction_d augment_CH_dim(Convex_hull_d& CH){
	int idim=0;
	std::vector<Point_d> hpoints;
	for (Vertex_iterator_d Vdim=CH.vertices_begin(); Vdim!=CH.vertices_end(); Vdim++){
		hpoints.push_back(CH.associated_point(Vdim));
	}
	Hyperplane_d hp(hpoints.begin(),hpoints.end(),hpoints[0]);
	return hp.orthogonal_direction();
}

