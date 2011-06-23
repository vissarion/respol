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
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>

std::ostream& operator<<(std::ostream& ost, const std::vector<size_t>& V);

// for fast determinant computation
#include <../include/fast_hashed_determinant.h>

// for indexed points
//#include <../include/indexed_point.h>

using namespace std;

// the kernel
#include <gmpxx.h>
typedef mpq_class Field;

// Cayley space kernel for the CH of the lifted points
typedef CGAL::Cartesian_d<Field> 									CK;
//typedef Indexed_Cartesian_d<Field> 								CK;
//typedef CK::IndexedPoint_d												IndexedPoint_d;

typedef CGAL::Triangulation<CK> 									CTriangulation;
typedef CTriangulation::Point_d 									CPoint_d;
typedef CK::Hyperplane_d													CHyperplane_d;
typedef CK::Vector_d															CVector_d;
typedef CTriangulation::Vertex_iterator						CVertex_iterator;
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
typedef FastHashedDeterminant<Field,CD-1>        HD;


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
	std::cout << d << D;
	//exit(0);
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

void init_normal_list(Triangulation& ppc, std::vector<PVector_d>& normal_list){
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
      PVector_d lft(PD,extreme_point.begin(),extreme_point.end());
			normal_list.push_back(lft);
		}
	}
	//normal_list.push_back(PVector_d(0,0,-1,0));
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
	  		//std::cout << simplex_vec << std::endl;
	  	  Field volume = abs(dets.homogeneous_determinant(simplex_vec));
	  		r[pit-proj.begin()] += volume;
			}
	  }
	}
	//std::cout << r << std::endl;
	return r;
}

vector<vector<Field> > lift(vector<vector<Field> >& points, PVector_d nli, std::vector<int> proj){
  vector<vector<Field> > lifted_points;
  int i=0,j=0;
  for (vector<vector<Field> >::iterator vit=points.begin(); vit!=points.end(); vit++){
  	vector<Field> lv = *vit;
  	if (proj[j] == i){
  		lv.push_back(nli[j]);
  		j++;
  	}	else { 
  		lv.push_back(0);
  	}
  	i++;
  	//std::cout << lv << "," << nli << ":" << proj << "\n";
  	lifted_points.push_back(lv);
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}

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
    
    // compute a hyperplane which has in its negative side the opposite point
    //CHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
              
    //CVector_d vec = hp.orthogonal_direction().vector();
    //std::cout << " vec=" << vec << std::endl;
    //bool is_upper = vec.cartesian(CD-1) > 0;
    //if (lower) is_upper=!is_upper;
    
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


vector<set<int> > LiftingTriangulation(vector<vector<Field> >& points, 
																			 PVector_d& nli, 
																			 std::vector<int>& proj, 
																			 map<vector<Field>,int>& points_index, 
																			 HD& dets,
																			 bool lower){
	// lift Cayley pointset
  vector<vector<Field> > P = lift(points,nli,proj);
  //std::cout << "points[" << i << j << "]" << points[i][j] << std::endl;
  //std::cout << points.no() << " " << points.rank() << std::endl;
  
  CTriangulation pc(CD);
  
  //std::cout << "lifted points:(dim=" << dim << "\n" << P << std::endl;
  for (int i=0; i<P.size(); i++){
		//std::cout << Point_d(CD,P[i].begin(),P[i].end()) << " point inserted" << std::endl;
  	pc.insert(CPoint_d(CD,P[i].begin(),P[i].end()));
  }
  
  // project CH i.e. triangulation
  int dcur = pc.current_dimension();
  vector<set<int> > t;
  //t = project_lower_upper_hull(pc,points_index,dcur,CD,lower);
  t = project_lower_upper_hull_fast(pc,points_index,dets,dcur,CD,lower);
  pc.clear();
  #ifdef PRINT_INFO  
    std::cout << t << std::endl;
  #endif
	return t;
}



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

std::vector<PVector_d> pc_normal_vectors_d(Triangulation& ppc){
	std::vector<PVector_d> normal_list;
	typedef std::vector<PSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  ppc.incident_full_cells(ppc.infinite_vertex(), out);
	//CPoint_d some_vertex_point = (*(inf_simplices.begin()))->vertex(1)->point();
	for( Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
    update_normal_list(ppc, normal_list, sit);
	}
	return normal_list;
}

void insert_new_Rvertex(Triangulation& ppc, std::vector<PVector_d>& normal_list, PPoint_d& new_point){
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


////////////////////////////////////////////////////////////
// the algorithm

void compute_res(std::vector<std::vector<Field> >& pointset, 
                 map<std::vector<Field>,int>& points_index, 
                 int m, 
                 std::vector<int>& mi, 
                 vector<int>& proj, 
                 HD& dets, 
                 int& numof_triangs, 
                 int& numof_init_Res_vertices, 
                 Triangulation& ppc){
									 
  // make a stack (stl vector) with normals vectors
  std::vector<PVector_d> normal_list_d;
  // initialize with some normal vectors
  init_normal_list(ppc, normal_list_d);
  int init_size = normal_list_d.size();
  
  while(!normal_list_d.empty()){
		#ifdef PRINT_INFO
		  std::cout << "dim=" << ppc.current_dimension() << std::endl;
			std::cout << "current normal vectors:" << " ";
			typedef std::vector<PVector_d>::iterator LDIterator_d;
			for (LDIterator_d it=normal_list_d.begin(); it!=normal_list_d.end(); ++it){
		  	std::cout << *it << " ";
		  }
		  std::cout << std::endl;
	  #endif
	  // take the first(last more efficient with vectors??) normal and remove it from normals stack
		PVector_d nli = *(normal_list_d.begin());
		normal_list_d.erase(normal_list_d.begin());
		
		// outer normal vector --> upper hull projection
	  vector<set<int> > t = LiftingTriangulation(pointset,nli,proj,points_index,dets,false);

	  if (t.size() != 0){//if t is indeed a triangulation!
		  // compute the new Res vertex from t
		  SRvertex new_vertex = compute_r_proj(t, pointset, m, mi, mi.size(),points_index,dets,proj);
		  #ifdef PRINT_INFO  
		    cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
		  #endif
		  // insert it in the complex ppc (if it is not already there) 
		  PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
		  insert_new_Rvertex(ppc,normal_list_d,new_point);
	  }
		numof_triangs++;
		if (numof_triangs == init_size)
			numof_init_Res_vertices = ppc.number_of_vertices();
	}
	print_res_vertices(ppc);
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// functions for new aglorithm


// lift to zero the points that are not to be projected
// there indices are not into vector<int>& proj

vector<vector<Field> > lift_to_zero(vector<vector<Field> >& points,
													  				vector<int>& proj){
															
  vector<vector<Field> > lifted_points;
  int i=0,j=0;
  vector<int>::iterator lpit=proj.begin();
  for (vector<vector<Field> >::iterator pit=points.begin(); pit!=points.end(); pit++){
  	if (pit-points.begin() == *lpit){
  		lpit++;
  	} else {
	  	vector<Field> lifted_point = *pit;
	  	lifted_point.push_back(0);
	  	lifted_points.push_back(lifted_point);
	  }
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute the Triangulation of the static lifted points

void StaticTriangulation(vector<vector<Field> >& points,  
																			 std::vector<int>& proj, 
																			 CTriangulation& CT){
	// lift the "static" points Cayley pointset 
	// (i.e. the points that are not going to be projected)
  vector<vector<Field> > P = lift_to_zero(points,proj);

  //std::cout << "lifted points:(dim=" << dim << "\n" << P << std::endl;
  for (int i=0; i<P.size(); i++){
		CPoint_d p(CD,P[i].begin(),P[i].end());
		//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
  	CT.insert(p);
  }
  //std::cout << "[";
	//for (CVertex_iterator vit = CT.vertices_begin(); vit != CT.vertices_end(); vit++)
	//	std::cout << vit->point() << "|"<<(vit->point()).index()<<"],[";
	//std::cout << std::endl;
	//exit(0);
}


// lift the points that are going to be projected using nli

vector<vector<Field> > lift_to_proj(vector<vector<Field> >& points, 
															      PVector_d nli, 
															      std::vector<int> proj){
																			
  vector<vector<Field> > lifted_points;
  int i=0,j=0;
  for (vector<vector<Field> >::iterator vit=points.begin(); vit!=points.end(); vit++){
  	vector<Field> lv = *vit;
  	if (proj[j] == i){
			//std::cout << nli << "|"<<proj<<std::endl;
  		lv.push_back(nli[j]);
  		lifted_points.push_back(lv);
  		j++;
  	}
  	i++;
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute lifting triangulation but not from scratch
// use a copy of CT and add new lifted points there
vector<set<int> > LiftingTriangulationDynamic(vector<vector<Field> >& points, 
																			 PVector_d& nli, 
																			 std::vector<int>& proj, 
																			 map<vector<Field>,int>& points_index, 
																			 CTriangulation& CT,
																			 HD& dets,
																			 bool lower){
	// make a copy of CT
	CTriangulation LCT(CT);
  
  // lift the points to be projected
  vector<vector<Field> > P = lift_to_proj(points,nli,proj);

  // add the new lifted points to the copy (LCT) 
  for (int i=0; i<P.size(); i++){
		//std::cout << Point_d(CD,P[i].begin(),P[i].end()) << " point inserted" << std::endl;
  	LCT.insert(CPoint_d(CD,P[i].begin(),P[i].end()));
  }
  
  // project LCT i.e. triangulation
  int dcur = LCT.current_dimension();
  vector<set<int> > t;
  //t = project_lower_upper_hull(LCT,points_index,dcur,CD,lower);
  t = project_lower_upper_hull_fast(LCT,points_index,dets,dcur,CD,lower);
  
  // TODO: destroy here!
  LCT.clear();
  
  #ifdef PRINT_INFO  
    std::cout << t << std::endl;
  #endif
	return t;
}

////////////////////////////////////////////////////////////
// a new (hopefully faster) algorithm
void compute_res_faster( std::vector<std::vector<Field> >& pointset, 
                 	 			 map<std::vector<Field>,int>& points_index, 
                  			 int m, 
				                 std::vector<int>& mi, 
				                 vector<int>& proj,
				                 HD& dets, 
				                 int& numof_triangs, 
				                 int& numof_init_Res_vertices, 
				                 Triangulation& T){
	
	CTriangulation CT(CD);
	StaticTriangulation(pointset,proj,CT);
												 
  // make a stack (stl vector) with normals vectors
  std::vector<PVector_d> normal_list_d;
  // initialize with some normal vectors
  init_normal_list(T, normal_list_d);
  int init_size = normal_list_d.size();
  
  while(!normal_list_d.empty()){
		#ifdef PRINT_INFO
		  std::cout << "dim=" << T.current_dimension() << std::endl;
			std::cout << "current normal vectors:" << " ";
			typedef std::vector<PVector_d>::iterator LDIterator_d;
			for (LDIterator_d it=normal_list_d.begin(); it!=normal_list_d.end(); ++it){
		  	std::cout << *it << " ";
		  }
		  std::cout << std::endl;
	  #endif
	  // take the first(last more efficient with vectors??) normal and remove it from normals stack
		PVector_d nli = *(normal_list_d.begin());
		normal_list_d.erase(normal_list_d.begin());
		
		// outer normal vector --> upper hull projection
	  vector<set<int> > t = LiftingTriangulationDynamic(pointset,nli,proj,points_index,CT,dets,false);
		
	  if (t.size() != 0){//if t is indeed a triangulation!
		  // compute the new Res vertex from t
		  SRvertex new_vertex = compute_r_proj(t, pointset, m, mi, mi.size(),points_index,dets,proj);
		  #ifdef PRINT_INFO  
		    cout << "\nnew Res vertex (up)= ( " << new_vertex << ")\n\n";
		  #endif
		  // insert it in the complex T (if it is not already there) 
		  PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
		  insert_new_Rvertex(T,normal_list_d,new_point);
	  }
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

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//UNUSED
/*
Polytope init_uniform(PointConfiguration& points, Chirotope& chiro, std::vector<int> r, map<vector<Field>, int>& points_index, int& numofextreme, int m, std::vector<int>& mi, int d){	
  Polytope Res;
  Field vec[]={-1,0,1};
  int D=d, base=sizeof(vec)/sizeof(Field);
   
  for (int i=0; i<pow(base,PD); ++i){
    vector<Field> extreme_point;
    vector<int> v=int2vectord(i,base,D);
    //copy(v.begin(),v.end(),ostream_iterator<int>(std::cout," "));
    for (vector<int>::iterator it=v.begin(); it!=v.end(); it++){
		  extreme_point.push_back(vec[*it]);
		  //std::cout << vec[*it] << " ";
	  }
	  //std::cout << std::endl;
    if (!is_zero(extreme_point)){
      PVector_d lft(PD,extreme_point.begin(),extreme_point.end());
      //CVector_d lft(3,1,1,2);
      //std::cout << lft << lft.dimension() << std::endl;
      SimplicialComplex t = LiftingTriangulation(points,lft,chiro,r,points_index,false);
      //std::cout << t << std::endl;
      if (!t.is_empty()){//if t is indeed a triangulation!  
        SRvertex new_vertex = compute_r(t, points, m, mi, d);
		    std::cout << new_vertex << std::endl;
		    Res.insert(new_vertex);
		  }
    }
  }
	return Res;
}
*/
/* FOR DEBUGGING!!!
// iterate through all facets and put all the points of upper (lower) hull facets to triang
// the upper hull is the default
SimplicialComplex project_lower_upper_hull(CTriangulation& pc,
                                           map<vector<Field>, int>& points_index,
                                           Chirotope& chiro,
                                           int dcur, 
                                           int dim, 
                                           bool lower){
	//std::cout << dim << "|" << pc.current_dimension() << std::endl;
	SimplicialComplex triang, triang2;
	typedef std::vector<CSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  pc.gather_incident_simplices(pc.infinite_vertex(), out);
	//CPoint_d some_vertex_point = (*(inf_simplices.begin()))->vertex(1)->point();
	for( Simplices::iterator sit = inf_simplices.begin(); 
	                         sit != inf_simplices.end(); ++sit ){
    CFacet ft(*sit, 0); // |ft| is a facet of the convex hull
    std::cout << "dim=" << pc.current_dimension() << "facet points:" << std::endl;
    std::vector<CPoint_d> facet_points;
    for (int i=1; i<=pc.current_dimension(); i++){
			facet_points.push_back((*sit)->vertex(i)->point());
			std::cout << (*sit)->vertex(i)->point() << " ";
		}
		std::cout << std::endl;
		
		std::vector<CPoint_d> facet_vectors;
		for (std::vector<CPoint_d>::iterator vit=facet_points.begin()+1; vit!=facet_points.end();
		     vit++){
			//std::cout << CPoint_d(vit->cartesian_begin(),vit->cartesian_end())-facet_points[0]
			//          << ";;" << CPoint_d(vit->cartesian_begin(),vit->cartesian_end());
		  facet_vectors.push_back(CPoint_d(vit->cartesian_begin(),vit->cartesian_end())-facet_points[0]);
		}
		std::cout << "det";
		std::vector<Field> cross_product;
		Matrix matr(facet_points.size()-1,facet_vectors.size());
		for (int k=0; k<facet_points.size(); k++){
			for (int i=0; i<facet_vectors.size(); i++){
	 			for (int j=0, jj=0; j<facet_points.size(); j++){
					if (j!=k){
					  matr(jj,i)=facet_vectors[i][j];
					  jj++;
					}
				}
			}
			if (k%2==0){  
			  std::cout << " " << det(matr);
				cross_product.push_back(det(matr));
			} else {
				std::cout << " " << -det(matr);
				cross_product.push_back(-det(matr));
			}
		}
		CVector_d det_vec(cross_product.begin(),cross_product.end());
		std::cout << "dot_products:: ";
		for (CVertex_iterator vit=pc.vertices_begin(); vit!=pc.vertices_end(); vit++){
			Field dot_product=0;
			for (int i=0; i<CD; i++)
				dot_product += vit->point().cartesian(i)*det_vec.cartesian(i);
			
			std::cout << "(" << vit->point() << ")" << dot_product << " ";
		} 
		std::cout << std::endl;
		CVector_d det_vec2 = -det_vec;
		std::cout << "dot_products2:: ";
		for (CVertex_iterator vit=pc.vertices_begin(); vit!=pc.vertices_end(); vit++){
			Field dot_product=0;
			for (int i=0; i<CD; i++)
				dot_product += vit->point().cartesian(i)*det_vec2.cartesian(i);
			
			std::cout << "(" << vit->point() << ")" << dot_product << " ";
		} 
		std::cout << std::endl;
		
		CPoint_d opposite_point = (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
    std::cout << "opposite_point=" << opposite_point;
    
    CHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
    CHyperplane_d hp2(facet_points[0],det_vec.direction());
    
    std::cout << " hp2=" << hp2.has_on_positive_side(opposite_point) << " hp=" << hp.has_on_positive_side(opposite_point);
    //std::cout << "dim=" << pc.current_dimension() << " " << hp.orthogonal_direction()<< "\n";
    CVector_d vec = hp.orthogonal_direction().vector();
    bool is_upper = vec.cartesian(pc.current_dimension()-1) > 0;
    //bool is_upper = det_vec.cartesian(pc.current_dimension()-1) > 0;
    //if (hp2.has_on_positive_side(opposite_point))
    //	is_upper =!is_upper;
    if (lower) is_upper=!is_upper;
    
    std::cout << " vec=" << vec << std::endl;
    std::cout << "is_upper=" << is_upper << std::endl;
    
    std::vector<Field> p_temp;
    for (int i=0; i<pc.current_dimension(); i++){
			p_temp.push_back(vec.cartesian(i)+facet_points[0].cartesian(i));
    }
    CPoint_d p(p_temp.begin(),p_temp.end());
    //std::cout << p;
    //std::cout << "##" << hp.has_on_positive_side(p) << "!!";
    
    CLocate_type loc_type;
  	CFace f(CD);
		CFacet ft2;
  	pc.locate(opposite_point,loc_type,f,ft2);
    //std::cout << "loc_type=" << loc_type << std::endl;
    //std::cout << vec << "|" << vec.cartesian(pc.current_dimension()-1) << "|" << pc.current_dimension();
    //CHyperplane_d hp2(facet_points[0],hp.orthogonal_direction());
    //std::cout << "$" << hp2.has_on_positive_side(opposite_point) << "%" << hp.orthogonal_direction() << "|" << hp.orthogonal_direction().vector() << "%";
    IntegerSet s;
		for (std::vector<CPoint_d>::iterator vit=facet_points.begin();
		      vit!=facet_points.end(); ++vit){
			std::vector<Field> ppoi;
			for (int j=0; j<dim-1; ++j) // d-1 cuts the last coordinate!!
				ppoi.push_back((*vit)[j]);			
			//std::cout << *vit << ",";
			s += points_index[ppoi];
		}
    if (is_upper && s.card()==dim){
			triang.insert(s);
		}
		else if (is_upper && s.card()<dim){
			#ifdef PRINT_INFO  
			//  std::cout << "found a low dimensional simplex. possible ERROR\n" << std::endl;
			#endif
		}
	}
	//std::cout << triang << std::endl;
  //std::cout << triang2 << std::endl;
  return triang;
}
*/
