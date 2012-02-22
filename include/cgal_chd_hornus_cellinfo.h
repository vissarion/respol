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

#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

// include d-dimensional things from CGAL
//#include <CGAL/Filtered_kernel_d.h>
//#include <CGAL/Pure_complex.h>
//#include <CGAL/Delaunay_complex.h>
//#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation.h>
#include <CGAL/Triangulation_vertex.h>
#include <CGAL/Triangulation_full_cell.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_d.h>

// for indexed points
//#include "indexed_point.h"

// the kernel
//#include <gmpxx.h>
//typedef mpq_class Field;
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq Field;


// Cayley space kernel; for the CH of the lifted points
typedef CGAL::Cartesian_d<Field>		CK;

//typedef CGAL::Filtered_kernel_d<Field> 						CK;
//typedef Indexed_Cartesian_d<Field> 							CK;
//typedef CK::IndexedPoint_d									IndexedPoint_d;

//typedef CGAL::Triangulation_vertex<CK,size_t>						tv;
//typedef CGAL::Triangulation_full_cell<CK>		tc;
//typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<PD>,tv,tc > tds;
// not working if we'll not pas all the parameters explicitly
//e.g. CGAL::Triangulation_data_structure<CGAL::Dimension_tag<CD> > tds; gives error
//typedef CGAL::Triangulation<CK,tds> 							CTriangulation;
typedef CGAL::Triangulation<CK>			 						CTriangulation;
typedef CTriangulation::Point_d 								CPoint_d;
typedef CK::Hyperplane_d										CHyperplane_d;
typedef CK::Vector_d											CVector_d;
typedef CTriangulation::Vertex_iterator							CVertex_iterator;
typedef CTriangulation::Vertex_handle							CVertex_handle;
typedef CTriangulation::Full_cell_iterator						CSimplex_iterator;
typedef CTriangulation::Full_cell_const_iterator				CSimplex_const_iterator;
typedef CTriangulation::Full_cell_handle						CSimplex_d;
typedef CTriangulation::Facet									CFacet;
typedef CTriangulation::Face 									CFace;
typedef CTriangulation::Facet 									CFacet;
typedef CTriangulation::Locate_type 							CLocate_type;

#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY
// extreme points d includes
#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>
//extreme points typedefs
typedef CGAL::Extreme_points_traits_d<CPoint_d>   EP_Traits_d;
#endif

// projection kernel; for the Resultant polytope

typedef CGAL::Cartesian_d<Field>									PK;

typedef CGAL::Triangulation_vertex<CK>								tv;
typedef CGAL::Triangulation_full_cell<CK,bool>						tc;
typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag,tv,tc > tds;
//typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<12>,tv,tc > tds;
typedef CGAL::Triangulation<PK,tds> 							Triangulation;
//typedef CGAL::Delaunay_triangulation<PK,tds> 							Triangulation;
//typedef CGAL::Triangulation<PK> 								Triangulation;

typedef PK::Direction_d											PVector_d;
typedef PK::Vector_d											  Vector_d;
typedef PK::Hyperplane_d										PHyperplane_d;
typedef Triangulation::Full_cell_handle							PSimplex_d;
typedef Triangulation::Finite_full_cell_iterator				PSimplex_finite_iterator;
typedef Triangulation::Finite_full_cell_const_iterator  		PSimplex_finite_const_iterator;
typedef Triangulation::Full_cell_iterator						PSimplex_iterator;
typedef Triangulation::Point_d 									PPoint_d;
typedef Triangulation::Face 									PFace;
typedef Triangulation::Facet 									PFacet;
typedef Triangulation::Facet_iterator							PFacet_iterator;
typedef Triangulation::Locate_type 								PLocate_type;
typedef Triangulation::Vertex_handle							PVertex_handle;
typedef Triangulation::Vertex_iterator							PVertex_iterator;

// some global variables for experiments 
// should be excluded
double conv_time = 0;
#ifdef RESTRICTED_RES
int restricted_num_Res;
#endif
//misc typedefs

typedef std::vector<Field>      SRvertex;
typedef std::vector<SRvertex>   Resvertex;
//typedef IntegerSet 			SRvertex;
typedef std::set<SRvertex>      Polytope;

// for fast determinant computation
#include "hashed_determinant.h"
typedef HashedDeterminant<Field>::Table                 HD;

// normal vectors data structure
//#include <../include/normal_vector_ds.h>
#include "normal_set_ds.h"

typedef Normal_Vector_ds<PVector_d,Field>					NV_ds;

// functions for printing and parsing
#include "print_functions.h"
#include "parse_functions.h"

// lrslib includes
#ifdef USE_LRSLIB
#include "lrs_cgal.h"
#endif

/////////////////////////////
//functions implementations

/////////////////////////////////////////////////////////////////
// math
//template <class Field>
//inline Field factorial(Field n)
//{
//  return (n == Field(1) || n == Field(0)) ? Field(1) : factorial(n - 1) * n;
//}

template <class Field>
inline Field factorial(Field n)
{
  return 1;
}

// apply cayley trick 
template <class NT_>
void cayley_trick(std::vector<std::vector<NT_> >& pointset,
									std::vector<int> mi){
  typedef NT_                                           Field;
	// construct the cayley vector, zeroes-ones
	std::vector<std::vector<Field> > cayley_vec;
	int i=0, j=-1, end=0;
	for (std::vector<int>::iterator vit=mi.begin(); vit!=mi.end(); vit++){
		end += *vit;
		while (i<end){
			std::vector<Field> cayley_point(D,0);
			if (j>=0)
				cayley_point[j] = 1;
			cayley_vec.push_back(cayley_point);
			i++;
		}
		//std::cout << cayley_vec << std::endl;
		j++;
	}
	//std::cout << cayley_vec << std::endl;
	// append cayley vectors to points of pointset
	//int p_index=0;
	typename std::vector<std::vector<Field> >::iterator cit=cayley_vec.begin();
  for(typename std::vector<std::vector<Field> >::iterator vit=pointset.begin();
      vit!=pointset.end();
      vit++,cit++){
		vit->insert(vit->end(),cit->begin(),cit->end());
		//points_index[*vit] = p_index++;
		//std::cout << *vit << "-->" << points_index[*vit] << std::endl;
	}
	// write the (homogenized) poinset to file for topcom test
	std::ofstream outfile;
  outfile.open("topcom_cayley.txt");
	print_vertices_hom(pointset,outfile);
} 

#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY
// compute the extreme points
template <class NT_>
void compute_extreme_points(std::vector<std::vector<NT_> >& pointset,
														std::vector<int>& mi,
														std::vector<int>& proj,
                            const ResPol::config &conf){
  typedef NT_                                           Field;
  std::vector<std::vector<Field> > extreme_pointset;
  std::vector<int> extreme_mi;
  // the nonspecialized-projected points will be first 
  // thus we have to change the projection
  std::vector<int> proj_updated;
	//std::cout<<pointset<<std::endl;
	//std::cout << mi << std::endl;
	int current_start=0, current_end=0;
	for(std::vector<int>::iterator imi=mi.begin(); imi!=mi.end(); ++imi){
		current_end+=(*imi);
		//std::cout << current_start << "-" << current_end << std::endl;
		std::vector<CPoint_d> points;
		// vector<vector<Field> >  -->  vector<CPoint_d>
		int specialized_points_position = extreme_pointset.size();
		int count_specialized_points=0;
		for (int i=current_start;i<current_end;++i){
			// if the point in spesialized, i.e not projected
			if (std::find(proj.begin(), proj.end(), i)==proj.end()){
				CPoint_d p(D,pointset[i].begin(),pointset[i].end());
				points.push_back(p);
        if(conf.verbose>1)
				  std::cout << p << "\n";
			} else {
				proj_updated.push_back(specialized_points_position);
				extreme_pointset.push_back(pointset[i]);
				++specialized_points_position;
				++count_specialized_points;
			}
		}
		//std::cout << std::endl;
		
		// compute the extreme points
		CGAL::Extreme_points_d<EP_Traits_d> ep(D);
		ep.insert(points.begin(), points.end());
		std::vector<CPoint_d> extreme_points;
		ep.get_extreme_points(std::back_inserter(extreme_points));
		
		//vector<CPoint_d>  -->  vector<vector<Field> >
		for (std::vector<CPoint_d>::iterator it=extreme_points.begin();
				 it!=extreme_points.end();
				 it++) {
			//std::cout << *it << std::endl;
			std::vector<Field> extreme_pointset_value;
			for (int i=0; i<D; ++i)
				extreme_pointset_value.push_back(it->cartesian(i));
			extreme_pointset.push_back(extreme_pointset_value);
		}
    if(conf.verbose>1){
      std::cout<<"\nExtreme points=" << extreme_points.size()
        << "(out of" << points.size()
        << ")\n"<<std::endl;
    }
		extreme_mi.push_back(extreme_points.size()+count_specialized_points);
		current_start=current_end;
	}
  if(conf.verbose>1){
    std::cout << "mi= " << mi << " extreme_mi=" << extreme_mi << std::endl;
    std::cout << "proj= " << proj << " proj_updated=" << proj_updated <<
      std::endl;
  }
	// update pointset, mi, proj
	pointset = extreme_pointset;
	mi = extreme_mi;
	proj = proj_updated;
}
#endif

#ifdef USE_EXTREME_SPECIALIZED_POINTS_ONLY
template <class Triang>
int count_extreme_vertices(const Triang &Res){
	typedef typename Triang::Vertex_const_iterator        VCI;
  typedef typename Triang::Point_d                      P;
  typedef typename P::Cartesian_const_iterator          PCCI;
  std::vector<PPoint_d> points;
  for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
    points.push_back(vit->point());
  }
	// compute the extreme points
	CGAL::Extreme_points_d<EP_Traits_d> ep(Res.current_dimension());
	ep.insert(points.begin(),points.end());
	std::vector<CPoint_d> extreme_points;
	ep.get_extreme_points(std::back_inserter(extreme_points));
	return extreme_points.size();
}
#endif

///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//TODO use this to have only CPoint_d and not vector<vector<Field> >
//transform vector<vector<Field> > to a vector of CPoint_d
template <class NT_>
void transform_pointset(const std::vector<std::vector<NT_> >& pointset,
												std::vector<CPoint_d>& cgal_pointset,
												const HD& dets){
    typedef NT_                                         Field;
		std::cout << pointset << std::endl;
		for (typename std::vector<std::vector<Field> >::const_iterator
         vit=pointset.begin(); vit!=pointset.end(); vit++){
			CPoint_d p(CD,vit->begin(),vit->end());
			//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
			//p.set_index(vit->second);
			//p.set_hash(&dets);
			p[0]=Field(-1);
			std::cout << p << "\n ";
		}
		CPoint_d p(CD);
		p[1]=1;
		std::cout << "*" << p<< std::endl;
		exit(1);
	}
	
/*
std::vector<int> proj_more_coord(const int d,
                                 int m,
                                 const std::vector<int>& mi){
	//project at the first coordinate of each mi
	std::vector<int> proj(d);
	proj[0]=0;
	proj[1]=1;
  int k=2;
	for (int i=1; i<3; i++){
		int mm=0;
		for (int j=0; j<i; j++)
			mm+=mi[j];
		proj[k++]=mm;
		proj[k++]=mm+1;
		//std::cout << mm << " ";
	}
	return proj;
}
*/

std::vector<int> full_proj(const int d, int m, const std::vector<int>& mi){
	std::vector<int> proj(d);
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

template <class NT_>
std::vector<std::pair<std::vector<NT_>,size_t> > lift_to_zero(
        const std::vector<std::vector<NT_> >& points,
        const std::vector<int>& proj){
  typedef NT_                                           Field;
  std::vector<std::pair<std::vector<Field>,size_t> > lifted_points;
  std::vector<int>::const_iterator lpit=proj.begin();
  for (typename std::vector<std::vector<Field> >::const_iterator pit=
         points.begin();
       pit!=points.end();
       pit++){
  	int point_index = pit-points.begin();
  	if (point_index == *lpit){
  		lpit++;
  	} else {
	  	std::vector<Field> lifted_point = *pit;
	  	lifted_point.push_back(Field(0));
	  	std::pair<std::vector<Field>,size_t>
        lifted_point_with_index(lifted_point,point_index);
	  	lifted_points.push_back(lifted_point_with_index);
	  }
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute the Triangulation of the static points (lifted to zero)
// the points to be lifted will not be present here

template <class NT_>
void StaticTriangulation(const std::vector<std::vector<NT_> >& points,
												 const std::vector<int>& proj,
												 CTriangulation& T,
												 const HD& dets){
  typedef NT_                                           Field;
	// lift the "static" points Cayley pointset
	// (i.e. the points that are not going to be projected)
  std::vector<std::pair<std::vector<Field>,size_t> > P =
    lift_to_zero(points,proj);

  //std::cout << "lifted points:(dim=" << dim << "\n" << P << std::endl;
  for (typename
         std::vector<std::pair<std::vector<Field>,size_t> >::const_iterator
         vit=P.begin();
       vit!=P.end();
       vit++){
		CPoint_d p(CD,vit->first.begin(),vit->first.end());
		//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
  	p.set_index(vit->second);
  	p.set_hash(&dets);
  	CVertex_handle v = T.insert(p);
  }
}


// lift the points that are going to be projected using nli

template <class NT_>
std::vector<std::pair<std::vector<NT_>,size_t> > lift_to_proj(
        const std::vector<std::vector<NT_> >& points,
        const PVector_d& nli,
        const std::vector<int>& proj){

  typedef NT_                                           Field;
  std::vector<std::pair<std::vector<Field>,size_t> > lifted_points;
  //int i=0;
  std::vector<int>::const_iterator proj_it=proj.begin();
  for (typename std::vector<std::vector<Field> >::const_iterator vit=
         points.begin();
       vit!=points.end() && proj_it!=proj.end();
       vit++){
  	std::vector<Field> lifted_point = *vit;
  	int pindex=vit-points.begin();
  	//std::cout << vit-points.begin() << "," << proj_it-proj.begin() << std::endl;
  	if (*proj_it == pindex){
			//std::cout << nli << "|"<<proj<<std::endl;
  		lifted_point.push_back(nli[proj_it-proj.begin()]);
  		std::pair<std::vector<Field>,size_t>
        lifted_point_with_index(lifted_point,pindex);
	  	lifted_points.push_back(lifted_point_with_index);
	  	proj_it++;
  	}
  }
  //std::cout << lifted_points << std::endl;
  return lifted_points;
}


// compute lifting triangulation but not from scratch
// use a copy of T and add new lifted points there
// TODO: not insert points with negative lift
// they will not affect upper hull! 
// EDIT: this is not true!! (an einai sthn akrh?) 
// be careful
// if all points have neg lift then Tl will be
// 2n dimensional so we have to project carefully
template <class NT_>
void LiftingTriangulationDynamic(
        const std::vector<std::vector<NT_> >& points,
        const PVector_d& nli,
        const std::vector<int>& proj,
        const CTriangulation& T,
        CTriangulation& Tl,
        const HD& dets,
        bool lower){

  typedef NT_                                           Field;
  // lift the points to be projected
  std::vector<std::pair<std::vector<Field>,size_t> > P =
    lift_to_proj(points,nli,proj);

  // add the new lifted points to the copy (Tl)
  for (typename
         std::vector<std::pair<std::vector<Field>,size_t> >::const_iterator
         vit=P.begin();
       vit!=P.end();
       vit++){
		//if (*(vit->first.end()-1) >= 0){
			CPoint_d p(CD,vit->first.begin(),vit->first.end());
	  	p.set_index(vit->second);
	  	p.set_hash(&dets);
	  	//std::cout << p <<"|"<< p.index() <<" point inserted" << std::endl;
	  	CVertex_handle v = Tl.insert(p);
		//}
  }
}


// iterate through all facets determine the upper (i.e. compute determinants),
// project and compute the r_vector (i.e. compute other determinants)
// TODO: the second determinants are minors of the first

std::vector<Field> project_upper_hull_r(const CTriangulation& pc,
                                        HD& dets,
                                        int dcur,
                                        int dim,
                                        const std::vector<int>& mi,
                                        const std::vector<int>& proj,
                                        bool lower,
                                        const ResPol::config &conf){
	// the result vector
	std::vector<Field> rho(PD,0);
	// compute the infinite simplices
	typedef std::vector<CSimplex_d> Simplices;
	Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  pc.incident_full_cells(pc.infinite_vertex(), out);

	//iterate through all infinite simplices
	for(Simplices::const_iterator sit = inf_simplices.begin();
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
			std::set<int> s;
			for (std::vector<CPoint_d>::iterator vit=facet_points.begin();
			      vit!=facet_points.end(); ++vit){
				s.insert(vit->index());
		  }
		 // std::cout<<s<<std::endl;
		  // if there are not full dimensional cells
		  if (s.size() < CD){
        if(conf.verbose>1){
          std::cout <<
            "NOT a triangulation. Abort current computation..." << std::endl;
        }
				/*std::vector<Field> rho_empty;
				//exit(-1);
				return rho_empty;*/
        return std::vector<Field>();
			}
      //secondary case
      if(conf.polytope_type==1){
        for(std::set<int>::const_iterator it2=s.begin(); it2!=s.end(); it2++){
          std::vector<int>::const_iterator pit=
            find(proj.begin(),proj.end(),*it2);
          if (pit != proj.end()){
            rho[pit-proj.begin()] += det;
          }
        }
      }else if(conf.polytope_type==0){
        //resultant case
        // check if the projection of the facet (i.e. a Minkowski cell)
        // is mixed (i.e. has exactly one vertex summand)
        int d = mi.size();
        std::vector<std::vector<int> > sets(d);
        int current_set=0, current_m=mi[0];
        for (std::set<int>::const_iterator it2=s.begin(); it2!=s.end(); it2++){
          if (*it2 >= current_m)
            current_m += mi[++current_set];
          sets[current_set].push_back(*it2);
        }
        std::vector<int> mixed_vertices;
        for (int i=0; i<d; i++){
          if (sets[i].size() == 1){
            mixed_vertices.push_back(i);
          }
        }
        if (mixed_vertices.size() == 1){
          // find where is the mixed vertex in the vector of projection
          // coordinates
          std::vector<int>::const_iterator pit =
                  find(proj.begin(),proj.end(),sets[mixed_vertices[0]][0]);
          if (pit != proj.end()){
            rho[pit-proj.begin()] += det;
          }
        }
      }else{
        std::cout<<"unknown polytope type"<<std::endl;
        exit(-1);
      }
    }
	}
	// to compute the real volume we have to divide by d!
	for (std::vector<Field>::iterator vit=rho.begin();
       vit!=rho.end();
       vit++){
		*vit=*vit/factorial(pc.current_dimension()-1);
	}
  return rho;
}

//////////////////////////////////////////////////////////////
// insert vertices in triangulation that holds the Resultant
// polytope

// this is for init
void insert_new_Rvertex(Triangulation& Res,
												const NV_ds& normal_list,
												const SRvertex& new_vertex,
												HD& Pdets,
                        const ResPol::config &conf){
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(Res.number_of_vertices());
	new_point.set_hash(&Pdets);
  if(conf.verbose>1){
    //std::cout << "one new R-vertex found !!! "<< std::endl;
  }
	//insert it to the triangulation
	PVertex_handle new_vert = Res.insert(new_point);
}


void update_cell_data(const Triangulation& Res,
											const PVertex_handle &vert){

	// construct the edge (new_vert,inf_vert)
	PFace edge(vert->full_cell());
	edge.set_index(0, vert->full_cell()->index(vert));
	edge.set_index(1, 0);

	typedef std::vector<PSimplex_d> Simplices;
	Simplices new_simplices;
	std::back_insert_iterator<Simplices> out(new_simplices);

	Res.incident_full_cells(edge, out);

	for(Simplices::iterator sit = new_simplices.begin();
													sit != new_simplices.end(); ++sit ){
		(*sit)->data() = false;
	}
}

// this is for augment
void insert_new_Rvertex2(Triangulation& Res,
												const SRvertex& new_vertex,
												HD& Pdets,
												PSimplex_d& near_cell,
                        const ResPol::config &conf){
	// insert the coordinates of the point as a column to the HashDeterminants matrix
	Pdets.add_column(new_vertex);
	// construct the new point
	PPoint_d new_point(PD,new_vertex.begin(),new_vertex.end());
	new_point.set_index(Res.number_of_vertices());
	new_point.set_hash(&Pdets);
  if(conf.verbose>1){
    //std::cout << "one new R-vertex found !!! "<< std::endl;
  }
	int prev_dim = Res.current_dimension();
	//insert it to the triangulation
	double tstartall = (double)clock()/(double)CLOCKS_PER_SEC;
	PVertex_handle new_vert = Res.insert(new_point,near_cell);
  double tstopall = (double)clock()/(double)CLOCKS_PER_SEC;
	conv_time += tstopall - tstartall;
  if(conf.verbose>1){
    std::cout << "insert in Res time= " <<  tstopall - tstartall << std::endl;
  }
	int cur_dim = Res.current_dimension();
	update_cell_data(Res,new_vert);
}

//////////////////////////////////////////////////////////

// TODO: make it more efficient, don't gather all infinite cells
int get_illegal_facet(Triangulation& Res,
											NV_ds& normal_list,
											PVector_d& current_vector,
											PSimplex_d& near_cell){
	//typedef std::vector<PSimplex_d> Simplices;
	//Simplices inf_simplices;
	//std::back_insert_iterator<Simplices> out(inf_simplices);
	//Res.incident_full_cells(Res.infinite_vertex(), out);
	PSimplex_iterator sit = Res.full_cells_begin();
	//bool is_neg;
	do{
		// find an illegal facet
		for( ; sit != Res.full_cells_end(); ++sit ){
			if (sit->data() == false && Res.is_infinite(sit)){
				sit->data() = true;
				break;
			}
		}
		near_cell = sit;
		//print_cells_data(Res);
		if (sit != Res.full_cells_end()){
			std::vector<PPoint_d> facet_points;
			for (int i=1; i<=Res.current_dimension(); i++){
				facet_points.push_back(sit->vertex(i)->point());
				//std::cout << (*sit)->vertex(i)->point() << " ";
			}
			//std::cout << std::endl;
			PPoint_d opposite_point = sit->neighbor(0)->vertex(sit->mirror_index(0))->point();
			// compute a hyperplane which has in its negative side the opposite point
			PHyperplane_d hp(facet_points.begin(),facet_points.end(),opposite_point,CGAL::ON_NEGATIVE_SIDE);
			//is_neg=true;
			current_vector = hp.orthogonal_direction();
			//for (PVector_d::Delta_const_iterator dit=current_vector.deltas_begin();
			//						dit!=current_vector.deltas_end(); dit++){
			//	if (*dit >=0 ){
			//		is_neg = false;
			//		break;
			//	}
			//}
			//std::cout << normal_list.put(current_vector) << std::endl;
			//return 1;
		}
	} while(normal_list.put(current_vector) == 0 //|| is_neg==true)
					 && sit != Res.full_cells_end());
	// check if we run out of normals
	if (sit != Res.full_cells_end())
		return 1;
	else
		return 0;
}


// compute a Res vertex by constructing a lifting triangulation
// project and compute the volumes of some cells

template <class NT_>
std::vector<NT_> compute_res_vertex(
        const std::vector<std::vector<NT_> >& pointset,
        const std::vector<int>& mi,
        int RD,
        const std::vector<int>& proj,
        HD& dets,
        const HD& Pdets,
        const Triangulation& Res,
        const CTriangulation& T,
        NV_ds& normal_list_d,
        const ResPol::config &conf){

  typedef NT_                                           Field;
	// take the first(last more efficient with vectors??) normal and remove it from normals stack
	//PVector_d nli = *(normal_list_d.begin());
	//normal_list_d.erase(normal_list_d.begin());
	PVector_d nli = normal_list_d.back();
	normal_list_d.pop_back();
	//std::cout << "normal vector:" << nli << std::endl;

	if (pointset.size() == proj.size()){
		// make a new T
		CTriangulation Tl(CD);
		// outer normal vector --> upper hull projection
	  LiftingTriangulationDynamic(pointset,nli,proj,T,Tl,dets,false);
		//print_res_vertices(Tl);
	
		// project Tl i.e. triangulation
	  int dcur = Tl.current_dimension();
	  std::vector<Field> new_vertex =
	    project_upper_hull_r(Tl,dets,dcur,CD,mi,proj,false,conf);
	
	  // TODO: destroy here!
	  Tl.clear();
	
    if(conf.verbose>1){
      std::cout << "new Res vertex (up)= ( " << new_vertex << ")" << std::endl;
    }
	  //std::cout << new_vertex << std::endl;
	  return new_vertex;	
	} else {
		// make a copy of T
		CTriangulation Tl(T);
		// outer normal vector --> upper hull projection
	  LiftingTriangulationDynamic(pointset,nli,proj,T,Tl,dets,false);
		//print_res_vertices(Tl);
	
		// project Tl i.e. triangulation
	  int dcur = Tl.current_dimension();
	  std::vector<Field> new_vertex =
	    project_upper_hull_r(Tl,dets,dcur,CD,mi,proj,false,conf);
	
	  // TODO: destroy here!
	  Tl.clear();
	
    if(conf.verbose>1){
      std::cout << "new Res vertex (up)= ( " << new_vertex << ")" << std::endl;
    }
	  //std::cout << new_vertex << std::endl;
	  return new_vertex;
	}
}

template <class NT_>
std::vector<NT_> compute_res_vertex2(
				                 const std::vector<std::vector<NT_> >& pointset,
				                 const std::vector<int>& mi,
				                 int RD,
				                 const std::vector<int>& proj,
				                 HD& dets,
				                 const HD& Pdets,
				                 const Triangulation& Res,
				                 const CTriangulation& T,
				                 const PVector_d &nli,
                         const ResPol::config &conf){

  typedef NT_                                           Field;
	// take the first(last more efficient with vectors??) normal and remove it from normals stack
	//PVector_d nli = *(normal_list_d.begin());
	//normal_list_d.erase(normal_list_d.begin());
	//PVector_d nli = normal_list_d.back();
	//normal_list_d.pop_back();
	//std::cout << "normal vector:" << nli << std::endl;

	if (pointset.size() == proj.size()){
		// make a new T
		CTriangulation Tl(CD);
		// outer normal vector --> upper hull projection
	  LiftingTriangulationDynamic(pointset,nli,proj,T,Tl,dets,false);
		//print_res_vertices(Tl);
	
		// project Tl i.e. triangulation
	  int dcur = Tl.current_dimension();
	  std::vector<Field> new_vertex =
	    project_upper_hull_r(Tl,dets,dcur,CD,mi,proj,false,conf);
	
	  // TODO: destroy here!
	  Tl.clear();
	
    if(conf.verbose>1){
      std::cout << "new Res vertex (up)= ( " << new_vertex << ")\n";
    }
	  //std::cout << new_vertex << std::endl;
	  return new_vertex;	
	} else {
		// make a copy of T
		CTriangulation Tl(T);
		// outer normal vector --> upper hull projection
	  LiftingTriangulationDynamic(pointset,nli,proj,T,Tl,dets,false);
		//print_res_vertices(Tl);
	
		// project Tl i.e. triangulation
	  int dcur = Tl.current_dimension();
	  std::vector<Field> new_vertex =
	    project_upper_hull_r(Tl,dets,dcur,CD,mi,proj,false,conf);
	
	  // TODO: destroy here!
	  Tl.clear();
	
    if(conf.verbose>1){
      std::cout << "new Res vertex (up)= ( " << new_vertex << ")\n";
    }
	  //std::cout << new_vertex << std::endl;
	  return new_vertex;
	}
	
}

// compute Res vertices until it builts a simplex

template <class NT_>
int initialize_Res(const std::vector<std::vector<NT_> >& pointset,
								 const std::vector<int>& mi,
								 int RD,
								 const std::vector<int>& proj,
								 HD& dets,
								 HD& Pdets,
								 Triangulation& Res,
								 const CTriangulation& T,
                 const ResPol::config &conf){

  typedef NT_                                           Field;
	int num_of_triangs=0;
  if(conf.verbose>1){
    std::cout << "dim=" << Res.current_dimension() << std::endl;
  }
  // make a stack (stl vector) with normals vectors and initialize
  NV_ds normal_list_d;
  //normal_list_d.simple_initialize();
  normal_list_d.initialize();
  
  // compute trinagulations using normals as liftings until we compute a simplex
  // or run out of normal vectors
  int minD = (PD>RD)?RD:PD;
  //std::cout << "minD:"<<minD<<"PD:"<<PD<<std::endl;
  while(Res.current_dimension()<minD && !normal_list_d.empty()){
    std::vector<Field> new_vertex=compute_res_vertex(
      pointset,mi,RD,proj,dets,Pdets,Res,T,normal_list_d,conf);
		// insert it in the complex Res (if it is not already there)
  	if (Pdets.find(new_vertex) == -1 && new_vertex.size() != 0)
  		insert_new_Rvertex(Res,normal_list_d,new_vertex,Pdets,conf);
		num_of_triangs++;
    if(conf.verbose>1){
      normal_list_d.print();
      std::cout<<"current number of Res vertices: "
        <<Res.number_of_vertices()
        <<" dim=" << Res.current_dimension()<<std::endl;
    }
	}
	return num_of_triangs;
}

// compute all Res vertices left by augmenting the simplex
// to the directions of the normal vectors of the facets

template <class NT_>
int augment_Res(const std::vector<std::vector<NT_> >& pointset,
								 const std::vector<int>& mi,
								 int RD,
								 const std::vector<int>& proj,
								 HD& dets,
								 HD& Pdets,
								 Triangulation& Res,
								 const CTriangulation& T,
                 const ResPol::config &conf){

  typedef NT_                                           Field;
	int num_of_triangs=0;
  if(conf.verbose>1){
    std::cout << "\n\nAUGMENTING RESULTANT POLYTOPE" << std::endl;
    //print_res_vertices(Res);
    std::cout << "dim=" << Res.current_dimension() << std::endl;
  }

  update_cell_data(Res,Res.infinite_vertex());
  //print_cells_data(Res);
  PVector_d current_vector;
  // make a stack (stl vector) with normals vectors used
  NV_ds normal_list_d;
	//hint cell
	PSimplex_d near_cell;
  clock_t t_start;
  while(get_illegal_facet(Res,normal_list_d,current_vector,near_cell)){
    if(conf.verbose>1){
      std::cout << "\nAUGmenting step " << num_of_triangs << std::endl;
      std::cout << "current normal" << current_vector << std::endl;
      t_start=clock();
    }
		std::vector<Field> new_vertex =
      compute_res_vertex2(
        pointset,mi,RD,proj,dets,Pdets,Res,T,current_vector,conf);
    if(conf.verbose>1){
      std::cout << "Compute triang time= "<<
        ((double)(clock()-t_start))/CLOCKS_PER_SEC<<std::endl;
    }
		
		// insert it in the complex Res (if it is not already there)
  	if (Pdets.find(new_vertex) == -1 && new_vertex.size() != 0){	
  		insert_new_Rvertex2(Res,new_vertex,Pdets,near_cell,conf);
		}
    if(conf.verbose>1){
      {
        int cells, triang_facets, facets, edges, vertices;
        f_vector(Res,cells,triang_facets,facets,edges,vertices,conf);
        std::cout<<"("<<cells<<","<<triang_facets<<","
          <<facets<<","<<edges<<","<<vertices<<"\n";
      }
    }
		#ifdef RESTRICTED_RES
			if (Res.number_of_vertices() == restricted_num_Res){
				return num_of_triangs;
			}
		#endif
		//if (Res.number_of_vertices() == 40){
		//		return num_of_triangs;
		//	}
		++num_of_triangs;
	}
	return num_of_triangs;
}

////////////////////////////////////////////////////////////
// the algorithm
template <class NT_>
std::pair<int,int> compute_res(
        const std::vector<std::vector<NT_> >& pointset,
        int m,
        const std::vector<int>& mi,
        int RD,
        const std::vector<int>& proj,
        HD& dets,
        HD& Pdets,
        Triangulation& Res,
        const ResPol::config &conf){

  typedef NT_                                           Field;
	//std::cout << "cayley dim:" << CayleyTriangulation(pointset) << std::endl;

  // construct an initial triangulation of the points that will not be projected
  CTriangulation T(CD);
  StaticTriangulation(pointset,proj,T,dets);
	//std::cout << "static dim:" << T.current_dimension() << std::endl;
	
  // start by computing a simplex
  int start_triangs=
    initialize_Res(pointset,mi,RD,proj,dets,Pdets,Res,T,conf);

  // augment simplex to compute the Res polytope
  int augment_triangs=
    augment_Res(pointset,mi,RD,proj,dets,Pdets,Res,T,conf);

  // number of triangulations computed
  std::pair<int,int> num_of_triangs(start_triangs,augment_triangs);

  return num_of_triangs;
}


////////////////////////////////////////////////////////////
// misc

template <class Triang>
double recompute_Res(const Triang &Res){
  std::vector<PPoint_d> Res_points;
  for (typename Triang::Vertex_const_iterator vit = 
       Res.vertices_begin(); 
       vit != Res.vertices_end(); vit++)
    Res_points.push_back(vit->point());
  Triangulation Res2(PD);
  double tstart, tstop;
	tstart = (double)clock()/(double)CLOCKS_PER_SEC;
  Res2.insert(Res_points.begin(),Res_points.end());
  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
  //std::cout << "CH Res offline time= " << tstop - tstart << std::endl;
  return tstop - tstart;
}

double compute_Res_offline(HD& Pdets,
                           Triangulation& Res){
  //Res.clear();
  Triangulation Res2(PD);
  std::vector<PPoint_d> Res_points;
  for (HD::iterator Pit=Pdets.begin(); 
                           Pit!=Pdets.end(); ++Pit){
		PPoint_d p(PD,Pit->begin(),Pit->end());
		p.set_hash(&Pdets);
		p.set_index(Pit-Pdets.begin());
		Res_points.push_back(p);
		//Res.insert(p);
  }
  double tstart, tstop;
	tstart = (double)clock()/(double)CLOCKS_PER_SEC;
  Res2.insert(Res_points.begin(),Res_points.end());
  tstop = (double)clock()/(double)CLOCKS_PER_SEC;
  return tstop - tstart;
}

Field volume(const Triangulation& Res, HD& Pdets){
#ifdef COMPUTE_VOL
	Field vol=0;
	for (PSimplex_finite_const_iterator cit=Res.finite_full_cells_begin();
             cit!=Res.finite_full_cells_end();
             cit++){
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
#else
  return Field(-1);
#endif
}

int num_of_vertices(CTriangulation& T){
	//PFace f;
	typedef std::vector<CFace> Faces;
	Faces edges;
	std::back_insert_iterator<Faces> out(edges);
  // incident_faces is not declared const, why?
	T.incident_faces(T.infinite_vertex(), 1, out);
	// Count the number of points on the convex hull
	std::cout << "[";
	for (Faces::const_iterator fit=edges.begin(); fit!=edges.end(); fit++){
		CPoint_d p = (*fit).vertex(1)->point();
		for (CPoint_d::Cartesian_const_iterator pit=p.cartesian_begin();
         pit!=p.cartesian_end();
         pit++){
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

template <class NT_>
int num_of_simplices(const CTriangulation& T,
                     const std::map<std::vector<NT_>, int>& points_index){

  typedef NT_                                           Field;
	// Count the number of points on the convex hull
	std::cout << "{";
	for (CSimplex_const_iterator fit=T.full_cells_begin();
       fit!=T.full_cells_end();
       fit++){
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
												PPoint_d& new_point,
                        const ResPol::config &conf){
	PLocate_type loc_type;
  PFace f(PD);
	PFacet ft;
  ppc.locate(new_point,loc_type,f,ft);
  if (loc_type==4 || loc_type==5){ //replace with OUTSIDE_AFFINE_HULL OUTSIDE_CONVEX_HULL
    if(conf.verbose>1){
      std::cout<<"one new R-vertex found !!! loc_type="<<loc_type<<std::endl;
    }
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
		;//std::cout << "no new R-vertex found"<< endl;
	}
}
*
*
*
void insert_new_Rvertex(Triangulation& ppc,
												NV_ds& normal_list,
												SRvertex& new_vertex,
												HD& Pdets,
                        const ResPol::config &conf){
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
    if(conf.verbose>1){
      std::cout<<"one new R-vertex found !!! loc_type="<<loc_type<<std::endl;
    }
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
		;//std::cout << "no new R-vertex found"<< endl;
	}
}
*
*

// iterate through all facets and put all the points of upper (lower) hull facets to triang the upper hull is the default
std::vector<std::set<int> > project_lower_upper_hull(CTriangulation& pc,
                                           std::map<std::vector<Field>, int>& points_index,
                                           int dcur,
                                           int dim,
                                           bool lower,
                                           const ResPol::config &conf){
	//std::cout << dim << "|" << pc.current_dimension() << std::endl;
	//SimplicialComplex triang, triang2;
	std::vector<std::set<int> > triang;
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
			std::set<int> s;
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
      //  if(conf.verbose>1){
      //    std::cout<<"found a low dimensional simplex. possible ERROR\n"<<
      //      std::endl;
      //  }
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
std::vector<Field> compute_r_proj(const std::vector<std::set<int> >& triang, std::vector<std::vector<Field> >& points, int m, std::vector<int>& mi, int d, std::map<std::vector<Field>,int>& points_index,HD& dets, std::vector<int> proj){
	double top1, top2;
	std::vector<Field> r(PD,0);
	for (std::vector<std::set<int> >::const_iterator it = triang.begin(); it != triang.end();it++){
		//std::vector<std::vector<Field> > simplex;
		std::vector<size_t> simplex_vec;
	  std::vector<std::vector<int> > sets(d);
	  int current_set=0, current_m=mi[0];
	  for (std::set<int>::const_iterator it2=it->begin(); it2!=it->end(); it2++){
			if (*it2 >= current_m)
				current_m += mi[++current_set];
			sets[current_set].push_back(*it2);
			//std::vector<Field> point = points[*it2];
			//point.push_back(1);
	    //simplex.push_back(point);
	    simplex_vec.push_back(*it2);
	  }
	  std::vector<int> mixed_vertices;
	  for (int i=0; i<d; i++){
	  	if (sets[i].size() == 1){
	  	  mixed_vertices.push_back(i);
	  	}
	  }
	  if (mixed_vertices.size() == 1){
	  	// find where is the mixed vertex in the vector of projection coordinates
	  	std::vector<int>::iterator pit = find(proj.begin(),proj.end(),sets[mixed_vertices[0]][0]);
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
std::vector<std::set<int> > project_lower_upper_hull_fast(CTriangulation& pc,
                                           std::map<std::vector<Field>, int>& points_index,
                                           HD& dets,
                                           int dcur,
                                           int dim,
                                           bool lower){
	//std::cout << dim << "|" << pc.current_dimension() << std::endl;
	//SimplicialComplex triang, triang2;
	std::vector<std::set<int> > triang;
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
			std::set<int> s;
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

// vim: ts=2
