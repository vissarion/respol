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

#ifndef PRINT_FUNCTIONS_H
#define PRINT_FUNCTIONS_H

#include <ostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <CGAL/Gmpz.h>
#include <normal_set_ds.h> // we need normals to compute f-vectors

/////////////////////////////////////////////////////////////////
// overload of << operators for various types

#ifndef __FFLAFLAS_print_utils_H
// this will be enough for SRvertex and Resvertex, which are respectively
// vector<Field> and vector<SRvertex>
template <class T>
std::ostream& operator<<(std::ostream& ost,const std::vector<T> &V){
  if(!V.size())
    return ost;
  typename std::vector<T>::const_iterator it=V.begin();
  std::cout<<(*it++);
  for(;it!=V.end();it++)
    ost<<","<<(*it);
  return ost;
}
#endif

// this will be enough for Polytope, which is set<SRvertex>
template <class T>
std::ostream& operator<<(std::ostream& ost,const std::set<T> &V){
  if(!V.size())
    return ost;
  typename std::set<T>::const_iterator it=V.begin();
  std::cout<<(*it++);
  for(;it!=V.end();it++)
    ost<<","<<(*it);
  return ost;
}

/////////////////////////////////////////////////////////////////
// functions to print some data structures

template <class T>
void print_vertices(std::vector<std::vector<T> >& Poly, std::ofstream& ofs){
  typedef typename std::vector<std::vector<T> >::const_iterator VVCI;
  typedef typename std::vector<T>::const_iterator               VCI;
  ofs << "[";
  for (VVCI Polyit=Poly.begin(); Polyit!=Poly.end(); Polyit++){
    ofs << "[";

    for (VCI it=Polyit->begin(); it!=Polyit->end(); it++){
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

// for input topcom file construction
template <class T>
void print_vertices_hom(const std::vector<std::vector<T> > &Poly,
                        std::ofstream& ofs){
  ofs << "[";
  for (typename std::vector<std::vector<T> >::const_iterator Polyit=
        Poly.begin();
       Polyit!=Poly.end();
       Polyit++){
    ofs << "[";

    for (typename std::vector<T>::const_iterator it=Polyit->begin();
         it!=Polyit->end();
         it++){
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

template <class Triang>
void print_res_vertices_with_index(const Triang &Res){
  // print the vertices of the res polytope
  std::cout << "[";
  for (typename Triang::Vertex_const_iterator vit = Res.vertices_begin();
       vit != Res.vertices_end();
       vit++)
    std::cout << vit->point() << " | " << vit->point().index() << "],[";
  std::cout << std::endl;
}

template <class Triang>
void print_res_vertices(const Triang &Res,
                        std::ostream& os){
  // print the vertices of the res polytope
  int number_of_vertices = 0;
  typedef typename Triang::Vertex_const_iterator        VCI;
  typedef typename Triang::Point_d                      P;
  typedef typename P::Cartesian_const_iterator          PCCI;
  //os << "dim=" << Res.current_dimension() << std::endl;
  for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
    os << "[";
    for (PCCI cit=vit->point().cartesian_begin();
         cit != vit->point().cartesian_end();
         cit++){
      os << (*cit).numerator();
      if (cit - vit->point().cartesian_begin() != vit->point().dimension()-1)
        os << ",";
    }
    //os << "|" << vit->point().index();
    os << "]";
    if (number_of_vertices++ != Res.number_of_vertices())
      os << ",";
  }
  os << std::endl;
}

// use maple to compute the number of extreme vertices of the
// Resultant polytope
template <class Triang>
int compute_extreme_res_vertices_maple(const Triang &Res){
#ifdef USE_MAPLE_CONVEX_HULL
  // write a file with maple commands that compute the extreme vertices
  std::ofstream outfile;
  outfile.open("maple_ch.mpl");
  outfile << "with(convex):" << std::endl;
  outfile << "points:=";
  print_res_vertices(Res,outfile);
  outfile << ":" << std::endl;
  outfile << "P1 := convhull(points);" << std::endl;
  outfile << "nops(vertices(P1));" << std::endl;
  outfile.close();

  // execute maple commands written in the previous file
  // and redirect the result in a new file
  std::cout << "Executing maple..." << std::endl;
  #define QUOTEME_(x) #x
  #define QUOTEME(x) QUOTEME_(x)
  std::string command=QUOTEME(MAPLE_EXECUTABLE);
  command+=" maple_ch.mpl > maple_ch_output.txt";
  #undef QUOTEME
  #undef QUOTEME_
  if (std::system(command.c_str())){
    std::cout << "Unable to execute command" << std::endl;
    exit(1);
  }
  // read the file of the results
  std::ifstream infile;
  infile.open("maple_ch_output.txt");
  std::string line;
  std::vector<std::string> lines;
  if (infile.is_open())
  {
    while ( infile.good() )
    {
      std::getline(infile,line);
      lines.push_back(line);
    }
    infile.close();
  }
  else std::cout << "Unable to open file" << std::endl;
  // the number of vertices is in the 5-th line from the end
  int n = atoi(lines[lines.end()-lines.begin()-5].c_str());

  return n;
#else
  return -1;
#endif
}

template <class Triang>
void print_res_facets_number(const Triang &Res){
  typedef typename Triang::Full_cell_handle                     Simplex;
  typedef typename Triang::Facet_const_iterator                 FCI;
  // print the vertices of the res polytope
  int number_of_facets = 0;
  for (FCI vit = Res.facets_begin(); vit != Res.facets_end(); vit++)
    number_of_facets++;
  std::cout << number_of_facets << std::endl;

  int number_of_inf_simplicies = 0;
  typedef std::vector<Simplex>                                  Simplices;
  Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  Res.incident_full_cells(Res.infinite_vertex(), out);
  std::cout << inf_simplices.size() << std::endl;
}

template <class Triang>
void f_vector(Triang &Res,
              int& cells,
	            int& triang_facets,
	            int& facets,
	            int& edges,
	            int& vertices){
  typedef typename Triang::Geom_traits                          GT;
  typedef typename GT::FT                                       FieldType;
  typedef typename Triang::Point_d                              PP_d;
  typedef typename GT::Direction_d                              PV_d;
  typedef typename GT::Hyperplane_d                             PH_d;
  typedef typename Triang::Full_cell_handle                     Simplex;
  typedef std::vector<Simplex>                                  Simplices;
  typedef typename Triang::Vertex_iterator                      VCI;
  typedef typename Triang::Finite_full_cell_const_iterator      PS_FCI;
  typedef Normal_Vector_ds<PV_d,FieldType>                      NV_ds;

  Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  Res.incident_full_cells(Res.infinite_vertex(), out);
  size_t finite_cells=0;
  for (PS_FCI cit=Res.finite_full_cells_begin();
             cit!=Res.finite_full_cells_end();
             cit++){
    ++finite_cells;
  }

  //for(typename Simplices::const_iterator sit = inf_simplices.begin();
  //    sit != inf_simplices.end();
  //    ++sit )
  //  std::cout << (*sit)->data();
  //std::cout << std::endl;
  int total_edges=0;
	for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
		typedef typename Triang::Face Face;
		typedef std::vector<Face> Faces;
		Faces edges;
		std::back_insert_iterator<Faces> out(edges);
		Res.incident_faces(vit, 1, out); // collect edges
		total_edges += edges.size();
		edges.clear();
	}
	//compute facets of Res polytope (not triangulated)
	int dim = Res.current_dimension();
	int num_facets = inf_simplices.end() - inf_simplices.begin();
	NV_ds Res_facets;
	for (typename std::vector<Simplex>::iterator sit=inf_simplices.begin();
	     sit!=inf_simplices.end();++sit){
		std::vector<PP_d> facet_points;
		for (int i=1; i<=Res.current_dimension(); i++){
			facet_points.push_back((*sit)->vertex(i)->point());
			//std::cout << (*sit)->vertex(i)->point() << " ";
		}
		//std::cout << std::endl;
    PP_d opposite_point=
      (*sit)->neighbor(0)->vertex((*sit)->mirror_index(0))->point();
    // compute a hyperplane which has in its negative side the opposite point
    PH_d hp(facet_points.begin(),
            facet_points.end(),
            opposite_point,
            CGAL::ON_NEGATIVE_SIDE);
    //is_neg=true;
    PV_d current_vector = hp.orthogonal_direction();
    Res_facets.insert(current_vector);
  }
#ifdef PRINT_INFO
	std::cout << "(cells" << ","
	          << "triang_facets" << ","
	          << "facets" << ","
	          << "edges" << ","
	          //<< "boundary edges" << ","
	          << "vertices)=";
#endif
  cells=finite_cells;
	triang_facets=num_facets;
	facets=Res_facets.size();
	edges=total_edges/2;
	vertices=Res.number_of_vertices();
}

template <class Triang>
void print_output(const Triang &Res,
									std::ostream& os){
  // print the vertices of the res polytope
  int number_of_vertices = 0;
  typedef typename Triang::Vertex_const_iterator        VCI;
  typedef typename Triang::Point_d                      P;
  typedef typename P::Cartesian_const_iterator          PCCI;
  for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
    for (PCCI cit=vit->point().cartesian_begin();
         cit != vit->point().cartesian_end();
         cit++){
      os << cit->numerator();
      if (cit->denominator() != CGAL::Gmpz(1)){
				std::cout << " NOT INTEGER coeff in RES pol";
				exit(-1);
			}
      if (cit - vit->point().cartesian_begin() != vit->point().dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << std::endl;
}

template <class Triang>
void print_polymake_testfile(const Triang &Res,
														 std::string ch_algo,
														 std::ostream& os){
  // print the vertices of the res polytope
  int number_of_vertices = 0;
  typedef typename Triang::Vertex_const_iterator        VCI;
  typedef typename Triang::Point_d                      P;
  typedef typename P::Cartesian_const_iterator          PCCI;
  //os << "dim=" << Res.current_dimension() << std::endl;
	os << "use Time::HiRes qw(gettimeofday tv_interval);\n";
	os << "use application 'polytope';\n";
	os << "prefer '"<< ch_algo <<"';\n";
	os << "my $p=new Polytope<Rational>;\n";
	os << "$p->POINTS=<<'.';\n";
  for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
    os << "1 ";
    for (PCCI cit=vit->point().cartesian_begin();
         cit != vit->point().cartesian_end();
         cit++){
      os << cit->numerator();
      if (cit->denominator() != CGAL::Gmpz(1)){
				std::cout << " NOT INTEGER coeff in RES pol";
				exit(-1);
			}
      if (cit - vit->point().cartesian_begin() != vit->point().dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << ".\n";
	os << "print ' ';\n";
	os << "print $p->N_POINTS;\n";
	os << "print ' ';\n";
	os << "print $p->N_VERTICES;\n";
	os << "print ' ';\n";
	os << "print $p->DIM;\n";
	os << "print ' ';\n";
	os << "my $t0 = [gettimeofday];\n";
	os << "my $f=$p->FACETS;\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print ' ';\n";
	os << "print $p->F_VECTOR;\n";
	//os << "print \"\\n\";\n";
  os << std::endl;
}

template <class Triang>
void generate_polymake_scripts(const Triang &Res){
	std::ofstream polymakefile1;
  polymakefile1.open("test_cdd.polymake");
  print_polymake_testfile(Res,"cdd",polymakefile1);

  std::ofstream polymakefile2;
  polymakefile2.open("test_lrs.polymake");
  print_polymake_testfile(Res,"lrs",polymakefile2);

  std::ofstream polymakefile3;
  polymakefile3.open("test_bb.polymake");
  print_polymake_testfile(Res,"beneath_beyond",polymakefile3);
}

template <class Field>
void print_polymake_volume_file(std::vector< std::vector<Field> >& Hrep,
                                std::ostream& os){
  typedef typename std::vector<Field> Ineq;
  typedef typename std::vector<Ineq> Ineqs;
  // print the vertices of the res polytope
  //int number_of_vertices = 0;
  //typedef typename Triang::Vertex_const_iterator        VCI;
  //typedef typename Triang::Point_d                      P;
  //typedef typename P::Cartesian_const_iterator          PCCI;
  //os << "dim=" << Res.current_dimension() << std::endl;
	os << "use application 'polytope';\n";
	os << "my $inequalities=new Matrix<Rational>([\n";
  for(typename Ineqs::iterator iit=Hrep.begin();
      iit!=Hrep.end();
      iit++){
    os << "[";
    for(typename Ineq::iterator it=iit->begin();
        it!=iit->end();
        it++){
      if (it+1 != iit->end())
        os << "\"" << *it << "\"" << ",";
      else os << "\"" << *it << "\"";
    }
    if (iit+1 != Hrep.end())
        os << "],\n";
      else os << "]\n";
  }
	os << "]);\n";
	os << "my $p=new Polytope<Rational>(INEQUALITIES=>$inequalities);\n";
	os << "my $vol = $p->VOLUME;\n";
	os << "print $vol;\n";
  os << std::endl;
}


/////////////////////////////////////////////////////////////////
// functions to print statistics

template <class Vol>
void print_statistics(int current_dim,
                      int numoftriangs,
                      int numoftriangs2,
                      int numofvertices,
                      int numofextremevertices,
                      double timeall,
                      double timedet,
                      const Vol &volume){
  std::cout << std::endl;
  std::cout << "Res dim   \t\t\t\t\t\t\t" << current_dim << std::endl;
  std::cout << "Num of triangs enumed (init+augment)\t\t\t\t"
            << numoftriangs+numoftriangs2 << " ("
            << numoftriangs << "+" << numoftriangs2
            << ")" << std::endl;
  std::cout << "Projected Res vertices (extreme; if enabled) \t\t\t"
            << numofvertices
            << "(" << numofextremevertices << ")"
            << std::endl;
  std::cout << "Time overall     \t\t\t\t\t\t" << timeall << std::endl;
  std::cout << "Determinant time \t\t\t\t\t\t" << timedet << std::endl;
  std::cout << "Volume   \t\t\t\t\t\t\t" << volume
            << " ~ " <<  CGAL::to_double(volume) << std::endl;
}

template <class Vol,class Triang>
void print_statistics_small(int Cdim,
                            int Pdim,
                            int current_dim,
                            int init_num_of_input_points,
                            int num_of_input_points,
                            int numoftriangs,
                            int numofvertices,
                            int numofextremevertices,
                            double timeall,
                            double timehull,
                            double timeofflinehull,
                            double timedet,
                            const Vol &volume,
                            Triang& Res){
  int cells, triang_facets, facets, edges, vertices;
  f_vector(Res,cells, triang_facets, facets, edges, vertices);
  std::cout << Cdim << " "
            << Pdim << " "
            << current_dim << " "
            << init_num_of_input_points << " "
            << num_of_input_points << " "
            << numoftriangs  << " "
            << numofvertices  << " "
            << numofextremevertices << " "
            << timeall << " "
            << timehull << " "
            << timeofflinehull  << " "
            << timedet << " "
            << volume << " "
            << cells << " "
            << triang_facets << " "
            << facets << " "
            << edges << 
            std::endl;  
  //if (!(numofvertices+facets >= numoftriangs)){
//		std::cout << "WRONG CONJECTURE\n";
//		std::cout << "WRONG CONJECTURE\n";
//		std::cout << "WRONG CONJECTURE\n";
//		std::cout << "WRONG CONJECTURE\n"; 
//  }
}


template <class Vol,class Triang>
void pretty_print_statistics(int Cdim,
                            int Pdim,
                            int current_dim,
                            int init_num_of_input_points,
                            int num_of_input_points,
                            int numoftriangs,
                            int numofvertices,
                            int numofextremevertices,
                            double timeall,
                            double timehull,
                            double timeofflinehull,
                            double timedet,
                            const Vol &volume,
                            Triang& Res){
  int cells, triang_facets, facets, edges, vertices;
  f_vector(Res,cells, triang_facets, facets, edges, vertices);
  std::cout << "\ndimension (Cayley(d),proj(m),Res): \t(" 
            << Cdim << ","
            << Pdim << ","
            << current_dim << ")\n"
            << "#(input points,prepreocessed):  \t(" 
            << init_num_of_input_points << ","
            << num_of_input_points << ")\n"
            << "#(triangs,points boundary of Res):\t(" 
            << numoftriangs << ","
            << numofvertices  << ")\n"
            << "f-vector Res (cells,\ntriang_facets,facets,edges,vertices):\t(" 
            << cells << ","
            << triang_facets << ","
            << facets << ","
            << edges << "," 
            << numofextremevertices << ")\n"
            << "time (total,CH,CH offline,dets):\t(" 
            << timeall << ","
            << timehull << ","
            << timeofflinehull  << ","
            << timedet << ")\n"
            << "volume:\t\t\t\t\t" << volume << "\n";  
}

template <class Vol>
void print_statistics_rand(int Cdim,
                            int Pdim,
                            int current_dim,
                            const Vol &total_vol,
                            int num_of_input_points,
                            int numoftriangs,
                            int numofvertices,
                            int numofextremevertices,
                            double timeall,
                            double timehull,
                            double timeofflinehull,
                            double timedet,
                            const Vol &volume){
  std::cout << Cdim << " "
            << Pdim << " "
            << current_dim << " "
            << total_vol << " "
            << num_of_input_points << " "
            << numoftriangs  << " "
            << numofvertices  << " "
            << numofextremevertices << " "
            << timeall << " "
            << timehull << " "
            << timeofflinehull  << " "
            << timedet << " "
            << volume
#ifndef RESTRICTED_RES
            << std::endl
#endif
            ;
}

#endif // PRINT_FUNCTIONS_H
// vim: ts=2:expandtab
