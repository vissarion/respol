#ifndef PRINT_FUNCTIONS_H
#define PRINT_FUNCTIONS_H

#include <ostream>
#include <fstream>
#include <vector>
#include <set>

/////////////////////////////////////////////////////////////////
// overload of << operators for various types

#ifndef __FFLAFLAS_print_utils_H
// this will be enough for SRvertex and Resvertex, which are respectively
// vector<Field> and vector<SRvertex>
template <class T>
std::ostream& operator<<(std::ostream& ost, const std::vector<T> &V) {
  for (typename std::vector<T>::const_iterator it=V.begin(); it!=V.end(); it++)
    ost << *it << ",";
  return ost;
}
#endif

// this will be enough for Polytope, which is set<SRvertex>
template <class T>
std::ostream& operator<<(std::ostream& ost, const std::set<T> &V) {
  for (typename std::set<T>::const_iterator it=V.begin(); it!=V.end(); it++)
    ost << *it << ",";
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
void print_res_vertices(const Triang &Res){
  // print the vertices of the res polytope
  int number_of_vertices = 0;
  typedef typename Triang::Vertex_const_iterator        VCI;
  typedef typename Triang::Point_d                      P;
  typedef typename P::Cartesian_const_iterator          PCCI;
  //std::cout << "dim=" << Res.current_dimension() << std::endl;
  for (VCI vit = Res.vertices_begin(); vit != Res.vertices_end(); vit++){
    std::cout << "[";
    for (PCCI cit=vit->point().cartesian_begin();
         cit != vit->point().cartesian_end();
         cit++){
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

template <class Triang>
void print_res_facets_number(const Triang &Res){
  typedef typename Triang::Full_cell_handle             Simplex;
  typedef typename Triang::Facet_const_iterator         FCI;
  // print the vertices of the res polytope
  int number_of_facets = 0;
  for (FCI vit = Res.facets_begin(); vit != Res.facets_end(); vit++)
    number_of_facets++;
  std::cout << number_of_facets << std::endl;

  int number_of_inf_simplicies = 0;
  typedef std::vector<Simplex>                          Simplices;
  Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  Res.incident_full_cells(Res.infinite_vertex(), out);
  std::cout << inf_simplices.size() << std::endl;
}

template <class Triang>
void print_cells_data(const Triang &Res){
  typedef typename Triang::Full_cell_handle             Simplex;
  typedef std::vector<Simplex> Simplices;

  Simplices inf_simplices;
  std::back_insert_iterator<Simplices> out(inf_simplices);
  Res.incident_full_cells(Res.infinite_vertex(), out);
  std::cout << "simplex data:";
  for(typename Simplices::const_iterator sit = inf_simplices.begin();
      sit != inf_simplices.end();
      ++sit )
    std::cout << (*sit)->data();
  std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////
// functions to print statistics

// this is the old version of print_statistics, see below for the newer one
void print_statistics(int numoftriangs,
                      int numofinitvertices,
                      int numofvertices,
                      double timeall){
  std::cout << std::endl;
  std::cout << "Number of triangs enumerated \t" << numoftriangs << std::endl;
  std::cout << "Projected Res vertices (INIT)\t" << numofinitvertices <<
    std::endl;
  std::cout << "Projected Res vertices \t\t" << numofvertices << std::endl;
  std::cout << "Time overall   \t\t\t" << timeall << std::endl;
}

template <class Vol>
void print_statistics(int numoftriangs,
                      int numoftriangs2,
                      int numofvertices,
                      double timeall,
                      const Vol &volume){
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

template <class Vol>
void print_statistics_small(int Cdim,
                            int Pdim,
                            int num_of_input_points,
                            int numoftriangs,
                            int numofvertices,
                            double timeall,
                            const Vol &volume){
  std::cout << Cdim << " "
            << Pdim << " "
            << num_of_input_points << " "
            << numoftriangs  << " "
            << numofvertices  << " "
            << timeall << " "
            << volume  << " "
            << std::endl;
}

#endif // PRINT_FUNCTIONS_H
// vim: ts=2:expandtab
