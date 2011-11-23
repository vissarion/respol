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

int D,CD,PD; // global variables needed by read_input
#include <parse_functions.h>
#include <CGAL/Gmpz.h>
#include <vector>
#include <cassert>

template <class NT_>
std::ostream& print_gfan_vectors(std::ostream &o,
                                 const std::vector<std::vector<NT_> > &points,
                                 const std::vector<int> &mi){
        assert(points.size());
        std::cout<<"{\n{";
        size_t m_idx=0,m_total=mi[0]-1;
        for(size_t i=0;i<points.size()-1;++i){
                o<<'('<<points[i]<<')';
                if(i==m_total){
                        o<<"},\n{";
                        m_total+=mi[++m_idx];
                }else{
                        o<<',';
                }
        }
        return o<<'('<<points[points.size()-1]<<")}\n}\n";
}

int main(){
        // NT is the number type of the input exp vectors
        typedef CGAL::Gmpz                              NT;
        // Pset is the set of the M_i
        typedef std::vector<std::vector<NT> >           Pset;
        // Cset is a coordinate set (for the Mi's and the projections)
        typedef std::vector<int>                        Cset;

        Pset points;
        Cset mi,proj;
        int m;

        switch(read_pointset(points,mi,proj,m)){
                case 1:
                        std::cerr<<"# implicitization"<<std::endl;
                        print_gfan_vectors(std::cout,points,mi);
                        break;
                case 2:
                        std::cout<<"# arbitrary projection"<<std::endl;
                        break;
                case 3:
                        std::cout<<"# generic polytope"<<std::endl;
                        break;
                default:
                        std::cerr<<"I do not know how to handle this case!"
                                <<std::endl;
                        exit(-1);
        }

        std::cout<<"---------------\n";
        std::cout<<"points = {("<<points[0]<<')';
        for(size_t i=1;i<points.size();++i)
                std::cout<<",("<<points[i]<<')';
        std::cout<<"}\nmi = ["<<mi<<
                "]\nproj = ["<<proj<<
                "]\nm = "<<m<<std::endl;;

        return 0;
}
