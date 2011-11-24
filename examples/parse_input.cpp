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
#include <algorithm>
#include <cstring>
#include <cassert>

template <class NT_>
std::ostream& print_gfan_vectors(std::ostream &o,
                                 const std::vector<std::vector<NT_> > &points,
                                 const std::vector<int> &mi){
        assert(points.size());
        o<<"{{";
        size_t m_idx=0,m_total=mi[0]-1;
        for(size_t i=0;i<points.size()-1;++i){
                o<<'('<<points[i]<<')';
                if(i==m_total){
                        o<<"},\n{";
                        m_total+=mi[++m_idx];
                }else{
                        o<<",";
                }
        }
        return o<<'('<<points[points.size()-1]<<")}}\n";
}

std::ostream& print_gfan_generic(std::ostream &o,int dimension){
        assert(dimension>0);
        o<<"(1";
        for(size_t i=1;i<dimension;++i)
                o<<",1";
        return o<<')'<<std::endl;
}

std::ostream& print_gfan_projections(std::ostream &o,
                                     const std::vector<int> &proj,
                                     int dimension){
        assert(proj.size());
        o<<'(';
        for(size_t i=0;i<dimension;++i){
                if(i!=0)
                        o<<',';
                if(std::search_n(proj.begin(),proj.end(),1,i)==proj.end())
                        o<<1;
                else
                        o<<0;
        }
        return o<<')'<<std::endl;
}

int main(int argc,char *argv[]){
        // NT is the number type of the input exp vectors
        typedef CGAL::Gmpz                              NT;
        // Pset is the set of the M_i
        typedef std::vector<std::vector<NT> >           Pset;
        // Cset is a coordinate set (for the Mi's and the projections)
        typedef std::vector<int>                        Cset;

        // process command-line
        bool verbose=false;
        bool valid;
        for(int i=1;i<argc;++i){
                valid=false;
                if(strcmp(argv[i],"-v")==0||strcmp(argv[i],"--verbose")==0){
                        valid=true;
                        verbose=true;
                }
                if(strcmp(argv[i],"-h")==0||strcmp(argv[i],"--help")==0){
                        valid=true;
                        std::cerr<<
                                "This program transform a ResPol input"<<
                                " read from stdin into a Gfan input, "<<
                                "written to stdout.\nValid options "<<
                                "are:\n-h, --help\tprint this help "<<
                                "message\n-v, --verbose\tbe verbose "<<
                                "(in stderr)"<<std::endl;
                                exit(-1);
                }
                if(!valid){
                        std::cerr<<"option "<<argv[i]<<" not recognized"<<
                                std::endl;
                        exit(-1);
                }
        }

        Pset points;
        Cset mi,proj;
        int m;

        switch(read_pointset(points,mi,proj,m)){
                case 1:
                        // the pipe symbol is omitted
                        // this is what we want to implement now
                        if(verbose){
                                std::cerr<<
                                "# implicitization; like in ex. (g), "<<
                                "run with:\n# (1) traversing tropical "<<
                                "resultant, gfan _resultantfan "<<
                                "--vectorinput --special\n# (2) normal "<<
                                "fan from stable intersection, gfan "<<
                                "_resultantfan --vectorinput --special "<<
                                "--projection\n# (3) normal fan from "<<
                                "tropical elimination, they say it is "<<
                                "'missing'"<<std::endl;
                        }
                        print_gfan_vectors(std::cout,points,mi);
                        print_gfan_projections(std::cout,proj,m);
                        break;
                case 2:
                        // projections specified after the pipe
                        if(verbose){
                                std::cerr<<
                                "# arbitrary projection; like in ex. "<<
                                "(d), run with:\n# (1) traversing "<<
                                "tropical resultant, gfan _resultantfan "<<
                                "--vectorinput --special\n# (2) normal "<<
                                "fan from stable intersection, gfan "<<
                                "_resultantfan --vectorinput --special "<<
                                "--projection\n# (3) normal fan from "<<
                                "tropical elimination, they say it is "<<
                                "'missing'"<<std::endl;
                        }
                        print_gfan_vectors(std::cout,points,mi);
                        print_gfan_projections(std::cout,proj,m);
                        std::cerr<<"not implemented"<<std::endl;
                        exit(-1);
                        break;
                case 3:
                        // pipe present, but no projections specified
                        if(verbose){
                                std::cerr<<
                                "# generic polytope; like in ex. "<<
                                "(a), run with:\n# (1) secondary fan, "<<
                                "gfan _secondaryfan\n"<<
                                "# (2) traversing tropical resultant, "<<
                                "gfan _resultantfan --vectorinput\n"<<
                                "# (3) normal fan from simple "<<
                                "description, gfan _resultantfan "<<
                                "--vectorinput --projection"<<std::endl;
                        }
                        std::cerr<<"not implemented"<<std::endl;
                        exit(-1);
                        break;
                default:
                        // it should never reach this point
                        std::cerr<<"I do not know how to handle this case!"
                                <<std::endl;
                        exit(-1);
        }

        /*std::cout<<"---------------\n";
        std::cout<<"points = {("<<points[0]<<')';
        for(size_t i=1;i<points.size();++i)
                std::cout<<",("<<points[i]<<')';
        std::cout<<"}\nmi = ["<<mi<<
                "]\nproj = ["<<proj<<
                "]\nm = "<<m<<std::endl;;
        */

        return 0;
}
