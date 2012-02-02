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

#include <cassert>
#include <ctime>
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/point_generators_d.h>
#include <../include/fast_hashed_determinant.h>

int main(const int argc, const char** argv){

        typedef CGAL::Gmpq                              Field;
        typedef CGAL::Cartesian_d<Field>                CK;
        typedef CGAL::Triangulation<CK>                 CTriangulation;
        typedef CTriangulation::Point_d                 CPoint_d;
        typedef FastHashedDeterminant<Field>            HD;

        if(argc<3){
                std::cout<<"usage: "<<argv[0]<<
                        " dimension numberofpoints"<<std::endl;
                exit(-1);
        }
        int D=atoi(argv[1]);
        int N=atoi(argv[2]);

        // 1. generate N random points in dimension D
        std::vector<CPoint_d> v;
        CGAL::Random_points_in_cube_d<CPoint_d> gen(D,100);
        for(size_t i=0;i<N;++i)
                v.push_back(*gen++);
        // show the points
        //for(size_t i=0;i<N;++i)
        //        std::cout<<" "<<v[i]<<std::endl;

        // 2. put everything in a FHD table
        HD hash_table;
        for(size_t i=0;i<N;++i){
                v[i].set_hash(&hash_table);
                std::vector<Field> coords;
                for(size_t j=0;j<D;++j)
                        coords.push_back(v[i][j]);
                hash_table.add_column(coords);
                v[i].set_index(i);
        }
        // show hash table
        //std::cout<<"---"<<std::endl;
        //hash_table.print_matrix(std::cout);

        // 3. compute a triangulation and count the time spent
        CTriangulation T(D);
        clock_t triangulation_time=clock();
        for(size_t k=0;k<N;++k)
                T.insert(v[k]);
        triangulation_time=clock()-triangulation_time;

        // 4. output size of the triangulation and spent time
        std::cout<<(double)triangulation_time/CLOCKS_PER_SEC<<std::flush;

        return 0;
}
