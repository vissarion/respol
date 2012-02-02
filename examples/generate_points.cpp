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

#include <iostream>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <CGAL/Gmpq.h>
#include <unistd.h>

typedef CGAL::Cartesian_d<double>                               Kd;
typedef Kd::Point_d                                             Point;

int main(int argc,char *argv[]){
        if(argc<3){
                std::cerr<<
                        "usage: "<<argv[0]<<" dimension nb_points"<<std::endl;
                exit(-1);
        }
        int dim=atoi(argv[1]);
        int nb_points=atoi(argv[2]);
        double size=100.0;
        std::cerr<<"Generating "<<nb_points<<" random points in a cube in "
                <<dim<<"D, coordinates from "<<-size<<" to "<<size<<std::endl;
        std::vector<Point> v;
        CGAL::Random r(getpid());
        CGAL::Random_points_in_cube_d<Point> gen(dim,size,r);
        //CGAL::Random_points_in_ball_d<Point> gen(dim,size,r);
        //CGAL::Random_points_on_sphere_d<Point> gen(dim,size,r);
        for(int i=0;i<nb_points;++i)
                v.push_back(*gen++);
        for(int i=0;i<nb_points;++i){
                std::cout<<dim;
                for(size_t j=0;j<dim;++j)
                        std::cout<<' '<<CGAL::Gmpq(v[i][j]);
                std::cout<<std::endl;
        }
        return 0;
}
