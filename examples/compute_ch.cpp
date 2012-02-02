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
#include <fstream>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Triangulation.h>
#include <ctime>
#ifdef USE_HASHED_DETERMINANTS
#warning using hashing determinants for ch example
#include "hashed_determinant.h"
#else
#warning NOT using hashing determinants for ch example
#endif

#ifdef USE_HASHED_DETERMINANTS
typedef HashedDeterminant<CGAL::Gmpq>                           HD;
#endif
typedef CGAL::Cartesian_d<CGAL::Gmpq>                           Kd;
typedef Kd::Point_d                                             Point;
typedef CGAL::Triangulation<Kd>                                 Triangulation;

int main(int argc,char *argv[]){
        if(argc<2){
                std::cerr<<
                        "usage: "<<argv[0]<<" input_file"<<std::endl;
                exit(-1);
        }
        std::fstream infile;
        Point p;
        size_t index=0;
        std::vector<Point> points;
#ifdef USE_HASHED_DETERMINANTS
        HD Pdets;
#endif
        infile.open(argv[1],std::fstream::in);
        int dim=0;
        while(infile>>p){
                dim=p.dimension();
#ifdef USE_HASHED_DETERMINANTS
                p.set_index(index++);
                p.set_hash(&Pdets);
#endif
                points.push_back(p);
#ifdef USE_HASHED_DETERMINANTS
                // now we have to add the point to the hash table
                std::vector<CGAL::Gmpq> column;
                for(size_t i=0;i<dim;++i)
                        column.push_back(p[i]);
                Pdets.add_column(column);
#endif
        }
        infile.close();
        clock_t elapsed_offline,elapsed_online;
        Triangulation triang_offline(dim),triang_online(dim);

        clock_t start=clock();
        //for(size_t i=0;i<points.size();++i)
        //        triang_offline.insert(points[i]);
        elapsed_offline=clock()-start;

        start=clock();
        triang_online.insert(points.begin(),points.end());
        elapsed_online=clock()-start;

        typedef Triangulation::Face Face;
        typedef std::vector<Face> Faces;
        Faces edges;
        std::back_insert_iterator<Faces> out(edges);
        triang_online.incident_faces(triang_online.infinite_vertex(),1,out);

        std::cerr<<dim<<'\t'<<
                points.size()<<'\t'<<
                //(double)elapsed_offline/CLOCKS_PER_SEC<<'\t'<<
                .0<<'\t'<<
                (double)elapsed_online/CLOCKS_PER_SEC<<'\t'<<
                edges.size()<<
                std::endl;
        //edges.clear();
        return 0;
}
