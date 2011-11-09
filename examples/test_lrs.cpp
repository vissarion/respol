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

#include <CGAL/Gmpz.h>
#include <lrs_cgal.h>

// Example:
// -consider the 2d points p(0,0), q(3,3), r(0,3) and s(1,2)
// -the convex hull are the hyperplanes x>=y, y<=3 and x>=0, with
// -normal vectors to the convex hull are (-1,-1), (0,-1) and (1,0),
// respectively
int main(){
        typedef CGAL::Gmpz                                      NT;
        typedef std::vector<NT>                                 Point;
        typedef std::vector<Point>                              Points;
        typedef std::set<std::vector<NT> >                      Normals;

        // construct the input, the points and the vector containing them
        Point p(2),q(2),r(2),s(2);
        p[0]=NT(0);p[1]=NT(0);
        q[0]=NT(3);q[1]=NT(3);
        r[0]=NT(0);r[1]=NT(3);
        s[0]=NT(1);s[1]=NT(2);
        Points input(4);
        input[0]=p;input[1]=q;input[2]=r;input[3]=s;
        // compute the convex hull
        Normals output=lrs_ch(input);
        std::cout<<"normal vectors:\n";
        for(typename Normals::const_iterator i=output.begin();
            i!=output.end();
            ++i){
                std::cout<<"[ ";
                for(typename std::vector<NT>::const_iterator j=i->begin();
                    j!=i->end();
                    ++j){
                        std::cout<<(*j)<<' ';
                }
                std::cout<<"]\n";
        }
        return 0;
}
