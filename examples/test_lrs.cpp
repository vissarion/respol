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

#undef NDEBUG
#include <cassert>

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <lrs_cgal.h>

// Example:
// -consider the 2d points p(0,0), q(3,3), r(0,3) and s(1,2)
// -the convex hull are the hyperplanes x<=y, y<=3 and x>=0, with
// -normal vectors to the convex hull are (1,-1), (0,1) and (-1,0),
// respectively
template <class NT_>
int test_lrs(){
        typedef NT_                                             NT;
        typedef std::vector<NT>                                 Point;
        typedef std::vector<Point>                              Points;
        typedef std::vector<NT>                                 Normal;
        typedef std::set<std::vector<NT> >                      Polytope;

        // construct the input, the points and the vector containing them
        Point p(2),q(2),r(2),s(2);
        p[0]=NT(0);p[1]=NT(0);
        q[0]=NT(3);q[1]=NT(3);
        r[0]=NT(0);r[1]=NT(3);
        s[0]=NT(1);s[1]=NT(2);
        Points input(4);
        input[0]=p;input[1]=q;input[2]=r;input[3]=s;
        // create a LRS_CH object
        LRS_CH<NT> ch_object(input);
        // output the H-representation
        Polytope output=ch_object.get_h_rep();
        // check that the output is three hyperplanes
        if(output.size()!=3)
                return -1;
        typename Polytope::const_iterator hi=output.begin();
        // check that the first element is [-1,1,0]
        if(hi->size()!=3)
                return -2;
        if((*hi)[0]!=-1||(*hi)[1]!=1||(*hi)[2]!=0)
                return -3;
        // check that the second element is [0,-1,3]
        ++hi;
        if(hi->size()!=3)
                return -4;
        if((*hi)[0]!=0||(*hi)[1]!=-1||(*hi)[2]!=3)
                return -5;
        // check that the third element is [1,0,0]
        ++hi;
        if(hi->size()!=3)
                return -6;
        if((*hi)[0]!=1||(*hi)[1]!=0||(*hi)[2]!=0)
                return -7;
        // now, check the normals
        std::set<Normal> normals=ch_object.get_normals();
        // check that there are three normals
        if(normals.size()!=3)
                return -8;
        typename std::set<Normal>::const_iterator ni=normals.begin();
        // check that the first element is [-1,0]
        if(ni->size()!=2)
                return -9;
        if((*ni)[0]!=-1||(*ni)[1]!=0)
                return -10;
        // check that the second element is [0,1]
        ++ni;
        if(ni->size()!=2)
                return -11;
        if((*ni)[0]!=0||(*ni)[1]!=1)
                return -12;
        // check that the third element is [1,-1]
        ++ni;
        if(ni->size()!=2)
                return -13;
        if((*ni)[0]!=1||(*ni)[1]!=-1)
                return -14;
        return 0;
}

int main(){
        assert(test_lrs<long>()==0);
        assert(test_lrs<CGAL::Gmpz>()==0);
        assert(test_lrs<CGAL::Gmpq>()==0);
        return 0;
}
