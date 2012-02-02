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
#include <fast_hashed_determinant.h>
#include <CGAL/Gmpz.h>
#include <vector>

#ifndef USE_ORIENTATION_DET
#error please set the flag USE_ORIENTATION_DET
#endif

int main(){
        typedef CGAL::Gmpz                              NT;
        typedef FastHashedDeterminant<NT>               HD;
        typedef std::vector<NT>                         Point;
        typedef std::vector<size_t>                     Index;
        Point p,q,r,s,unusedpoint;
        p.push_back(NT(1));
        p.push_back(NT(7));
        q.push_back(NT(2));
        q.push_back(NT(4));
        r.push_back(NT(3));
        r.push_back(NT(2));
        s.push_back(NT(3));
        s.push_back(NT(1));
        unusedpoint.push_back(8);
        unusedpoint.push_back(9);
        std::vector<Point> points;
        points.push_back(p);
        points.push_back(q);
        points.push_back(unusedpoint);
        points.push_back(r);
        points.push_back(s);
        HD dets(points.begin(),points.end());
        std::vector<NT> lift;
        lift.push_back(NT(4));
        lift.push_back(NT(10));
        lift.push_back(NT(2));
        lift.push_back(NT(7));
        Index idx;
        idx.push_back(0);
        idx.push_back(1);
        idx.push_back(3);
        idx.push_back(4);
        assert(dets.orientation(idx,lift)==NT(9));
        std::cout<<"correct execution"<<std::endl;
        return 0;
}
