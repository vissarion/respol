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
