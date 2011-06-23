#ifdef FAST_ORIENTATION_H
#define FAST_ORIENTATION_H

#include <CGAL/enum.h>
#include <CGAL/Real_embeddable_traits.h>
#include <cassert>

// This function object computes the orientation of d+1 points of dimension
// d. The parameter _FHD is the FastHashedDeterminant object, which defines
// the dimension of the points, the number type of its coordinates and the
// type of the index used to referenciate the points. The indices in the
// idx parameter must be already inserted in the hash object. This function
// replaces the object OrientationCd in
// Kernel_d/include/CGAL/function_objectsCd.h, line 230.
// TODO: check the sizes of idx and r; are they d+1 or d+2?
template <class _FHD>
struct FastOrientation{
        typedef _FHD                                    FHD;
        typedef FHD::NT                                 NT;
        typedef CGAL::Real_embeddable_traits<NT>        RET;
        typedef RET::Sgn                                Sgn;
        typedef FHD::Row                                Row;
        typedef FHD::Index                              Index;
        typedef CGAL::Orientation                       Orientation;
        Orientation operator()(FHD &hash,const Index &idx,const Row &r)const{
                //assert(idx.size()==1+hash.dimension());
                //assert(r.size()==1+hash.dimension());
                // CGAL::Orientation is the same as CGAL::Sign.
                return Sgn()(hash.homogeneous_determinant(idx,r));
        }
};

#endif // FAST_ORIENTATION_H
