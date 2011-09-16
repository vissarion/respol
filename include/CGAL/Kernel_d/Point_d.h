// Copyright (c) 2000,2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Seel

#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H
#define HACKED_CGAL_POINT_D_H

#warning Using hacked Point_d interface

#include <CGAL/Dimension.h>

// for fast determinant computation
#ifdef USE_CGAL_DET
  #include <../include/fast_hashed_determinant_non_homog.h>
#else
  #include <../include/fast_hashed_determinant.h>
#endif

namespace CGAL {

template <class pR>
class Point_d : public pR::Point_d_base
{ public:
  typedef typename pR::Point_d_base Base;
  typedef Point_d<pR>               Self;
  typedef pR R;

private:
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;
  typedef FastHashedDeterminant<FT>     FHD;

private:
  size_t _index;
  FHD *_hash;

public:

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dimension_tag<0>      Feature_dimension;
    template < typename Kernel2 >
        struct WithAnotherKernel
        {
            typedef Point_d<Kernel2>  Type;
        };

  Point_d(int d=0) : Base(d),_index(0),_hash(NULL) {}
  Point_d(int d, const Origin &o) : Base(d,o),_index(0),_hash(NULL) {}

  Point_d(int a, int b, int c = 1) :
    Base(RT(a),RT(b),RT(c)),_index(0),_hash(NULL) {}
  Point_d(const RT& a, const RT& b, const RT& c = 1) :
    Base(a,b,c),_index(0),_hash(NULL) {}
  Point_d(int a, int b, int c, int d) :
    Base(RT(a),RT(b),RT(c),RT(d)),_index(0),_hash(NULL) {}
  Point_d(const RT& a, const RT& b, const RT& c, const RT& d) :
    Base(a,b,c,d),_index(0),_hash(NULL) {}

  template <class InputIterator>
  Point_d (int d, InputIterator first, InputIterator last)
    : Base (d, first, last),_index(0),_hash(NULL) {}
  template <class InputIterator>
  Point_d(int d, InputIterator first, InputIterator last, const RT& D)
    : Base (d, first, last, D),_index(0),_hash(NULL) {}

  Point_d(const Self &p) : Base(p),_index(p.index()),_hash(p.hash()) {}
  Point_d(const Base& p) : Base(p),_index(0),_hash(NULL) {}

  size_t index()const{
        return _index;
  }

  size_t set_index(size_t newi){
        size_t oldi=_index;
        _index=newi;
        return oldi;
  }

  FHD* hash()const{
        return _hash;
  }

  void set_hash(const FHD *h){
        _hash=(FHD*)h;
  }

  // TODO: for the moment, we don't care about _index and _hash when doing
  // arithmetic operations
  Vector_d<R> operator-(const Origin& o) const
  { return Base::operator-(o); }
  Vector_d<R> operator-(const Self& q) const
  { return Base::operator-(q); }
  Self operator+(const Vector_d<R>& v) const
  { return Base::operator+(v); }
  Self operator-(const Vector_d<R>& v) const
  { return Base::operator-(v); }
  Self& operator+=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator+=(v)); }
  Self& operator-=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator-=(v)); }

};

} //namespace CGAL
#endif //CGAL_POINT_D_H
