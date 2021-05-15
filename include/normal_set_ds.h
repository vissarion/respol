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

#ifndef NORMAL_VECTOR_DS_H
#define NORMAL_VECTOR_DS_H

#ifdef USE_BOOST_FLAT_SET
#include <boost/version.hpp>
#if BOOST_VERSION < 104800
#error Flat set needs Boost 1.48 or newer
#endif
#include <boost/container/flat_set.hpp>
#else
#include <set>
#endif

#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_d.h>

// comparison between data, needed to keep them sorted!
template <class V>
struct lv_compare:
public std::binary_function<V,V,bool>
{
  bool operator()(const V& v1, const V& v2) const
  {
    return std::lexicographical_compare(v1.deltas_begin(),
                                        v1.deltas_end(),
                                        v2.deltas_begin(),
                                        v2.deltas_end());
  } 
};

template <class V, class DT>
class Normal_Vector_ds:
#ifdef USE_BOOST_FLAT_SET
public boost::container::flat_set<V,lv_compare<V> >
#else
public std::set<V,lv_compare<V> >
#endif
{
private:

  // typedefs
  typedef V                                             data;
#ifdef USE_BOOST_FLAT_SET
  typedef boost::container::flat_set<V,lv_compare<V> >  base;
#else
  typedef std::set<V,lv_compare<V> >                    base;
#endif
  typedef typename base::const_iterator                 base_const_iterator;
  typedef typename base::iterator                       base_iterator;

public:

  // initialization of ds
  void initialize(){
    DT vec[]={-1,0,1};
    int base=sizeof(vec)/sizeof(DT);

    std::cout << base << " " << PD << " " << pow(base,PD) << std::endl;
    //exit(0);

    for (int i=0; i<pow(base,PD); ++i){
      //std::cout << i << "/" << pow(base,PD) << std::endl;
      std::vector<DT> extreme_point;
      std::vector<int> v=int2vectord(i,base,PD);
      //copy(v.begin(),v.end(),ostream_iterator<int>(std::cout," "));
      for (std::vector<int>::iterator it=v.begin(); it!=v.end(); it++){
        extreme_point.push_back(vec[*it]);
        //std::cout << vec[*it] << " ";
      }
      //std::cout << std::endl;
      if (!is_zero(extreme_point)){
        data lft(PD,extreme_point.begin(),extreme_point.end());
        put(lft);
      }
    }
  }

  // RANDOM initialization of ds
  void random_initialize(){
      random_initialize(100);
  }

  void random_initialize(int k){
    // Instanciate a random point generator
    CGAL::Random rng;
    typedef CGAL::Random_points_in_ball_d<V> Random_points_iterator;
    Random_points_iterator rand_it(PD, 1.0, rng);

    // Generate 1 random point
    std::vector<V> points;
    CGAL::copy_n(rand_it, k, std::back_inserter(points));
    for (typename std::vector<V>::const_iterator it=points.begin();
         it!=points.end();++it){
      //std::cout << "random point=" << V(*it) << std::endl;
      put(V(*it));
    }
  }

  // initialization of ds
  void simple_initialize(){
    int vec[]={-1,1};
    for (int j=0; j<2; j++){
      for (int i=0; i<PD; i++){
        std::vector<int> extreme_point(PD,0);
        extreme_point[i]=vec[j];
        //for (int i=0; i<PD; i++){
        //  std::cout<<extreme_point[i]<<" ";
        //}
        //std::cout<<std::endl;
        data lft(PD,extreme_point.begin(),extreme_point.end());
        put(lft);
      }
    }
    //std::cout<<this->size()<<std::endl;
  }

  int put(data d){
    std::pair<base_iterator,bool> result = this->insert(d);
    return (result.second)?1:0;
  }

  V back(){
    base_iterator last = this->end();
    last--;
    return *last;
  }

  void pop_back(){
		this->erase(back());
  }

  // insert data
  //int put(const data &d){
  //  if (find(this->begin(),this->end(),d) == this->end()){
  //    this->push_back(d);
  //    return 1;
  //  }
  //  else
  //    return 0;
  //}

  // print data
  void print() const{
    std::cout<<"current normal data structure size: "<<this->size()<<std::endl;
  }

  void print_all() const{
    std::cout<<"current normal data structure: "<<this->size()<<std::endl;
    for (base_const_iterator it=this->begin(); it!=this->end(); it++)
      std::cout<<*it<<" ";
    std::cout<<std::endl;
  }

private:

  static std::vector<int> int2vectord(int k, int vash, int d){
    std::vector<int> b;
    while (k!=0){
      b.insert(b.begin(),k%vash);
      k=k/vash;
    }
    while (b.size() != d)
      b.insert(b.begin(),0);
    return b;
  }

  static bool is_zero(const std::vector<DT> &v){
    for (typename std::vector<DT>::const_iterator it=v.begin()+1;
         it!=v.end();
         it++){
      if (*it != DT(0))
        return false;
    }
    return true;
  }

};

#endif //NORMAL_VECTOR_DS_H
// vim: ts=2:expandtab
