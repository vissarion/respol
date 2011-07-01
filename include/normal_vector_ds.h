#ifndef NORMAL_VECTOR_DS_H
#define NORMAL_VECTOR_DS_H

#include <vector>
#include <set>

using namespace std;

// comparison between data, needed to keep them sorted!
template <class V>
struct lv_compare
{
	bool operator()(const V& v1, const V& v2) const
	{
		return lexicographical_compare(v1.cartesian_begin(),v1.cartesian_end(),v2.cartesian_begin(),v2.cartesian_end());
	}
};

template <class V, class DT>
//class Normal_Vector_ds : public set<V,lv_compare<V> > 
class Normal_Vector_ds : public vector<V> 
{
private:

	// typedefs
	typedef V																		data;
	//typedef set<V,lv_compare<V> >							base;	
	typedef vector<V>														base;
	typedef typename base::const_iterator				base_const_iterator;
	
private:
	base _normal_list;  
	
public:
	
	//constructors
	Normal_Vector_ds() : _normal_list() {
		simple_initialize();
	}
	
	// initialization of ds
	void initialize(){
		DT vec[]={-1,0,1};
	  int base=sizeof(vec)/sizeof(DT);
	   
	  for (int i=0; i<pow(base,PD); ++i){
	    vector<DT> extreme_point;
	    vector<int> v=int2vectord(i,base,PD);
	    //copy(v.begin(),v.end(),ostream_iterator<int>(std::cout," "));
	    for (vector<int>::iterator it=v.begin(); it!=v.end(); it++){
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
	
	// initialization of ds
	void simple_initialize(){
		int vec[]={-1,1};
		for (int j=0; j<2; j++){
			for (int i=0; i<PD; i++){
				vector<int> extreme_point(PD,0);
				extreme_point[i]=vec[j];
				for (int i=0; i<PD; i++){
					std::cout<<extreme_point[i]<<" ";
				}
				std::cout<<std::endl;
				data lft(PD,extreme_point.begin(),extreme_point.end());
				put(lft);
			}
		} 	
		std::cout<<this->size()<<std::endl;
	}
		
	//void put(data d){
	//	insert(d);
	//}
	
	// insert data
	void put(data d){
		if (find(this->begin(),this->end(),d) == this->end())
		  this->push_back(d);
	}
	
	// print data
	void print() const{
		std::cout<<"current normal data structure size: "<<this->size()<<std::endl;
	}	
	
	void print_all() const{
		std::cout<<"current normal data structure: "<<this->size();
		for (base_const_iterator it=this->begin(); it!=this->end(); it++)
			std::cout<<*it<<" ";
		std:cout<<std::endl;
	}

private:	

	static vector<int> int2vectord(int k, int vash, int d){
		vector<int> b;
		while (k!=0){
			b.insert(b.begin(),k%vash);
			k=k/vash;
		}
		while (b.size() != d)
			b.insert(b.begin(),0);
		return b;
	}
	
	static bool is_zero(vector<DT> v){
		for (typename vector<DT>::iterator it=v.begin()+1; it!=v.end(); it++){
			if (*it != DT(0))
				return false;
		}	
		return true;
	}	

};


#endif //NORMAL_VECTOR_DS_H
