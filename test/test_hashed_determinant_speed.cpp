#include <iostream>
#include <gmpxx.h>
#include "fast_hashed_determinant.h"
#include "hashed_determinant.h"
#include <ctime>

template <size_t size,class det>
clock_t time_det(){
        typedef mpq_class                                       NT;
        typedef HashedDeterminant<NT,size,det>                  HD;
        std::vector<std::vector<NT> > c(size);
        for(size_t i=0;i<size;++i)
                for(size_t j=0;j<size;++j)
                        c[i].push_back(NT(rand()));
        HD dets(c.begin(),c.end());

        std::vector<size_t> index;
        for(int i=size-1;i>=0;i--)
                index.push_back(i);
        clock_t start=clock();
        dets.determinant(index);
        clock_t end=clock();
        return end-start;
};

template <size_t size>
clock_t time_fast_det(){
        typedef mpq_class                                       NT;
        typedef FastHashedDeterminant<NT,size>                  HD;
        std::vector<std::vector<NT> > c(size);
        for(size_t i=0;i<size;++i)
                for(size_t j=0;j<size;++j)
                        c[i].push_back(NT(rand()));
        HD dets(c.begin(),c.end());

        std::vector<size_t> index;
        for(int i=size-1;i>=0;i--)
                index.push_back(i);
        clock_t start=clock();
        dets.determinant(index);
        clock_t end=clock();
        return end-start;
};

#define TEST_DIMENSION(_D,_N) \
        total_topcom=0; \
        total_linbox=0; \
        total_naive=0; \
        total_fast=0; \
        for(size_t i=0;i<_N;++i){ \
                total_topcom+=time_det<_D,det_topcom<_D,mpq_class> >(); \
                total_linbox+=time_det<_D,det_linbox<mpq_class> >(); \
                if((_D)<8) \
                        total_naive+=time_det<_D,det_naive<_D,mpq_class> >(); \
                total_fast+=time_fast_det<_D>(); \
        } \
        std::cout<<_D<<'\t'<<_N<<'\t'<< \
                ((double)total_topcom)/CLOCKS_PER_SEC<<'\t'<< \
                ((double)total_linbox)/CLOCKS_PER_SEC<<'\t'; \
        if((_D)<8) \
                std::cout<<((double)total_naive)/CLOCKS_PER_SEC; \
        else \
                std::cout<<"slow"; \
        std::cout<<'\t'<<((double)total_fast)/CLOCKS_PER_SEC<<std::endl;

int main(){
        clock_t total_topcom,total_linbox,total_naive,total_fast;
        std::cout<<"# size\tqty\ttopcom\tlinbox\tnaive\tfast"<<std::endl;
        TEST_DIMENSION(2,1000);
        TEST_DIMENSION(3,1000);
        TEST_DIMENSION(4,1000);
        TEST_DIMENSION(5,1000);
        TEST_DIMENSION(6,1000);
        TEST_DIMENSION(7,1000);
        TEST_DIMENSION(8,1000);
        TEST_DIMENSION(9,1000);
        TEST_DIMENSION(10,1000);
        TEST_DIMENSION(11,1000);
        TEST_DIMENSION(12,1000);
        TEST_DIMENSION(13,1000);
        TEST_DIMENSION(14,1000);
        TEST_DIMENSION(15,1000);
        TEST_DIMENSION(16,1000);
        TEST_DIMENSION(17,1000);
        TEST_DIMENSION(18,1000);
        TEST_DIMENSION(19,1000);
        TEST_DIMENSION(20,1000);
        return 0;
}
