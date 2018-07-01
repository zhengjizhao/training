#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <parallel/numeric>
#include <parallel/algorithm>
//#include <gsl/gsl_histogram.h>


void dump_to_file(std::string& path, std::vector<double>& data){
    std::filebuf fb;
	fb.open(path, std::ios::out);
	std::ostream os(&fb);
	for (int i = 0; i<data.size(); ++i) {
		os << data[i] << '\n';
	}
	fb.close();
};


double covariance(std::vector<double>& a, std::vector<double>& b, double meanA, double meanB){
    std::vector<double> temp(a.size(),0);
#pragma omp parallel for simd
    for(int i = 0;i<a.size();++i){
        temp[i] = (a[i]-meanA)*(b[i]-meanB);
    }
    return __gnu_parallel::accumulate(temp.begin(),temp.end(),0.0)/a.size();
};


double max_para(std::vector<double>& x){
    return x[std::distance(x.begin(),__gnu_parallel::max_element(x.begin(),x.end()))];
};

double min_para(std::vector<double>& x){
    return x[std::distance(x.begin(),__gnu_parallel::min_element(x.begin(),x.end()))];
};

//void sort_based_on_index(std::vector<double>& x,std::vector<unsigned int>& indx){
//    std::vector<double> temp(x);
//#pragma omp parallel for
//    for(unsigned int i = 0;i<x.size();++i){
//        x[i] = temp[indx[i]];
//    }
//}
void sort_based_on_index(std::vector<double>& x,std::vector<int>& indx){
    std::vector<double> temp(x);
#pragma omp parallel for
    for(unsigned int i = 0;i<x.size();++i){
        x[i] = temp[indx[i]];
    }
}
