// simply the declaration of some functions that 
// can be used with vectors and others to return some 
// statistical functions


#ifndef STAT_H
#define STAT_H

#include <vector>
// I always wonder about that bloody namespace std thingy..
// as I read somewhere that it was bad, but I can't quite 
// remember where..

float mean(std::vector<float>& v);
float std_dev(std::vector<float>& v);
float min(std::vector<float>& v);
float max(std::vector<float>& v);

float euclidean(std::vector<float>& v1, std::vector<float>& v2);
std::vector<float> z_score(std::vector<float>& v);
void zScore(std::vector<float>& v);                /// MODIFIES VECTOR V !!!!!
float maxMeanDeviation(std::vector<float>& v);        // returns the maximum deviation from the mean value calculated as abs((value-mean)/mean)

float median(std::vector<float> v);
float mad(std::vector<float>& v);
std::vector<float> m_score(std::vector<float>& v);

std::vector<float> d_series(std::vector<float>& v);  // see stat.cpp for description of this vague concept

std::vector<int> f_distribution(std::vector<float>& v, float minV, float maxV, int divs);
std::vector<int> l_distribution(std::vector<float>& v, float minV, float maxV, int divs);
std::vector<float> norm_median(std::vector<float> v);
std::vector<float> norm_mean(std::vector<float> v);   // these operate on and change a local copy, hence not passed by reference
std::vector<float> norm_min_median(std::vector<float> v); // value-min / median-min
std::vector<float> norm_min_mean(std::vector<float> v);   // value-min / mean-min.

float binomialProb(int N, int s, float p);

#endif




