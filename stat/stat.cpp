// some statistical functions like, hope this works..

#include <math.h>
#include <algorithm>
#include <vector>
#include "stat.h"
#include <iostream>
using namespace std;

float mean(vector<float>& v){
  float sum=0;
  for(int i=0; i < v.size(); i++){
    sum += v[i];
  }
  return(sum/v.size());
}
/*
float std_dev(vector<float>& v){
  float m = mean(v);
  float sqsum =0;
  for(int i=0; i < v.size(); i++){
    sqsum += (v[i]-m)*(v[i]-m);
  }
  return(sqrt(sqsum/(v.size()-1)));
}
*/
float std_dev(vector<float>& v){
  //
  // Note..      SS = sum of squared devieates (squared sum)
  // AND         SS = sum(X^2) - ((sum(X))^2 / N)
  // so          std = sqrt( SS / N-1)
  // which should be reasonably fast to calculate.. (rather than using the std function.. (which I could rewrite)
  if(v.size() > 1){
    float xsum = 0;
    float xsquaredSum = 0;
    for(int i=0; i < v.size(); i++){
      xsum += v[i];
      xsquaredSum += (v[i]*v[i]);
    }
    float SS = xsquaredSum - (xsum*xsum)/(float)v.size();
    //cout << "SS : " << SS << "\tstd : " << sqrt(SS/(float)(v.size()-1)) << endl;
    return(sqrt(SS/(float)(v.size()-1)));
    //    cout << "doing new standard deviation.. " <<sqrt( (xsquaredSum - ((xsum*xsum)/(float)v.size()) ) / (float)v.size()-1)
    // << endl;
    //return(sqrt( (xsquaredSum - ((xsum*xsum)/(float)v.size()) ) / (float)v.size()-1));
  }
  cout << "vector too small returning impossible std deviation of -1 " << endl;
  return(-1);        // hmmm,,,, hmmm
}


float min(vector<float>& v){
  float minimum;
  if(v.size() < 1){
    cout << "min function from stat.cpp.. -vector has size 0, returning an unintialized float." << endl
	 << "be aware that you should probably run screaming out of the room, and vanquish any" << endl
	 << "ambitions you ever had to get anything done" << endl;
    return(minimum);  // which is actually really bad, because we need to throw an exception...
  }
  minimum = v[0];
  for(int i=1; i < v.size(); i++){
    if(minimum > v[i]) { minimum = v[i]; }
  }
  return(minimum);
}

float max(vector<float>& v){
  float maximum;
  if(v.size() < 1){
    cout << "max function from stat.cpp.. -vector has size 0, returning an unintialized float." << endl
	 << "be aware that you should probably run screaming out of the room, and vanquish any" << endl
	 << "ambitions you ever had to get anything done" << endl;
    return(maximum);  // which is actually really bad, because we need to throw an exception...
  }  
  maximum = v[0];
  for(int i=1; i < v.size(); i++){
    if(maximum < v[i]) { maximum = v[i]; }
  }
  return(maximum);
}

float euclidean(vector<float>& v1, vector<float>& v2){
  // simply calculate the euclidean distance between two things.
  float d = 0;     // the distance..
  if(v1.size() != v2.size()){
    cout << "stat.cpp, euclidean function. vectors of different sizes. returning -1" << endl;
    d = -1;
    return(d);
  }
  for(int i=0; i < v1.size(); i++){
    d += pow(v1[i]-v2[i], 2);
  }
  d = sqrt(d);
  return(d/v1.size());
}

vector<float> z_score(vector<float>& v){
  vector<float> z(v.size());
  float m = mean(v);
  float s = std_dev(v);
  for(int i=0; i < z.size(); i++){
    z[i] = (v[i]-m)/s;
  }
  return(z);
}

void zScore(vector<float>& v){
  ///// MODIFIES THE COPY IT RECEIVES AS THIS IS FASTER..
  float m = mean(v);
  float s = std_dev(v);
  for(int i=0; i < v.size(); i++){
    v[i] = (v[i]-m)/s;
  }
}

float median(vector<float> v){
  // not particularly efficient, but what the hell
  if(v.size() < 1){
    return(0);
  }
  sort(v.begin(), v.end());
  int mid = v.size()/2;
  float median;
  if(v.size() %2 == 0){
    median = (v[mid] + v[mid-1])/2;
  }else{
    median = v[mid];
  }
  return(median);
}

float mad(vector<float>& v){
  // first get median..
  float m = median(v);
  // then go through each thingy and get the 
  // the difference
  vector<float> absdev(v.size());
  for(int i=0; i < v.size(); i++){
    absdev[i] = fabs(v[i] - m);
  }
  // then get the median of the deviations..
  float m_ad = median(absdev);
  return(m_ad);
}
  
vector<float> m_score(vector<float>& v){
  vector<float> ms(v.size());
  float m = median(v);
  float m_ad = mad(v);
  for(int i=0; i < ms.size(); i++){
    ms[i] = 0.6745 * ((v[i] - m)/m_ad);
  }
  return(ms);
}

vector<float> d_series(vector<float>& v){
  // a somewhat confused concept, but one that may have some useful parts to it
  // essentially returns a vector of v.size(-1) which represents the difference
  // between the first member of v, and the rest.. -may be useful in the case where
  // the first member is some sort of control, -then the resulting point in N-dimensional
  // space will tend to represent how different that 'gene' is from its control state..
  // or the reproducibility.. or something... -useful for 3 samples as well, as it allows 
  // a 2 d representation of the relationships.. in thingss...
  vector<float> d(v.size()-1);
  for(int i=1; i < v.size(); i++){
    d[i-1] = v[i]-v[0];
  }
  return(d);
}


vector<int> f_distribution(vector<float>& v, float minV, float maxV, int divs){
  // choice: 0 -all, 1, pm, 2 is mm;
  vector<int> dis(divs, 0);     // i.e. make a vector with divs members..
  // go through the vector members and initialise to 0. this may not be 
  // necessary, but I prefer to do it nevertheless.
  //for(int i=0; i < dis.size(); i++) { dis[i] = 0; }
  // count ..
  for(int i=0; i < v.size(); i++){
    if(v[i] >= minV && v[i] <= maxV){ 
      dis[(int)(((v[i]-minV)*(float)divs)/(maxV-minV))]++;
    }
  }
  return(dis);
}

vector<int> l_distribution(vector<float>& v, float minV, float maxV, int divs){
  // double log distribution,, not sure how to do this, but lets have a go.. 
  vector<int> dis(divs, 0);
  
  if(minV <= 0){
    cout << "returning a load of 0's as the minValue is less than 0 " << endl;
    return(dis);
  }
  float range = log(maxV)-log(minV);
  float lminV = log(minV);
  int index; 
  for(int i=0; i < v.size(); i++){
    if(v[i] >= minV && v[i] <= maxV){ 
      index = (int)(((log(v[i])-lminV)/range)*(float)divs);
      if(index > 0 && index < dis.size()){
	dis[index]++;
      }else{
	cout << "index out of range, index: " << index << "\tsize: " << dis.size() << endl;
	cout << "v[" << i << "] : " << v[i] << "\tminV: " << minV << "\trange: " << range << endl;
      }
      //      dis[(int)((log(v[i]-minV)/(range))*(float)divs)]++;
      //      dis[(int)log(((v[i]-minV)*(float)divs)/(maxV-minV))]++;
    }
  }
  return(dis);
} 

vector<float> norm_median(vector<float> v){
  // calculate the median value, and divide all members by it
  // simple..
  // operates on a local copy that gets modified..
  if(v.size() == 0) { return v; }
  sort(v.begin(), v.end());
  int mid = v.size()/2;
  float median;
  if(v.size() % 2 == 0){
    median = (v[mid] + v[mid-1])/2;
  }else{
    median = v[mid];
  }
  for(int i=0; i < v.size(); i++){
    v[i] = v[i]/median;
  }
  return v;
}

vector<float> norm_mean(vector<float> v){
  // just divides by the mean value..
  float sum = 0;
  for(int i=0; i < v.size(); i++) { sum += v[i]; }
  sum = sum/(float)v.size();  // convert to mean..
  for(int i=0; i < v.size(); i++) { v[i] = v[i]/sum; }
  return(v);
}

vector<float> norm_min_median(vector<float> v){
  if(v.size() == 0) { return(v); }
  if(v.size() == 1) {
    v[0] = 0;
    return(v);
  }
  sort(v.begin(), v.end());
  float median;
  float min = v[0];
  int mid = v.size()/2;
  if(v.size() % 2 == 0){
    median = (v[mid] + v[mid-1])/2;
  }else{
    median = v[mid];
  }
  // go through and normalise
  float min_med = median-min;
  for(int i=0; i < v.size(); i++){
    v[i] = (v[i]-min) / min_med;
  }
  return(v);
}

vector<float> norm_min_mean(vector<float> v){
  // find the min, and mean values..
  float mean = 0;
  float min = v[0];
  for(int i=0; i < v.size(); i++){
    mean += v[i];
    if(v[i] < min) { min = v[i]; }
  }
  mean = mean/v.size();
  float min_mean = mean-min;
  // go through and normalise
  for(int i=0; i < v.size(); i++){
    v[i] = (v[i]-min) / min_mean;
  }
  return(v);
}

float maxMeanDeviation(vector<float>& v){
  // returns the maximum deviation from the mean value in terms of abs((value-mean)/mean);
  float meanD = mean(v);
  float maxD = 0;
  for(int i=0; i < v.size(); i++){
    if( fabs((v[i] - meanD)/meanD) > maxD){
      maxD = fabs((v[i] - meanD)/meanD);
    }
  }
  return(maxD);
}

float binomialProb(int N, int s, float p){
  // calculates the binomial probability for seeing s number of successes
  // from N trials where ps is the probability of seeing a success
  // 
  //
  //   The binomial equation is..

  // 
  //                N!        
  //     P(s) =  --------- * (p)^s * (1-p)^(N-s)
  //             s! (N-s)!
  // 
  // which ofcourse becomes simpler if the probability of success and failure are the same,
  // but we can't really assume that. Never mind. 
  // factorials get tricky if they are large so we should probably reduce the equation in a managable manner
  // which shouldn't be too difficult, just follow my procedure in Hypgeo.cpp.. (~/cpp/stat/)
  // 
  // First work out the left hand ratio, using two vectors of ints..
  float ratio = 1;  // top over num 
  if(s != N){
    for(int i=1; i <= (N-s); i++){
      //cout << "multiplying ratio by : " << s+i << " / " << i << "\t";
      ratio = ratio * ((float)(s+i) / (float)i);
      //      ratio = ratio * ((float)(N-i) / (float)(s)); 
      // cout << "ratio: " << ratio << endl;
    }
  }
  // and now we have the ratio its very easy to do the calculation using the power function..
  return(ratio *  pow(p, s) * pow(1-p, N-s));
}
    

