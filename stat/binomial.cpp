#include "stat.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
  if(argc != 3){
    cerr << "usage : binomial N p(s)" << endl;
    exit(1);
  }
  int N = atoi(argv[1]);
  float p = atof(argv[2]);
  float psum = 0;
  for(int s=0; s <= N; s++){
    float ps = binomialProb(N, s, p);
    psum += ps;
    cout << N << "\t" << p << "\t" << s << "\t" << ps << endl;
  }
  cout << "Sum of probs is : " << psum << endl;
}
