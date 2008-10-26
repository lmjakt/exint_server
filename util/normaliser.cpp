//Copyright Notice
/*
    eXintegrator integrated expression analysis system
    Copyright (C) 2004  Martin Jakt & Okada Mitsuhiro
  
    This file is part of the eXintegrator integrated expression analysis system. 
    eXintegrator is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. 

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    PS. If you can think of a better name, please let me know...
*/
//End Copyright Notice

#include "normaliser.h"
#include "sorted_floats.h"

// don't use any stl containers, as these don't seem to work well with multiple processors

void Normaliser::zScore(float* v, unsigned int s){
  // find the mean and std deviation..
  // Note..      SS = sum of squared devieates (squared sum)
  // AND         SS = sum(X^2) - ((sum(X))^2 / N)
  // so          std = sqrt( SS / N-1)
  // which should be reasonably fast to calculate.. (rather than using the std function.. (which I could rewrite)
  float std;
  float mean;
  if(s > 1){
    float xsum = 0;
    float xsquaredSum = 0;
    for(int i=0; i < s; i++){
      xsum += v[i];
      xsquaredSum += (v[i]*v[i]);
    }
    mean = xsum/s;
    float SS = xsquaredSum - (xsum*xsum)/(float)s;
    std = sqrt(SS/(float)(s-1));
    // then go through and  modify
    for(int i=0; i < s; i++){
      v[i] = (v[i]-mean)/std;
    }
  }
}

void Normaliser::zScore(float** v, unsigned int ps, unsigned int es){
  for(int i=0; i < ps; i++){
    zScore(v[i], es);
  }
}

void Normaliser::mScore(float** v, unsigned int ps, unsigned int es){    // this algorithm, first tried with medians rather than means, but medieans probably too slow
  float std_sum = 0;
  float mean_sum = 0;
  if(es < 2){
    return;
  }
  float* means = new float[ps];
  float* std_devs = new float[ps];
  sorted_floats sorted_means(ps);
  sorted_floats sorted_stds(ps);
  for(int i=0; i < ps; i++){
    float xsum =0;
    float xsquaredSum = 0;
    for(int j=0; j < es; j++){
      xsum += v[i][j];
      xsquaredSum += (v[i][j]*v[i][j]);
    }
    means[i] = (xsum / (float)es);
    sorted_means.insert(means[i]);
    mean_sum += means[i];
    std_devs[i] = sqrt( (xsquaredSum - (xsum*xsum)/(float)es)/(es - 1) );
    sorted_stds.insert(std_devs[i]);
    std_sum += std_devs[i];
  }
  std_sum = std_sum / (float)ps;
  mean_sum = mean_sum / (float)ps;      // using median values is likely to be too expensive for this purpose.. 
  
  float median_mean = sorted_means.median();
  float median_std = sorted_stds.median();

  // and.. let's go. 
  for(int i=0; i < ps; i++){
    //float f = std_sum / std_devs[i];    //
    float f = median_std / std_devs[i];    //
    for(int j=0; j < es; j++){
      //v[i][j] = (f * (v[i][j] - means[i])) + mean_sum;
      v[i][j] = (f * (v[i][j] - means[i])) + median_mean;
    }
  }
  // and that is pretty much it as far as I can tell.. 
  delete []means;
  delete []std_devs;
}

float Euclid::sqEuclidean(float* v1, float* v2, unsigned int l){
  float d = 0;
  for(int i=0; i < l; i++){
    d += (v1[i] - v2[i])*(v1[i] - v2[i]);
  }
  return(d);
}

float Euclid::euclidean(float* v1, float* v2, unsigned int l){
  float d = 0;
  for(int i=0; i < l; i++){
    d +=  (v1[i] - v2[i])*(v1[i] - v2[i]);
  }
  return(sqrt(d));
}
