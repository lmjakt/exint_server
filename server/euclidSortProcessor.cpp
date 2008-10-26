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

#include "euclidSortProcessor.h"
#include <qthread.h>
#include <qmutex.h>
#include <qwaitcondition.h>
#include <math.h>
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "connectionobject.h"
#include <iostream>

QMutex euclidMutex;      // for pushing the scores thingy.. (reserve space in it!!).

EuclidSortProcessor::EuclidSortProcessor(ProbeSetSet2* ps, probe_set* p, uint st, uint sp, vector<dist_set>* sc, set<int>& chips, vector<uint>& expts){
  pSet = ps;
  ps1 = p;
  begin = st;
  end = sp;
  //indices = ind;
  scores = sc;
  clientChips = chips;
  experimentNumber = expts.size();
  // it seems that calling new with 0 doesn't cause a crash, so seems like we don't have to do much.
  experiments = new uint[experimentNumber];   
  for(uint i=0; i < expts.size(); i++){
    experiments[i] = expts[i];
  }
}

EuclidSortProcessor::~EuclidSortProcessor(){
  delete experiments;
  cout << "deleting the processor thread" << endl;
}

void EuclidSortProcessor::run(){
  float value;
  if(begin < end){      // otherwise this doesn't make much sense..
    for(int i=begin; i < end; i++){
      //cout << "i " << i << "\tindex " << pSet->data[i]->index << endl;
      if(pSet->data[i]->index && clientChips.count(pSet->probeData[i].chip)){
	value = heavyCompare(ps1, pSet->data[i], experiments, experimentNumber);
	//	value = heavyCompare(ps1, pSet->data[i], ps1->exptIndex, ps1->exptSize);
      }else{
	value = -1;
      }
      //euclidMutex.lock();     /////////// NOT NEEDED IF EITHER ASSIGNMENT IS ATOMIC.. -WHICH I GUESS IT SHOULD BE?? 
      (*scores)[i].index = pSet->data[i]->index;
      (*scores)[i].value = value;
      //euclidMutex.unlock();
    }
  }else{
    cout << "BIg Screup in euclidSortProcessor start is bigger than stop,, " << endl;
  }
}

float** EuclidSortProcessor::copyProbes(float** p, int ps, int es){
  float** cp = new float*[ps];
  for(int i =0; i < ps; i++){
    cp[i] = new float[es];
    for(int j=0; j < es; j++){
      cp[i][j] = p[i][j];
    }
  }
  return(cp);
}

void EuclidSortProcessor::delProbes(float** p, int ps){
  for(int i=0; i < ps; i++){
    delete []p[i];
  }
  delete []p;
}

void EuclidSortProcessor::zScore(float* v, int s){  //s for the size.. 
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

float EuclidSortProcessor::heavyCompare(probe_set* pset1, probe_set* pset2, uint* expts, uint es){
  //
  //   Calculate the euclidean distance between the probe set 1 and the probe set 2 using a heavy
  //   algorithm which calculates the mean euclidean distance of all against all probe sets.
  //   
  //   Hence passing the same probe set in one and 2 will not give a 0 distance, but will give a 
  //   value which is related to how closely the different probe sets correlate with each other.
  //   
  //   Obviously I can only make a comparison for the same experimental points, so I have 
  //   to take that into account. Use the slow approach of actually creating a new set of vectors
  //   and manually filling those where both have got values. Then z-score normalise these..
  //   
  // OK, here goes..

  float** v1 = new float*[pset1->probeNo];  
  float** v2 = new float*[pset2->probeNo];
  int v1size = pset1->probeNo;
  int v2size = pset2->probeNo;
  
  // and reserve enough memory for each one..
  for(int i=0; i < pset1->probeNo; i++){ v1[i] = new float[es]; }
  for(int i=0; i < pset2->probeNo; i++){ v2[i] = new float[es]; }
  int selEx = 0; 
  
  for(int i=0; i < es; i++){
    if(expts[i] < pset1->allExptNo){           // allExptNo should be the same for both !!
      if(pset1->exptLookup[expts[i]] != -1 && pset2->exptLookup[expts[i]] != -1){
	for(int j=0; j < pset1->probeNo; j++){
	  v1[j][selEx] = pset1->probes[j][pset1->exptLookup[expts[i]]];
	}
	for(int j=0; j < pset2->probeNo; j++){
	  v2[j][selEx] = pset2->probes[j][pset2->exptLookup[expts[i]]];
	}
	selEx++;
      }else{
	delProbes(v1, v1size);
	delProbes(v2, v2size);
	return(-1.0);
      }
    }else{
      delProbes(v1, v1size);
      delProbes(v2, v2size);
      return(-1.0);    // which looks rather a bit too ugly.. 
    }
  }
  // AND now just do an all against all comparison, and take the mean value of the comparisons.. 
  float distance = 0;
  //cout << "beginning the distance calculations" << endl;
  // do the zScore transformation for the two thingies... 
  
  for(int i=0; i < pset1->probeNo; i++){ zScore(v1[i], selEx);}
  for(int i=0; i < pset2->probeNo; i++){ zScore(v2[i], selEx);}
  
  for(int i=0; i < v1size; i++){
    for(int j=0; j < v2size; j++){
      distance += euclidean(v1[i], v2[j], selEx);  
    }
  }
  distance = distance / (float)(v1size * v2size);
  // and lets delete the two float**,, we can use the function..
  delProbes(v1, v1size);
  delProbes(v2, v2size);
  //  distance = distance / (float)(pset1->probes.size() * pset2->probes.size());
  return(distance);
}

float EuclidSortProcessor::euclidean(float* v1, float* v2, int s){
  //cout << "beginning of euclidean s is " << s << endl;
  float distance = 0;
  //  if(v1.size() != v2.size() || v1.size() == 0) { 
  //  cout << "bugger, BUGGER,, euclidean functions, the vectors are of different size, returning -1" << endl;
  //  return(-1); 
  //}
  if(s < 1){ return(-1);}
  for(int i=0; i < s; i++){
    distance += ((v1[i]-v2[i]) * (v1[i]-v2[i]));
    //    distance += pow((v1[i]-v2[i]), 2);
  }
  // NOTE PLEASE,, note that distance is now the square of the euclidean distance, and 
  // for this to be real, I should just return the square root,, this however is a little troublesome,
  // as the larger the number of dimensions, the larger the distance, and this would then screw things up
  // if I compare distances with many measurements and distances with few. So I just divide by the number of 
  // dimensions. Maybe this should be the number of degrees of freedome, or something, but I don't know!!
  //cout << "euclidean function returning square root of " << distance << " / " << s << endl;
  return((sqrt(distance))/s);
}
