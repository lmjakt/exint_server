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

#include "kClusterProcess.h"
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "connectionobject.h"        // just for dist_set, not so good hmm. 
#include <vector>
#include <set>
#include <qthread.h>
#include <iostream>
#include <qmutex.h>                 // for changing the variable thingy.. 
#include <time.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

KClusterProcess::KClusterProcess(uint cNo, uint* genes, uint genesSize, vector<probe_set*>* data, uint* expts, uint exptNo, bool lN, set<void*>* cs, QMutex* cM){
  k = cNo;
  N = exptNo;
  geneNo = 0;
  someEmpty = true;                     // we should probably set a flag to give the user an option as to allow empty clusters or not. 

  selectedExperiments = expts;
  probeIndices = new uint[genesSize];
  points = new float*[genesSize];
  centerDistances = new float[genesSize];
  membership = new uint[genesSize];
  localNorm = lN;
  clusterMutex = cM;
  parentClusters = cs;

  // k based parameters are initialised later on so that I can decrease k if it is too big.. 

  //// first go through and initialise the data in points. After that we can use the data to 
  ///  to create some random points..
 
  // first create a min and max array with N things in them
  min = new float[N];
  max = new float[N];
  for(int i=0; i < N; i++){ min[i] = max[i] = 0; }

  // also need some temporary values which we can normalise..
  float* temp;
  if(localNorm){
    temp = new float[N];
  }else{
    temp = new float[(*data)[0]->allExptNo];
  }
  // and then go through the genes and work which ones have good data we can use.. 
  
  probe_set* p;
  for(int i=0; i < genesSize; i++){
    if(genes[i] > data->size()){
      continue;
    }
    p = (*data)[genes[i]];
    if(!p->index){
      continue;
    }
    
    ////  then just check whether or not it has data for all of the experiments
    bool ok = true;
    for(int j=0; j < N; j++){
      if(expts[j] > p->allExptNo){
	cerr << "KClusterProcess initialiser.. experiment Index is out of range " << endl;
	ok = false;
	break;
      }
      if(p->exptLookup[expts[j]] == -1){
	ok = false;
	break;
      }
    }
    if(!ok){
      continue;
    }
    // OK, if I got here, I should be Ok, and I can start collecting the appropriate values..
    float* meanValues = new float[N];     // the mean values
    for(int j=0; j < N; j++){ meanValues[j] = 0; }
    
    if(localNorm){
      for(int j=0; j < p->probeNo; j++){
	for(int l=0; l < N; l++){
	  temp[l] = p->probes[j][p->exptLookup[expts[l]]];
	}
	normalise(temp, N);
	for(int l=0; l < N; l++){
	  meanValues[l] += temp[l];
	}
      }
    }else{
      for(int j=0; j < p->probeNo; j++){
	for(int l=0; l < p->exptSize; l++){
	  temp[l] = p->probes[j][l];
	}
	normalise(temp, p->exptSize);
	for(int l=0; l < N; l++){
	  meanValues[l] += temp[p->exptLookup[expts[l]]];
	}
      }
    }
    //// Now we have mean values, but it actually contains the sum values. But it doesn't matter as we will normalise shortly
    normalise(meanValues, N);
    points[geneNo] = meanValues;
    probeIndices[geneNo] = p->index;
    membership[geneNo] = k;
    
    // also work out whether or nor min and max need to be updated.
    for(int j=0; j < N; j++){
      if(min[j] > meanValues[j]){ min[j] = meanValues[j]; }
      if(max[j] < meanValues[j]){ max[j] = meanValues[j]; }
    }
    // INCREMENT THE GENE COUNTER
    geneNo++;
  }
  // if k is larger then geneNo, then we are really screwed, no good. so lets check and make sure..
  if(k > geneNo){
    k = geneNo;
  }
  maxDistances = new float[k];
  clusterSizes = new uint[k];
  centers = new float*[k];
  clusters = new uint*[k];
  for(int i=0; i < k; i++){
    clusterSizes[i] = 0;
    clusters[i] = 0;
  }

  // now the data is all initialised. hopefully it's OK. need to initialise the cluster Centers.
  initialiseCenters(min, max);
  // delete the temp values..
  delete []temp;
  delete []genes;
  minChanges = geneNo;
  noImprovementCount = 0;
}

KClusterProcess::~KClusterProcess(){
  cout << "Deleting a KClusterProcess " << endl;
  for(int i=0; i < geneNo; i++){
    delete []points[i];
  }
  delete []points;
  
  for(int i=0; i < k; i++){
    delete []clusters[i];
    delete []centers[i];
  }
  delete []clusters;
  delete []centers;
  delete []selectedExperiments;
  delete []probeIndices;
  delete []clusterSizes;
  delete []membership;
  delete []min;
  delete []max;
  delete []centerDistances;
  delete []maxDistances;
}

void KClusterProcess::initialiseCenters(float* min, float* max){
  // assume that the first dimension of the centers has already been allocated.
  srand(time(0));
  for(int i=0; i < k; i++){
    centers[i] = new float[N];
    for(int j=0; j < N; j++){
      centers[i][j] = min[j] + ( (max[j]-min[j]) * ((float)rand()/(float)RAND_MAX) );
    }
  }
}

void KClusterProcess::reallocateEmptyCenters(){
  // order clusters containing more than one member, i.e. not empty by their size,
  // allocate empty clusters to a random member within the largest clusters.. 
  // well, one for each.. 
  //cout << "beginning of reallocateEmptyCenters " << endl;
  // first we need to get a sorted index for the clusters based on the cluster max distances..
  int* tempClusterIndex = new int[k];
  for(int i=0; i < k; i++){ tempClusterIndex[i] = i; }
  //  cout << "just before sorting.. what is going on" << endl;
  quickSort(tempClusterIndex, maxDistances, 0, k-1);   // 0 now contains the smallest.. no good.. we want reverse order.. 
  //cout << "tempclusters index made and sorted " << endl;

  srand(time(0));
  int emptyCount = 0;
  int clusterIterator=0;
  for(int i=0; i < k; i++){
    if(clusterSizes[i] == 0){
      emptyCount++;
      // choose the appropriate cluster..
      while(maxDistances[tempClusterIndex[k-1-(clusterIterator % k)]] == 0){
	clusterIterator++;
      }
      //cout << "i : " << i << "  iterator " << clusterIterator << endl;
      // get a random integer in the range of the size of the chosen cluster..
      int clusterChoice = tempClusterIndex[k-1-(clusterIterator % k)];      // just for clarity, not for speed.. 
      clusterIterator++;
      //cout << "cluster Choice is " << clusterChoice << endl;
      int choice = clusterSizes[clusterChoice] * ((float)(rand()-1)/(float)RAND_MAX);   // will fail if,, rand() = RAND_MAX
      //cout << "And choice is     " << choice << endl;
      //cout << "which turns out to be index " << clusters[clusterChoice][choice] << endl;
      //      centers[i] = points[clusters[clusterChoice][choice]];
      for(int j=0; j < N; j++){
	centers[i][j] = points[clusters[clusterChoice][choice]][j];
	//	centers[i][j] = min[j] + ( (max[j]-min[j]) * ((float)rand()/(float)RAND_MAX) );
      }
    }
  }
  delete []tempClusterIndex;
  //cout << "reallocate Empty Centers : emptyCount is " << emptyCount << endl;
  if(emptyCount == 0){
    someEmpty = false;
  }
}

void KClusterProcess::normalise(float* v, int s){
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

float KClusterProcess::euclidean(float* v1, float* v2, int s){
  float e = 0;
  if(s <= 0){
    return(-1);
  }
  for(int i=0; i < s; i++){
    e += ((v1[i] - v2[i]) * (v1[i] - v2[i]));
  }
  return(sqrt(e)/s);
}

void KClusterProcess::run(){
  int* clusterCounts = new int[k];
  int maxIterations = 50;        // number of useless iterations.,, i.e. without an improvement.. 
  while(someEmpty){
    while(allocate()){
      calculateCenters();
      if(noImprovementCount > maxIterations)
	break;
    }
    // allocate the members in the clusters structure.. 
    // again first dimension already allocated. 
    //    cout << "finished calculating centers and am allocating other stuff" << endl;
    for(int i=0; i < k; i++){
      if(clusters[i] != 0){
	delete []clusters[i];
      }
      clusters[i] = new uint[clusterSizes[i]];
      clusterCounts[i] = 0;
      maxDistances[i] = 0;
    }
    //cout << "Set lots of values to 0" << endl;
    for(int i=0; i < geneNo; i++){
      clusters[membership[i]][clusterCounts[membership[i]]] = i;
      clusterCounts[membership[i]]++;
      if(centerDistances[i] > maxDistances[membership[i]]){
	maxDistances[membership[i]] = centerDistances[i];
      }
    }
    reallocateEmptyCenters();          // this way I can do some statistics on the thing.. 
  }
  // lock the clusterMutex, then insert our address into the thing..
  delete []clusterCounts;
  clusterMutex->lock();
  parentClusters->insert((void*)this);
  clusterMutex->unlock();
  // and that is then it.
}

int KClusterProcess::allocate(){
  int changes = 0;
  // set all clusterSizes to 0
  for(int i=0; i < k; i++){ clusterSizes[i] = 0; }
  
  for(int i=0; i < geneNo; i++){
    float minDistance = euclidean(points[i], centers[0], N);
    int minIndex = 0;
    for(int j=0; j < k; j++){
      float d = euclidean(points[i], centers[j], N);
      if(minDistance > d){
	minDistance = d;
	minIndex = j;
      }
    }
    if(minIndex != membership[i]){
      changes++;
      membership[i] = minIndex;
    }
    centerDistances[i] = minDistance;         // even if there is no change in the index, the center may have moved.. !! 
    clusterSizes[minIndex]++;
  }
  cout << "changes : " << changes << endl;
  if(changes < minChanges){
    minChanges = changes;
    noImprovementCount = 0;
  }else{
    noImprovementCount++;
  }
  return(changes);
}

void KClusterProcess::calculateCenters(){
  // set all centers to 0.
  for(int i=0; i < k; i++){
    if(clusterSizes[i]){
      for(int j=0; j < N; j++){
	centers[i][j] = 0;
      }
    }
  }
  // then go through all of the genes and move the centers around..
  for(int i=0; i < geneNo; i++){
    for(int j=0; j < N; j++){
      centers[membership[i]][j] += points[i][j];
    }
  }
  // and then go through and divide by the cluster sizes...
  for(int i=0; i < k; i++){
    if(clusterSizes[i]){
      for(int j=0; j < N; j++){
	centers[i][j] = centers[i][j]/(float)clusterSizes[i];
      }
    }
    //    cout << "cluster " << i << "   size : " << clusterSizes[i] << endl;
  }
}

void KClusterProcess::swapElements(int* a, int* b){
  int temp = *a;
  *a = *b;
  *b = temp;
}

int KClusterProcess::divideArray(int* a, float* v, int left, int right){
  int mid = left;
  while(1){
    while(right > left && v[a[right]] > v[a[mid]])
      right--;
    while(left < right && v[a[left]] <= v[a[mid]])
      left++;
    if(left == right)
      break;
    swapElements(&a[left], &a[right]);
  }
  if(v[a[mid]] > v[a[left]])
    swapElements(&a[mid], &a[left]);
  return(left);
}

void KClusterProcess::quickSort(int* indices, float* values, int left, int right){
  int split;
  if(left < right){
    split = divideArray(indices, values, left, right);
    quickSort(indices, values, left, split-1);
    quickSort(indices, values, split+1, right);
  }
  return;
}
