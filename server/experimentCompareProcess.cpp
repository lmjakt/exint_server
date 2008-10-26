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

#include "experimentCompareProcess.h"
#include "../raw/probe_set.h"
#include <vector>
#include <qthread.h>
#include <math.h>
#include <qmutex.h>
#include <set>

using namespace std;

//typedef unsigned int uint;

ExperimentCompareProcess::ExperimentCompareProcess(uint* gs, uint gsSze, vector<probe_set*>* data, uint* expts, uint exptN, set<void*>* processes, QMutex* mutex){
  genes = gs;         // not the actual genes, but rather the probe set indices..
  geneNo = gsSze;    // number of genes.. (unsigned int)..

  probes = data;      // just a pointer.. 

  experiments = expts;  // the experiment indices.. -unfortunately these have to be checked for every probe set.. 
  exptNo = exptN;       // this as a result of trying to reduce the dependancies in the program.. (i.e. no central index which might be wrong)

  compMutex = mutex;
  compProcesses = processes;

  cout << "ExperimentCompareProcess set some stuff " << endl;
  cout << "exptNo is " << exptNo << endl;
  // let's set up the distances struct.. and 0 it..
  distances = new float*[exptNo];
  for(uint i=0; i < exptNo; i++){
    cout << "\t" << i  << endl;
    distances[i] = new float[exptNo];  // which should all be set to 0..
    for(uint j=0; j < exptNo; j++){
      distances[i][j] = 0;
    }
  }
  // presumably there is a faster way of doing this, but who knows eh.. 
  // this is actually twice as big as it needs to be, but what can you do... eh.. 
  // at the end we'll fill in the blanks, maybe,, so it can be looked up in both directions..
  cout << "end of experimentCompareProcess constructor" << endl;
}

ExperimentCompareProcess::~ExperimentCompareProcess(){
  // remember to delete lots of things here.. hoo hoo yeahh.
  cout << "Experiment Compare Process destructor " << endl;
  for(uint i=0; i < exptNo; i++){
    delete []distances[i];
  }
  delete []distances;
  
  //////  delete pointers passed in the constructor.  -- this may not be very orthodox,, -but these structures should not be shared with other
  /////  threads (apart from the data pointer which we leave alone).. ------ This is less flexible,, but user beware..
  delete []genes;
  delete []experiments;
  cout << "Experiment compare process destructor done " << endl;
}

// a couple of useful functions.. 
void ExperimentCompareProcess::normalise(float** v, uint s1, uint s2){
  // take normalise each series of values.. (i.e. outer dimension is the probe number, inner dimension is the experiments.)
  if(s2 < 2){
    cerr << "ExperimentCompareProcess normalise function, s2 is smaller than 2, not good" << endl;
    return;
  }
  float std;
  float mean;
  float xsum;
  float xsquaredSum;
  float SS;
  for(uint i=0; i < s1; i++){       // for each probe pair profile
    // some notes about the algorithm..
    // Note..      SS = sum of squared devieates (squared sum)
    // AND         SS = sum(X^2) - ((sum(X))^2 / N)
    // so          std = sqrt( SS / N-1)
    // which should be reasonably fast to calculate.. (rather than using the std function.. (which I could rewrite)
    //float std;
    //float mean;
    xsum = 0;
    xsquaredSum = 0;
    for(int j=0; j < s2; j++){
      xsum += v[i][j];
      xsquaredSum += (v[i][j]*v[i][j]);
    }
    mean = xsum/s2;
    SS = xsquaredSum - (xsum*xsum)/(float)s2;
    std = sqrt(SS/(float)(s2-1));
    // then go through and  modify
    for(int j=0; j < s2; j++){
      v[i][j] = (v[i][j]-mean)/std;
    }
  }
  /// which should be everything needed here.. 
}

void ExperimentCompareProcess::mean(float** v, float* m, uint s1, uint s2){
  // first make sure all the values in m are 0..
  for(uint i=0; i < s2; i++){
    m[i] = 0;
  }
  if(s1 < 1){
    return;
  }
  // then simpley do some addition..
  for(uint i=0; i < s1; i++){
    for(uint j=0; j < s2; j++){
      m[j] += v[i][j];
    }
  }
  // WE do NEED to take the mean values,, rather than the sums, because we cannot assume that we will have the
  // the same number of probes for each probe set. Still, this doesn't matter too much..
  for(uint i=0; i < s2; i++){
    m[i] = m[i]/s1;
  }
}

float** ExperimentCompareProcess::growValueArray(float** v, uint oldSize, uint newSize){
  float** nv = new float*[newSize];
  for(uint i=0; i < newSize; i++){
    nv[i] = new float[exptNo];
  }
  // delete the old set of values..
  for(uint i=0; i < oldSize; i++){
    delete []v[i];
  }
  delete []v;
  return(nv);
}


// And I think that should be pretty much everything that we need to do at the moment..
// Let's do the comparison .... in the main run() function...

void ExperimentCompareProcess::run(){
  // alleluliah.. -- lets see what we can do..
  
  // basic process is as follows..
  //
  // 1. find probe set (i.e. iterate over the probe sets)
  // 2. for each probe set extract the relevant experimental data by looking up the values
  //    in the experiments array, and find the relevant indices in the lookup index..
  //    if any of the values are -1, then go to next.. 
  // 3. normalise the probe pair profiles extracted.. 
  // 4. obtain the mean for the normalised data..
  // 5. use the mean profile to increment the distances in the 2d distance matrix.

  // should be pretty straightforward..

  // first allocate some memory for some things we'll need..
  //  cout << "beginning of experimentcompare process run function" << endl;

  float* meanValues = new float[exptNo];       // anyway, this one is simple..
  int* exptIndex = new int[exptNo];   // found in the exptLookup array in the thingy.. 

  int probeAllocated = 20;                     // number of probes allocated for the thingy..
  float** values = new float*[probeAllocated]; 
  for(uint i=0; i < probeAllocated; i++){
    values[i] = new float[exptNo];             // where we extract the actual numbers to.. if to small then we can grow..
  }

  //cout << "allocated stuff in the experimentcompare process thingy" << endl;
  // ok, let's go through the actual data, and see how we get along...
  probeCounter = 0;    // number of probes used in the comparison.
  for(uint i=0; i < geneNo; i++){
    //cout << "checking for gene i: " << i << "  " << genes[i] << endl;
    // first check if we have a probe Set for this one..
    if(genes[i] >= probes->size()){
      cerr << "ExperimentCompareProcess run() received an probe index out of range, ignoring" << endl;
      continue;    // go to next one..
    }
    probe_set* p = (*probes)[genes[i]];    // for convenience in writing and reading..
    //    cout << "defined p" << endl;
    /// now go through the expts, and see if they have been well allocated or not..
    bool exptOk = true;      // ok unless we have some sort of error.. 
    for(uint j=0; j < exptNo; j++){
      // check that the value pointed to by experiments isn't out of bounds..
      if(experiments[j] >= p->allExptNo){
	exptOk = false;
	break;
      }
      exptIndex[j] = p->exptLookup[experiments[j]];
      if(exptIndex[j] == -1){
	exptOk = false;
	break;
      }
    }
    //cout << "defined the exptIndex thingy " << endl;
    if(!exptOk){
      cout << "exptOk is not ok " << endl;
      continue;    // go to top of loop.. 
    }
    // now check if we have to grow the values.. 
    if(p->probeNo > probeAllocated){
      cout << "trying to grow the values : " << endl;
      values = growValueArray(values, probeAllocated, p->probeNo);   // automatically deletes the old values.. ok.. 
      probeAllocated = p->probeNo;
    }
    //cout << "ok have checked the values " << endl;
    // should be ok.. ok.. let's collect numbers.. 
    for(uint j=0; j < p->probeNo; j++){
      for(uint k=0; k < exptNo; k++){
	// and this shouldn't require any more checking.. 
	values[j][k] = p->probes[j][exptIndex[k]];
      }
    }
    //cout << "collected the numbers " << endl;
    // Normalise the values..
    normalise(values, p->probeNo, exptNo);
    // Assign the mean values
    mean(values, meanValues, p->probeNo, exptNo);
    //cout << "normalised and meaned" << endl;
    /// ok, at this point all we have to do, is go through and increment all the distances..
    for(uint j=0; j < exptNo; j++){
      for(uint k=(j+1); k < exptNo; k++){        // this takes half the time, but note that we'll need a separate thingy to copy numbers.. later
	distances[j][k] += ((meanValues[j] - meanValues[k]) * (meanValues[j] - meanValues[k]));
	// could do the other way round here as well, but quicker to do a loop after going through everything..
      }
    }
    //cout << "done the comparison " << endl;
    // at this point we're pretty much done for this probe set..
    probeCounter++;
  }
  ///// Now copy the table to make it symmetrical. Not sure that it is really needed but for lookups, it is probably the easiest thing to do.
  /// and get the square root for each function.. 

  for(uint i=0; i < exptNo; i++){
    for(uint j= (i+1); j < exptNo; j++){
      //distances[i][j] = sqrt(distances[i][j]);
      distances[j][i] = distances[i][j];
    }
  }

  //// print out,, just temporary, until I work out how to do things.. 
//   cout << "printing distances between different experiments : " << endl;
//   for(int i=0; i < exptNo; i++){ cout << "\t" << experiments[i]; }
//   cout << endl;
//   for(uint i=0; i < exptNo; i++){
//     cout << i << " : " << experiments[i];
//     for(uint j=0; j < exptNo; j++){
//       cout << "\t" << distances [i][j];
//     }
//     cout << endl;
//   }
  //////////////  Alleluliah, at this point we'd probably need to alter some structure somewhere, and see how the thing happens, but for 
  /////////////   now, let's worry about the destructor element..
  delete []meanValues;
  delete []exptIndex;
  for(uint i=0; i < probeAllocated; i++){
    delete []values[i];
  }
  delete []values;
  
  // since we are done..
  compMutex->lock();
  compProcesses->insert((void*)this);
  compMutex->unlock();
}
