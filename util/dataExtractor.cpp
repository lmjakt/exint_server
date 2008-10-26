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

#include "dataExtractor.h"
#include "normaliser.h"

using namespace std;   // some vectors and maps used in probe_set.. 

float* ExData::mean(){
  if(exptNo == 0){
    return(0);
  }
  float* meanValues = new float[exptNo];
  for(int i=0; i < exptNo; i++){ meanValues[i] = 0; }
  float n = (float)probeNo;   // so I can just cast once.. 
  for(int i=0; i < probeNo; i++){
    for(int j=0; j < exptNo; j++){
      meanValues[j] += (values[i][j]/n);
    }
  }
  return(meanValues);     // User has to delete the structure or have a memory leak..
}

float** DataExtractor::copyProbes(float** p, int ps, int es){
  //cout << "copy Probes ps : " << ps << "  es : " << es << endl;
  float** cp = new float*[ps];
  for(int i =0; i < ps; i++){
    //cout << "\t\t\t\t\tcopy Probes i : " << i << endl;
    cp[i] = new float[es];
    for(int j=0; j < es; j++){
      cp[i][j] = p[i][j];
    }
  }
  return(cp);
}

ExData* DataExtractor::extract(probe_set* ps, uint* expts, uint es, bool requireAll){
  // make a new ExData, and fill this..
  ExData* ed = new ExData(ps->probeNo, es);   //

  // ok, find out the indices of the experiment indices that we want to use.
  uint* exIndex = new uint[es];
  uint selEx = 0;          // the number of selected experiments.
  for(int i=0; i < es; i++){
    //    cout << "looking for experiments : " << expts[i] << endl;
    if(expts[i] < ps->allExptNo){
      if(ps->exptLookup[expts[i]] != -1){
	exIndex[selEx] = ps->exptLookup[expts[i]];
	ed->expts[selEx] = expts[i];    // ok.. 
	selEx++;
	//ed->exptNo++;  -- leave this as 0, until we've filled in the data.. or return it.. 
      }
    }
  }
  // ok, so now we now how many experiments we are going to be able to get data from, so let's use that data and
  // extract the information..
  // if all required then return here to make things a little bit faster..
  if(requireAll && selEx < es){
    delete []exIndex;
    return(ed);    // which will be empty, I could delete it and return a 0, which would probably be faster to check for..
                   // but the c notation for checking this is difficult to read. -At least for me.. 
  }
  
  ed->exptNo = selEx;
  
  for(int i=0; i < ps->probeNo; i++){
    for(int j=0; j < selEx; j++){
      ed->values[i][j] = ps->probes[i][exIndex[j]];
    }
  }
  // ok, we are done,,
  // let's delete exIndex, and return ed..
  delete []exIndex;
  return(ed);
}

void deleteValues(float** v, uint s){
  for(int i=0; i < s; i++){
    delete []v[i];
  }
}

ExData* DataExtractor::globalNorm(probe_set* ps, uint* expts, uint es, int type, bool requireAll){
  float** values = copyProbes(ps->probes, ps->probeNo, ps->exptSize);
  // then do some normalisation..
  Normaliser norm;
  if(type == 1){
    norm.zScore(values, ps->probeNo, ps->exptSize);
  }
  if(type == 2){
    norm.mScore(values, ps->probeNo, ps->exptSize);
  }
  // and then do exactly as before, but replace ps->probes with values..
  ExData* ed = new ExData(ps->probeNo, es);   //

  // ok, find out the indices of the experiment indices that we want to use.
  uint* exIndex = new uint[es];
  uint selEx = 0;          // the number of selected experiments.
  for(int i=0; i < es; i++){
    if(expts[i] < ps->allExptNo){
      if(ps->exptLookup[expts[i]] != -1){
	exIndex[selEx] = ps->exptLookup[expts[i]];
	ed->expts[selEx] = expts[i];    // ok.. 
	selEx++;
	//ed->exptNo++;  -- leave this as 0, until we've filled in the data.. or return it.. 
      }
    }
  }
  // ok, so now we now how many experiments we are going to be able to get data from, so let's use that data and
  // extract the information..
  // if all required then return here to make things a little bit faster..
  if(requireAll && selEx < es){
    delete []exIndex;
    selEx = 0;     // hmm... 
    deleteValues(values, ps->probeNo);
    return(ed);    // which will be empty, I could delete it and return a 0, which would probably be faster to check for..
                   // but the c notation for checking this is difficult to read. -At least for me.. 
  }
  
  ed->exptNo = selEx;
  
  for(int i=0; i < ps->probeNo; i++){
    for(int j=0; j < selEx; j++){
      ed->values[i][j] = values[i][exIndex[j]];
    }
  }
  // ok, we are done,,
  // let's delete exIndex, and return ed..
  delete []exIndex;
  deleteValues(values, ps->probeNo);
  return(ed);
}
