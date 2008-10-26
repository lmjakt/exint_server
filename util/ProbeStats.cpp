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

#include "ProbeStats.h"

// the constructor is just empty so we just need to put in the functions here.

float* ProbeStats::devFromMean(ExData* data){
  // check we have some experiments in the data set..
  if(data->exptNo < 3 || data->probeNo < 3){    // not worthwhile to do anytyhing
    cerr << "Probe Stats dev From Mean returning null pointer due to no data or no expts.. " << endl;
    return(0);
  }
  float* devs = new float[data->probeNo];
  //for(int i=0; i < data->probeNo; i++){ devs[i] = 0; }
  float* mean = new float[data->exptNo];
  for(int i=0; i < data->exptNo; i++){ mean[i] = 0; }
  // could probably do the above with memset, but not sure exactly how 0 is encoded.. 
  float sigma = 0;    // the standard deviation.. or soemthing..
  // let's normalise the data..
  norm.zScore(data->values, data->probeNo, data->exptNo);   // modifies the original.
  // and let's make a mean profile
  for(int i=0; i < data->exptNo; i++){
    for(int j=0; j < data->probeNo; j++){
      //cout << data->values[j][i] << "\t";
      mean[i] += (data->values[j][i] / (float)data->probeNo);
    }
    //cout << endl;
  }
  norm.zScore(mean, data->exptNo);
  //cout << "mean : ";
  for(int i=0; i < data->exptNo; i++){
    //cout << mean[i] << "\t";
  }
  //  cout << endl;
  // we need to calculate euclidean distances so lets have an Euclid object called e
  Euclid e;  // don't know if this works or not.. 
  float n = float(data->probeNo - 1);
  for(int i=0; i < data->probeNo; i++){
    devs[i] = e.sqEuclidean(data->values[i], mean, data->exptNo);
    sigma += devs[i]/n;
    devs[i] = sqrt(devs[i]); // ok !!
  }
  //  cout << "sigma is : " << sigma << endl;
  sigma = sqrt(sigma);
  for(int i=0; i < data->probeNo; i++){
    //cout << i << "\t" << devs[i] << "----> ";
    devs[i] = devs[i]/sigma;
    //cout << devs[i] << endl;
  }
  return(devs);   // ok.. that's basically it.. hoho yeah.. 
}
