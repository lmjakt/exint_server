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

#include "experimentTracer.h"
#include "pathTracer.h"
#include "dataExtractor.h"
#include <iostream>

using namespace std;

// this only has a constructor and a destructor. The only purpose is to extract 
// appropriate data for the path tracer and start it. The idea is to make the pathtracer
// itself sufficiently abstract such that we can put any appropriate set of n-dimensional
// points into it, whilst at the same time trying to minimise the amount of code we change 
// in connectionobject.cpp (which is horribly bloated and stuff at the moment anyway -->

ExperimentTracer::ExperimentTracer(vector<probe_set*> ps, vector<uint> expt, std::set<PathTracer*>* tracers, QMutex* mutex, float sm){
  DataExtractor de;
  uint* expts = new uint[expt.size()];
  for(uint i=0; i < expt.size(); i++){ expts[i] = expt[i]; }
  
  // ok life is a bit complicated as we have to transpose the resulting array of values which we are going to obtain.
  uint dimNo = 0;      // this is actually going to be the same as the number of probe sets which we get to use..
  float** values = new float*[expt.size()];
  cout << "creating the values 2d array. expt size : " << expt.size() << endl;
  for(int i=0; i < expt.size(); i++){
    values[i] = new float[ps.size()];   // ok.. hmmm 
  }
  cout << "transposing values ps size : " << ps.size() << endl;
  for(uint i=0; i < ps.size(); i++){
    cout << i ;
    ExData* data = de.globalNorm(ps[i], expts, expt.size(), 2, true);
    cout << "-";
    if(!data->exptNo){
      delete data;
      continue;
    }
    float* mean = data->mean();
    // We are now at dimension number dimNo, so go through the different experiments (values) and copy the appropriate value
    for(uint j=0; j < expt.size(); j++){
      //cout << "- " << j << " " << expts[j] << " -";
      values[j][dimNo] = mean[j];
    }
    cout << "| ";
    delete []mean;
    delete data;       // this is a little expensive, but anyway, this is not important in the grand scheme of things.
    dimNo++;
    //cout << endl;
  }
  cout << endl;
  // ok the values have been set up.. we could do some normalisation on these before we hand over to the thingy,
  // but I think that I will leave it at this...
  // I would also suggest that in the future we provide some string with some information regarding what we are doing, so that it 
  // can be retrieved from the data. .. but..
  cout << "Making pathtracer .. " << endl;
  PathTracer* tracer = new PathTracer(values, expts, expt.size(), dimNo, tracers, mutex, sm);
  cout << "starting path tracer " << endl;
  tracer->start();
  cout << "and going on with life as it were .. " << endl;

  /// and I leave off here, and see what happens when we change stuff for the better... or something like that...
}

ExperimentTracer::~ExperimentTracer(){
  cout << "ExperimentTracer deleting itself, but not doing much else... " << endl;
}
