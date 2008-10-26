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

#ifndef DATAEXTRACTOR_H
#define DATAEXTRACTOR_H

#include "../raw/probe_set.h"
#include <iostream>

struct ExData {        // contains some some extracted expression data..
  float** values;
  uint* expts;
  uint exptNo;
  uint probeNo;
  uint targetNo;
  // and a constructor..
  ExData(uint pNo, uint exNo){  // the number of experiments and probe pairs contained
    exptNo = 0;                 // the actual number of experiments -- entered. set externally by the creator.
    probeNo = pNo;
    targetNo = exNo;   // the number of experiments we want to have..
    expts = new uint[targetNo];
    values = new float*[pNo];
    for(int i=0; i < probeNo; i++){
      values[i] = new float[targetNo];
    }
  }
  float* mean();       // return the mean value. Convenient to have here.. 
                       // NOTE, it is the users responsibility to delete this,, 
  // there isn't any easy way that I can keep track of how many experiments have been added without 
  // requiring that all the probes be added as an array, this would anyway, take too much time.. 
  ~ExData(){
    //cout << "calling delete on expts. " << endl;
    delete []expts;
    for(int i=0; i < probeNo; i++){
      //cout << "\tdelete on values [i] " << endl;
      delete []values[i];        // don't know if I need to do this, but seems reasonable... might crash on it though..
    }
    //cout << "done the inner loop" << endl;
    delete []values;
    //cout << "and now the outer as well " << endl;
  }
  // use with care.. manipulate the data structs externally.. as this will be faster, and anyway, we don't want to use this in too many places.
};
  

class DataExtractor {
  public :
    DataExtractor(){}
  ~DataExtractor(){}     // this class just has some functions.. declared as a class to make sure it occupies its
  // own memory, and thus should be re_entrant and thread safe. --i.e. each class has its own 
  // copy, -- there are probably other ways of doing this, but this is simple..
  
  ExData* extract(probe_set* ps, uint* expts, uint es, bool requireAll=false);    // return the data if OK.. 
  ExData* globalNorm(probe_set* ps, uint* expts, uint ex, int type, bool requireAll=false);  // return globally normalised data
  // 1 for zScore, 2 for mScore.. 
  // Note that it is the users responsibility to delete the pointer if not wanted. 
  // Note that the user can check if it worked or not by checking the value of selEx, if this, is less than the experiment number, then
  // it didn't work. But, -what's the point of the option then ??? can't remember, seems a little strange. 

  private :
    float** copyProbes(float** p, int ps, int es);   // ps- probe no,  es -experiment size, 

};

#endif
