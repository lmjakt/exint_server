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

#ifndef EUCLIDSORTPROCESSOR_H
#define EUCLIDSORTPROCESSOR_H

#include <qthread.h>
#include <qmutex.h>
#include <qwaitcondition.h>
//#include <iostream>
#include <vector>
#include <set>
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "connectionobject.h"        // just for dist_set, not so good hmm. 

using namespace std;


class EuclidSortProcessor : public QThread 
{
  public :
  // no signals and slots, so no Q_OBJECT. keep it simple..
  EuclidSortProcessor(ProbeSetSet2* ps, probe_set* p, uint st, uint sp, vector<dist_set>* sc, set<int>& chips, vector<uint>& expts);  // nothing for the moment..
  ~EuclidSortProcessor();          // delete the things that need it !!

  
  protected :
    void run();

  private :
    ProbeSetSet2* pSet;
  probe_set* ps1;      // the one we are comparing everything against. 
  //vector<int>* indices;
  vector<dist_set>* scores;
  uint begin;
  uint end;
  
  uint* experiments;    // the experiments that we are using
  uint experimentNumber; // the number of experiments   -- use instead of the vector since the functions have been written for this.
  set<int> clientChips;  // the arrays that we are interested in.


  float heavyCompare(probe_set* pset1, probe_set* pset2, uint* expts, uint es);    // can I wake it up like this ??
  void zScore(float* v, int s);
  float** copyProbes(float** p, int ps, int es);   // ps- probe no,  es -experiment size, 
  float euclidean(float* v1, float* v2, int s);
  void delProbes(float** p, int ps);

};

#endif
