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

#ifndef ANOVAPROCESSOR_H
#define ANOVAPROCESSOR_H

#include <qthread.h>
#include <qmutex.h>
#include <qwaitcondition.h>
//#include <iostream>
#include <vector>
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "connectionobject.h"        // just for dist_set, not so good hmm. 

using namespace std;


class AnovaProcessor : public QThread 
{
  public :
  // no signals and slots, so no Q_OBJECT. keep it simple..
  AnovaProcessor(ProbeSetSet2* p, uint st, uint sp, vector<uint>* ind, vector<dist_set>* sc);  // nothing for the moment..
  AnovaProcessor(ProbeSetSet2* p, uint st, uint sp, vector<uint>* ind, vector<dist_set>* sc, uint* exptS, uint exptSS);  // nothing for the moment..
  ~AnovaProcessor();          // delete the things that need it !!

  
  protected :
    void run();

  private :
    //QMutex aMutex;      // for pushing the scores thingy.. (reserve space in it!!). -- but I think it is wrong!! shouldn't be there. -- completely wrong. 
    ProbeSetSet2* pSet;
  vector<uint>* indices;
  vector<dist_set>* scores;
  uint* exptSelection;
  uint exptSelectionSize;
  bool userSelectedExperiments;       // did the user specify the experiments or not??
  int begin;
  int end;

  float anova(probe_set* pSet, uint* expts, uint es);    // can I wake it up like this ??
  void zScore(float* v, int s);
  float** copyProbes(float** p, int ps, int es);   // ps- probe no,  es -experiment size, 
  void delProbes(float** p, int ps);

};

#endif
