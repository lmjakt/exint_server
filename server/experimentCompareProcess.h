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

#ifndef EXPERIMENTCOMPAREPROCESS_H
#define EXPERIMENTCOMPAREPROCESS_H

#include <qthread.h>
#include "../raw/probe_set.h"
#include <vector>
#include <qmutex.h>
#include <set>

using namespace std;        // which apparently is a really bad thing to put in a header file.. 

// a thread that creates a structure that holds all agains all euclidean distances for
// a given set of experiments, and a given set of genes chosen by the user.

// -- given it's own thread as it will probably need to normalise all of the data and create a temporary
// -- data structure holding means of normalised data. -- i.e. it has to normalise all of the data.

// we can do this in two ways I suppose, normalise as we go along, or we can preallocate a huge amount of 
// data and then use that.. -- using less memory is probably good, and probably it is no more difficult..

class ExperimentCompareProcess : public QThread
{
  public :
    ExperimentCompareProcess(uint* gs, uint gsSze, vector<probe_set*>* data, uint* expts, uint exptN, set<void*>* processes, QMutex* mutex);   // we will probably need to pass some other data at some point as well. 
  ~ExperimentCompareProcess();
  //////// NOTE : that the destructor will destroy the pointers to gs and expts.. - so DO NOT SHARE these pointers 
  ////////        with any other processes, as this will cause memory violations and crashes. Use with care and read the documentation.


  /// private ..

  float** distances;           // a 2 d vector holding the values.. (there may be a more effective way of doing this, but I'm not sure..).
  uint* experiments;
  uint exptNo;                  // just the size of the experiments..
  int probeCounter;
  int probeNo(){
    return(geneNo);
  }
  /// terrible I know, but the overlying process can just read these and send to the client.. that's all, just for reading,, 

  private :
  uint* genes;    // the genes (gs.. ) -don't forget to delete.
  uint geneNo;
  vector<probe_set*>* probes;   // the probes.. access for stuff..



  QMutex* compMutex;
  set<void*>* compProcesses;

  void normalise(float** v, uint s1, uint s2);   // normalise the values in the 2d array (copied from the probes -- using the things..)..
  void mean(float** v, float* m,  uint s1, uint s2);  // mean value profile.. note s2 is the number of experiments, and m should be of s2 size.. -and is used for the thingy.  
  float** growValueArray(float** v, uint oldSize, uint newSize);  // if we should need to grow the value array to some new number.. 
  
  protected :
    void run();    // starts the new thread and everything.. 
  // think that might be all I need.. 
};

#endif
