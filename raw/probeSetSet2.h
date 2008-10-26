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

#ifndef PROBESETSET2_H
#define PROBESETSET2_H

#include "probe_set.h"
#include "dataStructs.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <qmutex.h>

using namespace std;    //??


class ProbeSetSet2
{
 public:
  ProbeSetSet2(const char* cinfo, const char* dataTable="data");            // initialiase from database..
  ProbeSetSet2();                         // do nothing??

  vector<probe_set*> data;                // lots of work to fix, but what the hell
  probe_set* emptySet;                    // a stupid kludge to overcome some trouble.. 
  vector<probe_data> probeData;           // can potentially change as a result of user interaction. may need mutexes. -be careful.
                                          // but it seems to be ok at the moment. 

  map<float, exInfo> experiments;
  //  map<int, string> chipDescriptions;
  map<int, chipInfo> chipDescriptions;
  QMutex* chipMutex;

  /// Ensembl genome data. Annotation for the genes..
  map<int, multimap<int, string> > ensemblAnnotation;
  map<int, string> ensemblFields;   // the fields in the ensemblAnnotation.. to save space.. 

  //// Linkage to ensembl gene indices..
  map<int, int> probeSetEnsemblIndex;   // probe Set index to ensembl Index
  multimap<int, int> ensemblProbeSetIndex;  // the other way around.. 
  
  // User annotation and sessions..
  map<int, set<int> > sessionLookup;      // map<probe_set-index, sessionIndex>
  map<int, sessionInformation> sessions;  // the index, followed by the information.
  QMutex* sessionMutex;
  
  // and annotation in the same vein.
  map<int, set<int> > annotationLookup;   // key is the gene index, and the set contains the sessions. 
  map<int, annotationInformation> userAnnotation;
  QMutex* annotationMutex;

  // user Information.. 
  map<int, userInformation> userTable;
  QMutex* userMutex;

  map<string, ifstream*> chromFiles;    // for reading sequences from files using chromosomal positions. be careful to protect this with mutexes. 
  QMutex* chromFileMutex;
  /// the two below are not mutable, so no need for any mutexes I hope.. 
  map<string, int> chromSizes;          // the chromosome sizes in the crhomFiles.. !!! 
  map<string, chromAnnotation*> chromosomeAnnotation;
  /// maybe we should add methods for changing this annotation..
  QMutex* chromAnnotationMutex;   // not terribly fine grained, but I doubt it will happen very often.. 

  /// in situ probe data.. 
  map<int, ishProbeData> ishProbes;     // contains some information about in situ probes, and pointers to genomic blast matches.. so the whole thing can be checked out..
  set<string> ishTextFields;
  set<string> ishFloatFields;
  set<string> ishClasses;     // so we can easily tell clients what these are,, just remember to update them when necessary... 
  QMutex* ishProbeDataMutex;             // lock if we want to change this information.. -- 

  const char* conninfo;       // 
 private:
  void setExpData(const char* conninfo);       // sets the data about the experiments
  void setData(const char* conninfo, const char* tableName="data");          // sets the actual data.. sounds painful. I'm not sure how the arrays will return..
  void guessGenes();                    // guess linkaged between the probe set matches and the genes. painful .. and inaccurate.. probably.. 
};

#endif
