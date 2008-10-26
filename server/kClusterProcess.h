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

#ifndef KCLUSTERPROCESS_H
#define KCLUSTERPROCESS_H

#include <qthread.h>
#include <qmutex.h>
#include <iostream>
#include <vector>           // only for initialisation of clusters.. which doesn't take much processing. 
#include <set>
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"

using namespace std;

class KClusterProcess : public QThread
{
  public :
    KClusterProcess(uint cNo, uint* genes, uint genesSize, vector<probe_set*>* data, uint* expts, uint exptNo, bool lN, set<void*>* cs, QMutex* cM);
  ~KClusterProcess();

  // some things that I may need to access. Ofcourse, I should make accessor functions, but that's just extra work. 

  bool localNorm;            // use local or global normalisation.. 
  int k;                     // the number of clusters, equal to cNo
  uint N;                     // the dimensionality of the space, equal to the exptNo
  uint geneNo;                // the actual number of genes clustered... 
  uint* selectedExperiments;  // the experimental indices, a total of N indices

  float** points;            // the data points,, first dimension has geneNo points, second has N points. 
  float** centers;           // the cluster centers.. first dimension has k points, second has N points.
  
  uint* clusterSizes;         // the sizes of the clusters,, 
  uint** clusters;             // the members of the clusters..  
  uint* probeIndices;         // the probe Indices clustered.. -so a total of geneNo probeIndices..

  protected :
    void run();

  private :
  void initialiseCenters(float* min, float* max);       // length is N anyway... 
  void reallocateEmptyCenters();                        // if we have empty centers, then reallocate these.. 
  void normalise(float* v, int s);
  int allocate();
  void calculateCenters();
  float euclidean(float* v1, float* v2, int s);

  // I am not sure yet. but I think I may need some methods for sorting.. so let's implement a quicksort..
  void swapElements(int* a, int* b);
  int divideArray(int* a, float* v, int left, int right);
  void quickSort(int* indices, float* values, int left, int right);

  // and then the data I need to work on this..

  bool someEmpty;            // true if some of the clusters are empty.. 
  uint* membership;           // the membership,, the indices are from probeIndices
  
  QMutex* clusterMutex;      // for changing the parentClusters.. 
  set<void*>* parentClusters; // to pass on myself in the afterlife..
  
  float* min;
  float* max;
  float* centerDistances;   // optional. maybe later I will store the distances from the centers.
  float* maxDistances;      // max distance for each cluster... 
  int minChanges;           // the minimum number of changes. set to geneNo initiallay
  int noImprovementCount;       // the number of iterations without an improvement.. 
};

#endif
