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

#ifndef PROBE_SET_H
#define PROBE_SET_H

#include <vector>
#include <map>
#include <string>
#include <qsocket.h>
#include "../cluster/cluster.h"
#include "../util/util.h"

struct probe_set {
  vector< vector<float> > probes;   // probes[i][j] ; i-probe pair number, j-experiment index
  vector<bool> excluded;         // true if filtered... 
  vector<int> exptIndex;         // the experimental indices.. 
  map<int, int> exptLookup;       // map for reverse looking up.. 
  vector<float> mean;               // mean values, initially set to the full mean, but recalculated after filtering.. -
  int index;                        // probe set index -- related to the affymetrix thingy
  float distance;                   // for keeping distances.. -temporary measurement.. 
  float e_quality;                  // mean euclidean distance of all pairwise comparisons.. (does both directions, so a bit redundant..)
  float anovaScore;               // set dynamically depending on scores. initialised to 0.. 
  probe_set();
  probe_set(int i);                 // constructor -i is the the index
  probe_set(int i, vector< vector<float> >& p);  // complete constructor
  probe_set(int i, vector< vector<float> >& p, vector<int>& eindices);
  float euclidean();                        // overloaded..  compares itself.. 
  float euclidean(probe_set& p);            // overloaded of the below..
  float euclidean(probe_set& p, vector<int>& comp_points);      // gets mean euclidean distance with other probe sets.. 
  int get_index();
  void euclidean_filter(float cut_off);
  void unFilter();      // make all excluded to false..
  float concordance();                  // returns the composite euclidean distance of a transposed set.. 
  float sum();                          // simply returns the sum of the all of the probe sets.. 
  float meanMaxMeanDeviation();         // calculates the mean of the maximum deviations from the mean for each probe pair.
  float anova(vector<int>& expts);
  float selSum(vector<int>& expts); 
  float meanStd();                      // the mean standard deviation of the individual probe pairs across the expt. series.
  float meanOverStd();                  // returns   (sum of means) / (sum of std deviations) across the experimental series. 
  float diff(int up, int down);
  float zCompare(vector<float>&, vector<int>);   // stupid name,, gives the euclidean distance between the points in the float vector for the experimental points given.. -- does a local z-score transformation of the values in the probe set to allow for small sets.. 
  float zMeanCompare(vector<float>&, vector<int>); // as above, but only compares agains the mean of the probe set.. 
  // does nothing if size of vector is 1.. 
  /// SOCKET write operation. Write data to socket.. if connected. Use the address of a QSocket as this seems to be quite 
  /// easy...
  bool writeDataToQSocket(QSocket* socket);  // returns true if it works, false otherwise.. 

};

struct uniGeneData {
  int index; // the mouse unigene index,, lets assume for now that we are dealing with mouse genes (allows me to use an int)
  string title;
  string gene;
  uniGeneData(int i, string t, string g);    // constructor.. 
};

struct celeraMatch {
  string celeraGene;      // the CG field in the celera entries..
  float expectation;      // blast expectation score..
  float match;            // the identity percentage.. sort of..
  string sf;
  string fn;
  string gn;
  string gs;
  string nd;             // I don't know the exact definition of these, but they are taken from celera entries..
  celeraMatch(string cg, float exp, float m, string SF, string FN, string GN, string GS, string ND);
};

struct probe_data {
  // I need to redesign this in view of the fact that one probe set can logically 
  // be linked to more than one unigene.. Sequence homology and other difficulties 
  // in assigning probe sets to single probe identities mean that we can't reasonably 
  // exclude some but include others.. 
  // simplest way is probably by creating a 'uid_data' struct that contains some stuff.
  int index;              // probe set index from p_sets
  string gbid;            // gb id from thingy
  string afid;            // af_id from thingy
  vector<uniGeneData> ugData;     // for the unigene fields.. 
  vector<celeraMatch> celeraMatches;  // celera genes with matches.. 
  //int ugid;               // unigene_id 
  //string uggene;          // unigene gene field
  //string ugtitle;         // unigene title field
  string afdes;           // affy description.. 
  string tigrDescription; // tigr annotation from table tigr_annotation.. 
  vector< vector<string> > go;    // go classification from dots..
  bool defined;            // set to false in empty constructor

  probe_data();       // empty constructor  -- just construct the thing manually.. -and hope it works..
  ~probe_data();      // empty destructor. probably not necessary..
};


struct probeSetSet {
  // a struct containing a set of probe sets. one raw, and one normalised..
  // well, sort of anyway... plus the relevant header information about the
  // probe sets. 
  vector<probe_set> raw;
  vector<probe_set> normalised;
  string headerInfo;      // expand this later to be more specific.. 
  vector<nadFileInfo> fileInformation;  // hope this is fine.. 
  vector<float> correctionFactors; // keep correction factors, initially set to 1, but can be dynamically changed in functions
  float maxCorrectionFactor;       // good for drawing use..
  int minIndex;
  int maxIndex;
  int chipIndex;      // but we should change this at some point.. 
  int fileNo;         // number of CEL files or number of experiments..
  vector<float> experimentIndices;    // one for each ...
  probeSetSet(string& filename);    // construct from file..
  probeSetSet(vector<string>& filenames); // construct from a list of files. Potentially dangerous..
  probeSetSet();                   // empty constructor..
  void setCorrectionValues(vector<probe_set*>& probes);  // sets the correction values based on the mean profile of a set of genes thought to be housekeepoing ones
  ~probeSetSet();                  // a destructor..
};

// Network creator...
probe_set readProbeDataFromQSocket(QSocket* socket, bool& ok);


vector<probe_set> probes_from_file(string& infile, string& type);
vector<float> get_values(string& line);
void probes_to_binary_file(vector<probe_set>& p_set, string& outfile);
vector<probe_set> probes_from_binary_file(string& name);
vector<int> eComparer(vector<probe_set>& probes, probe_set& probe, vector<int>& comp_points);
vector<int> anovaScorer(vector<probe_set*>& probes, vector<int>& comp_points);
vector<int> anovaScorer(vector<probe_set*>& probes, vector<int>& comp_points, vector<float>& rvalues); // overloaded.. 
// assigns the anova scores to the rvalues vector.. 
vector< vector<int> > meanPosNumberDeterminer(vector<probe_set*>& probes);
vector< vector<int> > probPosNumberDeterminer(vector<probe_set*>& probes);


vector<int> selSumScorer(vector<probe_set*>& probes, vector<int>& comp_points);
vector<int> meanOverStdSorter(vector<probe_set*>& probes);
vector<int> concordanceScorer(vector<probe_set>& probes);
vector<int> sumScorer(vector<probe_set*>& probes);
vector<int> diffScorer(vector<probe_set*>& probes, int up, int down); 
vector<int> zComparer(vector<probe_set*>& probes, vector<float>& values, vector<int>& tindex);
vector<int> zMeanComparer(vector<probe_set*>& probes, vector<float>& values, vector<int>& tindex);
vector<int> meanMaxMeanDeviationScorer(vector<probe_set*>& probes);
cluster_set kCluster(vector<probe_set*>& probes, int n);
cluster_set kCluster(vector<probe_set*>& probes, vector<int> expts, int n);  // overloaded.. clusters only on the experimental points in expts..
vector<float> meanMean(vector<probe_set*>& probes);      // just returns the mean of the values of a set of probes.. 

vector<probe_data> data_from_db();   // what the hell, eh.. 

#endif
