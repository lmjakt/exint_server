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

#include "dataStructs.h"
#include <vector>
#include <map>
#include <string>

using namespace std;



struct probe_set {
  float** probes;
  uint probeNo;       // the number of probes.
  uint exptSize;      // the numbers of experiments for this probe set (and the number of values for each probe) 
  uint exptCapacity;  // capacity of probes and exptIndex.

  uint* exptIndex;                 // the experimental indices. should have same number ..
  int* exptLookup;                // all experiments, with the value being the index of the -1 not defined.. 
  uint allExptNo;                // the size of all of the experiments. Tricky this !!... 

  float* mean; 

  int index;                        // probe set index -- related to the affymetrix thingy
  probe_set();
  probe_set(int i, uint allESize);
  probe_set(int i, vector< vector<float> >& p, vector<uint>& eindices, uint allESize);
  ~probe_set();
};

// allows the extension of the 
bool addDataTo_probe_set(probe_set* ps, vector<float>& p, uint expt);

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

struct blastAlignment {
  // hold the data for a single blast alignment
  blastAlignment();       // a constructor for nothing.. 
  blastAlignment(int as, int ae, int ms, int me, double e, int al, int m, int an);
  int affyStart;
  int affyEnd; 
  int mStart;
  int mEnd;
  double expectation;  // hopefully a really small number..
  int alignmentLength;  // account for gaps..
  int matches;          // number that match
  int affyN;                // number of N's in the affy sequence. 
  //  vector<int> affyGap; //the gap positions in the affy sequence..
  //  vector<int> mGap;    // gaps in the match sequence.. // -hmm, how to count these I wonder..
  //  vector<int> misMatch; // the positions of the mismatches. 
  //// whether to include the positions of n's I dont' know.. 
};

struct blastMatch {
  // essentially contains a series of blastAlignments and information about the
  // gene that is matchesd... each probe data may have a series of these..
  blastMatch();
  blastMatch(string mI, int l, int al);
  string matchId;
  int length;
  int af_length;      // redundant, but its easier to put here.. 
  // features ?? I'm not sure.. 
  multimap<string, string> description;     // as I don't know what I'll be putting in
  vector<blastAlignment> alignments;   // the alignments made between this and the other..
  // and I think that is enough for now.. // hmm. sounds trouble.. 
  // Putting the annotation here is inherently BAD,, -and I should resist it, but I need working version by 
  // the end of tomorrow where people can see the linkage, so I will put this in for now, but I am hoping that
  // I will get the chance of cleaning up the protocols and the data structures one of these days.. 
};

struct blastGenomeMatch {
  // information from a blast match. Gives the genomic context. -chromosome, positions along the chromosome
  // and the neighbouring gene. And a best guess. All taken from a tabe in thingy.. 
  blastGenomeMatch();
  blastGenomeMatch(string c, int l, int m, int h, int g, int af, bool o);
  string chromosome;          // string as we also have X and Un ??? 
  int lowerIndex;
  int midIndex;
  int higherIndex;
  int af_length;              // soo redundant. We have to fix this stuff up one of these days... 
  int guessGene;              // the index of the guessed one.. -- a little bit wrong, because it could theoretically be different for
  bool overlap;
  // the different alignments, though it is unlikely to be so.. 
  vector<blastAlignment> alignments;
  /// change this so that none of the gene descriptions are held here, but instead have a map<int, multimap<int, string> > for 
  /// the lookups, along with a map<int, string> for the fieldnames Keep this in the main probesetset data structure, and look up when we send the data. 
  /// When sending data look up what's available in these structures.. and just send it on.. 
};

struct probe_data {
  //////////// NOTE THAT THIS STRUCT SHOULD ONLY CONTAIN STATIC DATA SUCH THAT IT WILL NOT BE CHANGED BY THE 
  //////////// --- TO MAKE SURE WE DON'T NEED TO USE ANY MUTEXES OR OTHER UGLY THINGS... 
  int index;              // probe set index from p_sets
  int chip;
  string gbid;            // gb id from thingy
  string afid;            // af_id from thingy
  vector<uniGeneData> ugData;     // for the unigene fields.. 
  vector<probeSetMatchSet> probeSetMatches;    // pointers to genomic locations. maybe it should be a map of sorts, but I'm not sure yet. 
  int blastGuess;                      // the ensembl gene that seems to match the best given some silly checks.. 
  string sequence;                    // lets get the sequence as well. Era-san' been asking for it.. !!

  string afdes;           // affy description.. 
  string tigrDescription; // tigr annotation from table tigr_annotation.. 
  vector< vector<string> > go;    // go classification from dots..
  bool defined;            // set to false in empty constructor

  probe_data();       // empty constructor  -- just construct the thing manually.. -and hope it works.. // this is ugly !! 
  ~probe_data();      // empty destructor. probably not necessary..
  void addProbeSetMatch(probeSetMatch* psm);   // adds a probeSet match in the right place.. 
};

struct ish_annotation {
  int userId;
  int annotation_id; // important..  -- for performing updates.. -- use the actual database oid for this.. as we can easily get hold of this.
  string userName;   // it's redundant I know, but will probably make life more convenient.
  string fieldName;  // arghhe... 
  string annotation;
  float value;        // -- define either the float or the annotation string..
  bool numerical;     // -- did we define the value or the string.. 
  ish_annotation(int annId, int uid, string usr, string note, string fName){
    annotation_id = annId;
    userId = uid;
    userName = usr;
    annotation = note;
    numerical = false;
    fieldName = fName;
  }
  ish_annotation(int annId, int uid, string usr, float v, string fName){
    annotation_id = annId;
    userId = uid;
    userName = usr;
    value = v;
    numerical = true;
    fieldName = fName;
  }
  ish_annotation(){
    annotation_id = 0;
    userId = 0;
    userName = "null";
    annotation = "null";
    value = 0;
    numerical = false;
    fieldName = "null";
  }
};

struct ishProbeData {
  vector<ishProbeMatchSet*> probeMatches;   // the genomic matches..
  int probeId;
  string consensusSequence;
  string antisensePromoter;
  int afIdAssociation;     // the associated afId..
  string vectorName;
  int designLength;
  string probeName;        // maybe a blank string unless it has been defined. Can be changed by privileged users.. 
  string probeIdentifier;   // a unique string for the probe.. 
  int ensemblGuess;        // a guess for an ensembl gene.. --- I don't actually have this yet.. 
  /// some maps for annotation stuff..
  map<int, ish_annotation> textAnnotation;       // int is the oid, so that I can update these bastards at some point .. arghh.. oid should be unique.. I hope
  map<int, ish_annotation> numberAnnotation;    // anything that would measured by a number.. 
  map<int, ish_annotation> classification;      // -- does it fall into class A? -- and with what degree of confidence.
  //  multimap<string, ish_annotation> textAnnotation;
  //multimap<string, ish_annotation> numberAnnotation;    // anything that would measured by a number.. 
  //multimap<string, ish_annotation> classification;      // -- does it fall into class A? -- and with what degree of confidence.
                                          // I am not sure yet as to what access/modification policy we'll use for these
                                          // clearly as these are intended to be provided by users,, -there has to be some sort 
                                          // of policy in terms of who gets to change or modify the number.. 
                                          // we could set it up such that users get to keep their own set of numbers.. -with the default being 
                                          // the most common precedent, or some such criteria.. but this is likely to be rather difficult, and can
                                          // think about it later.. -- at the moment, just let the superuser change this.. (i.e. me). 
  /// it doesn't really make sense to have a multimap of the classifications, -- in the database the classifications .. and who made them and 
  /// and when will be detailed.. -the multimap at least let's us look at what the numbers are like, even if it doesn't let us check out users..
  /// --- which we need to do, but that will require a more complex data structure, -- which I can implement later on. However,, .. it shouldn't 
  //  be too bad to implement after having done this.. 
  ishProbeData();
  ishProbeData(int pId, string cons, string asp, int afid, string v, int dll, string pName, string p_ident, int eGuess){
    probeId = pId;
    consensusSequence = cons;
    antisensePromoter = asp;
    afIdAssociation = afid;
    vectorName = v;
    designLength = dll;
    probeName = pName;
    probeIdentifier = p_ident;
    ensemblGuess = eGuess; 
  }
  void addIshProbeMatch(ishProbeMatch* ipm);   // adds a ishProbeMatch to the ishProbeMatchSet vector.. if appropriate..
};

vector<probe_data> data_from_db(const char* conninfo);

#endif
