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

#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H


#include <string>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <qstring.h>

using namespace std;

struct statusMessage {
  unsigned int id;  // some identifier for a requester..
  int no;           // some internal identifier for a requester
  bool ok;          // if it worked ok..
  vector<string> errorMessages;  // error Messages..
  vector<int> errorCodes;        // not sure which one is better, but anyway, these are all going to be ok
  statusMessage(bool o=true){
    id = 0;
    no = 0;
    ok = o;
  }
  statusMessage(unsigned int i, bool o=true){
    id = i;
    no = 0;
    ok = o;
  }
  statusMessage(unsigned int i, int n, bool o=true){
    id = i;
    no = 0;
    ok = o;
  }
  statusMessage(unsigned int i, int n, vector<string> em, bool o=true){
    id = i;
    no = 0;
    ok = o;
    errorMessages = em;
  }
  statusMessage(unsigned int i, int n, vector<string> em, vector<int> ec, bool o=true){
    id = i;
    no = 0;
    ok = o;
    errorMessages = em;
    errorCodes = ec;
  }
  statusMessage(string message){
    id = no = 0;
    ok = true;
    errorMessages.push_back(message);
  }
};

struct Comment {
  int index;
  QString comment;
  Comment(int i, QString c){
    index = i;
    comment = c;
  }
  Comment(){
    index = -1;
    comment = "";  // not really necessary, but..
  }
};

struct RegionSpecification {
  string chromosome;
  uint begin;
  uint end;
  
  RegionSpecification(){
    chromosome = "NULL";
    begin = end = 0;
  }
  RegionSpecification(string c, uint b, uint e){
    chromosome = c;
    if(b < e){
      begin = b;
      end = e;
    }else{
      begin = e;
      end = b;
    }
  }
  bool merge(RegionSpecification& spec);     // if overlap then merge spec into this.. (don't change spec, just find overlap.. ).
};
      

struct Exon {
  // a general exon, with no coding information. (other exons should really inherit from this).
  ~Exon(){
    cout << "\t\t\t\tEXON BEING DESTROYED BY HEINOUS ACTION" << endl;    // though we shouldn't really see this.. 
  }
  Exon(){
    start = stop = strand = 0;
  }
  Exon(int st, int sp, int sd, int tb, int te){
    if(st < sp){
      start = st;
      stop = sp;
    }else{
      start = sp;
      stop = st;
    }
    if(sd == -1 || sd == 1){ 
      strand = sd;
    }else{
      strand = 0;
    }
    tbegin = tb;
    tend = te;    
  }
  int start;
  int stop;
  int strand;
  int tbegin;
  int tend;    // positions in the transcript.. 
};

struct Transcript {
  ~Transcript(){
    // should destroy the exons, implement later.. 
  }
  Transcript(){
    id = "null";
    chromosome = "null";
    source = "null";
    start = stop = strand = 0;
    exons = 0;
    exonMemSize = 0;
    exonNo = 0;
    start = stop = -1;
    strand = 0;
  }
  Transcript(string i, string s, string c, int l, int stnd){
    id = i;
    source = s;
    chromosome = c;
    length = l;
    start = stop = 0;
    exonMemSize = 5;
    exons = new Exon*[exonMemSize];
    exonNo = 0;
    start = stop = -1;
    strand = stnd;
  }
  void addExon(Exon* ex);    // put in cpp file as a little longer needs more checking, and would be kind of ugly here.. 
  void addExon(int st, int sp, int sd, int tb, int te);
  string id;
  string source;   // the source of the identification and the library
  string chromosome;
  int start;
  int stop;
  int strand;
  int length;    // the length of the actual transcript so we can work out how complete the alignment is.
  Exon** exons;

  unsigned int exonNo;
  unsigned int exonMemSize;
};


struct ensemblExon {
  ensemblExon();
  ensemblExon(string i, string c, int st, int sp, int cst, int csp, int sd);
  string id;     // the ensembl identifier, so I can work out unique things..
  string chromosome;
  int start;
  int stop;
  int codeStart;
  int codeStop;   // if non coding. Doesn't seem like there is an awful lot of things going on for this..
  int strand;     // -1 or 1 ,, -to use ensembl syntax.
};

struct ensemblTranscript {
  //ensemblTranscript();
  ensemblTranscript(int tIndex, string i, string c, int st, int sp, int sd);
  ~ensemblTranscript(){};
  void addExon(ensemblExon* ex);     // so I can grow it with some exons.. 
  int index;        // the db Index, should be unique perhaps.. (called transcript)
  string id;        // the ensembl ID.
  string chromosome;
  int start;
  int stop;
  int strand;
  ensemblExon** exons;    // array of pointers to exons, easier for memory management.. as I don't have to worry about copying actual exons when growing the array.
  int exonNo;            // the number of exons..
  private :
  int exonMemSize;       // the amount of memory allocated for the exons.. 
  // these last guys should really be private so that I don't screw them up in later code. Arghh. 
};

struct ensemblGene {
  ensemblGene();
  ensemblGene(int dbi, string ensid, string chr, int sd);
  //  ensemblGene(int dbi, string ensid, string extid, string des, string chr, int sd);
  int dbIndex;
  string ensemblId;
  //string externalId;
  //string description;

  // later - include a pointer for an annotation struct, and some means of inserting it, but for now.. 
  // leave it blank.. 
  string chromosome;
  ensemblTranscript** transcripts;
  void addTranscript(ensemblTranscript* t);
  void findLimits();            // find my start and stops.. 
  int transcriptNo;            // number of transcripts... -- should really be private.. hmm. 
  int strand;
  int start;
  int stop;    // check when adding transcripts. 
  private :
  int transcriptMemSize;       // the amount of memory allocated for transcripts... 
};

struct probeSetMatch {
  probeSetMatch();
  //  probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int afn, int m, int cs, int ce, double exp, string chr);
  //  probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int m, int cs, int ce, double exp, int strand, string chr, string af_match, string gen_match);
  probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int m, int cs, int ce, double exp, int strand, string chr);
  ~probeSetMatch(){};
  int dbIndex;    // the affy database Index
  int aIndex;     // the array index, at the moment it is always one less than dbIndex, so a little redundant, but..
  int afStart;
  int afEnd;
  int afLength;
  int alignLength;
  //int af_n_count;
  int match;
  int cStart;
  int cEnd;
  double expectation;
  int strand;
  string chromosome;    // should probably make this an int and have an index, but it's not that important.
  //string af_match_sequence;
  //string ensembl_match_sequence;   // wasteful, but at the most around 30 MB of sequence, which is not that bad
};

struct probeSetMatchSet {
  probeSetMatchSet(){
    minExpect = 1;
    minPos = maxPos = 0;
    chromosome = "";
  };    // empty constructor,, it's easier that way.. 
  probeSetMatchSet(probeSetMatch* psm){
    matches.insert(psm);
    dbIndex = psm->dbIndex;
    minPos = psm->cStart;
    maxPos = psm->cEnd;
    minExpect = psm->expectation;
    expectProduct = psm->expectation;
    chromosome = psm->chromosome;  // !! I hope.. 
    afLength = psm->afLength;
    matchSum = psm->match;                    // increment at later stages...
    mismatchSum = psm->alignLength-matchSum;  // and keep incrementing these later on as well.
    strand = psm->strand;
  }
  void insertMatch(probeSetMatch* psm);
  set<probeSetMatch*> matches;
  string chromosome;
  int dbIndex;      // the probe set id. -- we should not mix matches from different probe sets after all. 
  int minPos;
  int maxPos;   // the min and max positions..
  double minExpect;  // the minimumexpectation..
  double expectProduct;  // the product of the individual expectations.. sometimes I suppose this should be useful.. 
  int matchSum;
  int mismatchSum;
  int afLength;       // redundant I know, as it's in all of the probe_set_matches.. 
  int strand;
  friend bool operator<(const probeSetMatchSet& a, const probeSetMatchSet& b){
    //    return(a.minExpect < b.minExpect);
    return(a.matchSum > b.matchSum);
  }
};

struct ishProbeMatch {
  ishProbeMatch(){};
  ishProbeMatch(int pi, int pst, int psp, int cst, int csp, float pcnt, string c, int asid){
    probeIndex = pi;
    pStart = pst;
    pEnd = psp;
    //pLength = pl;
    //alignLength = all;
    //match = m;
    cStart = cst;
    cEnd = csp;
    //expectation = e;
    //strand = s;
    chromosome = c;
    //p_match_sequence = pms;
    //ensembl_match_sequence = ems;
    percent = pcnt;
    assemblyId = asid;
  }
  int probeIndex;
  int pStart;
  int pEnd;
  // int pLength;   known by the container.. 
  //  int alignLength;
  //  int match;
  int cStart;
  int cEnd;
  float percent;
  // double expectation;
  //  int strand;    // known by the container.. 
  int assemblyId;    // which assembly do I belong to.. 
  string chromosome;
  //string p_match_sequence;
  //string ensembl_match_sequence;
};

struct ishProbeMatchSet {
  ishProbeMatchSet(){
    //minExpect = 1;
    minPos = maxPos = 0;
    chromosome = "";
  };    // empty constructor,, it's easier that way.. 
  //  ishProbeMatchSet(ishProbeMatch* ipm, int strnd, int pl, int beg, int stp, int asid){
  ishProbeMatchSet(int dbi, string c, int strnd, int pl, int beg, int stp, int asid, float scr){
    //matches.insert(ipm);
    dbIndex = dbi;
    minPos = beg;
    maxPos = stp;
    //minExpect = ipm->expectation;
    //expectProduct = ipm->expectation;
    chromosome = c;  // !! I hope.. 
    pLength = pl;
    matchSum = 0;                    // increment at later stages...
    score = scr;
    //mismatchSum = ipm->alignLength-matchSum;  // and keep incrementing these later on as well.
    strand = strnd;
    assemblyId = asid;
  }
  void insertMatch(ishProbeMatch* psm);
  set<ishProbeMatch*> matches;
  string chromosome;
  int dbIndex;      // the probe set id. -- we should not mix matches from different probe sets after all. 
  int minPos;
  int maxPos;   // the min and max positions.. --- i.e. start and stop, we can get direclty from the thingy.. now.. 
  //  double minExpect;  // the minimumexpectation..
  //double expectProduct;  // the product of the individual expectations.. sometimes I suppose this should be useful.. 
  int matchSum;            // use the score.. from the new style table..
  float score;             // sum (lengths * percent) for all matches.. 
  //int mismatchSum;
  int pLength;       // redundant I know, as it's in all of the probe_set_matches.. 
  int strand;
  int assemblyId;
  friend bool operator<(const ishProbeMatchSet& a, const ishProbeMatchSet& b){
    return(a.score > b.score);
  }
};

struct genomicRegion {
  genomicRegion();
  genomicRegion(int s, int e, string chr);
  ~genomicRegion(){};             // don't delete the pMatches, as these could overlap with another region. 
  void addProbeSetMatch(probeSetMatch* psm);
  string chromosome;
  int start;
  int end;
  probeSetMatch** pMatches;       // array of pointers..
  int pMatchNo;
  int pMatchSize;                 // amount of memory allocated.. 
  //////////////////////  and the something for ensembl Gene predictions..
  ensemblGene** ensGenes;
  int ensGeneNo;
  int ensGeneSize;         // again, something for the amount of memory allocatoed.. 
  void addEnsGene(ensemblGene* ensg);
  ///////////////////      and something for the in situ probe matches..
  /////////////////        I should really virtualise this.. as they all take pointers 
  /////////////////        it would indeed be very simple to either use inheritance or void pointers of objects..
  ////////////////         and then in all of the functions just check the sizes.. -- only problem is knowing the 
  ///////////////          number of each type of object in the data.. but for now just keep adding it on until I 
  ////////////////         rationalise the whole damn thing.. .. !!
  ishProbeMatchSet** ishMatches;
  int ishMatchNo;
  int ishMatchSize;
  void addIshMatch(ishProbeMatchSet* ishpm);
  ////
  /// a universal transcript structure, which we may subclass in the future..
  Transcript** transcripts;
  unsigned int transcriptNo;
  unsigned int transcriptMemSize;
  void addTranscript(Transcript* transc);
};

struct chromAnnotation {         // just a series of genomicRegion for each chromosome.. 
  chromAnnotation();
  chromAnnotation(string chr, int size, int rsize);
  ~chromAnnotation();
  void addRegion(genomicRegion* gr); 
  string chromosome;
  genomicRegion** regions;    // array of pointers whaooooa.. 
  genomicRegion** regionsCovered(int st, int sp, int& n);   // returns a pointer and the number of regions covered by a range.. 
  int regionNo;
  int regionMemSize;          // amount of memory allocated.. 
  int regionSize;             // number of base pairs allocated.. per region.. 
  int chromosomeSize;
};

// and then we will add more stuff to it.. 

struct exInfo{
  exInfo();      
  int dbaseIndex; // just refers to the index in the database table. But we may make more use of this later.    
  int exptGroup; // what group does it refer to??
  int realDbaseIndex;    // a kludge,, using the dbase index as an index for an array, is obviously hopelessly inflexible.. -> but 
  ///                       call it dbaseIndex for now to allow me not to change to many things. Later change the name to something more reasonable. 
  float index;
  string shortName;
  string description;
  map<int, bool> chips;   // whether used or not
  map<string, string> parameters;
};  

// need some sort of structure to contain descriptions of the chips.. rather 
// than just having a string for it.. 
struct chipInfo{
  int index;
  string id;
  string description;
  set<int> equivs;    // chip equivalents, so that chip 6 can represent 4 and 5.. sort of..
  chipInfo(){
    index = -1;
    description = "null";
  }
  chipInfo(int i, string Id, string d){
    index = i;
    id = Id;
    description = d;
  }
  chipInfo(int i, string Id, string d, set<int> e){
    index = i;
    id = Id;
    description = d;
    equivs = e;
  }
  // and a function to check..  which is a very rudimentary check. 
  bool defined(int i){
    return((bool)equivs.count(i));
  }
};

struct sessionInformation {
  // holds information about one session. 
  int index;       // the index given in the expression database. (the sessions index!!) 
  int owner;       // the session owner. .. good to know
  string title;
  string description;
  vector<string> keywords;    // vector, cheaper than set.
  set<int> members;   // other members.. of the set.
  
  sessionInformation();
  sessionInformation(int i, int o, string t, string d, vector<string> kw, set<int> om);
  friend bool operator<(const sessionInformation& a, const sessionInformation& b){
    return(a.index < b.index);
  }
};

struct annotationInformation {
  // similar to the above, but holds information about an annotatoin. It's a bit simpler.
  int index;
  int owner;
  string annotation;
  set<int> members;      // not really members, but the best word I could think of at the moment.
  
  annotationInformation(int i, int o, string a, set<int> om);
  annotationInformation();
};

struct userInformation{
  int index;   // the user's index number
  string userName;
  string fullName;
  string labName;     // if these are defined.. 
  userInformation(int i, string u, string fn, string ln);
  userInformation();
};


#endif
