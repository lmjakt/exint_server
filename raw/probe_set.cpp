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

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <string>
#include <qstring.h>
//#include <libpq++.h>
#include <libpq-fe.h>
#include <algorithm>
#include <stdlib.h>
#include "probe_set.h"

using namespace std;

float* meanValues(vector< vector<float> >& p){
  // return the mean values..
  // assume length is the same for all. change this later if we need to.
  float* mean = new float[p[0].size()];
  for(uint i=0; i < p.size(); i++){
    for(uint j=0; j < p[i].size(); j++){
      mean[j] += p[i][j];
    }
  }
  for(uint i=0; i < p[0].size(); i++){
    mean[i] = mean[i]/p.size();
  }
  return(mean);
}


probe_set::probe_set(){
  index = 0;
  probeNo = 0;
  exptSize = 0;
  allExptNo = 0;
  exptCapacity = 0;
  exptIndex = 0;
  exptLookup = 0;
  probes = 0;
  exptIndex = 0;
}

probe_set::probe_set(int i, uint allESize){
  probeNo = 0;
  exptSize = 0;
  exptCapacity = 0;
  exptIndex = 0;

  index = i;
  allExptNo = allESize;
  exptLookup = new int[allExptNo];
  for(uint i=0; i < allExptNo; ++i)
    exptLookup[i] = -1;
  probes = 0;
  exptIndex = 0;

}

bool addDataTo_probe_set(probe_set* ps, vector<float>& p, uint expt)
{
  int expt_increment = 40;  // this should be settable.. 
  if(!ps->probeNo)
    ps->probeNo = p.size();
  if(ps->probeNo != p.size()){
    cerr << "addDataTo_probe_set vector p size " << p.size() << " != " << ps->probeNo << endl;
    exit(1);
  }
  // grow the data if necessary.. 
  uint exptCapacity = ps->exptCapacity;
  if(ps->exptSize == exptCapacity){
    exptCapacity += expt_increment;  // swap all values lower down.
    
    uint* exptIndex = new uint[exptCapacity];
    for(uint i=0; i < ps->exptSize; ++i)
      exptIndex[i] = ps->exptIndex[i];
    
    float** probes = new float*[ps->probeNo];
    for(uint i=0; i < ps->probeNo; ++i){
      probes[i] = new float[exptCapacity];
      for(uint j=0; j < ps->exptSize; ++j)
	probes[i][j] = ps->probes[i][j];
    }
    // This section should be atomic if we want to be able to grow
    // the data during operation. (Where to place mutexes is unclear,
    // but it should be possible.. )
    float** temp_p = ps->probes;
    uint* temp_ei = ps->exptIndex;

    ps->exptCapacity = exptCapacity;
    ps->exptIndex = exptIndex;
    ps->probes = probes;
    if(temp_p){
      for(uint i=0; i < ps->probeNo; ++i)
	delete []temp_p[i];
      delete []temp_p;
    }
    delete []temp_ei;
  }

  ps->exptIndex[ps->exptSize] = expt;
  ps->exptLookup[ expt ] = ps->exptSize;
  for(uint i=0; i < p.size(); ++i)
    ps->probes[i][ps->exptSize] = p[i];
  ps->exptSize++;
  return(true);
}

probe_set::probe_set(int i, vector< vector<float> >& p, vector<uint>& eindices, uint allESize){
  //  cout << "Creating probe set for index " << i << "  with address " << this << endl;
  index = i;

  probes = new float*[p.size()];
  probeNo = p.size();
  exptSize = p[0].size();
  allExptNo = allESize;
  exptCapacity = exptSize;
  for(uint i=0; i < p.size(); i++){
    probes[i] = new float[p[i].size()];
    for(uint j=0; j < p[i].size(); j++){
      probes[i][j] = p[i][j];
    }
  }
  if(eindices.size() != exptSize){
    cerr << "CRASH CRASH, exptSize is not the same as eindices.size " << endl;
    exit(1);
  }
  // set up the exptLookup, with all -1 values.. 
  exptLookup = new int[allESize];
  for(uint i=0; i < allESize; i++){ 
    exptLookup[i] = -1; 
  }     // default to -1;

  exptIndex = new uint[eindices.size()];
  for(uint i=0; i < eindices.size(); i++){
    exptIndex[i] = eindices[i];
    exptLookup[eindices[i]] = i;        
    // no error checking; initially, for speed,
    // but since clearly this is not limiting we might want to rethink this.
  }
  mean = meanValues(p);
}

probe_set::~probe_set(){
  delete []probes;
  delete []exptIndex;
  delete []exptLookup;
}

probe_data::probe_data(){   // is this really necessary ?
  defined = false;
  blastGuess = 0;
} 

void ishProbeData::addIshProbeMatch(ishProbeMatch* ipm){
  int maxRange = 100000;    // the size of an indiviidual thingy..
  for(uint i=0; i < probeMatches.size(); i++){
    ishProbeMatchSet* it = probeMatches[i];   // historical,, for shorthand.. 
    if(it->assemblyId == ipm->assemblyId){
      (*it).insertMatch(ipm);
      sort(probeMatches.begin(), probeMatches.end());    // which puts the one with the highest score first.. though it doesn't actually matter that much.. 
      return;
    }
  }
  //  probeMatches.push_back(ishProbeMatchSet(ipm));
  cerr << "Unassembled ishProbeMatch, bugger, maybe I should die here instead of bravely carrying on.. " << endl;
}

void probe_data::addProbeSetMatch(probeSetMatch* psm){
  // first go through all of the probeSetMatchSets and see if there is one that matches..
  // I'm going to make an arbitrary limit here of 100,000 bp.. i.e. anything within this gets included in the same probeSetMatchSet...
  int maxRange = 100000;
  vector<probeSetMatchSet>::iterator it;
  for(it = probeSetMatches.begin(); it != probeSetMatches.end(); it++){
    if((*it).chromosome == psm->chromosome && abs((*it).minPos - psm->cStart) < maxRange && abs((*it).maxPos - psm->cEnd) < maxRange && psm->dbIndex == (*it).dbIndex && (*it).strand == psm->strand){
      // OK, so that could give me a little bit larger range than maxRange, but what can you do eh?? 
      (*it).insertMatch(psm);
      sort(probeSetMatches.begin(), probeSetMatches.end());
      return;  // I'm done.. 
    }
  }
  // if I got here, I didn't find anywhere to stick this probe set match, so I better make a new probeSetMatchSet..
  probeSetMatches.push_back(probeSetMatchSet(psm));   // and that should be all I need to do in fact !!.
  sort(probeSetMatches.begin(), probeSetMatches.end());
}
  
  

probe_data::~probe_data(){   // and what about this??? 
}


vector<probe_data> data_from_db(const char* conninfo) {
  // open up a connection to a database backend,, -lets just use a normal ascii cursor.. 
  //  PgCursor data(conninfo, "portal");    // still don't know why the portal
  PGconn* conn = PQconnectdb(conninfo);

  // check the backend..
  if(!conn || (PQstatus(conn) != CONNECTION_OK)){
    cerr << "Unable to connect to database: " << conninfo << endl;
    exit(1);
  }
  PGresult* res = PQexec(conn, "select a.*, b.chip from probe_data a, p_sets b where a.index=b.index order by index");
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    cerr << "PQexec was not able to execute query\n"
	 << PQresultErrorMessage(res) << endl;
    PQclear(res);
    PQfinish(conn);
    exit(1);
  }
  
  int tuples = PQntuples(res); // data.Tuples();
  // create a vector of probe_data ..
  vector<probe_data> p_data(tuples);  // with the appropriate amount of space..
  // resize the vector afterwards. This is very ugly, but the best I can think of
  // at the time. I really need to sort this out, but it will have to do for now..
  // and go throuh and get all the thingies..
  int index;
  int p_index;
  int maxIndex = 0;
  for(int i=0; i < tuples; i++){
    index = atoi(PQgetvalue(res, i, 0)); // data.GetValue(i,0));
    // for this one we shall be defining the p_data[index-i]
    p_index = index-1;
    int ugid = atoi(PQgetvalue(res, i, 3)); // data.GetValue(i, 3));
    string uggene = PQgetvalue(res, i, 4); //data.GetValue(i, 4);
    string ugtitle = PQgetvalue(res, i, 5); // data.GetValue(i, 5);
    //uniGeneData tempdata(ugid, ugtitle, uggene);
    p_data[p_index].ugData.push_back(uniGeneData(ugid, ugtitle, uggene));
    if(!p_data[p_index].defined){
      p_data[p_index].index = index;
      p_data[p_index].defined = true;
      p_data[p_index].afid = PQgetvalue(res, i, 1); // data.GetValue(i, 1);
      p_data[p_index].gbid = PQgetvalue(res, i, 2); // data.GetValue(i, 2);
      p_data[p_index].afdes = PQgetvalue(res, i, 6); // data.GetValue(i, 6);
      p_data[p_index].tigrDescription = PQgetvalue(res, i, 7); // data.GetValue(i, 7);
      p_data[p_index].chip = atoi(PQgetvalue(res, i, 8)); // data.GetValue(i, 8));
      p_data[p_index].go.resize(0);         // don't know if that's good or not.. 
    }
    if(index > maxIndex){
      maxIndex = index;
    }
  }
  //close the cursor,, and lets leave it at that for the time being, it won't work anyway..
  p_data.resize(maxIndex);         // hopefully this will now work.. 
  PQclear(res);

  // Some GO data
  res = PQexec(conn, "select a.index, b.generation, c.description from p_sets a, af_go_gen b, go c where a.af_id=b.af_id and b.go=c.index and a.chip=1 order by a.index, b.generation desc");
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    cerr << "GO data query had some problem\n"
	 << PQresultErrorMessage(res) << endl;
    exit(1);
  }
  
  int go_gen;
  tuples = PQntuples(res);
  for(int i=0; i < tuples; i++){
    index = atoi(PQgetvalue(res, i, 0)); // data.GetValue(i, 0));
    p_index = index-1;
    go_gen = atoi(PQgetvalue(res, i, 1)); // data.GetValue(i, 1));
    if((int)p_data[p_index].go.size() < (go_gen)){
      p_data[p_index].go.resize(15);             // I DON'T UNDERSTAND THIS, SOMEONE PLS. EXPLAIN..
    }
    p_data[p_index].go[go_gen-1].push_back(PQgetvalue(res, i, 2)); // data.GetValue(i,2));
  }
  PQclear(res); 
  PQfinish(conn);

  return(p_data);
}

uniGeneData::uniGeneData(int i, string t, string g){
  index = i;
  title = t;
  gene = g;
}

celeraMatch::celeraMatch(string cg, float exp, float m, string SF, string FN, string GN, string GS, string ND){
  celeraGene = cg;
  expectation = exp;
  match = m;
  sf = SF;
  fn = FN;
  gn = GN;
  gs = GS;
  nd = ND;
  //cout << "Assigned values to celeraMatch and created a thingy, by, by.. " << endl;
}

blastAlignment::blastAlignment(){
  affyStart = -1;
  affyEnd = -1;
  mStart = -1;
  mEnd = -1;
  expectation = -1;
  alignmentLength = -1;
  matches = -1;
  affyN = -1;
  // assign all to impossible values so I can error check at some point
}

blastAlignment::blastAlignment(int as, int ae, int ms, int me, double e, int al, int m, int an){
  affyStart = as;
  affyEnd = ae;
  mStart = ms;
  mEnd = me;
  expectation = e;
  alignmentLength = al;
  matches = m;
  affyN = an;
}

blastMatch::blastMatch(){
  // do nothing it should be OK.
}

blastMatch::blastMatch(string mI, int l, int al){
  matchId = mI;
  length = l;
  af_length = al;
}

blastGenomeMatch::blastGenomeMatch(){
  // do nothing.. 
}

blastGenomeMatch::blastGenomeMatch(string c, int l, int m, int h, int g, int af, bool o){
  chromosome = c;
  lowerIndex = l;
  midIndex = m;
  higherIndex = h;
  guessGene = g;
  af_length = af;
  overlap = o;
}
