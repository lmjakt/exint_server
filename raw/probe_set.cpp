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
#include <libpq++.h>
#include <algorithm>
#include <stdlib.h>
#include "probe_set.h"

using namespace std;

float* meanValues(vector< vector<float> >& p){
  // return the mean values..
  // assume length is the same for all. change this later if we need to.
  //vector<float> mean;
  //if(p.size() == 0){
  //  return(mean);
  //}
  float* mean = new float[p[0].size()];
  //  mean.resize(p[0].size());
  for(int i=0; i < p.size(); i++){
    for(int j=0; j < p[i].size(); j++){
      mean[j] += p[i][j];
    }
  }
  for(int i=0; i < p[0].size(); i++){
    mean[i] = mean[i]/p.size();
  }
  return(mean);
}

//probe_set::probe_set(){
//}

probe_set::probe_set(){
  index = 0;
  probeNo = 0;
  exptSize = 0;
  allExptNo = 0;
  exptIndex = 0;
  exptLookup = 0;
}

probe_set::probe_set(int i, vector< vector<float> >& p, vector<uint>& eindices, uint allESize){
  //  cout << "Creating probe set for index " << i << "  with address " << this << endl;
  index = i;

  probes = new float*[p.size()];
  probeNo = p.size();
  exptSize = p[0].size();
  allExptNo = allESize;
  for(int i=0; i < p.size(); i++){
    probes[i] = new float[p[i].size()];
    for(int j=0; j < p[i].size(); j++){
      probes[i][j] = p[i][j];
    }
  }
  if(eindices.size() != exptSize){
    cerr << "CRASH CRASH, exptSize is not the same as eindices.size " << endl;
    exit(1);
  }
  // set up the exptLookup, with all -1 values.. 
  exptLookup = new int[allESize];
  for(int i=0; i < allESize; i++){ 
    exptLookup[i] = -1; 
  }     // default to -1;

  exptIndex = new uint[eindices.size()];
  for(int i=0; i < eindices.size(); i++){
    exptIndex[i] = eindices[i];
    exptLookup[eindices[i]] = i;        // and here we hope that everything is OK.. no error checking. it's faster this way.
  }
  //probes = p;
  mean = meanValues(p);
  //cout << "finished making probe set for index " << index << endl;
  //  for(int i=0; i < eindices.size(); i++){
  //  exptLookup.insert(make_pair(eindices[i], i));      // no checking so its up to the client.. bad, but faster.. 
  //}
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

// moved to dataStructs.cpp
// void probeSetMatchSet::insertMatch(probeSetMatch* psm){
//   matches.insert(psm);
//   // increment the sums.. if this region does not overlap with any of the previous matches..
//   ///// NOTE that we will here be screwed by multiple matches to the same sequence in the affymetrix probe set
//   ///// as may be the case if it mathces to repetitive sequences. However, I cannot at the present think of a good
//   ///// way of compansating for this. All ways I can think of are either too complicated, or don't work properly.
//   ///// It is also the case that perhaps if we have multiple mathces this does increase the likelihood of linkage 
//   ///// to a given gene, though, if it is only small part of the probe set then we will probably be screwed.
//   ///// For now, I will leave it, and see how well this simple stuff works, and if it seems to be a problem, then
//   ///// I will have to use one of the more complicated things... 
  
//   matchSum += (psm->pEnd - psm->pStart);
//   //score += ((float)(psm->pEnd - psm->pStart)) * psm->percent/(float)100; -- this is already set at the beginning so no need to worry about it.. 
//   //mismatchSum += (psm->alignLength - psm->match);  // ok.. !!
//   //expectProduct *= psm->expectation;

//   //if(psm->cStart < minPos){ minPos = psm->cStart; }
//   //if(psm->cEnd > maxPos){ maxPos = psm->cEnd; }
//   //if(psm->expectation < minExpect){ minExpect = psm->expectation; }
// }

// below moved to datastructs.cpp
// void ishProbeMatchSet::insertMatch(ishProbeMatch* ipm){
//   matches.insert(ipm);
//   matchSum += ipm->match;
//   mismatchSum += (ipm->alignLength - ipm->match);
//   expectProduct *= ipm->expectation;
//   if(ipm->cStart < minPos){ minPos = ipm->cStart; }
//   if(ipm->cEnd > maxPos){ maxPos = ipm->cEnd; }
//   if(ipm->expectation < minExpect){ minExpect = ipm->expectation; }
// }

void ishProbeData::addIshProbeMatch(ishProbeMatch* ipm){
  int maxRange = 100000;    // the size of an indiviidual thingy..
  //vector<ishProbeMatchSet>::iterator it;         .... -- gives me a parse error for some unknown reason. I really don't understand why.
  for(int i=0; i < probeMatches.size(); i++){
    //  for(it = probeMatches.begin(); it != probeMatches.end(); it++){
    ishProbeMatchSet* it = probeMatches[i];   // historical,, for shorthand.. 
    if(it->assemblyId == ipm->assemblyId){
    //  if((*it).chromosome == ipm->chromosome && abs((*it).minPos - ipm->cStart) < maxRange && abs((*it).maxPos - ipm->cEnd) < maxRange && ipm->probeIndex == (*it).dbIndex && (*it).strand == ipm->strand){
      //if((*it).chromosome == ipm->chromosome && abs((*it).minPos - ipm->cStart) < maxRange && abs((*it).maxPos - ipm->cEnd) < maxRange && ipm->probeIndex == (*it).dbIndex && (*it).strand == ipm->strand){
      (*it).insertMatch(ipm);
      sort(probeMatches.begin(), probeMatches.end());    // which puts the one with the highest score first.. though it doesn't actually matter that much.. 
      return;
    }
  }
  //  probeMatches.push_back(ishProbeMatchSet(ipm));
  //sort(probeMatches.begin(), probeMatches.end());
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
      //     (*it).matches.insert(psm);
      //if(psm->cStart < (*it).minPos){ (*it).minPos = psm->cStart; }
      //if(psm->cEnd > (*it).maxPos){ (*it).maxPos = psm->cEnd; }
      //if(psm->expectation < (*it).minExpect){ (*it).minExpect = psm->expectation; }
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
  //const char* dbname = "dbname=expression";
  PgCursor data(conninfo, "portal");    // still don't know why the portal
  // check the backend..
  if( data.ConnectionBad() ) {
    cerr << "Connection to database '" << conninfo << "' failed." << endl
	 << "Error returned: " << data.ErrorMessage() << endl;
    exit(1);
  }
  // lets send some commands to the database FOR THE MOMENT: limit to data for chip 1 to keep the thingy..
  //  if( !data.Declare("select a.index, a.af_id, a.id, a.suid, b.gene, b.title, a.description, c.description from p_sets a, uni_data b, tigr_annotation c where a.suid=b.index and a.af_id = c.af_id")){

  //// IMPORTANT NOTE : the probe_data table has a stupid structure and can have several lines for each probe set, even though it really shouldn't. It's based on a 
  //// join between a few different tables, (don't remember exactly, but .. ) afid_uid_map, p_sets, uni_data and tigr_annotation
  //// -- and I will need to do something to update this when adding new probes.. hmmm. hmmmm.


  if(!data.Declare("select a.*, b.chip from probe_data a, p_sets b where a.index=b.index order by index")){          // table created by the above joined, and stored in db for speed of access..
    cerr << "DECLARE CURSOR command failed" << endl;
    exit(1);
  }
  // then fetch the data... 
  if( !data.Fetch() ) {
    cerr << "FETCH ALL command didn't return tuples properly" << endl;
    exit(1);
  }
  int tuples = data.Tuples();
  // create a vector of probe_data ..
  vector<probe_data> p_data(tuples);  // with the appropriate amount of space..
  // resize the vector afterwards. This is very ugly, but the best I can think of
  // at the time. I really need to sort this out, but it will have to do for now..
  // and go throuh and get all the thingies..
  int index;
  int p_index;
  int maxIndex = 0;
  for(int i=0; i < tuples; i++){
    index = atoi(data.GetValue(i,0));
    // for this one we shall be defining the p_data[index-i]
    p_index = index-1;
    int ugid = atoi(data.GetValue(i, 3));
    string uggene = data.GetValue(i, 4);
    string ugtitle = data.GetValue(i, 5);
    //uniGeneData tempdata(ugid, ugtitle, uggene);
    p_data[p_index].ugData.push_back(uniGeneData(ugid, ugtitle, uggene));
    if(!p_data[p_index].defined){
      p_data[p_index].index = index;
      p_data[p_index].defined = true;
      p_data[p_index].afid = data.GetValue(i, 1);
      p_data[p_index].gbid = data.GetValue(i, 2);
      p_data[p_index].afdes = data.GetValue(i, 6);
      p_data[p_index].tigrDescription = data.GetValue(i, 7);
      p_data[p_index].chip = atoi(data.GetValue(i, 8));
      p_data[p_index].go.resize(0);         // don't know if that's good or not.. 
    }
    if(index > maxIndex){
      maxIndex = index;
    }
  }
  //close the cursor,, and lets leave it at that for the time being, it won't work anyway..
  p_data.resize(maxIndex);         // hopefully this will now work.. 
  data.Close();
  // then lets see if we can get ourselves some celera mapping,, I'm just worried that hmmm, but what the hell.. 
  // if(! data.Declare("select * from affy_cel_match_annot")){
//     cerr << "BUGGER, couldn't select from affy_cel_match_annot" << endl;
//   }
//   if(! data.Fetch() ){
//     cerr << "couldn't fetch the affy cel match data" << endl;
//   }
//   cout << "got " << data.Tuples() << " tuples for the affy_cel_match query " << endl;
//   int push_back_counter = 0;
//   for(int i=0; i < data.Tuples(); i++){
//     // let's hope we don't have a problem with memory allocation here.. hold on to your horses
//     // this will be a little difficult..
//     int p_index = atoi(data.GetValue(i, 0));
//     p_index--;
//     string CG = data.GetValue(i, 1);
//     float exp = atof(data.GetValue(i, 2));
//     float m = atof(data.GetValue(i, 3));
//     string SF = data.GetValue(i, 4);
//     string FN = data.GetValue(i, 5);
//     string GN = data.GetValue(i, 6);
//     string GS = data.GetValue(i, 7);
//     string ND = data.GetValue(i, 8);
//     // cout << "index: " << p_index << "\tCG : " << CG << "\tSF : " << SF << endl;
//     //cout << "internal index: " << p_data[p_index].index << endl;
//     p_data[p_index].celeraMatches.push_back(celeraMatch(CG, exp, m, SF, FN, GN, GS, ND));  // 
//     push_back_counter++;
//   }
//   cout << "pushed back a total of " << push_back_counter << endl;
//   data.Close();

  // then see if we can get ourselves some more data for the go..
  if( !data.Declare("select a.index, b.generation, c.description from p_sets a, af_go_gen b, go c where a.af_id=b.af_id and b.go=c.index and a.chip=1 order by a.index, b.generation desc")){
    cerr << "GO data declare didn't work so well.. " << endl;
    exit(1);
  }
  // fetch..
  if( !data.Fetch() ) {
    cerr << "Couldn't fetch all the GO data " << endl;
    exit(1);
  }
  int go_gen;
  tuples = data.Tuples();
  for(int i=0; i < tuples; i++){
    index = atoi(data.GetValue(i, 0));
    p_index = index-1;
    go_gen = atoi(data.GetValue(i, 1));
    if(p_data[p_index].go.size() < (go_gen)){
      p_data[p_index].go.resize(15);             // I DON'T UNDERSTAND THIS, SOMEONE PLS. EXPLAIN..
    }
    p_data[p_index].go[go_gen-1].push_back(data.GetValue(i,2));
  }
  data.Close();
  // And now finally, lets see if we can load the blastMatches,, just using the values in the thingy.. hmm.. 
  // the table ensembl_blast_matches has the following struture (at least at the moment..)..
  //
  //  afid                 int
  //  blast_match          int
  //  af_length            int
  //  ensembl_transcript   int
  //  ensembl_gid          int
  //  ensembl_length       int
  //  expectation          int
  //  af_start             int
  //  af_end               int
  //  ensembl_start        int
  //  ensembl_end          int
  //  align_length         int
  //  match                int
  //  af_n_count           int
  //
  //   -- ofcourse the ensembl_transcript and the ensemb_gid are ints referring to the tables 
  //  ensembl_genes and ensembl_transcripts,, -- 
  // 
  // to get information about the ensembl_transcripts, look in the table ensembl_annotation, -and to see 
  // what the fields are look in ensembl_fields.. 
  // 
  //        Now the question is this.. where to put the ensembl annotation. -We can just dump it, or parts of
  //        it into the the probeata_data structure, and this is sort of OK, -but it is likely to waste memory 
  //        and bandwidth, -and as such is not desirable. However the alternative is to put it into some maps
  //        in the probeSetSet2 structure, -and then just refer to it from the probe_data struct. -- However, this
  //        now becomes complicated, as I have to come up with a way of serialising this for sending down the line.
  //        -- I could have different kind of structures at either end of the pipe, -or I could get the client to request
  //        the data each time it comes across it. -- But all of these solutions seem a bit painful. The easiest way is 
  //        to stick the annotation into the blastMatch field,, -but I probably have to do this first, -as I can't really 
  //        do a sensible join on the blast match table.. Then I need data structs to keep track of where I put stuff. That's
  //        kind of painful, and seems as if it might be fragile, as I'll be creating lots of implicit dependancies. Hmm.
  //        have to give this some thought.. 
  
  // First do a join between the ensembl_blast_matches table and the annotation table and create
  // new blastMatch struct's for the appropriate probe_data records.. (in the p_data vector);
  //map<int, map<int, int> > blastMatchMap;          // key is the affy id, the inner map first is the blast match Id, and the final int is the vector position.. -- oh dear. sounds 
//   if(!data.Declare("select distinct b.afid, b.blast_match, a.name, c.field_name, d.annotation, b.ensembl_length, b.af_length from ensembl_genes a, ensembl_blast_matches b, ensembl_fields c, ensembl_annotation d where a.index=b.ensembl_gid and a.index=d.gene and c.index=d.field and (c.index = 1 or c.index = 3 or c.index = 15 or c.index = 21) order by b.blast_match")){
//     cerr << "Couldn't select the annotation for the unigene things. and so forth and o son" << endl;
//     cerr << data.ErrorMessage() << endl;
//   }
//   if(!data.Fetch()){            // this will probably use up a huge amount of memory, but what the hell, lets get a bigger machine
//     cerr << "Couldn't fetch the annotation for the ensembl blast matches" << endl;
//     cerr << data.ErrorMessage() << endl;
//   }
//   set<int> blastMatchIdSet;        // have I come across this one before.. 
//   map<int, int> blastMatchIdIndex; // map<blastMatchId, ensemblMatchIndex> --  as blastMatchId is unique for thingy
//   //int tuples = data.Tuples();        // this might be a case of causeing trouble again..
//   for(int i=0; i < data.Tuples(); i++){
//     int affyIndex = atoi(data.GetValue(i, 0));
//     affyIndex--;     // to make it fit with p_data counting.. hmmm... 
//     int blastMatchId = atoi(data.GetValue(i, 1));
//     string ensembl_gid = data.GetValue(i, 2);
//     string fieldName = data.GetValue(i, 3);
//     string annotation = data.GetValue(i, 4);
//     int ensembl_length = atoi(data.GetValue(i, 5));
//     int af_length = atoi(data.GetValue(i, 6));
//     //    cout << affyIndex << "\t" << blastMatchId << "\t" << ensembl_gid << "\t" << fieldName << "\t" << annotation << endl;
//     if(blastMatchIdSet.count(blastMatchId) == 0){
//       // its a new blast MatchId, and we need to push back the appropriate thingy
//       // and thingy.. and so on. 
//       blastMatchIdIndex.insert(make_pair(blastMatchId, p_data[affyIndex].ensemblMatches.size())); // hooooaaahhh. 
//       p_data[affyIndex].ensemblMatches.push_back(blastMatch(ensembl_gid, ensembl_length, af_length));
//       blastMatchIdSet.insert(blastMatchId);
//     }else{
//       /// and this is the painful bit. The blastmatch id has already been defined, in the given probe data.. -
//       /// because the ordering is by the blast_match, it must be the last one in the vector ensemblMatches
//       /// so we can assuem that..
//       int lastIndex = p_data[affyIndex].ensemblMatches.size() - 1;
//       p_data[affyIndex].ensemblMatches[lastIndex].description.insert(make_pair(fieldName, annotation));
//       // and that should be it. really..
//     }
//   }
//   /// and here is where it gets very messy. Because I don't know how many blastMatches for each 
//   /// probe_data, and I don't know how many alignments for each one of these.. but maybe I can do anyway
//   ///  -- use the blastMatchIdIndex,, to work out what goes where.. 
//   data.Close();
//   // and 
//   if(!data.Declare("select * from ensembl_blast_matches")){
//     cerr << "Couldn't declare the blast alignments query" << endl;
//     cerr << data.ErrorMessage() << endl;
//   }
//   if(!data.Fetch()){
//     cerr << "Couldn't fetch the blast alignments " << endl;
//     cerr << data.ErrorMessage() << endl;
//   }
//   for(int i=0; i < data.Tuples(); i++){
//     int afid = atoi(data.GetValue(i, 0));
//     afid--;
//     int blastMatchId = atoi(data.GetValue(i, 1));
//     int af_length = atoi(data.GetValue(i, 2));
//     int ensembl_length = atoi(data.GetValue(i, 5));
//     QString expectString(data.GetValue(i, 6));
//     bool ok; 
//     double expect = expectString.toDouble(&ok);
//     //float expectD = expectString.toDouble(&ok);
//     //    float expect = atof(
//     int af_start = atoi(data.GetValue(i, 7));
//     int af_end = atoi(data.GetValue(i, 8));
//     int ensembl_start = atoi(data.GetValue(i, 9));
//     int ensembl_end = atoi(data.GetValue(i, 10));
//     int align_length = atoi(data.GetValue(i, 11));
//     int match_length = atoi(data.GetValue(i, 12));
//     int af_n_count = atoi(data.GetValue(i, 13));
//     //cout << "blast Match Id: " << blastMatchId << "\tafid: " << afid << "\tstring: " << expectString << "\texpect: " 
//     //	 << expect << "\tdouble: " << expectD << "\tok: " << ok << endl;
//     // and lets make sure that we have defined the 
//     if(blastMatchIdIndex.count(blastMatchId) > 0){
//       int indexPosition = blastMatchIdIndex[blastMatchId];  // and this should be OK.. should check..
//       if(p_data[afid].ensemblMatches.size() > indexPosition){
// 	p_data[afid].ensemblMatches[indexPosition].alignments.push_back(blastAlignment(af_start, af_end, ensembl_start, ensembl_end, expect, align_length, match_length, af_n_count));
// 	// and that should be everything..
//       }else{
// 	cerr << "ensembl Index out index, for affy id " << afid << "\tblastMatchId " << blastMatchId << endl;
//       }
//     }else{
//       cerr << "no blastMatchIdIndex defined for " << blastMatchId << endl;
//     }
//   }
//   data.Close();


  ////////////////// I SHOULD REMOVE THE BELOW AS I DON'T ACTUALLY USE THE BELOW DATA STRUCT AT ALL.. --
  /////////////////// IT DOESN'T ACTUALLY CHANGE MUCH AS THE MEMORY USAGE IS NOT TOO BAD.. 
//   if(!data.Declare("select afid, af_length, chromosome, lower_gene, mid_gene, higher_gene, guess_gene, overlap, expect, af_start, af_end, ensembl_start, ensembl_end, align_length, match, af_n_count from ensembl_genome_matches")){
//     cerr << "Couldn't declare the ensembl_genome_matches" << endl
// 	 << data.ErrorMessage() << endl;
//   }
//   if(data.Fetch()){
//     cerr << "Couldn't fetch the ensemble_genome_matchs" << endl
// 	 << data.ErrorMessage() << endl;
//   }
//   for(int i=0; i < data.Tuples(); i++){        // this is pretty big one !!
//     // The current structure of the table is as follows..
//     int affyIndex = atoi(data.GetValue(i, 0));
//     affyIndex--;
//     if(affyIndex >= p_data.size()){ continue; } 
//     int lowerI = atoi(data.GetValue(i, 3));
//     // check to see if we have defined a blast Match for the thingy.. 
//     map<int, blastGenomeMatch>::iterator it;
//     it = p_data[affyIndex].genomeMatches.find(lowerI);
//     if(it == p_data[affyIndex].genomeMatches.end()){        // i.e. we haven't defined one for this one yet.. 
//       bool overlap = false;
//       if(data.GetValue(i, 7) == "t") { overlap = true; }
//       p_data[affyIndex].genomeMatches.insert(make_pair(lowerI, blastGenomeMatch(data.GetValue(i, 2), lowerI, atoi(data.GetValue(i, 4)), atoi(data.GetValue(i, 5)), atoi(data.GetValue(i, 6)), atoi(data.GetValue(i, 1)), overlap)));
//       // need to find again so we can still use it!!! 
//       it = p_data[affyIndex].genomeMatches.find(lowerI);        // which should now make more sense..
//     }
//     QString expectString = data.GetValue(i, 8);
//     bool ok;
//     double expect = expectString.toDouble(&ok);
//     (*it).second.alignments.push_back(blastAlignment(atoi(data.GetValue(i, 9)), atoi(data.GetValue(i, 10)), atoi(data.GetValue(i, 11)), atoi(data.GetValue(i, 12)), expect, atoi(data.GetValue(i, 13)), atoi(data.GetValue(i, 14)), atoi(data.GetValue(i, 15))));
//     cout << "pushing back the genome alignments for afid " << affyIndex << endl;
    
//   }
  // that should be all of the things we need for this purpose stored away in the database.. so it should be OK.. to proceed.. 
  
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
