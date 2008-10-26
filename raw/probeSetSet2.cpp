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

#include "probeSetSet2.h"
#include "../server/server.h"   // for experimental server definition.. 
#include "probe_set.h"
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <libpq++.h>
#include <qstring.h>
#include <stdlib.h>
#include <qmutex.h>

using namespace std;

vector<float> splitString(string numbers){
  // assumes structure is.. {"8.28282", "8.7837838", "2.383838"}
  vector<float> f;      // the real numbers!!!
  int first = 0;
  int next;             // for the positions..
  while(1){
    first = numbers.find_first_not_of(",{}", first);
    next = numbers.find_first_of(",}", first+1);
    //cout << "first: " << first << "\tnext: " << next << endl;
    if(first == numbers.npos || next == numbers.npos){
      break;
    }
    f.push_back(atof(numbers.substr(first, next-first).c_str()));
    first = next+1;
  }
  return(f);
}

  
ProbeSetSet2::ProbeSetSet2(){
  // do nothing;;
}   

ProbeSetSet2::ProbeSetSet2(const char* cinfo, const char* dataTable){
  // get the experimental data from the appropriate tables.
  conninfo = cinfo;
  probeData = data_from_db(conninfo);  // data for each probe set
  setExpData(conninfo);          // experiment and genomic data.. 
  cout << "setExpData done .. " << endl;
  //// guess gene linkage.. if nothing, using some stupid bloody algorithm..
  guessGenes();
  ishProbeDataMutex = new QMutex();
  sessionMutex = new QMutex();
  annotationMutex = new QMutex();
  userMutex = new QMutex();
  chromFileMutex = new QMutex;
  chromAnnotationMutex = new QMutex;
  chipMutex = new QMutex();

  cout << "making empty probe sets.. " << endl;
  // initialise the probe_set vector to have an entry for every probe_data entry..
  emptySet = new probe_set();  // ?
  data.resize(probeData.size());
  for(int i=0; i < probeData.size(); i++){
    data[i] = emptySet;   // so we start counting at the same point.. 
  }
  cout << "about to set data .. " << endl;
  setData(conninfo, dataTable);
}

void ProbeSetSet2::setData(const char* conninfo, const char* tableName){
  //cout << "at the beginning of setData" << endl;
  //string conninfo("dbname=");
  //conninfo += dbname;
  PgCursor cursor(conninfo, "portal");
  if(cursor.ConnectionBad()){
    cerr << "Connection to Database " << conninfo << " failed"
	 << "Error returned:        " << cursor.ErrorMessage() << endl;
    return;
  }
  ostringstream os;
  os << "select * from " << tableName << " order by probe, experiment";

  const char* query = os.str().c_str(); //"select * from data order by probe, experiment";
#ifdef EXPERIMENTAL_SERVER
  query = "select * from little_data order by probe, experiment";
  cout << "query should now be set to " << query << endl;
#endif
  cout << "query for the data base is : " << query << endl;
  //  if(!cursor.Declare("select * from data order by probe, experiment")){    // should make an if_def here.. 
  if(!cursor.Declare(query)){
    // if(!cursor.Declare("select * from data where probe < 1000 order by probe, experiment")){
    cerr << "Failed to declare data retrieval query" << endl;
    return;
  }
  int lastIndex =  -1; //atoi(data.GetValue(0, 0));
  int currentIndex;
  vector<float> pm;
  vector<float> mm;
  vector<uint> exptIndices;        // work it out..!!
  vector< vector<float> > delta;  // for the difference values.. 
  string pmString;
  string mmString;
  float exp;
  data.reserve(40000);
  int counter = 0;
  while(cursor.Fetch(1000)){
    if(cursor.Tuples() == 0){
      cerr << "fetched data, but didn't return any tuples, returning with nothing" << endl;
      return;
    }
    cout << "Getting " << cursor.Tuples() <<  " tuples   : " << counter  << endl;
    counter += 1000;
    for(int i=0; i < cursor.Tuples(); i++){
      //   cout << "\ti : " << i << endl;
      currentIndex = atoi(cursor.GetValue(i, 0));
      //cout << "current Index: " << currentIndex << endl;
      exp = atof(cursor.GetValue(i, 1));
      pmString = cursor.GetValue(i, 2);
      mmString = cursor.GetValue(i, 3);
      //cout << "pmString : :" << pmString << endl;
      pm = splitString(pmString);
      //cout << "Size of resulting pm vector " << pm.size() << endl;
      mm = splitString(mmString);
      //cout << "managed to split the strings!! " << endl;
      if(currentIndex != lastIndex){
	if(lastIndex != -1){
	  //cout << "current Index: " << currentIndex << "\tlast Index: " << lastIndex << endl;
	  int indexPos = lastIndex - 1;
	  if(indexPos < data.size() && probeData[indexPos].index == lastIndex){
	    data[indexPos] = new probe_set(lastIndex, delta, exptIndices, experiments.size());
	    //cout << "creating a probe set with index : " << lastIndex << "  at positions ; " << indexPos << endl;
	  }else{
	    cerr << "no index for probe set " << lastIndex << " don't know why .. " << endl;
	    cerr << "index pos is : " << indexPos << " probeData[" << indexPos << "].index is " << probeData[indexPos].index << endl; 
	    exit(1);
	  }
	}
	//	  data.push_back(new probe_set(lastIndex, delta, exptIndices, experiments.size()));
	// cout << "pushing back data for index: " << lastIndex << endl;
	delta.resize(pm.size());
	for(int i=0; i < delta.size(); i++){
	  delta[i].resize(0);
	}
	exptIndices.resize(0);
	lastIndex = currentIndex;
	// i.e. create a probeSet, and 
      }
      if(experiments.find(exp) == experiments.end()){
	// PUKE.. HAVE TO REWRITE COMPLETELY
	cerr << "PUKE,,, PUKE,,, PUKE !!! CAN'T USE STRUCTURES LIKE THESE, REWRITE APPLICATION YOU FOOOL!!" << endl;
	cerr << "couldn't find the experiment for experiment: " << exp << endl;
	return;
      }
      exptIndices.push_back(experiments[exp].dbaseIndex);
      if(pm.size() != delta.size() || mm.size() != delta.size()){
	cerr << "PUKE, going to segment fault if I don't sort this out.. vector sizes are different to each other.. " << endl;
	return;
      }
      for(int j=0; j < pm.size(); j++){
	delta[j].push_back(pm[j]-mm[j]);     // the actual thing!!!
      }
    }
    // and don't forget to do the last one..
  }
  cout << "Cursor doesn't seem to fetch anything anymore, what's going on " << endl;
  data.push_back(new probe_set(lastIndex, delta, exptIndices, experiments.size()));
  // and that should be it. Ugly as hell if you ask me.. and very sensitive..
}


void ProbeSetSet2::setExpData(const char* conninfo){
  // make a connection to the database..
  cout << "at the beginning of setExpData.. " << endl;
  //string conninfo("dbname=");
  //conninfo += dbname;
  PgCursor cursor(conninfo, "portal");
  if(cursor.ConnectionBad()){
    cerr << "Connection to Database " << conninfo << " failed"
	 << "Error returned: " << cursor.ErrorMessage() << endl;
    return;
  }
  // and lets do some queries.
  // get, experiments (floats), experimental short descriptions, experimental group,
  // and others from experiments and experiment_group including colours for the various things!!
  string query = "select index, experiment, expt_group, short_description, description from experiments order by experiment";
#ifdef EXPERIMENTAL_SERVER
  query = "select index, experiment, expt_group, short_description, description from little_experiments order by experiment";
#endif

  if(!cursor.Declare(query.c_str())){
    //  if(!cursor.Declare("select index, experiment, expt_group, short_description, description from experiments order by experiment")){
    cerr << "DECLARE CURSOR command failed" << endl;
    return;
  }
  if(!cursor.Fetch()){
    cerr << "FETCH command failed" << endl;
    return;
  }
  for(int i=0; i < cursor.Tuples(); i++){
    float experimentIndex = atof(cursor.GetValue(i, 1));
    experiments[experimentIndex].dbaseIndex = i;
    experiments[experimentIndex].realDbaseIndex = (atoi(cursor.GetValue(i, 0))-1);        // count from 0 in the application. 
    experiments[experimentIndex].exptGroup = atoi(cursor.GetValue(i, 2));
    experiments[experimentIndex].shortName = cursor.GetValue(i, 3);
    experiments[experimentIndex].description = cursor.GetValue(i, 4);
    experiments[experimentIndex].index = experimentIndex;
  }
  // close data, and let's work out which chips have been used..
  cursor.Close();
  // first lets get the chip descriptions available
  if(!cursor.Declare("select distinct index, id, description from chips")){
    cerr << "DECLARE CURSOR failed on chip table" << endl;
    return;
  }
  if(!cursor.Fetch()){
    cerr << "Fetch failed for the chips.. " << cursor.ErrorMessage() << endl;
    return;
  }
  for(int i=0; i < cursor.Tuples(); i++){
    int chipInt = atoi(cursor.GetValue(i, 0));
    string chipId = cursor.GetValue(i, 1);
    string chipDescription = cursor.GetValue(i, 2);
    chipDescriptions[chipInt] = chipInfo(chipInt, chipId, chipDescription);
    // set the bool used to false for the different experiments..
    //    map<float, exInfo>::iterator it;
    //for(it = experiments.begin(); it != experiments.end(); it++){
    //  (*it).second.chips[chipInt] = false;    // also inserts it.. 
    //}   // set to boolean later on if we find that they actually have been used.
  }
  cursor.Close();
  // and then lets look at the chipequivalents..
  if(!cursor.Declare("select * from chip_equivalents")){
    cerr << "Couldn't get chip equivalents, this is likely to cause real problems so lets die here " << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "Couldn't fetch the chip equivalents, lets die.. " << endl;
    exit(1);
  }
  for(uint i=0; i < cursor.Tuples(); i++){
    int chipId = atoi(cursor.GetValue(i, 0));
    int chipEquiv = atoi(cursor.GetValue(i, 1));
    map<int, chipInfo>::iterator it = chipDescriptions.find(chipId);
    if(it != chipDescriptions.end()){
      (*it).second.equivs.insert(chipEquiv);
    }
    map<float, exInfo>::iterator eit;
    for(eit = experiments.begin(); eit != experiments.end(); eit++){
      (*eit).second.chips[chipEquiv] = false;    // also inserts it.. 
    }   // set to boolean later on if we find that they actually have been used.

  }
  cursor.Close();
  // ok, then let's look at the cel_files table to see which of these chips have been 
  // used for the different experiments. This is not really a good way to do it, as it's
  // kind of accidental, but what the hell..

  query = "select distinct a.experiment, b.equivalent from cel_files a, chip_equivalents b where a.chip=b.chip";
  //  query = "select distinct experiment, chip from cel_files order by experiment";
#ifdef EXPERIMENTAL_SERVER
  query = "select distinct a.experiment, b.equivalent from little_cel_files a, chip_equivalents b where a.chip=b.chip";
  //  query = "select distinct experiment, chip from little_cel_files order by experiment";
#endif
  if(!cursor.Declare(query.c_str())){
    //  if(!cursor.Declare("select distinct experiment, chip from cel_files order by experiment")){
    cerr << "couldn't get the chip usage from the cel_files table" << endl;
    return;
  }
  if(!cursor.Fetch()){
    cerr << "couldnt' fetch the data. BUGGER.. " << endl;
    return;
  }
  for(int i=0; i < cursor.Tuples(); i++){
    float expIndex = atof(cursor.GetValue(i, 0));
    int chipIndex = atoi(cursor.GetValue(i, 1));
    if(experiments.find(expIndex) != experiments.end()){
      if(experiments[expIndex].chips.find(chipIndex) != experiments[expIndex].chips.end()){
	experiments[expIndex].chips[chipIndex] = true;
      }else{
	cerr << "couldn't find the chip index: " << endl
	     << "experiment Index: " << expIndex << endl
	     << "chip Index      : " << chipIndex << endl;
      }
    }else{
      cerr << "UNKNOWN EXPERIMENTAL INDEX, WHAT'S GOING ON HERE?? THIS IS JUST STRANGE" << endl;
    }
  }
  cursor.Close();
  // first lets get the fields from the ensembl_fields table..
  if(!cursor.Declare("select * from ensembl_fields_9")){
    cerr << "Couldn't declare on ensembl_fields " << endl
	 << cursor.ErrorMessage() << endl;
  }
  if(!cursor.Fetch()){
    cerr << "Couldn't fetch the ensembl fields" << endl
	 << cursor.ErrorMessage() << endl;
  }
  for(int i=0; i < cursor.Tuples(); i++){
    ensemblFields.insert(make_pair(atoi(cursor.GetValue(i, 0)), cursor.GetValue(i, 1)));
  }
  cursor.Close();
  // and lets get the actual annotation, just get all of it, we'll get more RAM soon, and 
  // anyway, it will pale in comparison to the amoutn we deal with later..
  if(!cursor.Declare("select * from ensembl_annotation_9")){
    cerr << "Couldnt declare cursor for the annotation" << endl
	 << cursor.ErrorMessage() << endl;
  }
  if(!cursor.Fetch()){
    cerr << "Couldn't fetch the annotation " << endl
	 << cursor.ErrorMessage() << endl;
  }
  map<int, multimap<int, string> >::iterator it;
  for(int i=0; i < cursor.Tuples(); i++){
    int enIndex = atoi(cursor.GetValue(i, 0));
    int enField = atoi(cursor.GetValue(i, 1));
    string annot = cursor.GetValue(i, 2);
    /// now here is a little tricky as we are dealing with a map and a multimap. 
    it = ensemblAnnotation.find(enIndex);
    if(it == ensemblAnnotation.end()){       // i.e. not defined,, then we just insert the pair of pairs.
      multimap<int, string> temp;
      temp.insert(make_pair(enField, annot));
      ensemblAnnotation.insert(make_pair(enIndex, temp));
      temp.erase(enField);     // just in case.. 
    }else{
      (*it).second.insert(make_pair(enField, annot));
    }
  }
  cursor.Close();
  /// Ok let's get the sessionInformation and populate the different maps and things..
  /// First the sessioninformation and session members.
  // ignore 0, Null session stuff. Hooah...
  // Unfortunately, I think I will actually get everything in 3 different queries,
  // one for the session information, one for the session keywords, and one for the associated
  // genes. That sucks badly, but other options are also not very nice. OK.. 
  if(!cursor.Declare("select index, user_id, title, description from sessions where index != 0")){
    cerr << "Couldn't declare the sessions query, dying." << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "couldn't fetch the session info" << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  set<int> tv;    // hmmm, not so nice.. but what the hell.
  int owner;
  int index;
  string title;
  string description;
  vector<string> nkw;   // null keywords, push later, this really sucks.. but can't be bothered to rewrite
  for(int i=0; i < cursor.Tuples(); i++){
    index = atoi(cursor.GetValue(i, 0));
    owner = atoi(cursor.GetValue(i, 1));
    title = cursor.GetValue(i, 2);
    description = cursor.GetValue(i, 3);
    sessions.insert(make_pair(index, sessionInformation(index, owner, title, description, nkw, tv)));
  }
  cursor.Close();
  // and now get the keywords.. 
  if(!cursor.Declare("select index, keyword from session_keywords")){
    cerr << "Couldn't select from keywords table" << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "Couldn't fetch the keywords..  bummer man" << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  map<int, sessionInformation>::iterator sit;
  string kw;
  for(int i=0; i < cursor.Tuples(); i++){
    index = atoi(cursor.GetValue(i, 0));
    kw = cursor.GetValue(i, 1);
    sit = sessions.find(index);
    if(sit != sessions.end()){
      (*sit).second.keywords.push_back(kw);
    }
  }
  cursor.Close();
  // and now do something a little more complicated. Bit wasteful but to get the actual 
  // session members, I need to do a 3 table join, and then do things. But otherwise pretty similar
  // to the above..
  if(!cursor.Declare("select distinct a.index, c.gene from sessions a, user_annotation b, user_annotation_genes c where a.index !=0 and a.index=b.session_id and b.index=c.annotation_id")){
    cerr << "couldn't get session members,, bugger " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "Couldn't fetch session members " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  int gIndex;
  for(int i=0; i < cursor.Tuples(); i++){
    index = atoi(cursor.GetValue(i, 0));
    gIndex = atoi(cursor.GetValue(i, 1));
    sit = sessions.find(index);
    if(sit != sessions.end()){
      (*sit).second.members.insert(gIndex);
      sessionLookup[gIndex].insert(index);
    }
  }
  cursor.Close();
  //////// and now that is essentially everything for the session information. Let's stop and see if it 
  ///////  compiles.. !!
  // OK let's get the user annotation. Use the rather crude filter,, -don't get empty string annotations,, -
  // these are most often just used for the insertion into sessions. Not so good, but can't think of something
  // perfect..
  if(!cursor.Declare("select index, user_id, annotation from user_annotation where annotation != ''")){
    cerr << "couldn't declare cursor for user_annotation " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "couldn't fetch user_annotation" << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  
  for(int i=0; i < cursor.Tuples(); i++){
    index = atoi(cursor.GetValue(i, 0));
    owner = atoi(cursor.GetValue(i, 1));
    description = cursor.GetValue(i, 2);
    userAnnotation.insert(make_pair(index, annotationInformation(index, owner, description, tv)));
  }
  // and lets get the genes.. .. hoo,, hooo yeahh.
  cursor.Close();
  if(!cursor.Declare("select annotation_id, gene from user_annotation_genes")){
    cerr << "couldn't declare cursor for user_annotation_genes" << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "couldn't fetch the genes for the user annotation, " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  map<int, annotationInformation>::iterator ait;
  for(int i=0; i < cursor.Tuples(); i++){
    index = atoi(cursor.GetValue(i, 0));
    gIndex = atoi(cursor.GetValue(i, 1));
    ait = userAnnotation.find(index);
    if(ait != userAnnotation.end()){
      (*ait).second.members.insert(gIndex);
      annotationLookup[gIndex].insert(index);
    }
  }
  cursor.Close();
  // and let's pick up the userInformation. Useful, so that users can see who actually thinks what..
  if(!cursor.Declare("select index, user_name, full_name, lab from users")){
    cerr << "Couldn't declare cursor for user table " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "couldn't fetch the user stuff " << endl
	 << cursor.ErrorMessage() << endl;
    exit(1);
  }
  for(int i=0; i < cursor.Tuples(); i++){
    int index = atoi(cursor.GetValue(i, 0));
    string user = cursor.GetValue(i, 1);
    string fullName = cursor.GetValue(i, 2);
    string labName = cursor.GetValue(i, 3);
    userTable.insert(make_pair(i, userInformation(index, user, fullName, labName)));
  }
  cursor.Close();   // hoo hhooo yeahhh....
  // and then lets pick up the files containing the chromosomal locations and create the map containing 
  // pointers to open ifstreams. -- these can then be used by the appropriate connectionobjects to read
  // sequence delineated by the user. -- Remember to protect these by the use of mutexes when seeking and 
  // reading from them, as different seeks could bugger everything up otherwise, and it might be tricky to 
  // work out what is going on.. -Obviously seekg is going to be a bit tricky otherwise. I think this may 
  // be kind of dangerous, but lets have a go anyway..
  vector<string> cFileNames;
  vector<string> chromNames; 
  int regionSize = 100000;    // divide chromosomal fragments into 100000 bp fragments.. 
  // although I don't need it at the moment, it might be an idea to pick up the chromosomal sizes here as well
  // so that I can create some sort of appropriate data structures filled with genomicRegions.. -- 
  // Anyway, I need this somewhere.. !!! .. 
  if(!cursor.Declare("select path, name, size from chrom_files_7")){
    cerr << "Couldn't declare cursor for chrom_files, might as well die. " << endl;
    exit(1);
  }
  if(!cursor.Fetch()){
    cerr << "Couldnt fetch cursor for chrom_files, dying" << endl;
    exit(1);
  }
  for(int i=0; i < cursor.Tuples(); i++){
    string temp = cursor.GetValue(i, 0);
    temp += "/";
    string temp2 = cursor.GetValue(i, 1);
    chromNames.push_back(temp2);
    int size = atoi(cursor.GetValue(i, 2));
    chromSizes.insert(make_pair(temp2, size));
    chromosomeAnnotation.insert(make_pair(temp2, new chromAnnotation(temp2, size, regionSize)));
    temp += temp2;
    cFileNames.push_back(temp);
  }
  // then just go and open all of the files..
  ifstream* tstream;   // hmm...
  for(int i=0; i < cFileNames.size(); i++){
    tstream = new ifstream(cFileNames[i].c_str(), ios::binary);
    if(!tstream){
      cerr << "Couldn't open file : " << cFileNames[i] << endl;
      exit(1);
    }
    chromFiles.insert(make_pair(chromNames[i], tstream));
  }
  /// and that then really should be everything. Lets close the cursor, and get on with it.
  cursor.Close();
  /// well as always there is more. Let's create chromosomal annotation. Base this first on the 
  /// probe set genome matches, -- so that we can look up a chromosomal region, and select probe 
  /// sets on the basis of the chromosomal map. Eventually tie the probe set identification more closely
  /// to individual gene loci, as far as we can. Later, I'll have to do some updating with new maps..
  /// but for now, well, we'll just include the information.. 
  string currentChromosome = "";
  map<string, chromAnnotation*>::iterator cait;
  cout << "OK I'm now at the beginning of setting data for the genomic alignments " << endl;
  /////////////////////////////// let's use the new set of tables containing information regarding alternative transcripts, 
  /////////////////////////////// and exon structure for this. .. All the chromosomal regions are already sorted so we don't need 
  ////////////////////////////// to make them, just to work out the location where to put the ensemble genes that we make. 
  //  if(!cursor.Declare("select a.chromosome, a.gene, a.id, a.strand, b.transcript, b.id, b.start, b.stop,  c.id, c.start, c.stop, c.codestart, c.codestop from ensembl_genes_9 a, ensembl_transcript_9 b, ensembl_exon_9 c where a.gene=b.gene and b.transcript=c.transcript order by a.chromosome, b.start, c.start")){
  if(!cursor.Declare("select a.chromosome, a.gene, a.id, a.strand, b.transcript, b.id, b.start, b.stop,  c.id, c.start, c.stop, c.codestart, c.codestop from ensembl_genes_9 a, ensembl_transcript_9 b, ensembl_exon_9 c where a.gene=b.gene and b.transcript=c.transcript order by a.gene, b.transcript")){
    cerr << "couldn't declare a cursor for gene structures " << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  //  if(!cursor.Fetch()){
  //  cerr << "couldn't fetch the gene structures " << endl << cursor.ErrorMessage() << endl;
  //  exit(1);
  //}
  //cout << "currentGeneId is now -1" << endl;
  int currentGeneId = -1;
  int currentTranscriptIndex = -1;  // so I know when to make a new transcript or gene thingy..
  ensemblTranscript* currentTranscript = 0;
  ensemblGene* currentGene = 0;     // use these to add things..
  int counter = 0;
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      counter++;
      //cout << "i is now " << i << endl;
      string chrom = cursor.GetValue(i, 0);
      int geneIndex = atoi(cursor.GetValue(i, 1));
      string geneId = cursor.GetValue(i, 2);
      int strand = atoi(cursor.GetValue(i, 3));
      int transcriptIndex = atoi(cursor.GetValue(i, 4));
      string transcriptId = cursor.GetValue(i, 5);
      int transcriptStart = atoi(cursor.GetValue(i, 6));
      int transcriptStop = atoi(cursor.GetValue(i, 7));
      string exonId = cursor.GetValue(i, 8);
      int exonStart = atoi(cursor.GetValue(i, 9));
      int exonStop = atoi(cursor.GetValue(i, 10));
      int exonCodeStart = atoi(cursor.GetValue(i, 11));
      int exonCodeStop = atoi(cursor.GetValue(i, 12));
      // if we have a new gene, then create it..
      if(geneIndex != currentGeneId){
	/// if the current gene is not 0, then add to the appropriate chromosomal region,, (make a function somewhere).. 
	//cout << "making new gene" << endl;
	if(currentGene){
	  currentGene->findLimits();
	  //cout << "found the limits in the current gene" << endl;
	  cait = chromosomeAnnotation.find(currentGene->chromosome);
	  //cout << "found something from chromosome annotation : " << currentGene->chromosome << endl;
	  if(cait != chromosomeAnnotation.end()){
	    int n=0;
	    genomicRegion** regions = (*cait).second->regionsCovered(currentGene->start, currentGene->stop, n);
	    //cout << "got the regions " << regions << "  covering n regions " << n << endl;
	    if(n){
	      for(int j=0; j < n; j++){
		//cout << "and n is " << n << "  and regions are " << regions[j] << endl;
		regions[j]->addEnsGene(currentGene);
	      }
	    }
	  }
	}
	//	cout << "just before making new gene " << endl;
	currentGene = new ensemblGene(geneIndex, geneId, chrom, strand);
	//cout << "new gene Index " << geneIndex << "chromosome " << chrom << "strand " << strand << endl;
	currentGeneId = geneIndex;
      }
      // check if we need to make a new transcript.. -- actually the way I'm doing things now is not great, -- as a transcript could
      // potentially be split into several genes, as I am ordering by transcript start position. hmm. 
      if(transcriptIndex != currentTranscriptIndex){
	//cout << "making new transcript transcriptIndex is :  " << transcriptIndex  << endl;
	currentTranscript = new ensemblTranscript(transcriptIndex, transcriptId, chrom, transcriptStart, transcriptStop, strand);
	currentTranscriptIndex = transcriptIndex;
	// and add to the current gene ..
	//cout << "and adding the new transcript to the currentGene " << endl;
	currentGene->addTranscript(currentTranscript);     // that is nice.. ;->
      }
      //      cout << "adding new exon to currentTranscript : " << currentTranscript << endl;
      //cout << " counter is                          : " << counter << endl;
      currentTranscript->addExon(new ensemblExon(exonId, chrom, exonStart, exonStop, exonCodeStart, exonCodeStop, strand));
      // and that is actually all that we have to do..
    }
  }
  cursor.Close();
  cout << "end of getting genes and transcripts.. " << endl;
  // I have a feeling that the above procedure will fail to add the last gene to the thingy.. hmm. bugger..
  // let's select the fantom transcripts from the database.. Shouldn't be too much of a problem..
  if(!cursor.Declare("select a.id, a.length, b.assembly, b.chromosome, b.strand, c.cstart, c.cend, c.fstart, c.fend from fantom_transcripts a, fantom_assemblies b, fantom_matches c where a.transcript=b.transcript and b.assembly=c.assembly order by a.transcript desc, b.rank, c.fstart")){
    cerr << "Couldn't select the fantom matches into the database, bugger, dooo suru, aifuru.. " << endl;
    exit(1);
  }
  cout << "looking for fantom matches .. " << endl;
  int currentAssembly = -1;
  Transcript* currTranscript = 0;
  Exon* exn = 0;
  // Note this is not very efficient, but, oh, well, who cares.. 
  // we make a new transcript for every new assembly,, as these are anyway I suppose on different parts of the chromosome if more than one.
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      // and get the values.. 
      string id = cursor.GetValue(i, 0);
      int length = atoi(cursor.GetValue(i, 1));
      int assemblyId = atoi(cursor.GetValue(i, 2));
      string chromosome = cursor.GetValue(i, 3);
      int strand = atoi(cursor.GetValue(i, 4));
      int cstart = atoi(cursor.GetValue(i, 5));
      int cend = atoi(cursor.GetValue(i, 6));
      int fstart = atoi(cursor.GetValue(i, 7));
      int fend = atoi(cursor.GetValue(i, 8));
      if(currentAssembly != assemblyId){         // we need to make a new assembly..
	if(currTranscript){                   // but if we have a made transcript then we need to add this to the appropriate chromosomal annotation
	  cait = chromosomeAnnotation.find(currTranscript->chromosome);
	  if(cait != chromosomeAnnotation.end()){
	    //	    cout << "adding fantom transcript to genomic region" << endl;
	    int n=0;
	    genomicRegion** regions = (*cait).second->regionsCovered(currTranscript->start, currTranscript->stop, n);
	    for(int j=0; j < n; j++){
	      regions[j]->addTranscript(currTranscript);
	    }
	  }else{
	    cout << "no chromosome region defined for fantom match on chromosome " << currTranscript->chromosome << endl;
	  }
	}
	// make a new transcript with the above data..
	//	cout << "making new fantom transcript" << endl;
	currTranscript = new Transcript(id, "Fantom", chromosome, length, strand);
	currentAssembly = assemblyId;
      }
      currTranscript->addExon(cstart, cend, strand, fstart, fend);
    }
  }
  //// note the last transcript will never get added to the genomic region.. so we better do manually... 
  if(currTranscript){
    cait = chromosomeAnnotation.find(currTranscript->chromosome);
    if(cait != chromosomeAnnotation.end()){
      int n=0;
      genomicRegion** regions = (*cait).second->regionsCovered(currTranscript->start, currTranscript->stop, n);
      for(int j=0; j < n; j++){
	regions[j]->addTranscript(currTranscript);
      }
    }
  } 
  cursor.Close();
  QString expString;   // need for the toDouble() function.. 
  //  set<string> doneChromes;    // ones that I have done already. Not so efficient but it's ok.. 

  /// OK, use the new table and the new probeSet match struct to get the relevant details. Note that the fact that I changed the probeSetMatch struct
  /// will probably screw up a whole load of things, but that can't really be helped at the moment. 
  cout << "am about to select from blast_genome_matches_7_b " <<  endl;
  if(!cursor.Declare("select * from blast_genome_matches_7_b")){
    cerr << "couldn't select from blast_genome_matches_7_b " << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    //    cout << "selectting from blast_genome matches 7_b and the cursorTuples are " << cursor.Tuples() << endl;
    for(int i=0; i < cursor.Tuples(); i++){
      //cout << "and i is now " << i << endl;
      int afid = atoi(cursor.GetValue(i, 0));
      int afLength = atoi(cursor.GetValue(i, 1));
      string chrom = cursor.GetValue(i, 2);
      expString = cursor.GetValue(i, 3);
      double expect = expString.toDouble();
      int afStart = atoi(cursor.GetValue(i, 4));
      int afEnd = atoi(cursor.GetValue(i, 5));
      int genStart = atoi(cursor.GetValue(i, 6));
      int genEnd = atoi(cursor.GetValue(i, 7));
      int alignLength = atoi(cursor.GetValue(i, 8));
      int matches = atoi(cursor.GetValue(i, 9));
      //string afMatchSeq = cursor.GetValue(i, 10);
      //string ensemblMatchSeq = cursor.GetValue(i, 11);
      int strand = atoi(cursor.GetValue(i, 10));
      // which should be ok.. 
      // try to find the chromosome..
      cait = chromosomeAnnotation.find(chrom);
      if(cait != chromosomeAnnotation.end()){
	int n=0;
	genomicRegion** regions = (*cait).second->regionsCovered(genStart, genEnd, n);
	//cout << "n is now " << n << endl;
	//	probeSetMatch* tempMatch = new probeSetMatch(afid, afid-1, afStart, afEnd, afLength, alignLength, matches, genStart, genEnd, expect, strand, chrom, afMatchSeq, ensemblMatchSeq);
	probeSetMatch* tempMatch = new probeSetMatch(afid, afid-1, afStart, afEnd, afLength, alignLength, matches, genStart, genEnd, expect, strand, chrom);
	if(afid-1 < probeData.size()){
	  probeData[afid-1].addProbeSetMatch(tempMatch);
	}else{
	  cerr << "Somehow we got an index that is not covered by the probeData structure, better die and clean it up" << endl
	       << "The index is " << afid << " - 1 " << endl;
	  exit(1);
	}
	for(int j=0; j < n; j++){
	  regions[j]->addProbeSetMatch(tempMatch);
	}
      }else{
	cerr << "No annotation for chromosome : " << chrom << endl;
      }
    }
  }
  cursor.Close();
  // first let's get some general information about the in situ probes and fill the map<int, ishProbeData> ishProbes with generalised inforamtion..
  if(!cursor.Declare("select probe, sequence_consensus, antisense_promoter, afid, vector, design_length, probe_name, ensembl_guess, probe_id from in_situ_probes")){
    cerr << "coulnd't select from in_situ_probes dying here" << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeId = atoi(cursor.GetValue(i, 0));
      string sequence = cursor.GetValue(i, 1);
      string antiprom = cursor.GetValue(i, 2);
      int afid = atoi(cursor.GetValue(i, 3));
      string v = cursor.GetValue(i, 4);
      int dLength = atoi(cursor.GetValue(i, 5));
      string pName = cursor.GetValue(i, 6);
      int ensGuess = atoi(cursor.GetValue(i, 7));
      string probeIdentifier = cursor.GetValue(i, 8);
      ishProbes.insert(make_pair(probeId, ishProbeData(probeId, sequence, antiprom, afid, v, dLength, pName, probeIdentifier, ensGuess)));
    }
  }
  cursor.Close();    // hey, that 's good, now we can just see if we can find the little buggers below.. 
  ////  each probe set may be associated with a set of annotation or comments.. - these are represented in three different multimaps in
  ///   the ishprobe data structure.. --need a three table joing or something like that..

  /// first text based annotation.. for the ish probes.. -- note that this is not the 
  if(!cursor.Declare("select a.index, a.user_name, b.field, c.probe_id, c.annotation, c.oid from users a, ish_annotation_text_fields b, ish_probe_text_annotation c where a.index=c.user_index and c.field=b.index")){
    cerr << "Coulnd't select from ish_probe_annotation will die here" << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  map<int, ishProbeData>::iterator mit;   // = ishProbes.find(probeId)
  while(cursor.Fetch(1000)){  // bit unnecessary as we are unlikely to have so many, but you never know..
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeIndex = atoi(cursor.GetValue(i, 3));
      mit = ishProbes.find(probeIndex);
      if(mit == ishProbes.end()){
	cerr << "Somehow we come across a disrepancy in the database. No ishProbe defined for index : " << probeIndex << endl;
	continue;
      }
      int usrIndex = atoi(cursor.GetValue(i, 0));
      string usrName = cursor.GetValue(i, 1);
      string fieldName = cursor.GetValue(i, 2);
      string annotation = cursor.GetValue(i, 4);
      int annotationId = atoi(cursor.GetValue(i, 5));
      (*mit).second.textAnnotation.insert(make_pair(annotationId, ish_annotation(annotationId, usrIndex, usrName, annotation, fieldName)));
      // and put the fields into the structure..
      ishTextFields.insert(fieldName);
    }
  }
  // and close the cursor..
  cursor.Close();
  if(!cursor.Declare("select a.index, a.user_name, b.field, c.probe_id, c.annotation, c.oid from users a, ish_probe_num_fields b, ish_probe_num_annotation c where a.index=c.user_index and c.field=b.index")){
    cerr << "couldn't select from ish_probe_num_annotation " << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeIndex = atoi(cursor.GetValue(i, 3));
      mit = ishProbes.find(probeIndex);
      if(mit == ishProbes.end()){
	cerr << "Somehow we come across a disrepancy in the database. No ishProbe defined for index : " << probeIndex << endl;
	continue;
      }
      int usrIndex = atoi(cursor.GetValue(i, 0));
      string usrName = cursor.GetValue(i, 1);
      string fieldName = cursor.GetValue(i, 2);
      float value = atof(cursor.GetValue(i, 4));
      int annotationId = atoi(cursor.GetValue(i, 5));
      (*mit).second.numberAnnotation.insert(make_pair(annotationId, ish_annotation(annotationId, usrIndex, usrName, value, fieldName)));
      ishFloatFields.insert(fieldName);
    }
  }
  cursor.Close();
  // and finally the classification table. Again very similar to the above. should be not too much problem ..
  if(!cursor.Declare("select a.index, a.user_name, b.class, c.probe_id, c.confidence, c.oid from users a, ish_probe_classes b, ish_probe_classification c where a.index=c.user_index and c.class=b.index")){
    cerr << "couldn't select from ish_probe_classification " << endl << cursor.ErrorMessage() << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeIndex = atoi(cursor.GetValue(i, 3));
      mit = ishProbes.find(probeIndex);
      if(mit == ishProbes.end()){
	cerr << "Somehow we come across a disrepancy in the database. No ishProbe defined for index : " << probeIndex << endl;
	continue;
      }
      int usrIndex = atoi(cursor.GetValue(i, 0));
      string usrName = cursor.GetValue(i, 1);
      string fieldName = cursor.GetValue(i, 2);
      float value = atof(cursor.GetValue(i, 4));   // in this case this is the confidence value.. (as long as the users actually use it properly and don't add their height or iq, or whatever ..
      int annotationId = atoi(cursor.GetValue(i, 5));
      (*mit).second.classification.insert(make_pair(annotationId, ish_annotation(annotationId, usrIndex, usrName, value, fieldName)));
      ishClasses.insert(fieldName);
    }
  }
  cursor.Close();

  ///// ISH PROBE MATCHES AND ASSEMBLIES IN TWO STEPS ... 
  ///// 1. get the assemblies from the assembly table, -insert into probedata things, and add to the appropriate genomic region..
  ///// 2. get the constituent matches and add these to the ish probe data thingy, -- as the region uses a pointer, everything is OK..
  if(!cursor.Declare("select a.probe, b.chromosome, b.strand, a.design_length, b.begin, b.stop, b.assembly, b.score from in_situ_probes a, ish_probe_assemblies b where a.probe = b.probe order by a.probe, b.score desc")){
    cerr << "Couldn't select from ish_probe_assemblies, bugger let's die" << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeId = atoi(cursor.GetValue(i, 0));
      string chromosome = cursor.GetValue(i, 1);
      int strand = atoi(cursor.GetValue(i, 2));
      int pLength = atoi(cursor.GetValue(i, 3));
      int cBegin = atoi(cursor.GetValue(i, 4));
      int cEnd = atoi(cursor.GetValue(i, 5));
      int assemblyId = atoi(cursor.GetValue(i, 6));
      float score = atof(cursor.GetValue(i, 7));
      // and let's make a probeSetMatchSet, and add it to the appropriate ishprobe data thing and add a pointer
      // to the appropriate genomic location.. 
      cait = chromosomeAnnotation.find(chromosome);
      if(cait != chromosomeAnnotation.end()){
	ishProbeMatchSet* tset = new ishProbeMatchSet(probeId, chromosome, strand, pLength, cBegin, cEnd, assemblyId, score);
	map<int, ishProbeData>::iterator iit = ishProbes.find(probeId);
	if(iit != ishProbes.end()){
	  (*iit).second.probeMatches.push_back(tset);
	}else{
	  cerr << "We got an ishProbeMatchSet for what appears to be an undefined index : " << probeId << "  will die to avoid further embaressment" << endl;
	  exit(1);
	}
	// add to the appropriate regions.. 
	int n = 0;
	genomicRegion** regions = (*cait).second->regionsCovered(cBegin, cEnd, n); // n is passed by reference.. 
	// First find the appropriate ishprobedata set..
	for(int j=0; j < n; j++){
	  regions[j]->addIshMatch(tset);
	}
      }
    }
  }
  cursor.Close();
  //// so now we should have all of the ish_probe_match_set regions specified, and it's time to fill in the data for the actual 
  ///  matches.   
  if(!cursor.Declare("select a.probe, b.fstart, b.fend, b.cstart, b.cend, b.percent, a.chromosome, a.assembly from ish_probe_assemblies a, ish_probe_matches b where a.assembly=b.assembly")){
    cerr << "Couldn't select from ish_probe_matches, dieing... " << endl;
    exit(1);
  }
  while(cursor.Fetch(1000)){
    if(!cursor.Tuples()){
      break;
    }
    for(int i=0; i < cursor.Tuples(); i++){
      int probeId = atoi(cursor.GetValue(i, 0));
      int pstart = atoi(cursor.GetValue(i, 1));
      int pend = atoi(cursor.GetValue(i, 2));
      int cstart = atoi(cursor.GetValue(i, 3));
      int cend = atoi(cursor.GetValue(i, 4));
      float percent = atof(cursor.GetValue(i, 5));
      string chromosome = cursor.GetValue(i, 6);
      int asid = atoi(cursor.GetValue(i, 7));
      // find the appropriate thingy.. 
      map<int, ishProbeData>::iterator iit = ishProbes.find(probeId);
      if(iit != ishProbes.end()){
	ishProbeMatch* tmatch = new ishProbeMatch(probeId, pstart, pend, cstart, cend, percent, chromosome, asid);
	(*iit).second.addIshProbeMatch(tmatch);
      }else{
	cerr << "We got an ishProbeMatchSet for what appears to be an undefined index : " << probeId << "  will die to avoid further embaressment" << endl;
	exit(1);
      }
    }
  }
  cursor.Close();

//   if(!cursor.Declare("select probe, p_start, p_end, probe_length, align_length, matches, ensembl_start, ensembl_end, expect, strand, chromosome, p_match_seq, ensembl_match_seq from ish_probe_blast")){
//     cerr << "Couldn't select from ish_probe_blast, dying here" << endl << cursor.ErrorMessage() << endl;
//     exit(1);
//   }
//   while(cursor.Fetch(1000)){
//     if(!cursor.Tuples()){
//       break;
//     }
//     cout << "Adding ish probe blast matches no: " << cursor.Tuples() << endl;
//     for(int i=0; i < cursor.Tuples(); i++){
//       int probeId = atoi(cursor.GetValue(i, 0));
//       int probeStart = atoi(cursor.GetValue(i, 1));
//       int probeEnd = atoi(cursor.GetValue(i, 2));
//       int probeLength = atoi(cursor.GetValue(i, 3));
//       int alignLength = atoi(cursor.GetValue(i, 4));
//       int matchLength = atoi(cursor.GetValue(i, 5));
//       int ensStart = atoi(cursor.GetValue(i, 6));
//       int ensEnd = atoi(cursor.GetValue(i, 7));
//       QString expectString = cursor.GetValue(i, 8);
//       double expect = expectString.toDouble();   // should have an OK in there and check if it actually works.. but there you go..
//       int strand = atoi(cursor.GetValue(i, 9));
//       string chromosome = cursor.GetValue(i, 10);
//       string p_match_seq = cursor.GetValue(i, 11);
//       string ensembl_match_seq = cursor.GetValue(i, 12);
//       // then make a new thing and insert it into the thingy.. 
//       // first find the appropriate regions to insert into..
//       cait = chromosomeAnnotation.find(chromosome);
//       if(cait != chromosomeAnnotation.end()){
// 	int n = 0;   // number of regions..
// 	genomicRegion** regions = (*cait).second->regionsCovered(ensStart, ensEnd, n); // n is passed by reference.. 
// 	ishProbeMatch* tempMatch = new ishProbeMatch(probeId, probeStart, probeEnd, probeLength, alignLength, matchLength, ensStart, ensEnd, expect, strand, chromosome, p_match_seq, ensembl_match_seq);
// 	for(int j=0; j < n; j++){
// 	  regions[j]->addIshMatch(tempMatch);   // hooo yeah..
// 	}
// 	map<int, ishProbeData>::iterator iit = ishProbes.find(probeId);
// 	if(iit != ishProbes.end()){
// 	  (*iit).second.addIshProbeMatch(tempMatch);
// 	}else{
// 	  cerr << "We got an ishProbeMatch for what appears to be an undefined index : " << probeId << "  will die to avoid further embaressment" << endl;
// 	  exit(1);
// 	}
	
//       }
//     }
//   }
//   cursor.Close();
  /// and that should be it. we should probably just comment out all of the following for now, and then later if everything works we can remove it completely.
  /// otherwise, I think we are done for now.. !!! .. 



 //  if(!cursor.Declare("select afid, af_start, af_end, af_length, align_length, af_n_count, match, ensembl_start, ensembl_end, expect, chromosome from ensembl_genome_matches order by chromosome, ensembl_start")){
//     cerr << "Couldn't declare cursor for genome alignments, in trying to make some chromosomal maps" << endl;
//     cerr << "Error message : " << cursor.ErrorMessage() << endl;
//     exit(1);
//   }
//   if(!cursor.Fetch()){
//     cerr << "Couldn't fetch the data for chrom maps, dying as usual" << endl;
//     cerr << "Error Message : " << cursor.ErrorMessage() << endl;
//     exit(1);
//   }
//   for(int i=0; i < cursor.Tuples(); i++){
//     string chr = cursor.GetValue(i, 10);
//     //cout << "i is now " << i << "  and the chromsome is " << chr << endl;
//     if(!chromosomeAnnotation.count(chr)){        // we haven't come across this chromosome before, need to make a new struct.. 
//       map<string, int>::iterator csit = chromSizes.find(chr);
//       if(csit == chromSizes.end()){    // bugger, it really shouldn't happen, and I'm not sure how to deal with it. Exit with error message
// 	cerr << "no size defined for chromsome : " << chr << endl
// 	     << "don't know what to do, will exit" << endl;
// 	exit(1);
//       }
//       //chromAnnotation tempannot(chr, (*csit).second, regionSize);
//       chromosomeAnnotation.insert(make_pair(chr, new chromAnnotation(chr, (*csit).second, regionSize)));
//     }          
//     if(currentChromosome != chr){
//       cout << "finding a new iterator.. " << endl;
//       cait = chromosomeAnnotation.find(chr);
//       currentChromosome = chr;   // ok I forgot this. !! 
//     }
//     /// so now (*cait).second is the chromAnnotation..thingy.. we could make a pointer to it but what the hell
//     /// get all the values from the table..
//     int afid = atoi(cursor.GetValue(i, 0));
//     int afstart = atoi(cursor.GetValue(i, 1));
//     int afend = atoi(cursor.GetValue(i, 2));
//     int aflength = atoi(cursor.GetValue(i, 3));
//     int align_length = atoi(cursor.GetValue(i, 4));
//     int af_n_count = atoi(cursor.GetValue(i, 5));
//     int match = atoi(cursor.GetValue(i, 6));
//     int cStart = atoi(cursor.GetValue(i, 7));
//     int cEnd = atoi(cursor.GetValue(i, 8));
//     expString = cursor.GetValue(i, 9);
//     double expect = expString.toDouble();
//     /// now we can just make a pointer to a probe Set,, and then see where we should put this one..
//     /// which we can work out from the different parameters. This sounds a bit dangerous.. !!!!!
//     probeSetMatch* tempMatch = new probeSetMatch(afid, afid-1, afstart, afend, aflength, align_length, af_n_count, match, cStart, cEnd, expect, chr);
//     // and then we just have to work out where to put it.. 
//     int startRegion = (cStart-1)/regionSize;
//     int endRegion = (cEnd-1)/regionSize;              // -1 to compensate for not counting from 0, i.e. 1 * regionSize should be in 0, not number 1. 
//     //    cout << "start Region is : " << startRegion << "   and endregion is " << endRegion << endl;
//     (*cait).second->regions[startRegion]->addProbeSetMatch(tempMatch);
//     if(endRegion != startRegion){                   // which is to say that the match overlaps two regions..
//       //cout << "adding a second time for the thingy " << endl; 
//       (*cait).second->regions[endRegion]->addProbeSetMatch(tempMatch);
//     }
//   }
//   cursor.Close();

//   // and let's get ourselves the ensembl genes and stick those into the chrom regions thingies. I'm not sure exactly how yet, buy maybe it will become
//   // obvious as I go through..
//   if(!cursor.Declare("select * from ensembl_chrom_genes order by chromosome, start")){
//     cerr << "Couldn't declare the cursor for the chrom_genes table, by bye" << endl;
//     cerr << cursor.ErrorMessage() << endl;
//     exit(1);
//   }
//   if(!cursor.Fetch()){
//     cerr << "couldn't fetch the data from chrom_genes table bybyb" << endl;
//     cerr << cursor.ErrorMessage() << endl;
//     exit(1);
//   }
//   for(int i=0; i < cursor.Tuples(); i++){
//     int gIndex = atoi(cursor.GetValue(i, 0));
//     string ensemblId = cursor.GetValue(i, 1);
//     string extId = cursor.GetValue(i, 2);
//     string description = cursor.GetValue(i, 3);
//     string chrom = cursor.GetValue(i, 4);
//     int start = atoi(cursor.GetValue(i, 5));
//     int stop = atoi(cursor.GetValue(i, 6));
//     int strand = atoi(cursor.GetValue(i, 7));
//     cait = chromosomeAnnotation.find(chrom);
//     if(cait == chromosomeAnnotation.end()){
//       cerr << "ensembl Gene creation,, no such chromosome found : " << chrom << endl;
//       continue;
//     }
//     if(start > 0 && stop > 0 && start < (*cait).second->chromosomeSize && stop < (*cait).second->chromosomeSize){
//       int startNo = start/(*cait).second->regionSize;
//       int stopNo = start/(*cait).second->regionSize;
//       /// and just make sure that both startNo and stopNo are smaller than regionNo..
//       if(startNo < (*cait).second->regionNo && stopNo < (*cait).second->regionNo){
// 	ensemblGene* tGene = new ensemblGene(gIndex, ensemblId, extId, description, chrom, start, stop, strand);
// 	(*cait).second->regions[startNo]->addEnsGene(tGene);
// 	cout << "inserting gene : " << gIndex << "\tdescription" << endl;
// 	if(startNo != stopNo){
// 	  (*cait).second->regions[stopNo]->addEnsGene(tGene);
// 	}
//       }
//     }
//   }
  /// and that should be enough.. 
}

void ProbeSetSet2::guessGenes(){
  // go through all of the probe set matches within the probe set data things and try to attach a best guess gene for each one..
  ofstream out("ensemblBlastGuesses");
  for(int i=0; i < probeData.size(); i++){
    //    cout << "guessing for probeData " << i << endl;
    int bestGuessIndex = 0;
    vector<int> downstreamDistances(probeData[i].probeSetMatches.size());
    vector<int> upstreamDistances(probeData[i].probeSetMatches.size());
    vector<ensemblGene*> downstreamCandidates(probeData[i].probeSetMatches.size());
    vector<ensemblGene*> upstreamCandidates(probeData[i].probeSetMatches.size());
    vector<float> downstreamScores(probeData[i].probeSetMatches.size());
    vector<float> upstreamScores(probeData[i].probeSetMatches.size());
    bool print = false;
    //bool print = (i == 20476 || i == 17953 || i == 21113);
    cout << "index : " << i << "   probe set matches size : " << probeData[i].probeSetMatches.size() << endl;
    for(int j=0; j < probeData[i].probeSetMatches.size(); j++){
      // get an iterator for the chromAnnotation..
      map<string, chromAnnotation*>::iterator it = chromosomeAnnotation.find(probeData[i].probeSetMatches[j].chromosome);
      if(it == chromosomeAnnotation.end()){
	continue;
      }
      
      // get the genomic regions covered by the current probeSetMatch.. 
      int n = 0;
      ensemblGene* upstream = 0;
      ensemblGene* downstream = 0;    // should be the two neigbouring genes.. if we find any, if not, then no guess I'm afraid... 
      
      int upstreamSmallest = 100000000;  // 100 million should be big enough..
      int downstreamSmallest = 100000000;  // 100 million should be big enough..
      bool overlap = false;
      double maxExpect = 1e-5;   /// minimum expectation product and minMatchSum length..           
      int minMatchLength = 38;   // --38, because I seem to remember that the smallest probe set might be something like 42 b.. 
      if(print){ cout << "guessing for : " << probeData[i].afid << endl; }
      if(probeData[i].probeSetMatches[j].matchSum < minMatchLength || probeData[i].probeSetMatches[j].expectProduct > maxExpect){
	if(print){ 
	  cout << "but ignoring as one of the parameters is not good enough" << endl
	       << "matchSum : " << probeData[i].probeSetMatches[j].matchSum << endl 
	       << "expectProduct : " << probeData[i].probeSetMatches[j].expectProduct << endl;
	}
	continue;
      }
      //// in order to not include imperfect repeats, go through all of the probe set matches in the set, and select those that
      //// either do not overlap (in the affymetrix positions), or where there is an overlap, include only the best scoring one..
      set<probeSetMatch*> selectedMatches;
      // cout << "selected Matches size : " << selectedMatches.size() << endl;
      int matchSum = 0;
      int mismatchSum = 0;
      //      int minPos = (*probeData[i].probeSetMatches[j].matches.begin())->cStart;
      //      int maxPos = (*probeData[i].probeSetMatches[j].matches.begin())->cEnd;
      set<probeSetMatch*>::iterator oit;   // outer iterator
      set<probeSetMatch*>::iterator iit;   // inner iterator.
      for(oit = probeData[i].probeSetMatches[j].matches.begin(); oit != probeData[i].probeSetMatches[j].matches.end(); oit++){
	bool keep = true;
	for(iit = probeData[i].probeSetMatches[j].matches.begin(); iit != probeData[i].probeSetMatches[j].matches.end(); iit++){
	  // do we have overlap, if we do,, then check if the outer one is worse, if it's worse than set keep to false.. 
	  //if(*oit == *iit){
	  //  continue;
	  //}    -- don't actually need to check, if they are the same, then the match is not smaller and it is kept anyway.. 
	  if( ((*oit)->afStart < (*iit)->afEnd) != ((*oit)->afEnd < (*iit)->afStart) ){
	    // we have overlap,, see if the quality of this one is better.. judge quality just by the number of matches... 
	    if((*oit)->match < (*iit)->match){
	      if(print){
		cout << "I would have been setting keep to false here, but I'm not bothering " << endl;
		cout << "oit afStart : " << (*oit)->afStart << "\toit afEnd " << (*oit)->afEnd << endl
		     << "iit afStart : " << (*iit)->afStart << "\tiit afEnd " << (*iit)->afEnd << endl
		     << "oit match   : " << (*oit)->match << "\tiit match " << (*iit)->match << endl;
	      }
	      keep = false;
	    }
	  }
	}
	if(keep){
	  selectedMatches.insert(*oit);
	}
      }
      /// and now determine minPos, maxPos, matchSum, misMatchSum for this set of alignments...
      /// and see how it goes...
      // make sure we keep some!!
      if(!selectedMatches.size()){
	if(probeData[i].probeSetMatches[j].matches.size()){
	  // something seriously wrong...
	  cerr << "selected Matches is 0, but is shouldn't be as we have some matches,, .. exit " << endl;
	  exit(1);
	}
	cerr << "selectedMatches Size is 0, suggests we may have some sort of stupid problem,, contine I suppose" << endl;
	continue;
      }
      int minPos = (*selectedMatches.begin())->cStart;
      int maxPos = (*selectedMatches.begin())->cEnd;
      for(oit = selectedMatches.begin(); oit != selectedMatches.end(); oit++){
	if((*oit)->cStart < minPos){ minPos = (*oit)->cStart; }
	if((*oit)->cEnd > maxPos){ maxPos = (*oit)->cEnd; }
	matchSum += (*oit)->match;
	mismatchSum += ((*oit)->alignLength - (*oit)->match);
      }
      if(print){
	cout << "maxPos :" << maxPos << "\toriginal maxpos : " << probeData[i].probeSetMatches[j].maxPos << endl
	     << "minPos : " << minPos << "\toriginal minpos : " << probeData[i].probeSetMatches[j].minPos << endl
	     << "matchSum: " << matchSum << "\toriginal matchsum: " << probeData[i].probeSetMatches[j].matchSum << endl
	     << "mismatchSum: " << mismatchSum << "\toriginal : " << probeData[i].probeSetMatches[j].mismatchSum << endl;
      }
      // and substitute these values below.. 
      genomicRegion** regions = (*it).second->regionsCovered(minPos-100000, maxPos+100000, n);
      //      genomicRegion** regions = (*it).second->regionsCovered(probeData[i].probeSetMatches[j].minPos-100000, probeData[i].probeSetMatches[j].maxPos+100000, n);
      // iterate through the ensemblGenes in the region and collect suitable candidates..
      for(int k=0; k < n; k++){
	for(int m=0; m < regions[k]->ensGeneNo; m++){
	  /// if the strand is correct, then collect, or at least check this one.. well, actually we only need one which is upstream, and one which is downstream
	  /// so in fact we can do something clever with it later on..
	  if(regions[k]->ensGenes[m]->strand == probeData[i].probeSetMatches[j].strand){
	    // check for overlap, if overlap, then set overlap to true and upstream and downstream to the current gene...
	    //	    if( (regions[k]->ensGenes[m]->start < probeData[i].probeSetMatches[j].maxPos) != (regions[k]->ensGenes[m]->stop < probeData[i].probeSetMatches[j].minPos)){
	    if( (regions[k]->ensGenes[m]->start < maxPos) != (regions[k]->ensGenes[m]->stop < minPos)){
	      if(print){ cout << "found an overlap with gene " << regions[k]->ensGenes[m]->dbIndex << endl; }
	      overlap = true;  //more than one gene can overlap.. unfortunately.. 
	      upstream = downstream = regions[k]->ensGenes[m];
	      upstreamSmallest = downstreamSmallest = 0;
	      //cout << "got an overlap, going to the next probeSetMatch set thingy " << endl;
	      break;   // no point continuning, as we cannot get a smaller distance than 0, and we'll just be replacing things otherwise.. 
	    }
	    /// else we have to check two things,, is it upstream, or is it downstream, and then is it closer than the closest thing..
	    if(regions[k]->ensGenes[m]->stop < minPos){   // the gene is upstream..
	      if(minPos - regions[k]->ensGenes[m]->stop < upstreamSmallest){
		upstream = regions[k]->ensGenes[m];
		upstreamSmallest = minPos - regions[k]->ensGenes[m]->stop;
		if(print){
		  cout << "upstream smallest: " << upstreamSmallest << "  for gene " << upstream->dbIndex << " ranging : " << upstream->start << " -- " << upstream->stop << endl;
		}
		//cout << "\tupstream smallest : " << upstreamSmallest << endl;
	      }
	    }else{     // so no overlap, then the ensGene must be downstream..
	      if(regions[k]->ensGenes[m]->start - maxPos < downstreamSmallest){
		downstream = regions[k]->ensGenes[m];
		downstreamSmallest = regions[k]->ensGenes[m]->start - maxPos;
		if(print){
		  cout << "downstream smallest: " << downstreamSmallest << "  for gene " << downstream->dbIndex << " ranging : " << downstream->start << " -- " << downstream->stop << endl;
		}
		//cout << "\tdownstream smallest : " << downstreamSmallest << endl;
	      }
	    }
	  }
	  if(overlap){ break; }
	}
	if(overlap){ break; }
      }
      ///// now we have either an overlap, (overlap is true, downstream and upstream are the same),,
      ///// or we have an upstream or downstream gene,, -or perhaps nothing.. for this probeSetSet..
      ///// we probably need to remember these, so that we can compare them to the other probeSetMatches later on and see how things are working
      ///// out.. 
      /// but now let's just see if the thing compiles, and doesn't crash... 
      downstreamDistances[j] = downstreamSmallest;
      upstreamDistances[j] = upstreamSmallest;
      downstreamCandidates[j] = downstream;
      upstreamCandidates[j] = upstream;
      int m = matchSum - mismatchSum;
      float alignmentScore = (float)(m*m)/(float)probeData[i].probeSetMatches[j].afLength;
      if(upstream){
	upstreamScores[j] = alignmentScore / (float)(upstreamSmallest + probeData[i].probeSetMatches[j].afLength);     // (the gene is upstream, so the match is downstream
      }else{
	upstreamScores[j] = 0;
      }
      if(downstream){
	downstreamScores[j] = alignmentScore / (float)(downstreamSmallest + probeData[i].probeSetMatches[j].afLength);     // (the gene is downstream, penalise ? ? for upstream or not..?
      }else{
	downstreamScores[j] = 0;
      }
      if(print){
	cout << "downstream smallest : " << downstreamSmallest << endl;
	if(downstream){
	  cout << "\tfor gene " << downstream->dbIndex << "  ranging  : " << downstream->start << " -- " << downstream->stop << endl;
	}
	cout << "\twith score : " << downstreamScores[j] << endl;
	cout << "upstream smallest : " << upstreamSmallest << endl;
	if(upstream){
	  cout << "\tfor gene " << upstream->dbIndex << "  ranging : " << upstream->start << " -- " << upstream->stop << endl;
	}
	cout << "\twith score : " << upstreamScores[j] << endl;
      }
    }
    //cout << "and we got past to here " << endl;
    /////// and now, we'll need to go through these guys, and choose one... I say, I say, usually it should be the first one as this is the one with the lowest
    //////  the scores are already calculated, so just go through each one of them, and find the biggest scoring one.. 
    /////   though I should really sort them. But I will have to declare another struct and stuff and then a comparison function.. 
    float maxScore = 0;
    int maxIndex = 0;
    int direction = 0;   // (-1 or +1 when defined. ).
    //cout << "and then here " << endl;
    for(int j=0; j < upstreamScores.size(); j++){
      if(upstreamScores[j] > maxScore){
	maxScore = upstreamScores[j];
	maxIndex = j;
	direction = 1;
      }
      if(downstreamScores[j] > maxScore){
	maxScore = downstreamScores[j];
	maxIndex = j;
	direction = -1;
      }
    }
    //cout << "but what about here, I think the crash is coming. " << endl;
    /// and now we have an index, and a direction, but actually it doesn't matter.. upstream or downstream, let's just find out the index, and tell the thingy that..
    int dbIndex = 0;
    if(upstreamCandidates.size()){
      if(direction == 1){
	if(upstreamCandidates[maxIndex]){
	  dbIndex = upstreamCandidates[maxIndex]->dbIndex;
	}
      }else{
	if(downstreamCandidates[maxIndex]){
	  dbIndex = downstreamCandidates[maxIndex]->dbIndex;
	}
      }
    }
    //cout << "I don't think I will get here " << endl;
    probeData[i].blastGuess = dbIndex;   // alleluliah.. but man that is really bad code up there.. have to make pretty one day.. 
    probeSetEnsemblIndex.insert(make_pair(i+1, dbIndex));
    ensemblProbeSetIndex.insert(make_pair(dbIndex, i+1));   // i.e. all the probe sets associated with a given ensembl gene.. 
    if(out){
      out << i+1 << "\t" << dbIndex << endl;     // somehow sometimes the thing isn't fixed.. arghhhh. bugger.. i.e. probeData[i].index is just wrong.. oh my god.. 
    }
  }
  cout << "end of guess genes.. " << endl;
}
	
  
