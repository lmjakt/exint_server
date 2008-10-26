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
#include <map>
#include <iostream>
#include <string>
#include <set>
#include "../require.h"
#include "probe_set.h"
#include <math.h>
#include <algo.h>
#include <libpq++.h>
#include <qsocket.h>
#include <qdatastream.h>
#include <qtextstream.h>

#include "../stat/stat.h"
#include "../cluster/cluster.h"
#include "../util/util.h"

using namespace std;

probe_set::probe_set(){
}

/*
probe_set::probe_set(int i){
  index = i;
  distance = 0;
  anovaScore = 0;
}

probe_set::probe_set(int i, vector< vector<float> >& p){
  index = i;
  probes = p;
  distance = 0;
  anovaScore = 0; 
  mean.resize(probes[0].size());
  exptIndex.resize(probes[0].size());
  for(int i=0; i < probes.size(); i++){
    excluded.push_back(false);
    for(int j=0; j < probes[i].size(); j++){
      mean[j] += probes[i][j];
    }
  }
  for(int i=0; i < mean.size(); i++){
    mean[i] = mean[i]/probes.size();
    exptIndex[i] = i;
    exptLookup.insert(make_pair(i, i));
  }
  e_quality = euclidean();
}

*/
probe_set::probe_set(int i, vector< vector<float> >& p, vector<int>& eindices){
  // eindices is the experimental indices, denoted by the values in p, 
  // this should map to something useful somewhere. Whatever the case is, the 
  // exptIndex[i] value will be used by member functions in comparisons between 
  // different probe Sets. This is probably the only reasonable way of doing this..
  index = i;
  probes = p;
  exptIndex = eindices;
  distance = 0;       /// ?????????
  mean.resize(probes[0].size());
  for(int i=0; i < probes.size(); i++){
    excluded.push_back(false);
    for(int j=0; j < probes[i].size(); j++){
      mean[j] += probes[i][j];
    }
  }
  for(int i=0; i < mean.size(); i++){
    mean[i] = mean[i]/probes.size();
    exptLookup.insert(make_pair(exptIndex[i], i));
  }
  e_quality = euclidean();
}

/*
int probe_set::get_index(){
  return(index);
}

void probe_set::euclidean_filter(float cut_off){
  // a simple filter. -finds the mean expression profile of the 
  // data, - then finds the euclidean distances from the mean of the 
  // the individual probes,, calculates the standard deviation, of the 
  // the distances, and marks those as falling outside of a certain number 
  // of std_deviations as being excluded.. -I may include the distance from the euclidean
  // in the file,, -I will have to learn to write .. -some functions..
  vector<float> mean(probes[0].size());   // a vector for the mean profile..
  vector<float> euclid_dist(probes.size());
  float mean_dist = 0;
  for(int i=0; i < mean.size(); i++){
    mean[i] = 0;
    for(int j=0; j < probes.size(); j++){
      mean[i] += probes[j][i];
    }
    mean[i] = mean[i]/probes.size();
  }
  // -- then get the euclidean.. values for each probe pair..
  for(int i=0; i < probes.size(); i++){
    euclid_dist[i] = 0;
    for(int j=0; j < probes[i].size(); j++){
      euclid_dist[i] += pow((probes[i][j]-mean[j]), 2);
    }
    euclid_dist[i] = euclid_dist[i]/probes[i].size();
    mean_dist += euclid_dist[i];
  }
  mean_dist = mean_dist/probes.size();
  // and lets get the standard deviation..
  float standard_dev = std_dev(euclid_dist);
  // and lets go through all the probe sets and set the flag to true if 
  // the euclidean distance is larger than the cut-off number of standard deviations
  // I have this feeling that there is something distinctly wrong with this, but what the hell
  for(int i=0; i < probes.size(); i++){
    if( ((euclid_dist[i]-mean_dist)/standard_dev) > cut_off){
      //if( (euclid_dist[i]/standard_dev) > cut_off){
      excluded[i] = true;
    }else{
      excluded[i] = false;
    }
  }
}

void probe_set::unFilter(){
  for(int i=0; i < excluded.size(); i++){
    excluded[i] = false;
  }
}

float probe_set::diff(int up, int down){
  // simple function. for a given probe set, calculates the mean and standard deviation
  // of the differences for each probe pair between to different experimental points.
  // returns  (mean diff) / (std dev of diffs)   which is essentially a t-score for 
  // a linked thingy..
  //
  
  float score = 0;
  //if(up >= probes[0].size() || down >= probes[0].size()) {   // not fool proof, but all the probes should have the same length
  //  return(score);
  //}
  if(probes.size() < 1) { return(score); }           // no point doing it if there's only one probe after all..
  // Have to find out if the points are defined or not, by going through the exptIndex so that I know if they've been 
  // defined or not.
  
  if(exptLookup.find(up) != exptLookup.end() && exptLookup.find(down) != exptLookup.end()){
    up = exptLookup[up];
    down = exptLookup[down];
  }else{
    return(score);
  }

  //  SS  =  (sum(X^2)) - ((sum(X))^2)/N 
  // and variance = SS/N-1              -- so use this to calculate the std deviation and mean.. using a single bypass
  float mean = 0;
  float sq_sum = 0;
  float sum_sq = 0;
  for(int i=0; i < probes.size(); i++){
    mean += (probes[i][up]-probes[i][down]);
    sq_sum += pow(probes[i][up]-probes[i][down], 2);
  }
  sum_sq = pow(mean, 2);
  mean = mean/probes.size();
  float SS = sq_sum - (sum_sq/probes.size());
  float std = sqrt(SS/(probes.size()-1));
  return(mean/std);
}

float probe_set::zCompare(vector<float>& values, vector<int> index){
  float distance = 0;
  // return if size is 0, or vector sizes are different..
  if(values.size() != index.size() || values.size() < 2) { return(500); }
  if(values.size() > probes[0].size()) { return(500); }  // stupid numbers if we cant compare
  int counter = 0;

  vector<float> temp(index.size());  // for keeping the temporary values in
  vector<float> tempz(index.size()); // for the temporary z -scores..
  //  vector<int> exIndex(index.size()); // for the translated experiment indices..
  for(int i=0; i < index.size(); i++){
    //cout << "i " << i << "\tindex " << index[i] << endl;
    if(exptLookup.find(index[i]) != exptLookup.end()){
      index[i] = exptLookup[index[i]];
      // cout << "index is now " << index[i] << endl; 
    }else{
      //      cout << "lookup failed for i: " << i << "  index: " << index[i] << endl;
      return(500);    // don't do the comparison, but return some stupid number..
    }
  }
  // go through the probe sets, ...
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < index.size(); j++){
      temp[j] = probes[i][index[j]];
    }
    tempz = z_score(temp);
    // and add up the score... thingy..
    for(int j=0; j < index.size(); j++){
      distance += pow((values[j]-tempz[j]), 2);
      counter ++;
    }
  }
  return(distance/(float)counter);
}

float probe_set::zMeanCompare(vector<float>& values, vector<int> index){
  // compares the values in values to the mean of the probe set values 
  // rather than using the individual experiment pattern. Normalises the 
  // selected mean values, but assumes that the comparison values have already been
  // normalised.
  float distance = 0; 
  if(values.size() != index.size() || values.size() < 2) { return(500); }
  int counter = 0;
  // find out what refers to what..
  vector<float> temp(index.size());      // for keeping the temporary mean values -- needed for z-comparison
  for(int i=0; i < index.size(); i++){
    if(exptLookup.find(index[i]) != exptLookup.end()){
      temp[i] = mean[exptLookup[index[i]]];
    }else{
      return(500);
    }
  }
  zScore(temp);      // normalises temp.. 
  for(int i=0; i < index.size(); i++){
    distance += pow((values[i]-temp[i]), 2);
    counter++;
  }
  return(distance/(float)counter);
}

*/

float probe_set::euclidean(){
  float e_distance = euclidean(*this);
  return(e_distance);
}



float probe_set::euclidean(probe_set& pset){
  // overloaded ,, -same as below, but compares all the points..
  vector<int> comp_points;
  for(int i=0; i < probes[0].size(); i++){
    comp_points.push_back(exptIndex[i]);        // i.e. the experimental point..
  }
  float e_distance = euclidean(pset, comp_points);
  return(e_distance);
}

float probe_set::euclidean(probe_set& pset, vector<int>& comp_points){
  float distance = 0;
  int counter = 0;
  int last = comp_points.size()-1;
  // use the exptLookup to make two vectors of points to be compared against each other
  // if you know what I mean...
  vector<int> thisIndex;
  vector<int> thatIndex;
  for(int i=0; i < comp_points.size(); i++){
    if(exptLookup.find(comp_points[i]) != exptLookup.end() && pset.exptLookup.find(comp_points[i]) != pset.exptLookup.end()){
      thisIndex.push_back(exptLookup[comp_points[i]]);
      thatIndex.push_back(pset.exptLookup[comp_points[i]]);
    }
  }

  // -- comp_points defines the points to be compared, must be less than thingy...
  // simple function that just compares two probe sets -probe against probe
  // on an all against all basis across a series.. cor blimey governor..
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < pset.probes.size(); j++){
      //      if(probes[i].size() > comp_points[last] && pset.probes[j].size() > comp_points[last]){
      for(int k=0; k < thisIndex.size(); k++){
	//	if(exptLookup.find(comp_points[k]) != exptLookup.end() && pset.exptLookup.find(comp_points[k]) != pset.exptLookup.end()){
	distance += (probes[i][thisIndex[k]]-pset.probes[j][thatIndex[k]]) * (probes[i][thisIndex[k]]-pset.probes[j][thatIndex[k]]);
	counter++;
	//	}
      }
      //      }
    }
  }
  if(counter < 2){
    return(500);   // ITS crap, but what the hell... 
  }
  distance = distance/counter;
  return(distance);
}

/*
float probe_set::sum(){
  // simply returns the sum of the whole set..
  // note that this is absolutely useless when looking at z-score transformations..
  float sum = 0;
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < probes.size(); j++){
      sum += probes[i][j];
    }
  }
  return(sum);
}

float probe_set::concordance(){
  // returns a score based on the similarity of the intrinsic pattern within each probe set 
  // at all time points. 
  // create a set of vectors each one containing the values of each one..
  vector< vector<float> > probe_values(probes[0].size());       // a transposed vector  probe_values[exp point][probe number]
  float sum=0;
  for(int i=0; i < probe_values.size(); i++){
    probe_values[i].resize(probes.size());
    for(int j=0; j < probes.size(); j++){
      probe_values[i][j] = probes[j][i];
      sum += fabs(probes[j][i]);
    }
    //probe_values[i] = z_score(probe_values[i]);  // z-score transform --damme, this isn't a good idea.. OK
  }
  // do an all against all comparison... 
  // or rather a pyramid one..
  float distance = 0;
  int counter = 0;
  for(int i=0; i < probe_values.size(); i++){
    for(int j=i; j < probe_values.size(); j++){
      // compare against each other..
      if(probe_values[i].size() == probe_values[j].size() ) { // or we'll get a segmentation error or worse
	for(int k=0; k < probe_values[i].size(); k++){
	  distance += pow((probe_values[i][k] - probe_values[j][k]), 2);
	  counter++;
	}
      }
    }
  }
  return(distance / ((float)counter*sum));
}

float probe_set::selSum(vector<int>& expts){
  // returns the sum of the values for all probe sets given in the vector expt. variable..
  // divided by the number of thingies..
  int counter = 0;
  float sum = 0;
  vector<int> exIndex;
  if(expts.size() < 1) { return(sum); }
  for(int i=0; i < expts.size(); i++){
    if(exptLookup.find(expts[i]) != exptLookup.end()){
      exIndex.push_back(exptLookup[expts[i]]);
    }
  }
  for(int i=0; i < exIndex.size(); i++){
    for(int j=0; j < probes.size(); j++){
	counter++;
	sum += probes[j][exIndex[i]];
    }
  }
  sum = sum/(float)counter;
  return(sum);
}

float probe_set::meanStd(){
  // returns the mean standard deviation across the experimental series
  // calculated individually for each probe pair.
  //
  float stdSum = 0;
  if(probes.size() < 1){
    return(stdSum);
  }
  for(int i=0; i < probes.size(); i++){
    stdSum += std_dev(probes[i]);
  }
  return(stdSum/(float)probes.size());
}

float probe_set::meanOverStd(){
  // as above, but divides by the mean of all of the values..
  float stdSum = 0;
  float allMean = 0;
  for(int i=0; i < probes.size(); i++){
    stdSum += std_dev(probes[i]);
    float localMean = 0;
    for(int j=0; j < probes[i].size(); j++){
      localMean += probes[i][j];
    }
    allMean += (localMean/probes[i].size());
  }
  if(stdSum > 0){
    return(allMean/stdSum);
  }
  return(0);
}
*/

float probe_set::anova(vector<int>& expts){
  // calculates the anova score across the experimental values given in the expts variable..
  // this is essentiall equal to the (variance between groups/degrees of freedom) / (variance within groups/degrees of freedom)..
  // equations and stuff obtained from http://davidmlabe.com/hyperstat/B918101
  float t_mean=0;       // the mean of all of the values..
  // groups are essentially experimental time points..
  // so get the mean.. for all ..
  int counter = 0;
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < expts.size(); j++){
      //cout << "expts : " << expts[j] << "   probes: " << i << endl;
      if(expts[j] < probes[i].size()){             // i.e. we are not overstepping ourselves..
	t_mean += probes[i][expts[j]];
	counter++;
      }
    }
  }
  t_mean = t_mean/counter;
  counter=0;          // then we can reuse it.. 
  // calculate the mean for each experimental time point... // SO THERE IS A MUCH BETTER WAY of doing this, but I'm not going to have time to check it in one hour..
  int probesUsed = 0;
  for(int j=0; j < probes.size(); j++){
    if(!excluded[j]){
      probesUsed++;
    }
  }

  vector<float> exp_means(expts.size());
  for(int i=0; i < expts.size(); i++){
    exp_means[i] = 0;
    for(int j=0; j < probes.size(); j++){
      if(!excluded[j]){
	exp_means[i] += probes[j][expts[i]];
      }
    }
    //    exp_means[i] = exp_means[i]/probes.size();
    exp_means[i] = exp_means[i]/probesUsed;
  }
  // and then.. for each experimental time point calculate the variance within the group... 
  int er_df =0;
  float er_variance = 0;
  float var = 0;
  for(int i=0; i < expts.size(); i++){
    er_variance = 0;
    for(int j=0; j < probes.size(); j++){
      if(!excluded[j]){
	var += pow((probes[j][expts[i]]-exp_means[i]), 2);
      }
    }
    //    er_variance += (var/probes.size());
    er_variance += (var/probesUsed);
    er_df += (probesUsed-1);
    //    er_df += (probes.size()-1);
  }
  // and then find the variance between the different experimental points..
  int b_df= expts.size()-1;           // -- between degrees of freedom//
  float b_variance = 0;
  for(int i=0; i < expts.size(); i++){
    b_variance += pow((exp_means[i]-t_mean), 2);
  }
  b_variance = b_variance/expts.size();
  // so I should be able to just calculate the anova scrore at this point....
  float a_score = (b_variance/(float)b_df)/(er_variance/(float)er_df);
  anovaScore = a_score;
  return(a_score);
}

/*
float probe_set::meanMaxMeanDeviation(){
  // returns mean of the maximum deviations from the mean for each probe pair
  float meanMax = 0;        // to hold the value as it gets calculated..
  for(int i = 0; i < probes.size(); i++){
    meanMax += maxMeanDeviation(probes[i]);
  }
  return(meanMax/probes.size());
}

// some network functions, provided for convenience.. 
// UNFORTUNATELY these at the moment rely on a qsocket, and a qdatastream,, and so are not 
// very portable. I could write them differently, but that would probably slow things down 
// quite quickly. Use QDataStream to avoid using atoi and atof and so forth.. 

bool probe_set::writeDataToQSocket(QSocket* socket){
  char* head = "<ProbeSet>";
  int size;     // use for various purposes.. 
  socket->writeBlock(head, 10);
  socket->writeBlock((const char*)&index, 4);
  size = exptIndex.size();
  socket->writeBlock((const char*)&size, 4);
  for(int i=0; i < exptIndex.size(); i++){
    socket->writeBlock((const char*)&exptIndex[i], 4);
  }
  size = probes.size();
  socket->writeBlock((const char*)&size, 4);
  for(int i=0; i < probes.size(); i++){
    size = probes[i].size();
    socket->writeBlock((const char*)&size, 4);
    for(int j=0; j < probes[i].size(); j++){
      socket->writeBlock((const char*)&probes[i][j], 4);
    }
  }
  // and let's write a terminator... oh lalala ... 
  char* terminator = "<ProbeSetEnd>";
  socket->writeBlock(terminator, 14);
  socket->flush();
  return(true);
}


// A network function that returns a probe set after reading the above from a network..
probe_set readProbeDataFromQSocket(QSocket* socket, bool& ok){
  int probeIndex;
  int n;    // general purpose for various reasons..
  char data[4]; // = new char[4]; // ??
  unsigned int bufSize = 4;
  probe_set pSet;
  
  socket->readBlock(data, bufSize);
  probeIndex = *(int*)data;
  //cout << "Reading data for probeIndex: \t" << probeIndex << endl;
  socket->readBlock(data, bufSize);
  n = (int)*data;
  vector<int> exIn(n);
  for(int i=0; i < exIn.size(); i++){
    socket->readBlock(data, bufSize);
    exIn[i] = *(int*)data;
  }
  socket->readBlock(data, bufSize);
  n = *(int*)data;
  vector< vector<float> > probeValues(n);
  if(socket->bytesAvailable() < 4*n*(1+exIn.size())){ //{ waitALittleBit(socket, 4*n*(1+exIn.size())); }
    cout << "should wait for 2 seconds.. " << endl;
    socket->waitForMore(2000);
    if(socket->bytesAvailable() < 4*n*(1+exIn.size())){
      ok = false;
      return(pSet);
    }
  }
  for(int i=0; i < n; i++){
    int num2;
    socket->readBlock(data, bufSize);
    num2 = *(int*)data;
    probeValues[i].resize(num2);
    for(int j=0; j < num2; j++){
      float* f; // = new float[1];
      socket->readBlock(data, bufSize);
      f = (float*)data;
      probeValues[i][j] = *f;
    }
  }
  // which should be it..
  char endToken[14];
  bufSize = 14;
  if(socket->bytesAvailable() < 14){ 
    cout << "should wait for 2 seconds.. for 13 bytes.. " << endl;
    socket->waitForMore(2000);
  }
  socket->readBlock(endToken, bufSize);
  QString endT(endToken);;
  if(endT.compare(QString("<ProbeSetEnd>")) == 0){
    pSet = probe_set(probeIndex, probeValues, exIn);
    ok = true;
  }else{
    cout << "bugger can't find a ProbeSetEnd token got something else: " << endToken << endl;
    ok = false;
  }
  return(pSet);
}
*/
// some convenience functions:
// from file -- reads in a file, and assigns the probe sets to a 
// an vector of probe_sets ..  -- I'm buggered
// returns just that -an array of vectors. Think I'm going to get memory 
// errors here..

/*
vector<probe_set> probes_from_file(string& infile, string& type){
  // open the file and read from line to line..
  // and eventually get verything..
  ifstream in(infile.c_str());
  assure(in, infile.c_str());
  vector<probe_set> probes;
  probes.reserve(15000);
  string line;
  string substring;
  //string type("pm");
  int start;
  int end;
  int index;
  int probe;
  string temp_type;
  vector< vector<float> > values;
  int last_index = 1;
  while(1){
    while(getline(in, line)){
      if(line.find("#") != line.npos){
	continue;
      }
      if(line.find(type) == line.npos){
	continue;
      }
      /// the first value is - the index (an int), the second value is the probe number
      // and the third is text that indicates if its a pm, mm, or thingy..
      start = line.find_first_of("\t", 0);
      substring = string(line, 0, start);
      index = atoi(substring.c_str());
      end = line.find_first_of("\t", start+1);
      substring = string(line, start+1, end-(start+1));
      probe = atoi(substring.c_str());    
      start = end+1;
      end = line.find_first_of("\t", start);
      temp_type = string(line, start, (end-start));
      if(index > last_index){
	probes.push_back(probe_set(last_index, values));
	values.resize(0);
	values.reserve(20);
      }
      if(temp_type == type){
	values.push_back(get_values(line));
      }
      last_index = index;
    }
    probes.push_back(probe_set(index, values));     // excludeds will be set in this method..
    break;
  }
  return(probes);
}

*/
vector<float> get_values(string& line){
  // takes a line, and extract the numbers containing the thingies..
  // first number is the index, second number is the 
  vector<float> values; // = new vector<float>;
  int start = line.find_first_of("\t", 0);
  start = line.find_first_of("\t", start+1);
  start = line.find_first_of("\t", start+1);   /// three columns of numbers... 
  int end;
  float value;
  string substring;
  while(1){
    end = line.find_first_of("\t\n", start+1);
    if(end == line.npos)
      break;
    substring = string(line, start+1, (end-(start+1)));
    value = atof(substring.c_str());
    //cout << "value: " << value << endl;
    values.push_back(value);
    start = end;
  }
  substring = string(line, start+1, (line.size()-(start+1)));
  value = atof(substring.c_str());
  values.push_back(value);
  return(values);
}

void f_write(ofstream& out, float& f){
  // writes f to out.. -assumes that f is of size 4.. -bad habit..
  out.write((unsigned char*)&f, 4);
}

void i_write(ofstream& out, int& n){
  // writes n to out .. -assumes that v is of size 4.. 
  out.write((unsigned char*)&n, 4);
  // this could probably be written as elem_type ..
}

void vf_write(ofstream& out, vector<float>& v){
  // writes the size of the vector, followed by the members of the vector
  // writes floats to out..
  int sz = v.size();
  i_write(out, sz);
  for(int i=0; i < v.size(); i++){
    out.write((unsigned char*)&v[i], 4);
  }
}

float f_read(ifstream& in){
  // returns the next float in a stream
  // if its there --assume it must be there.. 
  float* f = new float[1]; // I still don't know why, but I have to do this
  in.read(f, 4);
  return(*f);
}

int i_read(ifstream& in){
  // returns the next int in a stream..
  int* n = new int[1];
  in.read(n, 4);
  return(*n);
}

/*
void probes_to_binary_file(vector<probe_set>& p_set, string& outfile){
  // write the probe sets to a binary format file
  // first open a file
  int size;
  ofstream out(outfile.c_str(), ios::binary);
  // write the number of probe sets in the thingy..
  size = p_set.size();
  i_write(out, size);
  // go through each probe set and write out values..
  for(int i=0; i < p_set.size(); i++){
    // write out the index of the probe set
    i_write(out, p_set[i].index); //(the distance is a temp thing anyway)
    // write out the number of probes..
    size = p_set[i].probes.size();
    i_write(out, size);
    // then go through each probe pair and write out the thingies..
    for(int j=0; j < p_set[i].probes.size(); j++){
      // write out the number of experiments..
      size = p_set[i].probes[j].size();
      i_write(out, size);
      // and go through and write out each float..
      for(int k=0; k < p_set[i].probes[j].size(); k++){
	f_write(out, p_set[i].probes[j][k]);
      }
    }
  }
  // and close file..
  out.close();
}

vector<probe_set> probes_from_binary_file(string& name){
  // read probes from a binary file as written with the above function
  //  vector<int> comp_points;  // stupid, -- but for getting the distanc
  ifstream in(name.c_str());
  assure(in, name.c_str());
  int index;
  int p_num;
  int e_num;
  // get the first integer and create a probe_set file containing the thingies..
  int size = i_read(in);
  cout << "creating a vector of probe_sets with " << size << " members" << endl;
  vector<probe_set> probes(size);
  // and go throught the file and just read in everything..
  for(int i=0; i < probes.size(); i++){
    // get the probe index.. 
    index = i_read(in);
    probes[i].index = index;
    p_num = i_read(in);
    //cout << "index: " << index << "\tnumber of probes: " << p_num << endl;
    probes[i].probes.resize(p_num);
    probes[i].excluded.resize(p_num);
    for(int j=0; j < p_num; j++){
      e_num = i_read(in);
      if(probes[i].mean.size() < e_num ) { probes[i].mean.resize(e_num); }
      probes[i].excluded[j] = false;
      probes[i].probes[j].resize(e_num);                    // i.e. the number of experiments... 
      for(int k=0; k < e_num; k++){
	probes[i].probes[j][k] = f_read(in);
	probes[i].mean[k] += probes[i].probes[j][k];
      }
    }
    for(int j=0; j < e_num; j++){
      probes[i].mean[j] = (probes[i].mean[j]/p_num);
    }
    probes[i].e_quality = probes[i].euclidean();
  }
  return(probes);
}
    
probeSetSet::probeSetSet(){
  // empty constructor..
}

probeSetSet::probeSetSet(vector<string>& filenames){
  // Read in from a series of files containing nad type data.
  // parse the header and try to merge the files in a reasonable manner. This
  // has the potential to get me into lots of memory issues, as I'm still going to be making 
  // some assumptions about the structures of the datafiles. Argh..
  set<float> allExperiments;                   // all the experiment indices.. 
  vector<ifstream*> filestreams;
  fileInformation.resize(filenames.size());    // declared in header file.. 
  for(int i=0; i < filenames.size(); i++){
    filestreams.push_back(new ifstream(filenames[i].c_str()));
    assure(*filestreams[filestreams.size()-1], filenames[i].c_str());
    fileInformation[i].fileName = filenames[i];
  }
  // then go through the filestreams, and read in the header for each one,, 
  // later on do some parsing of that data... to see how it works.. maybe I'll just 
  // include the FileInfo vector in the probeSetSet struct. Seems like the easiest thing to 
  // do really..
  for(int i=0; i < filenames.size(); i++){
    int headerSize = i_read(*filestreams[i]);
    fileInformation[i].chipIndex = i_read(*filestreams[i]);
    fileInformation[i].minIndex = i_read(*filestreams[i]);
    fileInformation[i].maxIndex = i_read(*filestreams[i]);
    fileInformation[i].fileNo = i_read(*filestreams[i]);  // the number of files..
    for(int j=0; j < fileInformation[i].fileNo; j++){
      //    cout << "j : " << j << endl;
      fileInformation[i].experimentIndices.push_back(f_read(*filestreams[i]));
      allExperiments.insert(fileInformation[i].experimentIndices[fileInformation[i].experimentIndices.size()-1]);
    }
    int infoSize = i_read(*filestreams[i]);       // the size of the infostring..
    char* info = new char[infoSize];             // memory leak... ???
    (*filestreams[i]).read(info, infoSize);      // and buffer overflow??? I think its a distinct possibibility.. 
    string infotext = info;
    delete info;              // should be OK, but seems terribly wasteful, Think I might have a useful function in util.. 
    // the infotext fields are divided by tabs, so I need to split them by tabs..
    string delimit("\t");
    vector<string> terms = split_pat(infotext, delimit);
    // vector should now have :fileName description fileName description..
    // check size is OK..
    if(terms.size() != 2*fileInformation[i].fileNo){             // i.e. something is seriously screwed up..
      cout << "term.size doesn't equal 2* file No.  from probe_set.cpp, function probeSetSet(vector<string>) something " << endl;
      cout << "So we'll be dumping core now.. well probably anyway" << endl;
      cout << "terms.size : " << terms.size() << endl;
    }
    for(int j=0; j < terms.size(); j+=2){
      cout << terms[j] << "\t" << terms[j+1] << endl;
      fileInformation[i].fileNames.push_back(terms[j]);
      fileInformation[i].descriptions.push_back(terms[j+1]);
    }
    // OK this is great.. but.. but what now.. ... 
  }
  // OK. all of the information from the files have been read into a vector of 
  // of nadFileInfo (described in util.h ,, oh what a mess.. ) .. 
  // and now follows the difficult thing.. 
  // OK,, I need to find out exaclty how many different experimental points there are..
  // NOTE that I am taking the experimental point to be equivalent to the RNA source. That means that 
  // experimental replicates will tend to screw things up. Because I will then have 2 experiments with the 
  // same probe sets, and the same RNA source. I clearly need to find a way to treat replicates in a logical manner
  // but its not clear to me at the moment how to do this. 
  //
  // AS PER USUAL WHEN I CAN'T THINK OF AN ELEGANT SOLUTION.. I'M GOING TO IGNORE IT FOR NOW.
  // AND GET BACK TO IT LATER. ASSUME THAT THERE WILL NO DUPLICATES... OK LAHH.... ???
  // HOAHHH..
  
  // so go through each file given in the vector and assign the appropriated values to them...
  // note that I need a boolean vector to the thing in a bit.. 
  set<float>::iterator eit;     // to iterate through the allExperiments
  for(eit = allExperiments.begin(); eit != allExperiments.end(); eit++){
    experimentIndices.push_back(*eit);
    correctionFactors.push_back(1.0);
  }              /// And then use this for mapping experiment indices to stuff...

  for(int i=0; i < filenames.size(); i++){
    // first need to work out what my experiment indices should be like..
    vector<int> exptIndices(fileInformation[i].experimentIndices.size());  /// for creating the probe_sets 
    for(int j = 0; j < fileInformation[i].experimentIndices.size(); j++){
      for(int k = 0; k < experimentIndices.size(); k++){
	if(fileInformation[i].experimentIndices[j] == experimentIndices[k]){
	  exptIndices[j] = k;
	  continue;
	}
      }
    }
    
    int probeIndex = i_read(*filestreams[i]);
    while(probeIndex != -1){
      //   cout << probeIndex << endl;
      // need to make another probeSet constructor.. hmmmmmm 
      // that can cope with everything.
      // OK, but the probe_set constructor takes -- (int index, vector< vector<float> > p);
      int probeNo = i_read(*filestreams[i]);
      vector< vector<float> > mm(probeNo);
      vector< vector<float> > pm(probeNo);
      vector< vector<float> > de(probeNo);
      for(int j=0; j < probeNo; j++){
	int expNo = i_read(*filestreams[i]);
	mm[j].resize(expNo);
	pm[j].resize(expNo);
	de[j].resize(expNo);
	for(int k=0; k < expNo; k++){
	  mm[j][k] = f_read(*filestreams[i]);
	  pm[j][k] = f_read(*filestreams[i]);
	  de[j][k] = pm[j][k] - mm[j][k];
	}
      }
      raw.push_back(probe_set(probeIndex + fileInformation[i].minIndex -1, de, exptIndices));
      for(int j=0; j < de.size(); j++){
	zScore(de[j]);
      }
      normalised.push_back(probe_set(probeIndex + fileInformation[i].minIndex -1, de, exptIndices));
      probeIndex = i_read(*filestreams[i]);
    }
  }
}

void probeSetSet::setCorrectionValues(vector<probe_set*>& probes){
  // set correction values to 1/(value/highest mean value) of a set of 
  // supposed housekeeping genes
  vector<float> meanValues = meanMean(probes);
  if(meanValues.size() != correctionFactors.size()){
    cout << "Hey that didn't work, meanValues is different size to correctionValues.." << endl;
    cout << "mean values " << meanValues.size() << "\tcorrectionValues " << correctionFactors.size() << endl;
    return;
  }
  // find the max value..
  float maxValue = 0;
  for(int i=0; i < meanValues.size(); i++){
    if(meanValues[i] > maxValue){
      maxValue = meanValues[i];
    }
  }
  if(maxValue == 0){
    cout << "Max value is 0, this is no good, .. exiting setCorrectionValues in probe_set.cpp " << endl;
    return;
  }
  // and go through and set the correction values. Leave as 1 if mean is 0, as this 
  // indicates that the mean wasn't set for some silly reason..
  maxCorrectionFactor = 0;
  for(int i=0; i < meanValues.size(); i++){
    if(meanValues[i] > 0){
      correctionFactors[i] = 1 / (meanValues[i]/maxValue);
      if(correctionFactors[i] > maxCorrectionFactor){
	maxCorrectionFactor = correctionFactors[i];
      }
      cout << "correction Factor for " << i << "\t: " << correctionFactors[i] << endl;
    }
  }
}

probeSetSet::~probeSetSet(){
  //cout << "destroying probe set .." << endl;
}

probeSetSet::probeSetSet(string& filename){
  // DON'T USE THIS FUNCTION, IT'S GOT INCOMPLETE FUNCTIONALITY, BUT I'M KEEPING IT HERE FOR HISTORICAL PURPOSES
  // reads in a set of probes from a nad format file, 
  // and assigns these to a probeSetSet structure
  // reads in raw data, and performs a z-score normalisation.. 
  // Don't read from a normalised file.. or we'll just waste a whole load of time.. 

  // Calculate difference values, by default, but I may give some options later on..
  // for normalisation, or other things, but for now... just do the difference values 
  // and make sure that everything works..
  ifstream in(filename.c_str());
  assure(in, filename.c_str());
  int headerSize = i_read(in);     // reads an integer in..
  chipIndex  = i_read(in);
  minIndex   = i_read(in);
  maxIndex   = i_read(in);
  fileNo     = i_read(in);
  for(int i=0; i < fileNo; i++){
    experimentIndices.push_back(f_read(in));
  }
  int infoSize = i_read(in);         // the size of the string of characters..
  char* info = new char[infoSize];
  in.read(info, infoSize);
  headerInfo = info;          // seems a bit wasteful.. but hmmmm...
  // Now to read in the actual data.. woe be me if I can remember how to do this..
  int probeIndex = i_read(in);     // the index of the probe.. (it may be that this should be added to minIndex, 
  while(probeIndex != -1){         // end character is -1
    // the probe_set constructor takes -- (int index, vector< vector<float> > p);
    int probeNo = i_read(in);
    vector< vector<float> > mm(probeNo);
    vector< vector<float> > pm(probeNo);
    vector< vector<float> > de(probeNo);     // for the delta.. 
    for(int i=0; i < probeNo; i++){
      int expNo = i_read(in);
      mm[i].resize(expNo);
      pm[i].resize(expNo);
      de[i].resize(expNo); 
      for(int j=0; j < expNo; j++){
	mm[i][j] = f_read(in);
	pm[i][j] = f_read(in);
	de[i][j] = pm[i][j] - mm[i][j];
      }
    }
    raw.push_back(probe_set(probeIndex, de));       // make the raw one from de/.. 
    for(int i=0; i < de.size(); i++){
       zScore(de[i]);
    }
    normalised.push_back(probe_set(probeIndex, de));    // and the normalised one.. 
    probeIndex = i_read(in);
  }
}

struct dist_set {
  int index;
  float distance;
};
    
struct comp_pset : public binary_function<float, float, bool> {
  bool operator()(dist_set x, dist_set y) { return x.distance < y.distance; }
};

struct r_comp_pset : public binary_function<float, float, bool> {
  bool operator()(dist_set x, dist_set y) { return x.distance > y.distance; }
};

*/
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



probe_data::probe_data(){   // is this really necessary ?
  defined = false;
} 

probe_data::~probe_data(){   // and what about this??? 
}


vector<probe_data> data_from_db() {
  // open up a connection to a database backend,, -lets just use a normal ascii cursor.. 
  const char* dbname = "dbname=expression";
  PgCursor data(dbname, "portal");    // still don't know why the portal
  // check the backend..
  if( data.ConnectionBad() ) {
    cerr << "Connection to database '" << dbname << "' failed." << endl
	 << "Error returned: " << data.ErrorMessage() << endl;
    exit(1);
  }
  // lets send some commands to the database FOR THE MOMENT: limit to data for chip 1 to keep the thingy..
  //  if( !data.Declare("select a.index, a.af_id, a.id, a.suid, b.gene, b.title, a.description, c.description from p_sets a, uni_data b, tigr_annotation c where a.suid=b.index and a.af_id = c.af_id")){
  if(!data.Declare("select * from probe_data order by index")){          // table created by the above joined, and stored in db for speed of access..
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
  if(! data.Declare("select * from affy_cel_match_annot")){
    cerr << "BUGGER, couldn't select from affy_cel_match_annot" << endl;
  }
  if(! data.Fetch() ){
    cerr << "couldn't fetch the affy cel match data" << endl;
  }
  cout << "got " << data.Tuples() << " tuples for the affy_cel_match query " << endl;
  int push_back_counter = 0;
  for(int i=0; i < data.Tuples(); i++){
    // let's hope we don't have a problem with memory allocation here.. hold on to your horses
    // this will be a little difficult..
    int p_index = atoi(data.GetValue(i, 0));
    p_index--;
    string CG = data.GetValue(i, 1);
    float exp = atof(data.GetValue(i, 2));
    float m = atof(data.GetValue(i, 3));
    string SF = data.GetValue(i, 4);
    string FN = data.GetValue(i, 5);
    string GN = data.GetValue(i, 6);
    string GS = data.GetValue(i, 7);
    string ND = data.GetValue(i, 8);
    // cout << "index: " << p_index << "\tCG : " << CG << "\tSF : " << SF << endl;
    //cout << "internal index: " << p_data[p_index].index << endl;
    p_data[p_index].celeraMatches.push_back(celeraMatch(CG, exp, m, SF, FN, GN, GS, ND));  // 
    push_back_counter++;
  }
  cout << "pushed back a total of " << push_back_counter << endl;
  data.Close();

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
  return(p_data);
}

/*  
// Another convenience function.. -take a vector of a probe set and a probe set reference..
// compare all the thingies against the given index using the euclidean thingy
// and then return a sorted index with the most similar at the beginning..

vector<int> eComparer(vector<probe_set>& probes, probe_set& probe, vector<int>& comp_points){
  // create a dist set for the sort function
  vector<dist_set> p_keys(probes.size());
  // and just go through..
  for(int i=0; i < probes.size(); i++){
    probes[i].distance = probe.euclidean(probes[i], comp_points);
    p_keys[i].index = i;
    p_keys[i].distance = probes[i].distance;
  }
  sort(p_keys.begin(), p_keys.end(), comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> anovaScorer(vector<probe_set*>& probes, vector<int>& comp_points){
  //create a dist set.. for the sort function..
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    // unfortunately we have to remap the comp_points
    vector<int> ecomp;
    for(int j=0; j < comp_points.size(); j++){
      if(probes[i]->exptLookup.find(comp_points[j]) != probes[i]->exptLookup.end()){
	ecomp.push_back(probes[i]->exptLookup[comp_points[j]]);
      }
    }
    p_keys[i].distance = probes[i]->anova(ecomp);
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> anovaScorer(vector<probe_set*>& probes, vector<int>& comp_points, vector<float>& rvalues){
  // create a dist set.. for the sort function..
  // also modifies the rvalues vector so that it contains the individual anova values
  // in order to enable the drawing of distribution curves and the like.. 
  vector<dist_set> p_keys(probes.size());
  rvalues.resize(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    // unfortunately we have to remap the comp_points
    vector<int> ecomp;
    for(int j=0; j < comp_points.size(); j++){
      if(probes[i]->exptLookup.find(comp_points[j]) != probes[i]->exptLookup.end()){
	ecomp.push_back(probes[i]->exptLookup[comp_points[j]]);
      }
    }
    p_keys[i].distance = probes[i]->anova(ecomp);
    rvalues[i] = p_keys[i].distance;
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector< vector<int> > meanPosNumberDeterminer(vector<probe_set*>& probes){
  //      Very simply,, attempts to determine the number of transcripts present in a given sample
  //      based on the assumption that the mean value of a given raw value point is related to the 
  //      probability of the transcript being present. However, as we have no idea as to where to draw
  //      a reasonable threshold, and this is likely to be different depending on the sample, provide 
  //      a vector with the counts between 0 - 100 percent of the maximum value. Completely ignore values
  //      below 0,, i.e. no counting..
  
  //      FIRST find out the maximum value and divide the area between it and 0 into a 100 points..
  //      WARNING:::::
  //      THE WAY THIS IS WRITTEN, IT WILL NOT WORK WHERE THE EXPERIMENTAL NUMBERS DON'T ADD UP. -- WHICH SUCKS.. 
  //      FIX IT LATER, WANT A QUICK FIX FOR NOW..
  
  float maxValue = 0;
  int meanMaxLength = 0;      // though this will give the incorrect results if we use probe sets with 
  for(int i=0; i < probes.size(); i++){
    if(meanMaxLength < probes[i]->mean.size()){
      meanMaxLength = probes[i]->mean.size();
    }
    for(int j=0; j < probes[i]->mean.size(); j++){
      if(probes[i]->mean[j] > maxValue){
	maxValue = probes[i]->mean[j];
      }
    }
  }
  // then create our vector of thingies.. 
  vector< vector<int> > counts(100);    // go through and resize the vectors with the values being set up to 0..
  for(int i=0; i < counts.size(); i++){
    counts[i].resize(meanMaxLength);    // initialised to 0 anyway, I think.. but check later..
    for(int j=0; j < counts[i].size(); j++){
      counts[i][j] = 0;
    }
  }
  vector<float> sums(meanMaxLength);
  // then go through the values and do some math's to see how it fits together..
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < probes[i]->mean.size(); j++){
      // ahh now for the maths' find the starting point..
      int end = (int)((probes[i]->mean[j]/maxValue)* (float)(counts.size()));
      sums[j] += probes[i]->mean[j];
      //if(end == counts.size()) { end--; }
      //if(end >= 0){
      //	counts[end][j]++;
      //}
      for(int h = 0; h < end; h++){
      	counts[h][j]++;
      }
    }
  }
  for(int i=0; i < sums.size(); i++){
    cout << i << "\t" << sums[i] << endl;
  }
  return(counts);
}

vector< vector<int> > probPosNumberDeterminer(vector<probe_set*>& probes){
  // tries to determine number of expressed genes by calculating the probability of seeing 
  // the number of positive PM-MM differences for a given experimental condition
  // uses a very simple positive / negative discrimination and uses a simple non-optimised
  // binomial distribution curve.
  int divNumber = 20;     // Number of divisions to divide the distribution into..
  // first find out the longest experimental length..
  int maxELength = 0; 
  for(int i=0; i < probes.size(); i++){
    if(probes[i]->mean.size() > maxELength) { maxELength = probes[i]->mean.size(); }
  }
  // min prob = 0, max prob =1 .. calculate accumulative probabilities to get prob of seeing more than.. 
  vector< vector<int> > probCount(divNumber+1);
  for(int i=0; i < probCount.size(); i++){
    probCount[i].resize(maxELength);
  }
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < probes[i]->mean.size(); j++){     // and then go through and work out where to put stuff
      int posCount = 0;
      int totCount = 0;
      for(int k=0; k < probes[i]->probes.size(); k++){
	if(probes[i]->probes[k][j] > 0){
	  posCount++;
	}
	totCount++;
      }
      // and now calculate the cumulative probability..
      float probSum = 0;
      for(int k=0; k < posCount; k++){
	probSum += binomialProb(totCount, k, 0.5);      // horribly inefficient, but lets see how it works.
      }
      // so now the probability of seeing more than posCount-1 is simply 1-probSum..
      int interval = (int)((1-probSum)*divNumber);
      if(interval < 0){
      cout << "some problem with calculating the cumulative binomial stuff. work it out.. " << endl;
      	cout << "crashing... !" << endl;
      }
      //cout << "probSum: " << probSum << "\tinterval: " << interval << "\tpos: " << posCount << "\ttot: " << totCount << endl;
      probCount[interval][j]++;
    }
  }
  return(probCount);
}

vector<int> selSumScorer(vector<probe_set*>& probes, vector<int>& comp_points){
  // create a dist set for the sort function..
  //for(int i=0; i < comp_points.size(); i++) { cout << "sel point: " << comp_points[i] << endl; }
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->selSum(comp_points);
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> concordanceScorer(vector<probe_set>& probes){
  //create a dist set for the sort function.. 
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = i;
    p_keys[i].distance = probes[i].concordance();
  }
  sort(p_keys.begin(), p_keys.end(), comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> sumScorer(vector<probe_set*>& probes){
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->sum();
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> meanOverStdSorter(vector<probe_set*>& probes){
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->meanOverStd();
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> meanMaxMeanDeviationScorer(vector<probe_set*>& probes){
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->meanMaxMeanDeviation();
  }
  sort(p_keys.begin(), p_keys.end(), comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < indices.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> diffScorer(vector<probe_set*>& probes, int up, int down){
  // create a dist set for the sort function
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->diff(up, down);
  }
  sort(p_keys.begin(), p_keys.end(), r_comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

cluster_set kCluster(vector<probe_set*>& probes, int n){
  // simply collect all the mean expression profiles (make sure all are of the same length)
  // and then pass them to the function, and hope for the best..
  vector< vector<float> > points(probes.size());   // for keeping the means in..
  for(int i=0; i < probes.size(); i++){
    points[i] = probes[i]->mean;
  }
  cluster_set cSet(points, n);
  return cSet;
}

cluster_set kCluster(vector<probe_set*>& probes, vector<int> expts, int n){
  // as above, but only collect values from the experiments in the expts variable..
  vector< vector<float> > points(probes.size());
  for(int i=0; i < probes.size(); i++){
    points[i].resize(expts.size());
    for(int j=0; j < expts.size(); j++){
      if(expts[j] < probes[i]->mean.size()){
	points[i][j] = probes[i]->mean[expts[j]];
      }else{
	cout << "kcluster function, probe_set.cpp, experimental point out of bounds, this could cause trouble\n";
      }
    }
    zScore(points[i]);       // modifies the vector. This is actually quite important if we are taking a subset, as we are interested in the change..
  }
  cluster_set cSet(points, n);
  return cSet;
}
  

vector<int> zComparer(vector<probe_set*>& probes, vector<float>& values, vector<int>& tindex){
  // create a dist set for the sort function..
  vector<dist_set> p_keys(probes.size());
  //for(int i=0; i < tindex.size(); i++){
    // cout << "index value for " << i << "\t is " << tindex[i] << endl;
  //}
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;
    p_keys[i].distance = probes[i]->zCompare(values, tindex);
  }
  sort(p_keys.begin(), p_keys.end(), comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<int> zMeanComparer(vector<probe_set*>& probes, vector<float>& values, vector<int>& tindex){
  vector<dist_set> p_keys(probes.size());
  for(int i=0; i < probes.size(); i++){
    p_keys[i].index = probes[i]->index-1;     //hmmmmmm, this is one really bad habit that I'm going to have to sort
    p_keys[i].distance = probes[i]->zMeanCompare(values, tindex);
  }
  sort(p_keys.begin(), p_keys.end(), comp_pset());
  vector<int> indices(p_keys.size());
  for(int i=0; i < p_keys.size(); i++){
    indices[i] = p_keys[i].index;
  }
  return(indices);
}

vector<float> meanMean(vector<probe_set*>& probes){
  // just returns the mean of a set of vectors.
  // OK, I really don't have to worry too much about the performance of this,
  // so I should just go through, and go through the set of probes first to 
  // work out the appropriate length of the thingy..
  int vectorSize = 0;
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < probes[i]->exptIndex.size(); j++){
      if(vectorSize < probes[i]->exptIndex[j]){
      vectorSize = probes[i]->exptIndex[j];
      }
    }
  }
  // ok. now create..
  vectorSize++;                                  // remember that the above refers to the last index, not the ...
  vector<float> meanValues(vectorSize, 0);       // not sure about this, but hell might as well try..
  vector<int> counter(vectorSize, 0);            // 0 is a bit dangerous. just have to remember to check..
  // now go through the things and make sure that we use the appropriate things..
  for(int i=0; i < probes.size(); i++){
    for(int j=0; j < probes[i]->mean.size(); j++){
      // use the lookup table to make sure that we get the appropriate index,,,
      meanValues[probes[i]->exptIndex[j]] += probes[i]->mean[j];
      counter[j]++;
    }
  } 
  // OK, now just go through and divide by the counter, and we'll see how it goes..
  for(int i=0; i < vectorSize; i++){
    if(counter[i] > 0){
      meanValues[i] = meanValues[i]/counter[i];  // otherwise its just 0 anyway... OK..
    }
  }
  return(meanValues);
}

*/
/*
  int main(int argc, char* argv[]){
  string file(argv[1]);
  vector<probe_set> probes = probes_from_binary_file(file);
  
  
  cout << "lets see about writing the probes to a binary file format.. " << endl;
  probes_to_binary_file(probes);
  cout << "have a look for a file called binary_test " << endl;
  string binfile("binary_test");
  cout << "starting to read in the probes from the binary file " << endl;
  vector<probe_set> probes2 = probes_from_binary_file(binfile);
  cout << "finished" << endl;


  for(int i=12468; i < probes.size(); i++){
  cout << "index: " << probes[i].index << endl;
  for(int j=0; j < probes[i].probes.size(); j++){
  cout << "\t" << j;
  for(int k=0; k < probes[i].probes[j].size(); k++){
  cout << "\t" << probes[i].probes[j][k];
	//cout << "\t" << probes2[i].probes[j][k];
       }
       cout << endl;
     }
  }

  cout << "seems that things have worked ok" << endl
       << "size of probes: " << probes.size() << endl;
  // try comparing by euclidean to all others.. 
  int tester = 6007;
  string pause;
  //vector<float> distances(probes.size());
  vector<dist_set> p_keys(probes.size());
  while(1){
    cout << "well do something stupid  : ";
    cin >> pause;
    tester = atoi(pause.c_str());
    cout << "comparing " << probes[tester].index << " vs all others" << endl;
    for(int i=0; i < probes.size(); i++){
      probes[i].distance = probes[tester].euclidean(probes[i]);
      p_keys[i].index = i;
      p_keys[i].distance = probes[i].distance;
      //cout << "distance " << tester << ":\t" << distances[i] << endl;
    }
    cout << "ALL THE DISTANCE DONE, HOW DID THAT GO??? " << endl;
    sort(p_keys.begin(), p_keys.end(), comp_pset());
    for(int i=0; i < probes.size(); i++){
      cout << i << "\t" << probes[p_keys[i].index].index << "\t" << probes[p_keys[i].index].distance << endl;
    }
  }
}

*/













