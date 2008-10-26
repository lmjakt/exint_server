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

#include "anovaProcessor.h"
#include <qthread.h>
#include <qmutex.h>
#include <qwaitcondition.h>
#include <math.h>
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "connectionobject.h"
#include <iostream>

AnovaProcessor::AnovaProcessor(ProbeSetSet2* p, uint st, uint sp, vector<uint>* ind, vector<dist_set>* sc){
  pSet = p;
  begin = st;
  end = sp;
  indices = ind;
  scores = sc;
  userSelectedExperiments = false;
}

AnovaProcessor::AnovaProcessor(ProbeSetSet2* p, uint st, uint sp, vector<uint>* ind, vector<dist_set>* sc, uint* exptS, uint exptSS){
  cout << "processor constructor" << endl;
  pSet = p;
  begin = st;
  end = sp;
  indices = ind;
  scores = sc;
  userSelectedExperiments = true;
  exptSelection = exptS;
  exptSelectionSize = exptSS;
  cout << "end of processor constructor" << endl;
  //  cout << "finished creating anova processor with handle " << currentThread() << endl
  //   << "begin : " << begin << "  end : " << endl;
}
AnovaProcessor::~AnovaProcessor(){
  cout << "deleting the processor thread" << endl;
}

void AnovaProcessor::run(){
  float value;
  if(begin < end){      // otherwise this doesn't make much sense..
    //cout << "anova processor : " << "  begin: " << begin << "  end: " << end << endl;
    for(int i=begin; i < end; i++){
      //      cout << "\t\tanova processor : begin: " << begin << "  end: " << end << "  i: " << i  << endl;
      //cout << "\t\tfor index : " << (*indices)[i] << endl;
      if(userSelectedExperiments){
	value = anova(pSet->data[(*indices)[i]], exptSelection, exptSelectionSize);
      }else{
	value = anova(pSet->data[(*indices)[i]], pSet->data[(*indices)[i]]->exptIndex, pSet->data[(*indices)[i]]->exptSize);
      }
      //      aMutex.lock();  // this is guaranteed not to be changed by the other thread(s). -- should not need a mutex here.
      //cout << "value is " << value << endl;
      (*scores)[i].index = pSet->data[(*indices)[i]]->index;   // though this can be 0..  
      (*scores)[i].value = value;
      
      // aMutex.unlock();
    }
  }else{
    //    cout << "BIg Screup in anovaProcessor anova function, start is bigger than stop,, " << endl;
  }
}

float** AnovaProcessor::copyProbes(float** p, int ps, int es){
  //cout << "copy Probes ps : " << ps << "  es : " << es << endl;
  if(!ps){
    return(0);
  }
  float** cp = new float*[ps];
  for(int i =0; i < ps; i++){
    //cout << "\t\t\t\t\tcopy Probes i : " << i << endl;
    cp[i] = new float[es];
    for(int j=0; j < es; j++){
      cp[i][j] = p[i][j];
    }
  }
  return(cp);
}

void AnovaProcessor::delProbes(float** p, int ps){
  for(int i=0; i < ps; i++){
    delete []p[i];
  }
  delete []p;
}

void AnovaProcessor::zScore(float* v, int s){  //s for the size.. 
  // find the mean and std deviation..
  // Note..      SS = sum of squared devieates (squared sum)
  // AND         SS = sum(X^2) - ((sum(X))^2 / N)
  // so          std = sqrt( SS / N-1)
  // which should be reasonably fast to calculate.. (rather than using the std function.. (which I could rewrite)
  float std;
  float mean;
  if(s > 1){
    float xsum = 0;
    float xsquaredSum = 0;
    for(int i=0; i < s; i++){
      xsum += v[i];
      xsquaredSum += (v[i]*v[i]);
    }
    mean = xsum/s;
    float SS = xsquaredSum - (xsum*xsum)/(float)s;
    std = sqrt(SS/(float)(s-1));
    // then go through and  modify
    for(int i=0; i < s; i++){
      v[i] = (v[i]-mean)/std;
    }
  }
}

float AnovaProcessor::anova(probe_set* ps, uint* expts, uint es){
  //  cout << " 1 ";
  if(!ps->index){
    return(-1.0);
  }
  //  cout << " 2 ";
  // calculates the anova score across th  experimental values given in the expts variable..
  // this is essentiall equal to the (variance between groups/degrees of freedom) / (variance within groups/degrees of freedom)..
  // equations and stuff obtained from http://davidmlabe.com/hyperstat/B918101
  float t_mean=0;       // the mean of all of the values..
  // groups are essentially experimental time points..
  // so get the mean.. for all ..
  int counter = 0;
  //cout << " 3 ";
  float** probes = copyProbes(ps->probes, ps->probeNo, ps->exptSize);
  if(!probes){
    return(-1.0);
  }
  //cout << "4";
  //  vector< vector<float> > probes = pset->probes; // then normalise it..
  //cout << "\t\t\t\tcalling zScore for anova with begin : " << begin << endl;
  for(int i=0; i < ps->probeNo; i++){
    zScore(probes[i], ps->exptSize);
  }
  //cout << " 5";
  //cout << "\t\t\t\tzScore done for : " << begin << endl;

  // OK slight complication.. the int* expts contain the experiment indices that we want to 
  // look at. These do not directly map to the indices of the probe set vectors, and indeed may 
  // not be represented at all. This means that we have to create a second index that is specifice
  // to this one by using the exptLookup map

  //   OK, this isn't entirely true. the other option is just to check for presence, i.e. doesn't equal -1
  //   in the int* exptLookup, but we would have to do this for each probe set, i.e. 16 times most of the time
  //   so it seems that it might be faster to make a temporary thingy, and check properly in that one.. 

  uint* eindex = new uint[es];
  int selEx = 0;      // the selected experiments.. 
  //cout << " 6";
  for(int i=0; i < es; i++){
    //if(expts[i] < ps->exptSize){
    if(expts[i] < ps->allExptNo){
      if(ps->exptLookup[expts[i]] != -1){
	eindex[selEx] = ps->exptLookup[expts[i]];
	selEx++;
      }
    }
  }
  //cout << " 7";
  if(!selEx){
    return(-1.0);
  }
  //cout << "\t\t\t\tset the eindex for " << begin << endl;
  for(int i=0; i < ps->probeNo; i++){
    for(int j=0; j < selEx; j++){
      t_mean += probes[i][eindex[j]];
      counter++;
    }
  }
  //cout << " 8";
  //cout << "\t\t\t\tworked out the mean for " << begin << endl;
  t_mean = t_mean/counter;
  counter=0;          // then we can reuse it.. 
  // calculate the mean for each experimental time point... // SO THERE IS A MUCH BETTER WAY of doing this, but I'm not going to have time to check it in one hour..
  //  int probesUsed = probes.size();  // we don't have any excluded at the moment.
  uint probesUsed = ps->probeNo;  // we don't have any excluded at the moment.

  //vector<float> exp_means(expts.size());
  float* exp_means = new float[selEx];
  for(int i=0; i < selEx; i++){
    exp_means[i] = 0;
    for(int j=0; j < ps->probeNo; j++){
      exp_means[i] += probes[j][eindex[i]];
    }
    //    exp_means[i] = exp_means[i]/probes.size();
    exp_means[i] = exp_means[i]/probesUsed;
  }
  //cout << " 9";
  //cout << "\t\t\t\tdone the exp_means for " << begin << endl;
  // and then.. for each experimental time point calculate the variance within the group... 
  int er_df =0;
  float er_variance = 0;
  float var = 0;
  for(int i=0; i < selEx; i++){
    er_variance = 0;
    for(int j=0; j < ps->probeNo; j++){
      var += ((probes[j][eindex[i]]-exp_means[i]) * (probes[j][eindex[i]]-exp_means[i]));
      //      var += pow((probes[j][eindex[i]]-exp_means[i]), 2);
    }
    //    er_variance += (var/probes.size());
    er_variance += (var/probesUsed);
    er_df += (probesUsed-1);
    //    er_df += (probes.size()-1);
  }
  // cout << " 10";
  //cout << "\t\t\t\tdone the er_variance for thread with begin : " << begin << endl;
  // and then find the variance between the different experimental points..
  int b_df= es-1;           // -- between degrees of freedom//
  //  int b_df= expts.size()-1;           // -- between degrees of freedom//
  float b_variance = 0;
  for(int i=0; i < selEx; i++){
    b_variance += ((exp_means[i]-t_mean) * (exp_means[i]-t_mean));
    //    b_variance += pow((exp_means[i]-t_mean), 2);
  }
  //cout << " 11";
  //cout << "and done something more for thread " << begin << endl;
  b_variance = b_variance/es;
  //  b_variance = b_variance/expts.size();
  // so I should be able to just calculate the anova scrore at this point....
  float a_score = (b_variance/(float)b_df)/(er_variance/(float)er_df);
  delProbes(probes, ps->probeNo);
  delete []eindex;
  delete []exp_means;
  //cout << " 12" << endl;
  //cout << "\t\t\t\treturning the anova score ..  " << begin << endl;
  return(a_score);
}


