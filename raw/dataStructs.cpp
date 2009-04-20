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

#include "dataStructs.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <stdlib.h>

using namespace std;

bool RegionSpecification::merge(RegionSpecification& spec){
  if(spec.chromosome != chromosome){
    return(false);
  }
  // check for overlap...
  if(begin < spec.end != end < spec.begin){   // we have overlap,, swell, this one to include both..
    if(spec.end > end){ end = spec.end; }
    if(spec.begin < begin){ begin = spec.begin; }
    return(true);
  }
  return(false);
}

exInfo::exInfo(){
  dbaseIndex = 0;
  exptGroup = 0;
  realDbaseIndex = 0;
  index = 0;
  shortName = "Null";
  description = "Null";
}


sessionInformation::sessionInformation(){
  // do nothing;
}

sessionInformation::sessionInformation(int i, int o, string t, string d, vector<string> kw, set<int> om){
  index = i;
  owner = o;
  title = t;
  description = d;
  keywords = kw;
  members = om;
}

annotationInformation::annotationInformation(){
  // do nothing/
}

annotationInformation::annotationInformation(int i, int o, string a, set<int> om){
  index = i;
  owner = o;
  annotation = a;
  members = om;
}

userInformation::userInformation(){
  // empty..
}

userInformation::userInformation(int i, string u, string fn, string ln){
  index = i;
  userName = u;
  fullName = fn;
  labName = ln;
}

probeSetMatch::probeSetMatch(){
  dbIndex = -1;
  aIndex = -1;   // not possible..
}

//probeSetMatch::probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int afn, int m, int cs, int ce, double exp, string chr){
//probeSetMatch::probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int m, int cs, int ce, double exp, int sd, string chr, string af_match, string gen_match){
probeSetMatch::probeSetMatch(int dbi, int ai, int afs, int afe, int afl, int all, int m, int cs, int ce, double exp, int sd, string chr){
  dbIndex = dbi;
  aIndex = ai;
  afStart = afs;
  afEnd = afe;
  afLength = afl;
  alignLength = all;
  //af_n_count = afn;
  match = m;
  cStart = cs;
  cEnd = ce;
  expectation = exp;
  strand = sd; 
  chromosome = chr;
  //  af_match_sequence = af_match;
  //ensembl_match_sequence = gen_match;
}

void probeSetMatchSet::insertMatch(probeSetMatch* psm){
  matches.insert(psm);
  // increment the sums.. if this region does not overlap with any of the previous matches..
  ///// NOTE that we will here be screwed by multiple matches to the same sequence in the affymetrix probe set
  ///// as may be the case if it mathces to repetitive sequences. However, I cannot at the present think of a good
  ///// way of compansating for this. All ways I can think of are either too complicated, or don't work properly.
  ///// It is also the case that perhaps if we have multiple mathces this does increase the likelihood of linkage 
  ///// to a given gene, though, if it is only small part of the probe set then we will probably be screwed.
  ///// For now, I will leave it, and see how well this simple stuff works, and if it seems to be a problem, then
  ///// I will have to use one of the more complicated things... 
  
  matchSum += psm->match;
  mismatchSum += (psm->alignLength - psm->match);  // ok.. !!
  expectProduct *= psm->expectation;

  if(psm->cStart < minPos){ minPos = psm->cStart; }
  if(psm->cEnd > maxPos){ maxPos = psm->cEnd; }
  if(psm->expectation < minExpect){ minExpect = psm->expectation; }
}


void ishProbeMatchSet::insertMatch(ishProbeMatch* ipm){
  matches.insert(ipm);

  matchSum += (ipm->pEnd - ipm->pStart);
//   mismatchSum += (ipm->alignLength - ipm->match);
//   expectProduct *= ipm->expectation;
//   if(ipm->cStart < minPos){ minPos = ipm->cStart; }
//   if(ipm->cEnd > maxPos){ maxPos = ipm->cEnd; }
//   if(ipm->expectation < minExpect){ minExpect = ipm->expectation; }
}

void Transcript::addExon(Exon* ex){
  if(!exonNo){
    start = ex->start;
    stop = ex->stop;
  }else{
    if(start > ex->start){ start = ex->start; }
    if(stop < ex->start){ stop = ex->stop; }
  }
  if(ex->strand != strand){
    cerr << "Wrong strand, ignore.. " << endl;
    delete ex;
    return;
  }
  if(exonMemSize <= exonNo){
    exonMemSize = exonMemSize * 2 + 1;    // if we start with 0, we can not just multiply by two..
    Exon** newExons = new Exon*[exonMemSize];
    // copy over old data..
    for(uint i=0; i < exonNo; i++){
      newExons[i] = exons[i];
    }
    if(exons){
      delete exons;    // 
    }
    exons = newExons;
  }
  // and finally ...
  exons[exonNo] = ex;
  exonNo++;
}

void Transcript::addExon(int st, int sp, int sd, int tb, int te){
  if(sd != strand){
    cerr << "Exon has wrong strand should be : " << strand << "   but is : " << sd << endl;
    return;
  }
  Exon* ex = new Exon(st, sp, sd, tb, te);
  addExon(ex);
}

ensemblGene::ensemblGene(){
  dbIndex = -1;
  ensemblId = "";
  //  externalId = "";
  //description = "";
  chromosome = "";
  start = 0;
  stop = 0;
  strand = 0;
}

ensemblExon::ensemblExon(){
  id = -1;
}

ensemblExon::ensemblExon(string i, string c, int st, int sp, int cst, int csp, int sd){
  id = i;
  chromosome = c;
  start = st;
  stop = sp;
  codeStart = cst;
  codeStop = csp;
  strand = sd;
}

//ensemblTranscript::ensemblTranscript(){
//  id = "";
//}

ensemblTranscript::ensemblTranscript(int tIndex, string i, string c, int st, int sp, int sd){
  index = tIndex;
  id = i;
  chromosome = c;
  start = st;
  stop = sp;
  strand = sd;
  exonNo = 0;
  exonMemSize = 10;      // perhaps a bit wasteful, but what the hell..
  exons = new ensemblExon*[exonMemSize];
}

void ensemblTranscript::addExon(ensemblExon* ex){
  //cout << "adding exon to transcript current exonNo : " << exonNo << endl
  //     << "exon Mem size is                         : " << exonMemSize << endl;
  if(exonNo < exonMemSize){
    exons[exonNo] = ex;
    exonNo++;
    return;
  }
  exonMemSize = exonMemSize * 2 + 1;    // argh.. maybe very wasteful..
  ensemblExon** tempExons = new ensemblExon*[exonMemSize];
  for(int i=0; i < exonNo; i++){
    tempExons[i] = exons[i];
  }
  tempExons[exonNo] = ex;
  exonNo++;
  delete []exons;
  exons = tempExons;
}

//ensemblGene::ensemblGene(int dbi, string ensid, string extid, string des, string chr, int sd){
ensemblGene::ensemblGene(int dbi, string ensid, string chr, int sd){
  dbIndex = dbi;
  ensemblId = ensid;
  //externalId = extid;
  //description = des;
  chromosome = chr;
  strand = sd;
  transcriptNo = 0;
  transcriptMemSize = 5;
  transcripts = new ensemblTranscript*[transcriptMemSize];
}

void ensemblGene::findLimits(){
  // go through my pointers and stuff, and the work out what my limits are..
  if(transcriptNo < 1){
    start =0;
    stop = 0;
    return;
  }
  start = transcripts[0]->start;
  stop = transcripts[0]->stop;
  for(int i=0; i < transcriptNo; i++){
    if(start > transcripts[i]->start){ start = transcripts[i]->start; }
    if(stop < transcripts[i]->stop){ stop = transcripts[i]->stop; }     // probably excessive, but can't be sure.. hmm. 
    for(int j=0; j < transcripts[i]->exonNo; j++){
      if(transcripts[i]->exons[j]->start < start){ start = transcripts[i]->exons[j]->start; }
      if(transcripts[i]->exons[j]->stop > stop){ stop = transcripts[i]->exons[j]->stop; }
    }
  }
}

void ensemblGene::addTranscript(ensemblTranscript* t){
  if(transcriptNo < transcriptMemSize){
    transcripts[transcriptNo] = t;
    transcriptNo++;
    return;
  }
  transcriptMemSize = transcriptMemSize * 2 + 1;     // + 1 to make sure that everything is OK.. 
  ensemblTranscript** tempTranscripts = new ensemblTranscript*[transcriptMemSize];
  for(int i=0; i < transcriptNo; i++){
    tempTranscripts[i] = transcripts[i];
  }
  tempTranscripts[transcriptNo] = t;
  transcriptNo++;
  delete []transcripts;
  transcripts = tempTranscripts;
}

genomicRegion::genomicRegion(){
  start = 0;
  end = 0;
  chromosome = "0";
  pMatchSize = 0;
  pMatchNo = 0;
  pMatches = 0;
  //// and the genomic regions..
  ensGeneNo = 0;
  ensGeneSize = 0;
  ensGenes = 0;
  ishMatches = 0;
  ishMatchNo = 0;
  ishMatchSize = 0;
  // and transcripts..
  transcripts = 0;
  transcriptNo = transcriptMemSize = 0;
}

genomicRegion::genomicRegion(int s, int e, string chr){
  start = s;
  end = e;
  chromosome = chr;
  pMatchNo = 0;
  pMatchSize = 1 + abs((end - start)/25000);          // generous.. 1/25000,, 
  pMatches = new probeSetMatch*[pMatchSize];      // so an array of pointers.. not too big.. 
  ///// and for the ensembl genes.. again assume something like 1/25000..
  ensGeneNo = 0;
  ensGeneSize = 1 + abs((end - start)/25000);          // generous.. 1/25000,, 
  ensGenes = new ensemblGene*[ensGeneSize];        // hoo la..  +1 so always more than 0.. 
  // and the same for the ishMatches..
  ishMatchNo = 0;
  ishMatchSize = 1 + abs((end - start)/25000);          // generous.. 1/25000,, 
  ishMatches = new ishProbeMatchSet*[ishMatchSize];
  // ok la. 
  transcriptNo = 0;
  transcriptMemSize = ensGeneSize;
  transcripts = new Transcript*[transcriptMemSize];
}

void genomicRegion::addIshMatch(ishProbeMatchSet* ishpm){
  if(ishMatchNo < ishMatchSize){
    ishMatches[ishMatchNo] = ishpm;
    ishMatchNo++;
    return;
  }
  ishMatchSize = ishMatchSize * 2 + 1;  // double plus add one in case we are starting with 0.. 
  ishProbeMatchSet** tempMatches = new ishProbeMatchSet*[ishMatchSize];
  for(int i=0; i < ishMatchNo; i++){
    tempMatches[i] = ishMatches[i];
  }
  tempMatches[ishMatchNo] = ishpm;
  ishMatchNo++;
  if(ishMatches){
    delete []ishMatches;
  }
  ishMatches = tempMatches;
}

void genomicRegion::addProbeSetMatch(probeSetMatch* psm){
  if(pMatchNo < pMatchSize){
    pMatches[pMatchNo] = psm;
    pMatchNo++;
    return;
  }
  // if not we have to make a new set of Pmatches..
  pMatchSize = pMatchSize * 2 +1;   // double the size..
  probeSetMatch** tempMatches = new probeSetMatch*[pMatchSize];
  for(int i=0; i < pMatchNo; i++){
    tempMatches[i] = pMatches[i];
  }
  tempMatches[pMatchNo] = psm;
  pMatchNo++;
  if(pMatches){           // it could be 0, in which case deleting it is probably not a good idea. 
    delete []pMatches;      // this should delete the pointers to the pointer to the object,, and should be safe. 
  }
  pMatches = tempMatches;
}

void genomicRegion::addEnsGene(ensemblGene* ensg){
  if(ensGeneNo < ensGeneSize){
    ensGenes[ensGeneNo] = ensg;
    ensGeneNo++;
    return;
  }
  // have to increase the size..
  ensGeneSize = ensGeneSize*2+1;
  ensemblGene** tempGenes = new ensemblGene*[ensGeneSize];
  for(int i=0; i < ensGeneNo; i++){
    tempGenes[i] = ensGenes[i];
  }
  tempGenes[ensGeneNo] = ensg;
  ensGeneNo++;
  if(ensGenes){
    delete ensGenes;
  }
  ensGenes = tempGenes;
}

void genomicRegion::addTranscript(Transcript* transc){
  if(transcriptNo >= transcriptMemSize){
    transcriptMemSize = transcriptMemSize * 2 + 1;
    Transcript** newTranscripts = new Transcript*[transcriptMemSize];
    for(uint i=0; i < transcriptNo; i++){
      newTranscripts[i] = transcripts[i];
    }
    if(transcripts){
      delete transcripts;
    }
    transcripts = newTranscripts;
  }
  transcripts[transcriptNo] = transc;
  transcriptNo++;
}

genomicRegion** chromAnnotation::regionsCovered(int st, int sp, int& n){
  n = 0;
  if(st > sp){ return(0); }
  //cout << "start is " << st << "  and stop is " << sp << endl;
  if(st < 1){ st = 1; }
  if(sp > chromosomeSize){ sp = chromosomeSize; }
  int startRegion = (st - 1)/regionSize;
  int endRegion = (sp-1)/regionSize;
  //cout << "startRegion " << startRegion << endl
  //     << "endRegion  " << endRegion << endl 
  //     << "regionSize " << regionSize << endl;
  if(endRegion < regionNo){
    n = endRegion-startRegion + 1;
    return(regions + startRegion);   // should be ok.. 
  }
  return(0);
}


chromAnnotation::chromAnnotation(){
  chromosome = "0";
  regions = 0;
  regionNo = 0;
  regionMemSize = 0;
  regionSize = 0;
  chromosomeSize = 0;
}

chromAnnotation::chromAnnotation(string chr, int size, int rsize){
  //cout << "chrom Annotation constructor" << endl;
  chromosome = chr;
  chromosomeSize = size;
  regionSize = rsize;
  int regionMemSize = (size/rsize)+1;   // which will have enough space.. for the sequence.. 
  //cout << "regionMemSize is " << regionMemSize << endl;
  regionNo = 0;
  //cout << "creating the new genomicRegion Array" << endl;
  regions = new genomicRegion*[regionMemSize];
  int i=0;
  int start;
  int end;
  while(i < regionMemSize-1){
    //cout << "and making more things, i is " << i << endl;
    start = i*rsize + 1;   // counting from one
    end = (i+1)*rsize;     // which actually means that the
    //cout << "calling new genomicRegion, with start : " << start << "  end, " << end << " and chrom " << chr << endl;
    regions[i] = new genomicRegion(start, end, chr);
    i++;
    regionNo++;
  }
  start = i*rsize + 1;
  //cout << " and calling for the last one size is " << size  << "  and i " << i << endl;
  regions[i] = new genomicRegion(start, size, chr);
  regionNo++;
  //  cout << "end of chromAnnotation constructor" << endl;
}


chromAnnotation::~chromAnnotation(){
  cout << "deleteing a chromAnnotation " << endl;
  for(int i=0; i < regionNo; i++){
    cout << "and i is " << i << endl;
    delete regions[i];     // remember regions[i] is a pointer itself
  }
  delete []regions;
}

void chromAnnotation::addRegion(genomicRegion* gr){
  if(regionNo < regionMemSize){
    regions[regionNo] = gr;
    regionNo++;
    return;
  }
  regionMemSize = regionMemSize * 2;
  genomicRegion** tempRegions = new genomicRegion*[regionMemSize];
  for(int i=0; i < regionNo; i++){
    tempRegions[i] = regions[i];
  }
  tempRegions[regionNo] = gr;
  regionNo++;
  if(regions){
    delete []regions;
  }
  regions = tempRegions;
}

