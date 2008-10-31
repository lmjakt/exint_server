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

#include "connectionobject.h"
#include "anovaProcessor.h"
#include "euclidSortProcessor.h"
#include "kClusterProcess.h"
#include "experimentCompareProcess.h"
#include "flatExptCompare.h"
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "../raw/dataStructs.h"
#include "../netArray/netArray.h"
#include "../protocol/protocol.h"
#include "../experiment/experiment.h"
#include "../util/dataExtractor.h"
#include "../util/normaliser.h"
#include "../util/ProbeStats.h"
#include "../util/experimentTracer.h"
#include "blastClient.h"
//#include "../../stat/stat.h"
#include <qapplication.h>
#include <qthread.h>
//#include <qsocket.h>
#include <qsocketdevice.h>
#include <qhostaddress.h>
#include <qstring.h>
#include <qtextstream.h>
#include <qcstring.h>
#include <qdatastream.h>
#include <qmutex.h>
#include <qevent.h>
#include <vector>
#include <set>
#include <map>
#include <libpq++.h>
#include <unistd.h>
#include <crypt.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>      // for sort();
#include <stdlib.h>     // for rand()
#include <sstream>
#include <netinet/in.h>
#include <errno.h>
#include <signal.h>
#include <string.h>

using namespace std;


typedef unsigned int uint;   

struct comp_set : public binary_function<float, float, bool> {
  bool operator()(dist_set x, dist_set y) { return x.value < y.value; }
};

struct r_comp_set : public binary_function<float, float, bool> {
  bool operator()(dist_set x, dist_set y) { return x.value > y.value; }
};

ConnectionObject::ConnectionObject(QSocketDevice* s, ProbeSetSet2* p, QWidget* sp){
  //dataWait = QWaitCondition();
  socket = s;
  socketNumber = socket->socket();
  addressString = socket->peerAddress().toString();
  pSet = p;
  server = sp;
  commandCounter = 0;

  /// SIGNAL HANDLER /////
  sa.sa_handler = SIG_IGN;
  //sa.sa_mask = 0;
  sigemptyset(&sa.sa_mask);           // allows other signals to interrupt the signal handling or something like that. 
  sa.sa_flags = 0;
  sigaction(SIGPIPE, &sa, (struct sigaction *)NULL);        /// ???? can I do that??? 
  
  //  void (ConnectionObject::*funcPtr)(int) = &(this->catchAlarm);
  //alarmAction.sa_handler = funcPtr;
  alarmAction.sa_handler = ConnectionObject::catchAlarm;  // but this might give problems with several threads ?? 
  sigfillset(&alarmAction.sa_mask);
  alarmAction.sa_flags = 0;
  
  sigaction(SIGALRM, &alarmAction, (struct sigaction *)NULL);        /// ?? should be ok.. both just ignore.. do nothing.. 
 
  ///////////////////////////////////////////////////////
  ///// CHANGE this by using a qmutex owned by pSet -- should be more effective..  
  // But that is silly, we are still running in the parent thread.. 
  // so there is no possibility of a conflict.. and anyway, this is kind of atomic.. 
  updateUsers = false;

  //  data = d;
  bufSize = 500;
  
  authenticated=false;
  userId = 0;   // impossible id.. so in case we get some problems
  logInTime = 0; 
  currentUserName = "null";
  commandState = 0;                // reset states and messages.
  currentMessage = "";
  lastMessage = "";
  currentCommand = "";
  lastCommand = "";
  conninfo = pSet->conninfo;
  sendData = new NetArray(1000000, socketNumber, "sendBuffer");      // a lot, but I currently allow users to take 500, 000 bp of data,,  
  inData = new NetArray(10000, socketNumber, "receiveBuffer");          // only reads the amount buffered by network,,i.e. one packet, so 10, 000 should be much more than it will ever use. 
  // set up the client chips. default these to being all the chips, and then have some method for changing these..
  for(map<int, chipInfo>::iterator it = pSet->chipDescriptions.begin(); it != pSet->chipDescriptions.end(); it++){
    clientChips.insert((*it).first);   // which is the chip number..
  }
  //  sendData.resize(0);
  //sendData.reserve((pSet->data.size() * 4) + 1000);      // should be more data than we need. 
  //socketNumber = socket->socket();
  stayConnected = true;
  

  //connect(socket, SIGNAL(readyRead()), this, SLOT(readyRead()));     // read from client.. --- which perhaps could crash,, anyway, we don't really want this, so let's remove.
  ////////
  //connect(socket, SIGNAL(delayedCloseFinished()), this, SLOT(sendDeleteMe()) );
  QThread::start();
  cout << "end of Connection Object Constructor socket number: " << socketNumber << endl;
}

void ConnectionObject::run(){
  cout << "beginning of run function, stayConnected is : " << stayConnected << endl;
  while(stayConnected){
    readClient();
  }
  cout << "emitting deleteMe just before closing the eventLoop of the thread" << endl;
  // rather than emit a sequence, let's send an event.. 
  QCustomEvent* connectionEnded = new QCustomEvent((QEvent::Type)9999, (void*)socket);
  QApplication::postEvent(server, connectionEnded);
}


ConnectionObject::~ConnectionObject(){
  QTextStream ts(socket);
  logOut();
  ts << "<droppingConnection><>";
  cout << "delete thingy sent.. " << endl;
  delete sendData;
  delete inData;
  //socket->flush();  DON'T FLUSH, DON'T CLOSE AS THIS MAY FORCE THE PROGRAM TO WRITE ON A CLOSED SOCKET,
  //socket->close();  GIVING RISE TO A BROKEN PIPE AND A PROGRAM CRASH.. 
  //delete socket;
}

void ConnectionObject::checkFlags(){
  if(updateUsers){
    sendUserInformation();
    qApp->lock();
    updateUsers = false;
    qApp->unlock();
  }    // perhaps use the mutex when calling the size option.. hmm. 
  clusterMutex.lock();
  if(clusterProcesses.size()){
    clusterMutex.unlock();  // this seems very messy, but the sendKClusters process locks the mutex again. -- hmm, this maybe really not a very good idea.. .. 
    sendKClusters();  // potentially more than one..  lock the mutex again when accessing.. but maybe it is not so good. ... 
  }
  clusterMutex.unlock();  // that's a double unlock I believe.. --not good.. have to think about this.. 
  
  if(compareExperimentProcesses.size()){
    sendExperimentDistances();
  }
  if(flatExptComparers.size()){
    sendFlatExptDistances();
  }
  
  exptTracerMutex.lock();
  while(exptTracers.size()){
    set<PathTracer*>::iterator it = exptTracers.begin();
    (*it)->wait();   // make sure..
    PathTracer* temp = (*it);
    // need a function to send the data to the thingy..
    sendExperimentTrace(temp);
    exptTracers.erase(it);
    delete temp;
  }
  exptTracerMutex.unlock();

  // maybe this is a simpler model..
  set<BlastClient*>::iterator it;
  while(blastClients.size()){
    it = blastClients.begin();
    //  for(it = blastClients.begin(); it != blastClients.end(); it++){
    if((*it)->finished()){
      BlastClient* tempclient = (*it);
      blastClients.erase(it);
      // should really do something with the client.. like map it and send the information back, or something like that.. 
      delete tempclient;
    }
  }
}

vector<QString> ConnectionObject::splitString(QString& word, char delimiter){
  vector<QString> words(0);  // bad..
  int start = 0;
  int end = word.find(delimiter, start);
  while(end != -1){
    words.push_back(word.mid(start, (end-start)));
    start = end+1;
    end = word.find(delimiter, start);
  }
  return(words);
}

//void ConnectionObject::readyRead(){
  // don't read, just enter the thread..
  //cout << "received input from socket: " << socketNumber << endl;
  //dataWait.wakeOne();
  //  readClient();
//}

void ConnectionObject::readClient(){
  //  cout << "at the top of readClient" << endl;
  if(!authenticated && !authenticate()){
    stayConnected = false;
    cout << "Trying to authenticate : " << endl;
    sApp(string("<Message>"));
    sApp(string("Couldn't authenticate user name/ passwd. Don't sweat it could be my fault."));
    sApp(string("<MessageEnd>"));
    writeArray();
    cout << "Disconnecting the client" << endl;
    return;                                // don't read any more.
  }
  char c;
  int numberRead;  
  // don't read in one byte at a time from the socket as this is supposed to be slow, but rather read into the inData structure
  // and then parse from this..
  //alarm(2);      // set the alarm,, 
  while( (numberRead = inData->appendFromSocket(socketNumber, inData->capacity())) > 0 ){    // returns 0 if connection closed.. 
    //alarm(0);   // cancel the alarm.. 
    //cout << "just after numberRead sorted out,, numberRead is now : " << numberRead << endl;
    //cout << "and size of inData is now : " << inData->size() << endl;
    while(inData->hasMoreData()){          // just can I read one more char from this..
      c = *inData->readChar();
      switch(commandState){
      case 0:                                  // looking for a '<' beginning of command..
	if(c == '<'){
	  commandState = 1;
	}
	break;
      case 1:
	if(c != '>'){
	  currentCommand.append(c);
	}else{
	  commandState = 2;
	}
	break;
      case 2:
	// read the message parameters until we get to the '>' delimiter
	if(c != '>'){
	  currentMessage.append(c);
	}else{
	  commandState = 0;
	  lastCommand = currentCommand;
	  currentCommand = "";
	  lastMessage = currentMessage;
	  currentMessage = "";
	  cout << "command fully read in socket: " << socketNumber << endl; //": " << lastCommand << endl;
	  checkFlags();                     // in case something needs to be done.. not ideal as we cannot guarantee that the user sees the update if no activity.
	  parseCommand();                     // looks at lastCommand, and lastMessage and tries to do something
	  commandCounter++;                   // regardless of whether a correct command or not.. 
	  // if something goes wrong,, return..
	  if(!stayConnected){
	    return;
	  }
	}
      }
    }
    //cout << "inData has no more data " << endl;
    inData->empty();    
  }
  //////////////// if the append from Socket function is interrupted it will return an error number if this is EINTR then we should just keep reading unless
  ///////////////  so just call readClient again..
  if(numberRead == -EINTR){
    //cout << "received an EINTR check flags and other stuff see how things are going.. " << endl;
    inData->empty();  // this should not happen,, -as nothing should be put into the inData if we get an error,.. but .. 
    checkFlags();    // check if we need to send anything.. 
    return;
    //    readClient();    // recursion !!! -- and let's set an alarm to see how that goes... 
  }
  //cout << "numberRead at this point should be 0 is it? " << numberRead << endl;
  if(numberRead < 0){
    reportError("readClient encountered an Error ", -numberRead);
    //    cerr << "Error description is : " << strerror(-numberRead) << endl;
  }
  stayConnected = false;
}

void ConnectionObject::parseCommand(){
  /// need to evaluate the contents of the QString, as I can't switch on them,, 
  /// could be done more elegantly with a map of functions, which evaluate to a number
  /// followed by a switch, but this is a long way around as well. 
  //cout << "in the Parse Command Function. Command is " << lastCommand << endl;
  if(lastCommand == QString("getProbeSet")){
    cout << "calling sendProbeSet" << endl;
    sendProbeSet();
    //cout << "sendProbeSet returned" << endl;
    return;
  }
  if(lastCommand == "sendRegion"){
    cout << "Calling send Region : " << endl;
    sendRegion();
    return;
  }
  if(lastCommand == QString("euclidCompare")){
    cout << "calling doEuclidCompare" << endl;
    doEuclidCompare();
    //cout << "euclidCompare returned" << endl;
    return;
  }
  if(lastCommand == QString("dbLookup")){
    cout << "calling doDBLookup" << endl;
    doDBLookup();
    //cout << "doDBLookup returned" << endl;
    return;
  }
  if(lastCommand == "genDBLookup"){
    cout << "Calling doGenDBLookup" << endl;
    doGenDBLookup();
    return;
  }
    
  if(lastCommand == QString("closingConnection")){
    cout << "got a closingConnection call" << endl;
    stayConnected = false;
    return;
  }
  if(lastCommand == QString("getExptInformation")){
    cout << "calling writeExptInfo() " << endl;
    writeExptInfo();
    //cout << "writeExptInfo returned" << endl;
    return;
  }
  if(lastCommand == "sendChipInfo"){
    cout << "writing chip information to the client .. " << endl;
    writeChipInfo();
    return;
  }
  if(lastCommand == QString("getDBChoices")){
    cout << "calling writeDBChoices() " << endl;
    writeDBChoices();
    //cout << "writeDBChoices returned" << endl;
    return;
  }
  if(lastCommand == QString("setClientIndex")){
    cout << "calling readClientIndex" << endl;
    readClientIndex();
    //cout << "readClientIndex returned" << endl;
    return;
  }
  if(lastCommand == "setClientChips"){
    cout << "setting the client chips for stuff " << endl;
    setClientChips();
    return;
  }
  if(lastCommand == QString("doAnovaSort")){
    cout << "calling Anova Sort" << endl;
    anovaSort();
    //cout << "anovaSort returned " << endl;
    return;
  }
  if(lastCommand == "doAnovaSelectSort"){
    cout << "calling anovaSelectSort " << endl;
    anovaSelectSort();
    //cout << "anovaSelectSort returned " << endl;
    return;
  }
  if(lastCommand == QString("doRawComparison")){
    cout << "calling rawComparisonSort " << endl;
    rawComparisonSort();
    //cout << "rawComparisonSort returned " << endl;
    return;
  }
  if(lastCommand == QString("doDiff")){
    cout << "calling doDiffSort" << endl;
    doDiffSort();
    //cout << "doDiffSort returned" << endl;
    return;
  }
  if(lastCommand == QString("doMeanComparison")){
    cout << "calling meanComparisonSort" << endl;
    meanComparisonSort();
    //cout << "meanComparisonSort returned" << endl;
    return;
  }
  if(lastCommand == QString("pleaseFlush")){
    cout << "calling doFlush" << endl;
    doFlush();
    //cout << "doFlush returned " << endl;
    return;
  }
  if(lastCommand == QString("changePassword")){
    cout << "calling doChangePassword" << endl;
    doChangePassword();
    //cout << "doChangePassword returned" << endl;
    return;
  }
  if(lastCommand == "createNewUser"){
    cout << "calling doCreateNewUser" << endl;
    doCreateNewUser();
    //cout << "doCreateNewUser returned" << endl;
    return;
  }
  if(lastCommand == QString("getStatCollection")){
    cout << "calling collectStats" << endl;
    collectStats();
    //cout << "collectStats returned" << endl;
    return;
  }
  if(lastCommand == "newSession"){
    cout << "calling createSession " << endl;
    createSession();
    //cout << "createSession returned" << endl;
    return;
  }
  if(lastCommand == "newAnnotation"){
    cout << "calling newAnnotation" << endl;
    newAnnotation();
    //cout << "newAnnotation returned" << endl;
    return;
  }
  if(lastCommand == "updateSessionDescription"){
    cout << "calling updateSessionInformation " << endl;
    updateSessionInformation();
    //cout << "updateSessionInformation returned " << endl;
    return;
  }
  if(lastCommand == "updateUserAnnotation"){
    cout << "calling updateUserAnnotation " << endl;
    updateAnnotation();
    //cout << "updateAnnotation returned " << endl;
    return;
  }
  if(lastCommand == "sendSessionInfo"){
    cout << "calling sendSessionInformation " << endl;
    sendSessionInformation();
    //cout << "sendSessionInformation returned " << endl;
    return;
  }
  if(lastCommand == "getUserInformation"){
    cout << "calling sendUserInformation " << endl;
    sendUserInformation();
    //cout << "sendUserInformation returned " << endl;
    return;
  }
  if(lastCommand == "getGenomeSequence"){
    cout << "calling sendGenomeSequence " << endl;
    sendGenomeSequence(); // have to work out how to do this.. 
    //cout << "sendGenomeSequence returned " << endl;
    return;
  }
  if(lastCommand == "getProbeSetSequence"){
    cout << "calling sendProbeSetSequence " << endl;
    sendProbeSetSequence();
    //cout << "sendProbeSetSequence returned" << endl;
    return;
  }
  if(lastCommand == "getEnsemblTranscript"){
    cout << "calling sendEnsemblTranscriptSequence " << endl;
    sendEnsemblTranscriptSequence();
    //cout << "sendEnsemblTranscriptSequence returned " << endl;
    return;
  }
  if(lastCommand == "getEnsemblPeptide"){
    cout << "calling sendEnsemblPeptide " << endl;
    sendEnsemblPeptide();
    //cout << "sendEnsemblPeptide returned " << endl; 
    return;
  }
  if(lastCommand == "doKCluster"){
    cout << "calling doKCluster" << endl;
    doKCluster();
    //cout << "doKCluster returned" << endl;
    return;
  }
  if(lastCommand == "compareExperiments"){
    cout << "calling doExperimentCompare" << endl;
    doExperimentCompare();
    //cout << "doExperimentCompare returned " << endl;
    return;
  }
  if(lastCommand == "flatCompareExperiments"){
    cout << "calling doFlatExptCompare" << endl;
    doFlatExptCompare();
    //cout << "doFlatExptCompare returned" << endl;
    return;
  }
  if(lastCommand == "traceExperiments"){
    cout << "Calling traceExperiments" << endl;
    traceExperiments();
    return;
  }
  if(lastCommand == "sendGenomicRegion"){
    cout << "calling parseGenomicRegionRequest" << endl;
    parseGenomicRegionRequest();
    //cout << "parseGenomicRegionRequest returned " << endl;
    return;
  }
  if(lastCommand == "getIshThumbnails"){
    cout << "calling sendIshThumbnails " << endl;
    sendIshThumbnails();   // do the parsing and everything.. it should'nt be that tricky.. maybe.. 
    //cout << "sendIshThumbnails returned" << endl;
    return;
  }
  if(lastCommand == "getIshImage"){
    cout << "calling sendIshImage" << endl;
    sendIshImage();
    //cout << "sendIshImage returned" << endl;
    return;
  }
  if(lastCommand == "getRegionForIshProbe"){
    cout << "calling parseIshProbeReqionRequest" << endl;
    parseIshProbeRegionRequest();
    //cout << "parseIshProbeRegionRequest returned " << endl;
    return;
  }
  if(lastCommand == "associateIshProbeWithEnsemblGene"){
    cout << "calling associateIshProbeWithEnsemblGene" << endl;
    associateIshProbeWithEnsemblGene();
    //cout << "associateIshProbeWithEnsemblGene returned " << endl;
    return;
  }
  if(lastCommand == "setIshProbeName"){
    cout << "calling setIshProbeName " << endl;
    setIshProbeName();
    //cout << "setIshProbeName returned" << endl;
    return;
  }
  if(lastCommand == "insertIshTextAnnotation"){
    cout << "calling insertIshProbeTextAnnotation" << endl;
    insertIshProbeTextAnnotation();
    //cout << "insertIshProbeTextAnnotation returned" << endl;
    return;
  }
  if(lastCommand == "updateIshProbeTextAnnotation"){
    cout << "calling updateIshProbeTextAnnotation" << endl;
    updateIshProbeTextAnnotation();
    //cout << "uppdateIshProbeTextAnnotation returned " << endl;
    return;
  }
  if(lastCommand == "insertIshFloatAnnotation"){
    cout << "calling insertIshProbeFloatAnnotation " << endl;
    insertIshProbeFloatAnnotation();
    //cout << "insertIshProbeFloatAnnotation returned " << endl;
    return;
  }
  if(lastCommand == "insertIshClassification"){
    cout << "calling insertIshProbeClassification " << endl;
    insertIshProbeClassification();
    //cout << "insertIshProbeClassification returned " << endl;
    return;
  }
  if(lastCommand == "getIshTextFields"){
    cout << "calling sendIshTextFields " << endl;
    sendIshTextFields();
    //cout << "sendIshTextFields returned " << endl;
    return;
  }
  if(lastCommand == "getIshFloatFields"){
    cout << "calling sendIshFloatFields " << endl;
    sendIshFloatFields();
    //cout << "sendIshFloatFields returned " << endl;
    return;
  }
  if(lastCommand == "getIshClasses"){
    cout << "calling sendIshClasses" << endl;
    sendIshClasses();
    return;
    //cout << "sendIshClasses returned" << endl;
  }
  if(lastCommand == "getProtocols"){
    cout << "calling sendProtocols" << endl;
    sendProtocols();
    return;
  }
  if(lastCommand == "getProtocol"){
    cout << "cannling sendProtocol" << endl;
    sendProtocol();
    return;
  }
  if(lastCommand == "commitProtocolToDB"){
    cout << "Receiving new protocol : " << endl;
    receiveNewProtocol();
    return;
  }
  if(lastCommand == "commitExperimentToDB"){
    cout << "Committing experiment to DB" << endl;
    commitExperimentToDB();
    return;
  }
  if(lastCommand == "SendExperiments"){
    cout << "sending experiments from database " << endl;
    sendExperiments();
    return;
  }
  if(lastCommand == "SendTissues"){
    sendTissues();
    return;
  }
  if(lastCommand == "SendIshAnnotationFields"){
    sendIshAnnotationFields();
    return;
  }
  if(lastCommand == "MakeIshTissue"){
    makeIshTissue();
    return;
  }
  if(lastCommand == "MakeIshAnnotationField"){
    makeIshAnnotationField();
    return;
  }
  if(lastCommand == "CommitIshImageToDB"){
    commitIshImageToDB();
    return;
  }
  if(lastCommand == "doBlast"){
    cout << "calling Do Blast, let's see what happens " << endl; 
    doBlast();
    return;
  }
  if(lastCommand == "compareCellGroups"){
    cout << "comparing cell groups : " << endl;
    compareCellGroups();
    return;
  }
  if(lastCommand == QString("Received")){
    //int id = lastMessage.toInt();
    //map<int, vector<QByteArray> >::iterator it;
    //it = clientRequests.find(id);
    //if(it != clientRequests.end()){
      //      QByteArray* tArray = &(*it).second;
      //clientRequests.erase(id);
      //delete tArray;
    //}
    cout << "Got a received command" << endl;
    return;
  }
  if(lastCommand == "DevsFromMean"){           // calculate number of standard deviations from the mean.. for probe pair data..
    devsFromMean();
    return;
  }
  if(lastCommand == "ExportMeans"){
    exportMeans();             // export means of current index, or something like that using some sort of normalisation. 
    return;
  }
  cout << "Got an Unknown command am ignoring : " << lastCommand << endl;
       
}

void ConnectionObject::sendUserInformation(){
  // farily simple, just send the info in the userTable..
  // before everything else send the currentId.. so the user knows who he is..
  map<int, userInformation>::iterator it;
  sApp("<userInformation>");
  qiApp(userId);
  pSet->userMutex->lock();
  qiApp(pSet->userTable.size());
  for(it = pSet->userTable.begin(); it != pSet->userTable.end(); it++){
    qiApp((*it).second.index);
    sApp((*it).second.userName);
    sApp((*it).second.fullName);
    sApp((*it).second.labName);
  }
  pSet->userMutex->unlock();
  sApp("<userInformationEnd>");
  writeArray();
  //  doFlush();
  //  doFlush();
  //doFlush();
}


void ConnectionObject::sendProbeSet(){
  // get the index from the message part, and then just use the function in the 
  // probe_set.cpp file..
  bool ok;
  int regionSize = 200000;         // the amount of genomic region to send to the client .. 
  //  double maxExpect = 1e-20;        // and the largest random expectation to send..
  vector<QString> words = splitString(lastMessage, '|');
  //cout << "\tsendProbeSet : words size is " << words.size() << endl;
  if(words.size() < 1){
    cerr << "\tsendProbeSet : words size is smaller than 1 " << endl;
    return;
  }
  if(words.size() > 1){
    uint rSize = words[1].toInt(&ok);
    if(ok){
      regionSize = rSize;   // if not ok, just leave it at 200000,, ok.. 
    }
  }
  //cout << "socket number " << socketNumber << "  SendProbeSet : regionSize is : " << regionSize << endl;
  uint i = words[0].toInt(&ok);
  //cout << "\t\tSEND PROBE SET INDEX IS : (before decrement) " << i << endl;
  if(ok){
    i--;                         // db counts from 1, index counts from 0 
    if(i < pSet->data.size()){
      cout << "\tcalling writeProbeSet i: " << i << endl;
      //if(pSet->data[i]->index){
      writeProbeSet(pSet->data[i]);
      //}
      //cout << "\twriteProbeSet returned" << endl;
      //pSet->normalised[i].writeDataToQSocket(socket);
    }else{
      cerr << "i is too large no probe Set defined i: " << i << endl;
      return;     
    }
    if(i < pSet->probeData.size()){
      /// and go through the blastGenomeMatches and check if we want to send some genomic regions to the client..
      // we will set a maximum expect, -this should be user changeable, and we can use some other values as well, later on, but for 
      // now we'll go with the minExpect.. 
      //cout << "calling writeProbeData i: " << i << endl;
      writeProbeData(&pSet->probeData[i]);
      //cout << "writeProbeData returned" << endl;
      double maxExpect = 1e-6;
      // send a notification..
      sApp("<sendingRegions>");
      qiApp(1);      // the target -tells the client that these regions are associated with a probe set data rather than an in situ probe, or.. 
      qiApp(i+1);    // the index, but maybe we should rewrite that to read the index from the probe data structure.. 
      sApp("<sendingRegionsEnd>");
      writeArray();
      //cout << "Size of probeData matches is : " << pSet->probeData[i].probeSetMatches.size() << endl;
      for(uint j=0; j < pSet->probeData[i].probeSetMatches.size(); j++){
	if(pSet->probeData[i].probeSetMatches[j].minExpect < maxExpect){
	  int range = pSet->probeData[i].probeSetMatches[j].maxPos - pSet->probeData[i].probeSetMatches[j].minPos;
	  int margin = 100;
	  if(range < regionSize){
	    margin = (regionSize-range)/2;
	  }
	  // and pass the command to send the genomic region..
	  //  cout << "socketNumber : " << socketNumber <<  "  calling sendGenomicRegion i: " << i << "  j: " << j << endl
	  //     << "margin is                  : " << margin << endl;
	  sendGenomicRegion(pSet->probeData[i].probeSetMatches[j].chromosome, pSet->probeData[i].probeSetMatches[j].minPos-margin, pSet->probeData[i].probeSetMatches[j].maxPos+margin, i+1, 1);
	  //cout << "sendGenomic Region returned : " << endl;
	}
      }
      sApp("<finishedSendingRegions>");
      qiApp(1);      // the target
      qiApp(i+1);    // the index
      sApp("<finishedSendingRegionsEnd>");
      writeArray();
    }else{
      cerr << "can't find the probeData struct you fool " << endl;
    }
  }else{
    cerr << "couldn't get an integer from the lastMessage can't send probe set" << endl;
  }
  //  cout << "end of sendProbeSet " << endl;
}

void ConnectionObject::sendRegion(){
  /// ok, this one should include the paramters to use,, surely.. bit extra, but easier to put here.
  //  to use in choosing probe sets anyway, if not in which ones to send.. 
  // arghh.
  //cout << "This is Send Region how are you ?"  << endl;

  vector<QString> words = splitString(lastMessage, '|');
  //for(int i=0; i < words.size(); i++){
  //  cout << words[i] << "\t";
  //}
  cout << endl;
  // we need a total of .. 5 words
  if(words.size() != 5){
    cerr << "sendRegion function words size is not 5 but : " << words.size() << endl;
    return;
  }
  bool ok;      // and now for the damn boring ugly code, one day..
  uint index = words[0].toInt(&ok);
  if(!ok || index >= clientRegions.size()){
    cerr << "index may be too large or otherwise screwed up.. index : " << index << "   size of clientRegions : " << clientRegions.size() << endl;
    return;
  }
  int maxExpect = words[1].toInt(&ok);
  if(!ok){
    cerr << "Coulnd't get MaxIndex in sendRegion function.. " << endl;
    return;
  }
  int minLength = words[2].toInt(&ok);
  if(!ok){
    cerr << "Coulnd't get minLength in sendRegion Fuction string is : " << words[2] << endl;
    return;
  }
  float minMatch = words[3].toFloat(&ok);
  if(!ok){
    cerr << "Couldn't get minMatch from words[3] .. which is " << words[3] << endl;
    return;
  }
  int target = words[4].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get target from " << words[4] << "  in sendRegion Function" << endl;
    return;
  }
  // ok, we should be ok, let's set the parameters..
  pSetThresholds.setThresholds(maxExpect, minLength, minMatch);
  // and then we just send the sequence..
  vector<int> includedProbeSets;
  const RegionSpecification* spec = &clientRegions[index];   // for shorthand..
  //  cout << "Calling sendGenomic Region with chrom " << spec->chromosome << "  begin : " << spec->begin << "  end : " << spec->end << endl;
  sendGenomicRegion(spec->chromosome, spec->begin, spec->end, 0, target, includedProbeSets, true);
  // and then set the index... 
  writeIndex(includedProbeSets, "Region Select");
  // which will automatically look for some data.. hmm.
}



void ConnectionObject::doEuclidCompare(){
  // get the index from the message part and then do the euclidSort function..
  bool ok;
  bool notOk = false;
  vector<QString> words = splitString(lastMessage, '|');
  vector<uint> experiments;
  // the first int is the probe set id, the following are the experiments for which to do the thingy.. 
  // words should have at least 3 words..
  
  // the probe set index ..
  uint i;

  if(words.size() < 3){
    cerr << "doEuclidCompare words is too small : " << words.size() << endl;
    // but rather than return and fail, let's select all experiments which are active for the given
    i = lastMessage.toUInt(&ok) - 1;
    // and then 
    if(ok && i < pSet->data.size() && pSet->data[i]->index){
      for(uint j=0; j < pSet->data[i]->exptSize; j++){
	experiments.push_back(pSet->data[i]->exptIndex[j]);
      }
    }else{
      cerr << "doEuclidCompare unable to get probe set info with index " << i << endl;
      return;
    }
  }else{
    i= words[0].toInt(&ok) - 1;
    //  uint i = lastMessage.toInt(&ok);
    notOk = !ok;
    for(uint i=1; i < words.size(); i++){
      experiments.push_back(words[i].toUInt(&ok));
      if(!ok){ notOk = true; }
    }
    if(notOk){
      cerr << "Some problem with getting an experiment in the doEuclidCompare function" << endl;
      return;
    }
  }

  if(i < pSet->data.size() && pSet->data[i]->index){
    euclidSort(pSet->data[i], experiments);
  }else{
    cerr << "Couldn't do euclid Compare, index out of range.. " << i << endl;
  }

}
    

bool ConnectionObject::authenticate(){
  // I used to use the crypt() function, but since I was unable to find it on Windows
  // I gave up and wrote my own stupid little function that converts the password into
  // three integer values using a one-way function. Probably not very secure, but..

  cout << "trying to authenticate.. " << endl;
  char* data = inData->readBinaryFromSocket(socketNumber, 4);   // gets created by function.. 
  if(!data){
    cerr << "Connection Broken or other error return false " << endl;
    return(false);
  }
  int userLength = *(int*)data;
  userLength = ntohl(userLength);  // change the order..
  delete data;
  // then get the string..
  char* cuname = inData->readBinaryFromSocket(socketNumber, userLength, true);  // zero terminated string
  if(!cuname){
    cerr << "ConnectionObject::autheticate Unable to get user name returning false" << endl;
    return(false);
  }
  currentUserName = cuname;
  delete cuname;

  // then I should get 3 keys  
  // not sure if should be sizeof(int) or just f. In any case doing the wrong thing will screw things up
  data = inData->readBinaryFromSocket(socketNumber, sizeof(int) * 3);  
  if(!data){
    cerr << "Connection Broken or other error return false " << endl;
    return(false);
  }
  int key1 = ntohl(*(int*)data);
  int key2 = ntohl(*(int*)(data + sizeof(int)) );
  int key3 = ntohl(*(int*)(data + sizeof(int) * 2) );
  delete data;

  cout << "User : " << currentUserName << endl;
  
  // look up the key.. in the users table.. 
  ostringstream query;       // thread safe???? maybe not.. hmm. 
  query << "select index from users where user_name = '" << currentUserName << "' and key1 = " << key1 << " and key2 = " << key2 << " and key3 = " << key3;
  // and then we just have to see if we return anything.. 

  // I don't quite remember why I'm using the PgCursor interface. I have a feeling that it is broken
  PgCursor conn(conninfo, "portal");
  if(conn.ConnectionBad()){
    cerr << "ConnectionObject::authenticate db connection not good" << endl;
    return(false);
  }
  if(!conn.Declare(query.str().c_str())){
    cerr << "ConnectionObject::authenticate conn.declare command didn't work: " << query << endl;
    return(false);
  }
  if(!conn.Fetch()){
    cerr << "ConnectionObject::authenticate conn.fetch didn't fetch anything" << endl;
    return(false);
  }
  int tuples = conn.Tuples();
  if(tuples > 0){
    authenticated = true;
    userId = atoi(conn.GetValue(0, 0));
    logInTime = time(&logInTime);
    logIn();
    return(true);
  }
  return(false);
}

void ConnectionObject::logIn(){
  if(!logInTime)
    return;
  ostringstream cmd;
  cmd << "insert into user_log values (" << userId << ", " << logInTime << ", 1, current_time, current_date, '" 
      << addressString.latin1() << "', " << socketNumber  << ")";
  int rows = dbCommand(cmd.str());
  if(rows != 1){
    cerr << "ConnectionObject::logOut cmd Rows is " << rows << endl;
  }
}

void ConnectionObject::logOut(){
  if(!logInTime)
    return;
  ostringstream cmd;
  cmd << "insert into user_log values (" << userId << ", " << logInTime << ", 2, current_time, current_date, '" 
      << addressString.latin1() << "', " << socketNumber << ")";
  int rows = dbCommand(cmd.str());
  if(rows != 1){
    cerr << "ConnectionObject::logOut cmd Rows is " << rows << endl;
  }
}


int ConnectionObject::dbCommand(string cmd){
  PGconn* conn = PQconnectdb(conninfo);
  if(PQstatus(conn) != CONNECTION_OK){
    cerr << "ConnectionObject::dbCommand database connection failed" << endl;
    PQfinish(conn);
    return(-1);
  }
  //char* escapedString = new char[cmd.length()*2 + 1]; 
  //int charsWritten = PQescapeString(escapedString, cmd.c_str(), cmd.length());
  //PGresult* res = PQexec(conn, escapedString);
  PGresult* res = PQexec(conn, cmd.c_str());
  int cmdTuples = atoi(PQcmdTuples(res));
  PQfinish(conn);
  //delete escapedString;
  return(cmdTuples);
}

bool ConnectionObject::createSession(){
  //cout << "hello there lastMessage is " << lastMessage << endl;
  vector<QString> words = splitString(lastMessage, '|');        // use a pipe for the delimiter, I don't think the user will use this!
  if(words.size() == 0){ return(false); }
  // otherwise.. lets do some stuff.
  ostringstream sessionInsert;
  /// the problem is quite simple. I dont know what index I will get when I insert into the sessions table, but it is imperative
  /// that I can keep this index during the work period. I can look up this index from the sequence table, or
  //  Ok, get the nextvalue from the session_list, then use this for the both the thing and the other thingyy. using work and others..

  PgDatabase conn(conninfo);
  conn.Exec("select nextval('session_list')");
  if(conn.Tuples() != 1){
    return(false);
  }
  vector<string> keywords;    // for the sessionInformation creator
  string title = words[0].latin1();
  string description;
  set<int> om;       // othermembers, just needed for the constructor.. ugly.. but there you go. 
  int sessionId = atoi(conn.GetValue(0, 0));   // should be OK!!
  bool transactionOk = false;     // set to false if any problems.. 
  // and now begin work.. !!
  sessionInsert << "insert into sessions values (" << sessionId << "," << userId << ", timestamp 'now', timestamp 'now', '"
		<< words[0] << "', ";
  if(words.size() > 1){
    sessionInsert << "'" << words[1] << "'";
    description = words[1].latin1();
  }else{
    sessionInsert << "NULL";
    description = "";
  }
  sessionInsert << ")";
  if(conn.Exec("begin work")){ transactionOk = true; }
  if(transactionOk){
    if(conn.Exec(sessionInsert.str().c_str())){ 
      for(int i=2; i < words.size(); i++){
	ostringstream keywordInsert;
	keywordInsert << "insert into session_keywords values (" << sessionId << ", '" << words[i] << "')";
	if(!conn.Exec(keywordInsert.str().c_str())){ transactionOk = false; }
	keywords.push_back(words[i].latin1());
      }
    }
  }
  int checkValue;
  string errorMessage = "Create Session Failed \n";
  if(transactionOk){
    checkValue = conn.Exec("commit work");
    if(checkValue == PGRES_COMMAND_OK){
      //cout << "everything OK, " << endl;
      sApp("<sessionCreated>");
      qiApp(sessionId);
      sApp("<sessionCreatedEnd>");
      writeArray();
      //////////////  more to the point.. update the pSet->sessions data structure..
      pSet->sessionMutex->lock();
      pSet->sessions.insert(make_pair(sessionId, sessionInformation(sessionId, userId, title, description, keywords, om)));
      pSet->sessionMutex->unlock();
      return(true);
    }else{
      cerr << "We have a problem. PG_RES_COMMAND_OK is not so good. " << endl;
      cerr << conn.ErrorMessage() << endl;
      errorMessage += conn.ErrorMessage();
      //writeString(errorMessage);
      return(false);
    }
  }else{
    checkValue = conn.Exec("rollback work");
    cerr << "Rolling back the session insert for some unfathomable reason" << endl;
    cerr << "conn.Exec check value for rollback is " << checkValue << endl;
  }
  return(false);
  /// and do some error checking.. 
}

bool ConnectionObject::newAnnotation(){
  //cout << "Annotation thing. lastMessage is " << lastMessage << endl;
  vector<QString> words = splitString(lastMessage, '|');
  // hmmm
  // words should be at least 2 words long. Which is to say, -- Hmm, I'm not sure how it counts
  // so lets just print to console first and do some experiments.
  if(words.size() < 2){ return(false); }
  bool ok;
  int sessionId = words[0].toInt(&ok); 
  int temp;
  if(!ok) { return(false); }
  vector<int> geneIndices;
  //cout << "words size is " << words.size() << endl;
  for(int i=2; i < words.size(); i++){
    temp = words[i].toInt(&ok);
    if(ok){
      geneIndices.push_back(temp);
    }else{
      cerr << "Couldn't convert string to number. Aborting annotation " << endl;
      return(false);
    }
    //cout << "word " << i << ":\t" << words[i] << endl;
  }
  //  if(geneIndices.size() == 0){
  //  geneIndices.push_back(0);
  //}
  ostringstream annotationInsert;
  PgDatabase conn(conninfo);
  conn.Exec("select nextval('annotation_count')");
  if(conn.Tuples() != 1){
    return(false);
  }
  int annoteId = atoi(conn.GetValue(0, 0));
  if(sessionId == 0){
    annotationInsert << "insert into user_annotation (index, session_id, annotation, user_id) values (" << annoteId << ", "
		     << sessionId << ", '" << words[1] << "', " << userId << " )";
  }else{
    annotationInsert << "insert into user_annotation (index, session_id, annotation, user_id) "
		     << "select " << annoteId << ", a.index, '" << words[1] << "', a.user_id from sessions a "
		     << "where a.index=" << sessionId << " and a.user_id= " << userId;
  }
  bool insertOk = false;
  if(conn.Exec("begin work")) { insertOk = true; }
  //cout << annotationInsert.str() << endl;
  set<int> associatedGenes;
  conn.Exec(annotationInsert.str().c_str());
  if(conn.CmdTuples() > 0){
    for(int i=0; i < geneIndices.size(); i++){
      ostringstream geneInsert;
      geneInsert << "insert into user_annotation_genes values (" << annoteId << ", " << geneIndices[i] << ")";
      if(!conn.Exec(geneInsert.str().c_str())){ insertOk = false; }
      associatedGenes.insert(geneIndices[i]);
      //cout << geneInsert.str() << endl;
    }
  }else{
    insertOk = false;
  }
  if(insertOk){
    if(conn.Exec("commit work")){
      ////// OK, now everything is Ok.. 
      cout << "Annotation added OK " << endl;
      sApp("<annotationAdded>");
      qiApp(annoteId);
      sApp("<annotationAddedEnd>");
      writeArray();
      //// and more to the point update the various appropriate datastructures.. the lookups and anothers..
      if(sessionId > 0){          // i.e. the comment and the genes are tied to a session,, update sessionLookup.. 
	pSet->sessionMutex->lock();
	for(int i=0; i < geneIndices.size(); i++){
	  pSet->sessionLookup[geneIndices[i]].insert(sessionId);   // hope this is thread safe.. hmm..
	}
	pSet->sessionMutex->unlock();
      }
      if(words[1].length() > 0){    // i.e. there's actually some desription given in the thingy..
	pSet->annotationMutex->lock();
	pSet->userAnnotation.insert(make_pair(annoteId, annotationInformation(annoteId, userId, words[1].latin1(), associatedGenes)));
	for(int i=0; i < geneIndices.size(); i++){
	  pSet->annotationLookup[geneIndices[i]].insert(annoteId);
	}
	pSet->annotationMutex->unlock();
      }
      return(true); 
    }else{
      cout << "Annotation add failed for some reason or other " << endl; 
      string errorMessage = "Annotation add failed ";
      errorMessage += conn.ErrorMessage();
      //writeString(errorMessage);
      return(false);
    }
  }else{
    conn.Exec("rollback work"); 
    return(false);
  }
}

bool ConnectionObject::updateSessionInformation(){
  // tricky one this, as we have to make sure that the session is owned by the thingy. Oh well never mind..
  vector<QString> words = splitString(lastMessage, '|');
  // message should have at least index, title description and then an optional number of keywords which we will 
  if(words.size() < 3){ return(false); }
  // have to deal with separately. 
  bool ok;
  int sessionId = words[0].toInt(&ok);
  vector<string> keywords;
  string title = words[1].latin1();
  string description = words[2].latin1();
  if(!ok){ return(false); } 
  ostringstream update;
  update << "update sessions set u_time = abstime(timestamp 'now'), title = '" << words[1].latin1()
	 << "', description='" << words[2].latin1() << "' where index=" << sessionId
	 << " and user_id=" << userId;
  PgDatabase conn(conninfo);
  conn.Exec(update.str().c_str());
  if(conn.CmdTuples() == 1){  // everything OK, lets remove the old keywords, and insert the new ones..
    conn.Exec("begin");
    ostringstream delString;
    delString << "delete from session_keywords where index=" << sessionId;
    conn.Exec(delString.str().c_str());
    for(int i=3; i < words.size(); i++){
      ostringstream insString;
      insString << "insert into session_keywords values (" << sessionId << ", '"
		<< words[i].latin1() << "')";
      conn.Exec(insString.str().c_str());
      keywords.push_back(words[i].latin1());
    }
    conn.Exec("commit");
    ///////  assume that that worked, and now update the appropriate thingy..
    pSet->sessionMutex->lock();
    map<int, sessionInformation>::iterator it = pSet->sessions.find(sessionId);
    if(it != pSet->sessions.end()){
      (*it).second.title = title;
      (*it).second.description = description;
      (*it).second.keywords = keywords;
    }
    pSet->sessionMutex->unlock();
  }
  
  cout << "updated. tuples affected : " << conn.CmdTuples() << endl;
}

bool ConnectionObject::updateAnnotation(){
  vector<QString> words = splitString(lastMessage, '|');
  // message should have only two words, the index, and the description.
  if(words.size() != 2){
    cerr << "update Annotation, words size is not 2 " << endl
	 << "                         size is    : " << words.size();
    return(false);
  }
  bool ok; 
  int index = words[0].toInt(&ok);
  if(!ok){
    cerr << "couldn't get a number from the QString : " << endl; //<< words[0] << endl;
    return(false);
  }
  // ok, so now we have a description, and we have something else as well, we just have to do an update
  // on the right table.. connect to database..
  ostringstream updateStatement;
  updateStatement << "update user_annotation set annotation = '" << words[1]
		  << "' where index = " << index << " and user_id = " << userId;
  PgDatabase conn(conninfo);
  if(!conn.Exec(updateStatement.str().c_str())){
    cerr << "couldn't update with statement : " << endl;
      // << updateStatement.str() << endl;
    return(false);
  }
  if(conn.CmdTuples() == 1){    // everything Ok.. 
    cout << "updated annotation succsefully " << endl;
    // and let's change the map<int, annotationInformation>
    pSet->annotationMutex->lock();
    map<int, annotationInformation>::iterator it = pSet->userAnnotation.find(index);
    if(it != pSet->userAnnotation.end()){
      (*it).second.annotation = words[1].latin1();
    }else{
      cerr << "couldn't find an annotationinformation for index : " << index << endl
	   << "Can't update annotations data structure though it seems database is updated" << endl;
    }
    pSet->annotationMutex->unlock();
    return(true);  // even if the memory structure is not so good.. 
  }else{
    cerr << "couldn't update database or updated more than one row"
	 << "                                        rows affected : " << conn.CmdTuples() << endl;
    return(false);
  }
  // I shouldn't be able to get here, but just in case..
  return(false);
}
  
void ConnectionObject::associateIshProbeWithEnsemblGene(){
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 2){
    cerr << "associateIshProbeWithEnsemblGene, received the wrong number of words, expected 2, got : " << words.size() << endl;
    return;
  }
  bool ok1, ok2;
  int probeIndex = words[0].toInt(&ok1);
  int ensemblGeneIndex = words[1].toInt(&ok2);
  if(!ok1 || !ok2){
    cerr << "couldn't get ints from probeIndex or ensemblGeneIndex, returning with no update" << endl;
    return;
  }
  /// then check that this ishprobe actually exists .. 
  map<int, ishProbeData>::iterator it;
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(probeIndex);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    cerr << "there is No ish probe data for index : " << probeIndex << endl;
    return;
  }
  /// unlock,,, the mutex, so others can check things. The database interaction might take 
  /// some time. Unfortunately we will have to check that the things exist again,, but never mind..
  pSet->ishProbeDataMutex->unlock();
  ostringstream uquery;
  ostringstream insert;
  // then check do we actually have the privilege of associating genes with thingies.. not many people do
  // if we don't we really should send a message to the the guy in charge..
  uquery << "select a.index from users a, user_privileges b where a.index = " << userId
	 << " and a.index=b.user_index and b.privilege='modify_in_situ_probe_name'";
  PgDatabase conn(conninfo);
  if(!conn.Exec(uquery.str().c_str())){
    cerr << "Couldn't execute query for checking privileges whilst trying to associate an ish probe with a gene" << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(!conn.Tuples()){
    cerr << "it would appear that the user doesn't have privilegese to change stuff " << endl;
    return;
  }
  // if we get here, I suppose that things are looking up...
  insert << "update in_situ_probes set ensembl_guess=" << ensemblGeneIndex << " where probe = " << probeIndex;
  if(!conn.Exec(insert.str().c_str())){
    cerr << "Couldn't update the table, am not sure why this should be the case, but returning nevertheless" << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(!conn.CmdTuples()){
    cerr << "Although the update statement would have seemed to be ok, no tuples affected, return with no action taken" << endl;
    return;
  }
  /// and now, we want to change the internal data here, so we better do stuff..
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(probeIndex);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    cerr << "there is No ish probe data for index : " << probeIndex << endl;
    return;
  }
  (*it).second.ensemblGuess = ensemblGeneIndex;   // hmmm, can I do this.. might be I'm not allowed..
  pSet->ishProbeDataMutex->unlock();
}

void ConnectionObject::setIshProbeName(){
  // parse the request.. should contain,, int probe id, string name.. very simple..
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 2){
    cerr << "setIsh probe name failed, words size is not 2, but " << words.size() << endl;
    return;
  }
  bool ok;
  int probeId = words[0].toInt(&ok);
  if(!ok){
    cerr << "couldn't get an int from words[0] " << endl; //words[0] << endl;
    return;
  }
  map<int, ishProbeData>::iterator it;
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(probeId);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    cerr << "Couldn't find an ishprobe with an index " << probeId << endl;
    return;
  }
  pSet->ishProbeDataMutex->unlock();
  /// check if we are allowed to change the names of the ish probes.. 
  ostringstream uquery;
  uquery << "select a.index from users a, user_privileges b where a.index = " << userId
	 << " and a.index=b.user_index and b.privilege='modify_in_situ_probe_name'";
  PgDatabase conn(conninfo);
  if(!conn.Exec(uquery.str().c_str())){
    cerr << "Couldn't execute query for checking privileges whilst trying to name an ish probe" << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(!conn.Tuples()){
    cerr << "it would appear that the user doesn't have privilegese to name ish probes, by by .. " << endl;
    return;
  }
  // and lets make a string for updating the table..
  ostringstream update;
  ////////////////////////////////  try to use the cleaning function from the pqsql.. have to check if it is available..
  //
  char* escapedString = new char[words[1].length()*2 + 1]; 
  int charsWritten = PQescapeString(escapedString, words[1].latin1(), words[1].length());

  update << "update in_situ_probes set probe_name = '" << escapedString << "' where probe = " << probeId;
  // and delete the string
  delete escapedString;
  if(!conn.Exec(update.str().c_str())){
    cerr << "couldn't update the name columne in the in situ probes table.. " << endl
	 << conn.ErrorMessage() << endl;
    return; 
  }
  if(!conn.CmdTuples()){
    cerr << "Although the update statement was accepted, no row was affected hmmm " << endl;
    return;
  }
  // now just update the internal data stuff here...  --
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(probeId);
  if(it == pSet->ishProbes.end()){
    cerr << "very bizarre, the ishprobe data seems to have disappeared while I updated the database table, man that sucks.." << endl;
    return;
  }
  (*it).second.probeName = words[1].latin1();
  pSet->ishProbeDataMutex->unlock();
  // done.. maybe .... 
}

void ConnectionObject::insertIshProbeTextAnnotation(){
  // client sends.. -- data id, ish probe id, bool fieldNamenew, fieldName, annotation.. 
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 5){
    cerr << "insertIshProbeTextAnnotation words size is not 5, but " << words.size() << endl;
    return;
  }
  bool ok1, ok2, ok3;
  int id = words[0].toInt(&ok1);
  int ishProbeId = words[1].toInt(&ok2);
  int fieldNameNew = words[2].toInt(&ok3);
  if(!ok1 || !ok2 || !ok3){
    cerr << "insertIshProbeTextAnnotation, failed to obtain one of the integer paramters" << endl;
    return;
  }
  // first check if we have an appropriate ish probe .. otherwise no point in doing this eh.. 
  map<int, ishProbeData>::iterator it;
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    // no good, get out of here..
    cerr << "No such ish probe id : " << ishProbeId << endl;
    return;
  }
  pSet->ishProbeDataMutex->unlock();
  // it should remain a pointer to the appropriate data structure, but we cannot really assume that noone else gets rid of that
  // structure in the meanwhile, so we will have to check it out later on for updating stuff...
  // clean the fieldName and the annotation name.
  char* escapedFName = new char[words[3].length() * 2 + 1];
  char* escapedAnnot = new char[words[4].length() * 2 + 1];
  PQescapeString(escapedFName, words[3].latin1(), words[3].length());
  PQescapeString(escapedAnnot, words[4].latin1(), words[4].length());
  // and let's see if we how to do this properly.. --
  // 1. we have to check if the fieldname exists.. - i.e. look up in the ish_annotation_text_fields
  //    if it exists, then use the integer given and do an insert into the thingy..
  //    if it doesn't exist, then insert, then check the name.. -- ahh, this is a real bummer..
  //    and then insert.. -- and hopefully everything will be fine..
  ostringstream fieldCheck;
  fieldCheck << "select index from ish_annotation_text_fields where field = '" << escapedFName << "'";
  PgDatabase conn(conninfo);
  if(!conn.Exec(fieldCheck.str().c_str())){
    cerr << "couldn't check the text annotation field name from the database return without inserting " << endl
	 << conn.ErrorMessage() << endl;
    /// I should call some function that sends back the error message as well, but I will make that later..
    delete escapedFName;
    delete escapedAnnot;
    return;
  }
  if(!conn.Tuples()){
    // we need to insert the field name..
    ostringstream fInsert;
    fInsert << "insert into ish_annotation_text_fields (field) values ('" << escapedFName << "')";
    if(!conn.Exec(fInsert.str().c_str())){
      cerr << "Couldn't insert the field name into the text annotation field value table.. bummer " << endl
	   << conn.ErrorMessage() << endl;
      delete escapedFName;
      delete escapedAnnot;
      return;
    }
    if(conn.CmdTuples() != 1){
      cerr << "Although the insert would have appeared to been successful, it seems that it didnt' affect any rows.. very strange " << endl
	   << "Number of rows affected : " << conn.CmdTuples() << endl;
      delete escapedFName;
      delete escapedAnnot;
      return;
    }
    // unfortuanetely, I know have to try and see if I can get the index again..
    if(!conn.Exec(fieldCheck.str().c_str())){
      cerr << "-- after insertion.. couldn't check the text annotation field name from the database return without inserting " << endl
	   << conn.ErrorMessage() << endl;
      /// I should call some function that sends back the error message as well, but I will make that later..
      delete escapedFName;
      delete escapedAnnot;
      return;
    }
    if(!conn.Tuples()){
      cerr << "even I tried to insert it and it seems to have worked, something screwed up somewhere, return without doing anything " << endl;
      delete escapedFName;
      delete escapedAnnot;
      return;
    }
  }
  /// if I get here, conn.Tuples is something,, that's good.. and it would appear that I've inserted it into the database structure..
  /// that is nice.. maybe I need to change the pSet data set to indicate the changes... hmm. 
  // no need to update the textAnnotation map until we insert into the database..

  // it is possible that we have a new field name, so we better try to insert it into the thingy..
  pSet->ishProbeDataMutex->lock();
  pSet->ishTextFields.insert(words[3].latin1());
  pSet->ishProbeDataMutex->unlock();
  // -- easier than checking and things.. --- in the long run should check, and .. emit some signal if this happens..
  // -- but for now, leave it at this.. 
  int fieldIndex = atoi(conn.GetValue(0, 0));
  //cout << "The field Index for : " << escapedFName << "  is " << fieldIndex << endl;
  // we don't need the escapedFName any more as we will just use original string.. thingy..
  delete escapedFName;
  // make an annotInsert statement..
  ostringstream annotInsert;
  annotInsert << "insert into ish_probe_text_annotation (user_index, probe_id, field, annotation) values (" << userId << ", "
	      << ishProbeId << ", " << fieldIndex << ", '" << escapedAnnot << "')";
  delete escapedAnnot;
  if(!conn.Exec(annotInsert.str().c_str())){
    cerr << "Bloody hell, couldn't insert the annot insert thingy,, byby, " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(conn.CmdTuples() != 1){
    cerr << "Conn Cmd Tuples is not 1,, but : " << conn.CmdTuples() << endl
	 << "what is going on " << endl;
    cerr << "the insert statement was : " << endl; // annotInsert.str() << endl;
    return;
  }
  // if here everything ok,, now just have to update the internal data structure.. bit of a pain that..
  // need the oid from the thingy.. let's guess that we can do.. 
  //int annotationId = (int)conn.OidValue();

  //cout << "and the OidStatus is : " << conn.OidStatus() << endl;
  int annotationId = atoi(conn.OidStatus());
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    cerr << "What the hell, the ish probe stuff disappeared during the inserts,, tonight, I will cry in my beer maybe " << endl;
    pSet->ishProbeDataMutex->unlock();
    return;
  }
  /////////////// fix this,,, 
  (*it).second.textAnnotation.insert(make_pair(annotationId, ish_annotation(annotationId, userId, currentUserName, words[4].latin1(), words[3].latin1())));
  pSet->ishProbeDataMutex->unlock();
  // and we are done... 
}

void ConnectionObject::updateIshProbeTextAnnotation(){
  vector<QString> words = splitString(lastMessage, '|');
  // message size should only be 2..
  if(words.size() != 3){
    cerr << "updateIshProbeTextAnnotation words size is not 3 but : " << words.size() << endl;
    return;
  }
  bool ok;
  bool ok2;
  int oid = words[0].toInt(&ok);
  int probeId = words[1].toInt(&ok2);
  if(!ok || !oid || !ok2){
    cerr << "Trouble with the oid identifier, not sure what's up here.. word is :"  << endl; // words[0] << endl;
    return;
  }
  /// ish_annotation now contains the field name as well.. and is organised by the oid, so that we can find it..
  /// however, we still need to know the index of the probe.. 
  char* escapedNote = new char[words[2].length()*2 + 1];
  PQescapeString(escapedNote, words[2].latin1(), words[2].length());  // ok..
  ostringstream update;
  update << "update ish_probe_text_annotation set annotation ='" << escapedNote << "' where oid=" << oid 
	 << " and user_index=" << userId;
  delete escapedNote;    // we don't actually need it anymore.. 
  PgDatabase conn(conninfo);
  if(!conn.Exec(update.str().c_str())){
    cerr << "Coulnd't update the database table " << endl
	 << conn.ErrorMessage() << endl;
    //    delete escapedNote;
    return;
  }
  if(conn.CmdTuples() != 1){
    cerr << "Command Tuples is : " << conn.CmdTuples() << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  map<int, ishProbeData>::iterator it;
  map<int, ish_annotation>::iterator iit;   
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(probeId);
  if(it != pSet->ishProbes.end()){
    iit = (*it).second.textAnnotation.find(oid);
    if(iit != (*it).second.textAnnotation.end()){
      (*iit).second.annotation = words[2].latin1();
    }
  }
  pSet->ishProbeDataMutex->unlock();
  
}
    

void ConnectionObject::insertIshProbeFloatAnnotation(){
  // 
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 5){
    cerr << "insertIshProbeFloatAnnotation words of wrong size : " << words.size() << "should be 5 " <<  endl;
    return;
  }
  bool ok1, ok2, ok3, ok4;
  int id = words[0].toInt(&ok1);
  int ishProbeId = words[1].toInt(&ok2);
  int fieldNameNew = words[2].toInt(&ok3);   // although we don't actually use this ..
  float value = words[4].toFloat(&ok4); 
  if(!ok1 || !ok2 || !ok3 || !ok4){
    cerr << "Couldn't get one of the parameters so returning without insertion " << endl;
    return;
  }
  // first check if we have the ish probe, otherwise doesn't make any sense.
  map<int, ishProbeData>::iterator it;
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    cerr << "No such ish probe id : " << ishProbeId << endl;
    return;
  }
  pSet->ishProbeDataMutex->unlock();
  
  /// at this point we need to escape the field name so that we can use it in db queries..
  // try not to forget to delete it.. !!
  char* escapedFName = new char[words[3].length() * 2 +1];
  PQescapeString(escapedFName, words[3].latin1(), words[3].length());
  // then use to check if the fieldname is available..
  ostringstream fieldCheck;
  fieldCheck << "select index from ish_probe_num_fields where field = '" << escapedFName << "'";
  PgDatabase conn(conninfo);
  if(!conn.Exec(fieldCheck.str().c_str())){
    cerr << "insertIshProbeFloatAnnotation couldn't select from annotation fields bugger " << endl;
    delete escapedFName;
    return;
  }
  if(!conn.Tuples()){
    cerr << "new field. let's try to insert the field into the database.. ho yeaah " << endl;
    ostringstream insert;
    insert << "insert into ish_probe_num_fields (field) values ('" << escapedFName << "')";
    if(!conn.Exec(insert.str().c_str())){
      cerr << "couldn't insert the field name into the ish_probe_num_fields " << endl
	   << conn.ErrorMessage() << endl;
      // I really don't know what is going on here, but there you go... 
      delete escapedFName;
      return;
    }
    /// just try to select it using the original thingy..
    if(!conn.Exec(fieldCheck.str().c_str())){
      cerr << "insertIshProbeFloatAnnotation, can't select from the field table on second occasion.. " << endl
	   << conn.ErrorMessage() << endl;
      delete escapedFName;
      return;
    }
    if(!conn.Tuples()){
      cerr << "although we tried to insert the field name we seem to have failed, not sure why error message : " << endl
	   << conn.ErrorMessage() << endl;
      delete escapedFName;
      return;
    }
  }
  // if we get here, then we really should be able to get the field index from the bugger, and then we can delete the field name
  // and try to insert the actual classification thingy.. bugger...
  int fieldIndex = atoi(conn.GetValue(0, 0));
  delete escapedFName;
  ostringstream insert;   // new context..
  insert << "insert into ish_probe_num_annotation (user_index, probe_id, field, annotation) values (" << userId << ", " 
	 << ishProbeId << ", " << fieldIndex << ", " << value << ")";
  if(!conn.Exec(insert.str().c_str())){
    cerr << "couldn't insert into ish_probe_num_annotation " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(conn.CmdTuples() != 1){
    cerr << "Although the insert was allowed nothing inserted due to some error " << endl
	 << conn.ErrorMessage();
    return;
  }
  int annotationId = atoi(conn.OidStatus());
  // and at this point we can update the rest of the things..
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    cerr << "The ish probe annotation seems to have disapperaed while I was looking the other way, trouble man.. " << endl;
    pSet->ishProbeDataMutex->unlock();
    return;
  }
  (*it).second.numberAnnotation.insert(make_pair(annotationId, ish_annotation(annotationId, userId, currentUserName, value, words[3].latin1())));
  // and might as well just do an insert into the fields.. for the field name
  pSet->ishFloatFields.insert(words[3].latin1());   // it shares the same mutex for convenience..
  pSet->ishProbeDataMutex->unlock();
  // and we are done.. 
}

void ConnectionObject::insertIshProbeClassification(){
  //
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 5){
    cerr << "insertIshClassification words of wrong size : " << words.size() << "should be 5 " <<  endl;
    return;
  }
  bool ok1, ok2, ok3, ok4;
  int id = words[0].toInt(&ok1);
  int ishProbeId = words[1].toInt(&ok2);
  int fieldNameNew = words[2].toInt(&ok3);   // although we don't actually use this ..
  float value = words[4].toFloat(&ok4); 
  if(!ok1 || !ok2 || !ok3 || !ok4){
    cerr << "Couldn't get one of the parameters so returning without insertion " << endl;
    //cout << words[0] << ", " << words[1] << ", " << words[2] << ", " << words[3] << ", " << words[4] << endl;
    return;
  }
  map<int, ishProbeData>::iterator it;
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    pSet->ishProbeDataMutex->unlock();
    cerr << "No such ish probe id : " << ishProbeId << endl;
    return;
  }
  pSet->ishProbeDataMutex->unlock();
  
  /// at this point we need to escape the field name so that we can use it in db queries..
  // try not to forget to delete it.. !!
  char* escapedFName = new char[words[3].length() * 2 +1];
  PQescapeString(escapedFName, words[3].latin1(), words[3].length());
  // then use to check if the fieldname is available..
  ostringstream fieldCheck;
  fieldCheck << "select index from ish_probe_classes where class = '" << escapedFName << "'";
  PgDatabase conn(conninfo);
  if(!conn.Exec(fieldCheck.str().c_str())){
    cerr << "insertIshProbeClassification couldn't select from class table..  bugger " << endl;
    delete escapedFName;
    return;
  }
  if(!conn.Tuples()){
    cerr << "new class. let's try to insert the class into the database.. ho yeaah " << endl;
    ostringstream insert;
    insert << "insert into ish_probe_classes (class) values ('" << escapedFName << "')";
    if(!conn.Exec(insert.str().c_str())){
      cerr << "couldn't insert the class name into the ish_probe_classes " << endl
	   << conn.ErrorMessage() << endl;
      // I really don't know what is going on here, but there you go... 
      delete escapedFName;
      return;
    }
    /// just try to select it using the original thingy..
    if(!conn.Exec(fieldCheck.str().c_str())){
      cerr << "insertIshProbeClassification, can't select from the field table on second occasion.. " << endl
	   << conn.ErrorMessage() << endl;
      delete escapedFName;
      return;
    }
    if(!conn.Tuples()){
      cerr << "although we tried to insert the class name we seem to have failed, not sure why error message : " << endl
	   << conn.ErrorMessage() << endl;
      delete escapedFName;
      return;
    }
  }
  /// everything should be OK.. 
  int fieldIndex = atoi(conn.GetValue(0, 0));
  delete escapedFName;
  ostringstream insert;   // new context..
  insert << "insert into ish_probe_classification (user_index, probe_id, class, confidence) values (" << userId << ", " 
	 << ishProbeId << ", " << fieldIndex << ", " << value << ")";
  if(!conn.Exec(insert.str().c_str())){
    cerr << "couldn't insert into ish_probe_num_annotation " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(conn.CmdTuples() != 1){
    cerr << "Although the insert was allowed nothing inserted due to some error " << endl
	 << conn.ErrorMessage();
    return;
  }
  int annotationId = atoi(conn.OidStatus());
  // ok evertying should be fine.. 
  pSet->ishProbeDataMutex->lock();
  it = pSet->ishProbes.find(ishProbeId);
  if(it == pSet->ishProbes.end()){
    cerr << "The ish probe annotation seems to have disapperaed while I was looking the other way, trouble man.. " << endl;
    pSet->ishProbeDataMutex->unlock();
    return;
  }
  (*it).second.classification.insert(make_pair(annotationId, ish_annotation(annotationId, userId, currentUserName, value, words[3].latin1())));
  // and might as well just do an insert into the fields.. for the field name
  pSet->ishClasses.insert(words[3].latin1());   // it shares the same mutex for convenience..
  pSet->ishProbeDataMutex->unlock();
  // and we are done.. 
}

void ConnectionObject::sendSessionInformation(){
  // do in three steps..
  // 1. get the sessions.. with user id, the title, and the description.. 
  // 2. get the keywords..
  // 3. get the associated genes.
  // 
  // I need some kind of datastructures to hold the data before I send it. hhmmm. 
  cout << "sendSessionInformation function. " << endl;

  map<int, sessionData> sData;    // hope this will work ok with threads.. Hmm. who knows..
  PgDatabase conn(conninfo);
  if(!conn.Exec("begin")){ 
    cerr << "couldn't begin " << endl;
    return; 
  }
  if(!conn.Exec("set transaction isolation level serializable")) { 
    cerr << "Couldn't set transcation isolation level " << endl;
    return; 
  }
  if(!conn.Exec("select index, user_id, title, description from sessions")){
    cerr << "Couldn't select from sessions " << endl;
    cerr << conn.ErrorMessage() << endl;
    return; 
  }
  for(int i=0; i < conn.Tuples(); i++){
    sessionData tData;
    tData.sessionIndex = atoi(conn.GetValue(i, 0));
    tData.userId = atoi(conn.GetValue(i, 1));
    tData.sessionTitle = conn.GetValue(i, 2);
    tData.description = conn.GetValue(i, 3);
    sData.insert(make_pair(tData.sessionIndex, tData));
  }
  // and lets get the keywords..
  if(!conn.Exec("select a.index, b.keyword from sessions a, session_keywords b where a.index=b.index order by a.index")){
    cerr << "couldn't select from keywords table " << endl;
    return; 
  }
  map<int, sessionData>::iterator it;
  int oldSession = -1; // should be impossible.. 
  for(int i=0; i < conn.Tuples(); i++){
    int session = atoi(conn.GetValue(i, 0));
    if(session != oldSession){
      oldSession = session;
      it = sData.find(session);
    }
    if(it != sData.end()){
      (*it).second.keyWords.insert(conn.GetValue(i, 1));
    }
  }
  // and then lets get the words,, in a a very similar way..
  if(!conn.Exec("select distinct a.session_id, b.gene from user_annotation a, user_annotation_genes b where a.index=b.annotation_id and a.session_id > 0")){
    cerr << "couldn't select from user_annotation_genes table.. " << endl;
    return;
  }
  oldSession = -1;   // again..
  for(int i=0; i < conn.Tuples(); i++){
    int session = atoi(conn.GetValue(i, 0));
    if(session != oldSession){
      oldSession = session;
      it = sData.find(session);
    }
    if(it != sData.end()){
      (*it).second.geneIds.insert(atoi(conn.GetValue(i, 1)));
    }
  }
  // at which point we've got all of the data and we can commit the cursor and send it to the local client.. 
  if(!conn.Exec("commit")){ 
    cerr << "couldn't commit for some obscure reason " << endl;
    return; 
  }
  // and then lets do some appending.. !!
  sApp("<sessionDescription>");
  qiApp(sData.size());
  for(it = sData.begin(); it != sData.end(); it++){
    qiApp((*it).first); // well it's the index.. !
    if((*it).second.userId == userId){
      qiApp(1);
    }else{
      qiApp(0);
    }
    //  qiApp((*it).second.userId);
    sApp((*it).second.sessionTitle);
    sApp((*it).second.description);
    qiApp((*it).second.keyWords.size());
    set<string>::iterator sit;
    for(sit = (*it).second.keyWords.begin(); sit != (*it).second.keyWords.end(); sit++){
      sApp((*sit));
    }
    qiApp((*it).second.geneIds.size());
    set<int>::iterator iit;
    for(iit= (*it).second.geneIds.begin(); iit != (*it).second.geneIds.end(); iit++){
      qiApp(*iit);
    }
  }
  sApp("<sessionDescriptionEnd>");
  writeArray();  // and there we go. no error reporting at all, but what can you do!! 
}

void ConnectionObject::sendProbeSetSequence(){
  // get the request id and stuff..
  vector<QString> words = splitString(lastMessage, '|');
  // first word is the affyId, and second word is the request Id
  if(words.size() != 2){
    cerr << "sendProbeSetSequence : wrong number of words, have " << words.size() << "  instead of the expected 2 words" << endl;
    return;
  }
  bool ok;
  uint requestId = words[1].toInt(&ok);
  if(!ok){
    cerr << "sendProbeSetSequence : couldn't get an int from the second word " << endl;
    return;
  }
  uint afIndex = words[0].toInt(&ok);
  if(!ok){
    cerr << "sendProbeSetSequence : couldn't get an int from the second word " << endl;
    return;
  }
  /// words[0] is now the int, i.e. the probe set request, so we have to do a joing.. 
  // format the string.. 
  string request = string("select a.sequence, a.af_id from seq a, p_sets b where a.af_id=b.af_id and b.index =") + words[0].latin1();  // shopuld beOK
  // do a database lookup to get the sequence and then write this along with things..
  //const char* conninfo = "dbname=expression";
  PgCursor conn(conninfo, "portal");
  if(!conn.Declare(request.c_str())){
    cerr << "couldn't declare the probe set sequence request " << endl;
    return;
  }
  if(!conn.Fetch()){
    cerr << "couldn't fetch the data,, don't know what's going on " << endl;
    return;
  }
  if(conn.Tuples() != 1){
    cerr << "conn.tuples isn't 1. returning without doing anything, but should do this differently" << endl;
    return;
  }
  string sequence = conn.GetValue(0, 0);  // there is only one !!
  string afid = conn.GetValue(0, 1);
  // and write the message and get out of here..
  sApp("<dnaSequence>");
  qiApp(requestId);
  qiApp(0); // the type of sequence..
  sApp(afid);  // hmm, what an ugly mix..
  qiApp(afIndex);
  sApp(sequence);
  sApp("<dnaSequenceEnd>");
  writeArray();
  // as usual, no error handling, so it kind of sucks.. 
}

void ConnectionObject::sendEnsemblPeptide(){
  // get the words.. 
  vector<QString> words = splitString(lastMessage, '|');
  bool ok;
  // and lets get some numbers from the thingy..
  if(words.size() < 2){
    cerr << "sendEnsemblPeptide no words " << endl;
    return;
  }
  int transcriptIndex = words[0].toInt(&ok);
  if(!ok){ 
    cerr << "coulndn't get a gene index" << endl;
    return;
  }
  int requestId = words[1].toInt(&ok);
  if(!ok){
    cerr << "coulnd't get a request Id " << endl;
    return;
  }
  ostringstream queryString;
  //  queryString << "select a.name, a.gene, b.sequence from ensembl_peptides a, ensembl_peptide_sequence b "
  //	      << "where a.index = b.index and a.gene=" << geneIndex;
  queryString << "select a.id, b.id, a.sequence from ensembl_peptide_9 a, ensembl_transcript_9 b "
	      << "where a.transcript=b.transcript and a.transcript = " << transcriptIndex;
  PgCursor conn(conninfo, "portal");
  if(!conn.Declare(queryString.str().c_str())){
    cerr << "couldn't declare the query for the peptide sequence" << endl;
    return;
  }
  if(!conn.Fetch()){
    cerr << "coulnd't fetch peptide sequences from database " << endl;
    return;
  }
  sApp("<dnaSequence>");
  qiApp(requestId);
  qiApp(4);  // peptide sequence 
  int tuples = conn.Tuples();
  qiApp(transcriptIndex);
  qiApp(tuples);  // so I know how many..
  for(int i=0; i < tuples; i++){
    sApp(conn.GetValue(i, 0));
    qiApp(atoi(conn.GetValue(i, 1)));
    sApp(conn.GetValue(i, 2));
  }
  sApp("<dnaSequenceEnd>");
  writeArray();
  // and that should be pretty much everything.. 
  
}

void ConnectionObject::sendEnsemblTranscriptSequence(){
  vector<QString> words = splitString(lastMessage, '|');
  // first word is the transcript index, second one is the request index.. not so troublesome..
  // but..
  if(words.size() != 2){
    cerr << "sendEnsemblTranscriptSequence words size is not 2, but " << words.size() << endl;
    return;
  }
  bool ok;
  int index = words[0].toInt(&ok);
  if(!ok){
    cerr << "couldn't get the ensembl Transcript index from word 0. word 0 is " << endl; //words[0] << endl;
    return;
  }
  int requestId = words[1].toInt(&ok);
  if(!ok){
    cerr << "couldn't get the request id (sendEnsemblTranscriptSequence) from words[1] " << endl; // words[1] << endl;
    return;
  }
  // hmm, now how to do this..
  // remember that the index is actually the gene index, -not the transcript index. Hence we need to get the 
  // all transcript sequences associated with a given gene index.. -which looks like a 3 table join. 
  ostringstream queryString;  // -8 is the strand.. 
  queryString << "select b.id, b.chromosome, b.start, a.sequence, c.annotation from ensembl_transcript_sequence_9 a, ensembl_transcript_9 b, ensembl_annotation_9 c where a.transcript = b.transcript and c.gene = b.gene and c.field = 8 and b.transcript = " << index;
  
  // Ok, that should be it, now get the usual suspect connection to the database..
  // maybe one day keep one connection instead of making new ones each time.. oh well.. when I'm old and wise.
  //const char* conninfo = "dbname=expression";
  PgCursor conn(conninfo, "portal");
  if(!conn.Declare(queryString.str().c_str())){
    cerr << "Couldn't declare the cursor for the transcript query thing.. never mind " << endl;
    return;
  }
  if(!conn.Fetch()){
    cerr << "Couldn't fetcht the data for transcript query thing, bugger " << endl;
    return;
  }
  int tuples = conn.Tuples();
  if(tuples == 0){
    cerr << "query didn't retrieve anythying. That sucks, " << endl;
    //return;
  }
  // ok do this the dangerous way..
  sApp("<dnaSequence>");
  qiApp(requestId);
  qiApp(1);    // sequence type..
  qiApp(index);  // just to remind the client, what it might be looking for. 
  qiApp(tuples);  // number of sequences.. 
  for(int i=0; i < tuples; i++){
    sApp(conn.GetValue(i, 0));
    sApp(conn.GetValue(i, 1));
    qiApp(atoi(conn.GetValue(i, 2)));
    sApp(conn.GetValue(i, 4));  // the strand
    sApp(conn.GetValue(i, 3)); // the sequence .. 
  }
  sApp("<dnaSequenceEnd>");
  writeArray();
}
	  
    
    
  

void ConnectionObject::sendGenomeSequence(){
  // read the sequence from the appropriate file using the genome coordinates and seekg.. to get to the right place in the
  // file. Then get the sequence and send to the owner. am not entirely sure how to convert to string, but we'll have a look
  // and see what we can do. 
  // make sure to protect the thing with a mutex, so that we don't alter the seek position between reads and stuff like that.
  // 
  // obviously can't read in stupid amount of information, as this is likely to end up with us running out of memory. Make maximum size
  // of 10, 0000 bp or something like that. We could send stuff in chunks, and thus send longer stretches of dna. Hmm, might be worthwhile, though, 
  // that would require something smarter down the other end, and it might be better for the client to keep requesting chunks, of the appropriate size
  // we can ofcourse change this in the future, but for now lets keep the thing simple. Ok..
  // still don't know how to display the sequence at the other end, but will work on it.
  //
  cout << "At the beginning of sendGenomeSequence, have to work on this for a bit, many things to sort out" << endl;

  int maxSize = 5000000;     // the maximum size to send..  be generous, 100, 000 is ok,, only 100 kB.. hoowah. 
  // first parse the message and see if it makes any sense. 
  vector<QString> words = splitString(lastMessage, '|');     // pipe probably the most common..
  if(words.size() != 4){
    cerr << "sendGenomeSequence : words size is not 4, but : " << words.size() << endl;
    return;
  }
  // first word is the chromosome name
  // second word is the start position
  // last word is the stop position.. - 
  
  // first check if we have chromosome by that name, --
  // -- seems that sgi consider concurrent reads (& finds to be thread safe??)
  // so I don't need to lock anything here.. 
  map<string, ifstream*>::iterator it = pSet->chromFiles.find(words[0].latin1());
  if(it == pSet->chromFiles.end()){
    cerr << "sendGenomeSequence : no such chromosome : " << endl; //words[0] << endl;
    return;
  }
  // now get the start and stop positions..
  bool ok;
  uint start = words[1].toInt(&ok);
  if(!ok){
    cerr << "sendGenomeSequence : couldn't get a number from the second word : " << endl; //words[1] << endl;
    return;
  }
  uint stop = words[2].toInt(&ok);
  if(!ok){
    cerr << "sendGenomeSequence : couldn't get a number from the third word : " << endl; //words[2] << endl;
    return;
  }
  int requestId = words[3].toInt(&ok);
  if(!ok){
    cerr << "sendGenomeSequence : couldn't get a number from the id, what is going on?  " << endl; //  words[3] << endl;
    return;
  }
  //cout << "send Genome Sequence : chromsome " << words[0] << "   start : " << start << "   end : " << stop << "   id : " << requestId << endl;
  // make sure both are above 0..
  if(!(start > 0 && stop >= start)){             // the ensembl referencing counts sequence positions from 1, not 0, like we do in the file. 
    cerr << "sendGenomeSequence : not very happy with start and stop numbers try again" << endl
	 << "           start   : " << start << "\tstop : " << stop << endl;
    return;
  }
  // ok, now everything should be pretty much hunky dory and all that, now we should be able to do more stuff. 
  if((stop - start) > maxSize){
    cerr << "requested range too large, setting stop to : " << start+maxSize << endl;
    stop = start + maxSize;
  }
  // and now we need to to the delicate locking procedure and stuff. Let's look for a QMutex..
  int length = stop - start + 1;     // inclusive reading,, -hence need an extra 1...
  char* sequence = new char[length];
  pSet->chromFileMutex->lock();                     // this is probably the wrong way to do it.. 
  // seek to the right location..
  (*it).second->seekg(start);
  if((*it).second->eof()){
    cerr << "tried seeking past end of file, that's not particularly good return no sequence " << endl;
    delete sequence;
    (*it).second->clear();
    pSet->chromFileMutex->unlock();
    return;
  }
  // then read the sequence into the thingy..
  (*it).second->read(sequence, length);
  // and determine how much was read..
  int readLength = (*it).second->gcount();
  // and make sure that it's ok
  if((*it).second->fail()){
    (*it).second->clear();
  }
  // as the glength tells us what happened, I think we can assume everything is OK.. if it's 0, it's not so good, but what the hell
  // Unlock the mutex.
  pSet->chromFileMutex->unlock();
  // and then let's write a nice message to the client as usual..
  sApp("<dnaSequence>");
  // lets fill in the chromosome and the numbers as well..
  qiApp(requestId);  // make sure the idRequest is OK.. !!! .. 
  qiApp(2);          // the type of sequence.. 
  sApp(words[0].latin1());  // what an unholy mix 
  qiApp(start);
  qiApp(readLength);  // i.e. not the stop, but how much read. 
  sApp(sequence, readLength);      // just hope that we don't have stupid things happening.
  // ok, now lets delete the sequence buffer
  delete sequence;
  // and send the terminator, and then write the array..
  sApp("<dnaSequenceEnd>");
  writeArray();
  // and voila, the genome brought to you.. courtesy of a couple of open files. Let's see if that actually works
}  
  

void ConnectionObject::doDBLookup(){
  // could probably be faster if I set up a connection to begin with, rather than do each time.. 
  // first, parse the message..
  char delimiter = '#';
  int delimPosition = lastMessage.find(delimiter);
  int mode;
  //const char* conninfo = "dbname=expression";
  QString query;
  if(delimPosition != -1){
    query = lastMessage.left((unsigned int)delimPosition);
    int rightMost = lastMessage.length()-delimPosition-1;
    if(rightMost > 0){
      QString modeString = lastMessage.right((unsigned int)rightMost);
      bool ok;
      mode = modeString.toInt(&ok);
      if(!ok){
	cerr << "Couldn't get a database mode (int) variable, assuming 0" << endl;
      }
    }else{
      mode = 0;
    }
  }else{
    cerr << "can't parse the query string: doDBLookup" << endl; //<< lastMessage << endl;
    return;
  }
  // and now format the query string...
  // but we should use the PqEscapeString function here.. so that we don't fail for all the wrong reasons.
  char* escapedQuery = new char[query.length()*2 + 1];
  PQescapeString(escapedQuery, query.latin1(), query.length());
  QString fullQuery;
  QTextStream ts(fullQuery, 2);     // 2 is the io read mode.. I don't know what it means.. 
  switch(mode){
  case 0:
    ts << "select distinct a.index from probe_data a, p_sets b where (a.title ~* '" << escapedQuery
       << "' or a.gene ~* '" << escapedQuery << "' or a.description ~* '" << escapedQuery
       << "') and a.index=b.index and b.type!='Control Sequence' order by index";
    break;
  case 1:
    ts << "select distinct gene from ensembl_annotation_9 where annotation ~* '" << escapedQuery << "'";
    break;
  case 2:
    ts << "select index from p_sets where af_id ~* '" << escapedQuery << "'";
    break;
  case 3:
    ts << "select distinct a.index from p_sets a, af_go_gen d, go e where a.af_id = d.af_id  and d.go=e.index and e.description ~* '"
       << escapedQuery << "'";
    break;
  case 4:
    ts << "select distinct afid from in_situ_probes where probe_id ~* '" << escapedQuery
       << "'" << " and afid > 0";
    break;
  case 5:
    ts << "select distinct a.index from p_sets a, era_pdgf_ug d where a.suid=d.uid and d.description ~* '" << escapedQuery << "'";
    break;
  case 6:
    ts << "select distinct a.index from p_sets a, body_map_genes b, body_map_counts c where c.bmid = b.index and b.uid !=0 and b.uid=a.suid and c.tissue="
       << escapedQuery << "' and c.count > 0";
    break;
  default :
    ts << "select distinct index from probe_data where title ~* '" << escapedQuery
       << "' or gene ~* '" << escapedQuery << "' or description ~* '" << escapedQuery
       << "' order by index";
  }
  delete escapedQuery;
  //cout << "Full Query is : " << fullQuery.latin1() << endl;
  // and now get a database connection..
  //PgCursor conn(conninfo, "portal");
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "connection not good" << endl;
    return;
  }
  if(!conn.Exec(fullQuery.latin1())){
    cerr << "command didn't work: " << fullQuery << endl;
    return;
  }
  //if(!conn.Fetch()){
  // cerr << "couldn't fetch the data" << endl;
  // return;
  //}
  int tuples = conn.Tuples();
  int pIndex;  
  pair< multimap<int, int>::iterator, multimap<int, int>::iterator > range;
  vector<int> index;
  index.reserve(tuples);
  for(int i=0; i < tuples; i++){
    pIndex = atoi(conn.GetValue(i, 0));
    if(mode !=1){  
      if(pIndex <= pSet->data.size()){
	index.push_back(pIndex);   
      }
    }else{
      range = pSet->ensemblProbeSetIndex.equal_range(pIndex);         // this is probably dangerous. should set up a simpler structure where we can just do a lookup.. 
      for(multimap<int, int>::iterator it=range.first; it != range.second; it++){
	if((*it).second <= pSet->data.size()){
	  index.push_back((*it).second);
	}
      }
    }
  }
  // and write to the socket..
  ostringstream description;
  description << "db (" << mode << ") : " << query;   
  writeIndex(index, description.str());
  
}

void ConnectionObject::doGenDBLookup(){
  // first parse the message and work out what we want to do... 
  vector<QString> words = splitString(lastMessage, '#');
  if(words.size() != 7){
    cerr << "doGenDBLookup words size is not 7 but : " << words.size() << "  returning void.. without doing anything.. " << endl;
    return;
  }
  // then try to get some numbers out of there.
  QString term = words[0];
  bool ok;    // painful conversions.... I would prefer a straight.. if
  int table = words[1].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get a table id : " << endl;
    return;
  }
  int upMargin = words[2].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get a upMargin.. " << endl;
    return;
  }
  int downMargin = words[3].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get a down margin " << endl;
    return;
  }
  int maxExpect = words[4].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get a max expect : " << endl;
    return;
  }
  int minLength = words[5].toInt(&ok);
  if(!ok){
    cerr << "Couldn't get a minLength " << endl;
    return;
  }
  float minMatch = words[6].toFloat(&ok);
  if(!ok){
    cerr << "Couldn't get a minMatch " << endl;
    return;
  }
  pSetThresholds.setThresholds(maxExpect, minLength, minMatch);    // ok.. 
  //////////// there has to be a prettier way of doing this.. maybe a wrapper or something, but that will probably not be that good either..
  /// anyway,, now we need to make a query string.. maybe using something like an ostringstream or something.. 
  char* escapedTerm = new char[term.length()*2 + 1];
  PQescapeString(escapedTerm, term.latin1(), term.length());    // and then use it when we need it..
  QString query;    // 
  QTextStream ts(query, 2);
  switch(table){
  case 0:
    ts << "select a.chromosome, a.start, a.stop, a.strand from ensembl_transcript_9 a, ensembl_genes_9 b, ensembl_annotation_9 c "
       << "where a.gene=b.gene and c.gene=a.gene and c.annotation ~* '" << escapedTerm << "' order by a.chromosome, a.start";
    break;
  case 1:
    ts << "select a.chromosome, a.begin, a.stop, a.strand from fantom_assemblies a, fantom_transcripts b, fantom_annotation c "
       << "where a.transcript=b.transcript and a.transcript=c.transcript and c.annotation ~* '" << escapedTerm << "' order by a.chromosome, a.begin";
    break;
  default :
    cerr << "No suc table in doGenDBLookup .. thingy.... " << endl;
  }
  // ok. we need to make a list of RegionSpecification  ,, and a vector.. 
  // the list first includes all genomic regions, whereas the vectore includes only the merged ones.. or something..
  // Ok..
  list<RegionSpecification> specList;
  vector<RegionSpecification> specVector;
  list<RegionSpecification>::iterator it;
  list<RegionSpecification>::iterator current;     // the one we are checking.. 
  list<RegionSpecification>::iterator pit;    // need two so we can safely erase elements while traversing the thingy.. complicated.. 
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "connection not good" << endl;
    return;
  }
  //cout << "Query " << endl << query << endl;
  if(!conn.Exec(query.latin1())){
    cerr << "command didn't work: " << query << endl;
    return;
  }
  int tuples = conn.Tuples();
  for(int i=0; i < tuples; i++){
    int begin = atoi(conn.GetValue(i, 1));
    int end = atoi(conn.GetValue(i, 2));
    int strand = atoi(conn.GetValue(i, 3));
    if(strand == 1){
      begin -= upMargin;
      end += downMargin;
    }else{
      begin -= downMargin;
      end += upMargin;
    }
    if(begin < 0){ begin = 1; }   // no negative numbers... 
    specList.push_back(RegionSpecification(conn.GetValue(i, 0), begin, end));
    //cout << "Making speclist with chromosome : " << conn.GetValue(i, 0) << "  start : " << begin << "\tend : " << end << endl;
  }
  /// and then for the clever bit.. well, maybe not that clever..
  while(specList.size()){
    current = specList.begin();
    it = specList.begin();
    it++;      // we don't do anything.. with the second one.. 
    while(it != specList.end()){
      if((*current).merge(*it)){      // current has been expanded to include it, which now needs to be removed..
	//cout << "\t\tMerged and inflated.. " << endl;
	pit = it;
	it++;
	specList.erase(pit);          // should be safe to erase after we have incremented the iterator.
	continue;
      }
      it++;
    }
    // and here, we push the list..
    specVector.push_back(*specList.begin());
    specList.erase(current);
  }
  /// and all regions should now be singletons of some sort...
  /// let's test this perhaps...
  //for(int i=0; i < specVector.size(); i++){
  //  cout << i << "\tchrom : " << specVector[i].chromosome << "\tbegin : " << specVector[i].begin << "\t end : " << specVector[i].end << endl;
  //}
  /// and now we do some magic..
  clientRegions = specVector;
  sApp("<clientRegionsChanged>");
  qiApp(specVector.size());
  sApp("<clientRegionsChangedEnd>");
  writeArray();       // which is very simple.. but now I need a function to send these regions. Well, I sort of have, don't I.. 
}
    
void ConnectionObject::writeIndex(vector<uint> v, string term, bool setClientIndex){      // ugly hack,, to allow certain things..
  vector<int> v2(v.size());
  for(uint i=0; i < v.size(); i++){
    v2[i] = (int)v[i];
  }
  writeIndex(v2, term, setClientIndex);
}

void ConnectionObject::writeIndex(vector<int> v, string term, bool setClientIndex){
  sApp(string("<ProbeIndex>"));
  vector<int> tv;
  tv.reserve(v.size());
  /////////// UGLY HACK WITH SOME GOOD AND SOME BAD PROPERTIES COMING UP. READ THE CODE AND WEEP
  /// Lots of ugly stuff coming from the fact that clientIndex is uint, and the other stuff is all vector<int>
  // this really should be sorted out, but would involve rewriting a lot of stuff.. 
  vector<uint> v2(v.size());
  
  for(uint i=0; i < v.size(); i++){
    if(v[i] > 0 && v[i] <= pSet->data.size() && pSet->data[v[i]-1]->index && clientChips.count(pSet->probeData[v[i]-1].chip)){
      tv.push_back(v[i]);
    }
    if(setClientIndex){   // which is default behaviour.. 
      v2[i] = (uint)v[i];
    }
  }
  if(setClientIndex){
    clientIndex = v2;
  }
  qiApp(tv.size());
  for(int i=0; i < tv.size(); i++){
    qiApp(tv[i]);
  }
  sApp(term);             // the reason for changing the index.. 
  sApp("<ProbeIndexEnd>");
  writeArray();    // packets up with parseable header..
}

void ConnectionObject::readClientIndex(){
  // not sure if it's good behaviour to check the clientChips when setting this, or nt to bother..
  // as I'm not sure I'm not going to implment the check.. 

  //cout << "readClientIndex() lastCommand length: " << lastCommand.length() << endl;
  if(lastMessage.length() == 0) { return; }
  bool ok;
  uint n;
  clientIndex.resize(0);
  int start=0;
  int end=lastMessage.find(":", start);   // for parsing purposes.
  //cout << "\tend : " << end;
  while(end != -1 && end != 0){
    QString num = lastMessage.mid(start, (end-start));
    n = num.toInt(&ok);
    n--;
    if(ok && n < pSet->data.size()){
      clientIndex.push_back(n);        // subtract one to get the proper indexing. (client uses 1 starting, vector uses 0);
      //cout << "pushing " << n << " into the client Index " << endl;
    }else{
      //cout << "problem parsing the clientIndex num: " << num << endl;
    }
    start = end+1;            // remove the : character..
    end = lastMessage.find(":", start);
  }
}

void ConnectionObject::setClientChips(){
  vector<QString> words = splitString(lastMessage, '|');
  clientChips.erase(clientChips.begin(), clientChips.end());
  bool ok;
  pSet->chipMutex->lock();
  for(uint i=0; i < words.size(); i++){
    int n = words[i].toInt(&ok);
    if(ok && pSet->chipDescriptions.count(n)){
      clientChips.insert(n);
    }else{
      cerr << "Some problem in assigning a chip with id : " << n << endl;
    }
  }
  pSet->chipMutex->unlock();
  // and write the resulting new index to the user..
  writeIndex(clientIndex, "clientChipsChanged", false);        // note this doesn't change the clientIndex. 
  // this may seem counterintuitive, but it means that if we remove a chip, and then we add it back, 
  // we should recover the lost ones. Note that the clientIndex has to be set each time that the client 
  // wants to do something on a set of indices. (we might want to check the chipInfo for this as well,
  // but it is not completely clear that we want to do that.
}
  
  

void ConnectionObject::writeProbeSet(probe_set* pset){  /// am ignoring the buffer size for now.. 

  sApp(string("<ProbeSet>"));
  qiApp(pset->index);
  qiApp(pset->exptSize);
  for(int i=0; i < pset->exptSize; i++){
    qiApp(pset->exptIndex[i]);
  }
  qiApp(pset->probeNo);
  for(int i=0; i < pset->probeNo; i++){
    qiApp(pset->exptSize);
    for(int j=0; j < pset->exptSize; j++){
      fApp(pset->probes[i][j]);
    }
  }
  // and let's write a terminator... oh lalala ... 
  sApp(string("<ProbeSetEnd>"));
  writeArray();
}


void ConnectionObject::doFlush(){
  sApp(string("<padding><paddingEnd>"));
  writeArray();
}

void ConnectionObject::writeFileInfo(){            //??????????????????????
}

// void dIWrite(QDataStream& ds, int n){
//   ds.writeRawBytes((const char*)&n, 4);
// }

// void dFWrite(QDataStream& ds, float f){
//   ds.writeRawBytes((const char*)&f, 4);
// }

// void dSWrite(QDataStream& ds, string s){
//   dIWrite(ds, s.size()+1);
//   ds.writeRawBytes(s.c_str(), s.size()+1);
// }

void ConnectionObject::writeExptInfo(){
  // serialise the map<float, exinfo> experiments      to a socket
  sApp(string("<ExptInfo>"));
  int size = pSet->experiments.size();
  qiApp(size);     // the casting may not be so good.. signed and unsigned.. 
  //map<float, exInfo>::iterator it;
  for(map<float, exInfo>::iterator it = pSet->experiments.begin(); it != pSet->experiments.end(); it++){
    fApp((*it).first);
    qiApp((*it).second.dbaseIndex);
    qiApp((*it).second.exptGroup);
    fApp((*it).second.index);                  /// This looks like it might be broken,,, bugger.. check it out.. later.. if I remember.. 
    sApp((*it).second.shortName);
    sApp((*it).second.description);
    qiApp((*it).second.chips.size());
    cout << (*it).first << "\t" << (*it).second.dbaseIndex << "\t" << (*it).second.exptGroup << "\t" << (*it).second.description << endl;
    map<int, bool>::iterator cit;
    for(cit = (*it).second.chips.begin(); cit != (*it).second.chips.end(); cit++){
      qiApp((*cit).first);
      sendData->app(*(char*)&((*cit).second));
    }
  }
  sApp(string("<ExptInfoEnd>"));
  writeArray();
}

void ConnectionObject::writeChipInfo(){
  sApp("<ChipInfo>");
  pSet->chipMutex->lock();
  qiApp(pSet->chipDescriptions.size());
  for(map<int, chipInfo>::iterator it=pSet->chipDescriptions.begin(); it != pSet->chipDescriptions.end(); it++){
    qiApp((*it).second.index);
    sApp((*it).second.id);
    sApp((*it).second.description);
    qiApp((*it).second.equivs.size());
    for(set<int>::iterator iit=(*it).second.equivs.begin(); iit != (*it).second.equivs.end(); iit++){
      qiApp(*iit);
    }
  }
  pSet->chipMutex->unlock();
  sApp("<ChipInfoEnd>");
  writeArray();
}


void ConnectionObject::writeArray(){
  //cout << "beginning of writeArray " << endl;
  int n = 0;
  int size = sendData->size();
  ostringstream header; 
  header << "#" << n << ":" << size << "#";
  string head = header.str();

  if(writen(head.c_str(), head.size()) < 0){
    cerr << "Failed to write header : writen returned less than 0, abort command " << endl;
    // stayConnected = false;
    return;
  }
  //cout << "just before writing the data : " << endl;
  if(writen(sendData->data, sendData->size()) < 0){
    cerr << "Failed to write data : writen returned less than 0, abort command " << endl;
    //    stayConnected = false;
    return;
  }
  //socket->writeBlock(sendData->data, sendData->size());
  //cout << "just after sending the data, and just before trying to flush" << endl;
  // if(connected()){
//     cout << "no error returned by qsocket device " << endl;
//     socket->flush();    // but this is never good ,,, 
//     cout << "after flushing" << endl;
//   }else{
//     stayConnected = false;
//   }
  // no flushing implemented on a qsocket device .. 
  //qApp->unlock();
  sendData->empty();       // just sets the counter back to 0.. the memory is still there. 
  //cout << "after sendData empty .. " << endl; 
 //cout << "Address of sendData " << (int)sendData << endl;
  //void (NetArray::*emptyAddress)() = &NetArray::empty;
  //cout << "Address of sendData empty function " << (int)&(sendData->*emptyAddress) << endl;
  // the question is if the address of the empty function is the same between different threads
  // which I suppose we can approach by printing it out.. 
  // cout << "\tAddress of sendData->empty() function: " << (int)&(sendData->empty) << endl;
  
}
    
ssize_t ConnectionObject::writen(const char* dptr, size_t n){
  size_t nleft = n;
  ssize_t nwritten;
  const char* ptr = dptr;
  //cout << "writen nleft is :" << nleft << endl;
  
  while(nleft > 0){
    //cout << "\tcalling write, nleft is " << nleft << endl;
    if( (nwritten = write(socketNumber, ptr, nleft)) <= 0) {
      // cout << "some error occured with errno " << errno << endl;
      if( errno == EINTR ){
	nwritten = 0;
	//cout << "\t\twriten received an interrupt, socket number is " << socketNumber << endl;
      }else{
	reportError("writen Received an error. Message :", errno);
	return(-1);
      }
    }
    nleft -= nwritten;
    ptr += nwritten;
    //cout << "\tnleft    is now : " << nleft << endl
    //	 << "\tnwritten is now : " << nwritten << endl;
  }
  return(n);
}

void ConnectionObject::reportError(string errorMessage, int errorNumber){
  // basically do a switch on the error number and then do cerr to give 
  // some message...
  size_t eLength = 100;   // ahh, typical ugly C string handling
  char stringError[eLength];
  if( strerror_r(errorNumber, stringError, eLength) ){
    cerr << errorMessage << "\n\t" << stringError << endl;
    return;
  }
  cerr << errorMessage << "\n\tUnknown error number" << endl;
}

void ConnectionObject::writeStatus(statusMessage message){
  sApp("<StatusMessage>");

 
  qiApp(message.id);   // hmmm, I wonder what happens if number to low, this should be fine though, it just writes the bytes.
  qiApp(message.no);
  int okInt = (int)message.ok;
  qiApp(okInt);    // anyway, it will go to an int.. 
  qiApp(message.errorMessages.size());
  for(int i=0; i < message.errorMessages.size(); i++){
    cerr << "error message is : " << message.errorMessages[i] << endl;
    sApp(message.errorMessages[i]);
  }
  qiApp(message.errorCodes.size());
  for(int i=0; i < message.errorCodes.size(); i++){
    qiApp(message.errorCodes[i]);
  }
  sApp("<StatusMessageEnd>");
  writeArray();
}

void ConnectionObject::writeDBChoices(){
  // just write the DB choices to the socket..
  vector<string> dbChoices; // very simple at the moment
  dbChoices.push_back("Annotation");
  dbChoices.push_back("Ensembl Annotation");
  dbChoices.push_back("Affymetrix Id");
  dbChoices.push_back("Gene Ontology");
  dbChoices.push_back("Ish Probe Id");
  // and send this to the client..
  //sendData.resize(0);
  sApp(string("<DBChoices>"));
  qiApp(1);
  qiApp(dbChoices.size());
  for(int i=0; i < dbChoices.size(); i++){
    sApp(dbChoices[i]);
  }
  sApp(string("<DBChoicesEnd>"));
  writeArray();
  //// and let's do region choises as well
  vector<string> regChoices;
  regChoices.push_back("Ensembl Genes");
  regChoices.push_back("Fantom");
  sApp(string("<DBChoices>"));
  qiApp(2);
  qiApp(regChoices.size());
  for(int i=0; i < regChoices.size(); i++){
    sApp(regChoices[i]);
  }
  sApp(string("<DBChoicesEnd>"));
  writeArray();
  // and maybe that's all that needs be done.. 
  // maybe this is so small that it doesn't get fed into the system, flush a couple of times..
  //doFlush();
  //doFlush();
  //doFlush();
}

void ConnectionObject::writeProbeData(probe_data* pdata){
  //cout << "writing Probe Data to the socket" << endl;
  // write the probe_data to the socket,,
  //if(!connected()) { return; }
  // write everything in the order given in the header file.
  // header..
  //sendData.resize(0);
  sApp(string("<ProbeData>"));
  qiApp(pdata->index);
  sApp(pdata->gbid);
  sApp(pdata->afid);
  qiApp(pdata->blastGuess);
  qiApp(pdata->ugData.size());
  
  for(int i=0; i < pdata->ugData.size(); i++){
    qiApp(pdata->ugData[i].index);
    sApp(pdata->ugData[i].title);
    sApp(pdata->ugData[i].gene);
  }

  ///////// THIS OBVIOUSLY COULD BENEFIT FROM SOME SORT OF READER/WRITE LOCKS, BUT FOR NOW KEEP THE CODE SIMPLE AND 
  ///////// SEE HOW IT GOES... --- 
  pSet->sessionMutex->lock();
  map<int, set<int> >::iterator slit=pSet->sessionLookup.find(pdata->index);
  map<int, sessionInformation>::iterator msit;
  if(slit != pSet->sessionLookup.end()){
    //    cout << "got a session information iterator and am going with it" << endl;
    //cout << "(*slit).first is " << (*slit).first << endl;
    //cout  << "Size is " << (*slit).second.size() << endl;
    qiApp((*slit).second.size());
    set<int>::iterator sit;
    for(sit = (*slit).second.begin(); sit != (*slit).second.end(); sit++){
      msit = pSet->sessions.find(*sit);
      qiApp((*msit).second.index);
      qiApp((*msit).second.owner);
      sApp((*msit).second.title);
      sApp((*msit).second.description);
      qiApp((*msit).second.keywords.size());
      for(vector<string>::iterator kwit = (*msit).second.keywords.begin(); kwit != (*msit).second.keywords.end(); kwit++){
	sApp(*kwit);
      }
      qiApp((*msit).second.members.size());
      for(set<int>::iterator memit = (*msit).second.members.begin(); memit != (*msit).second.members.end(); memit++){
	qiApp(*memit);
      }
    }
  }else{
    qiApp(0);
  }
  pSet->sessionMutex->unlock();
  /// let's see if that compiles.. !! 
  /// and lets do the same with the annotationInformation.. !! hwoo heayy...
  //cout << "looking up annotation information" << endl;
  pSet->annotationMutex->lock();
  map<int, set<int> >::iterator anit = pSet->annotationLookup.find(pdata->index);
  map<int, annotationInformation>::iterator infit;
  if(anit != pSet->annotationLookup.end()){
    //cout << "got an iterator not equal to the end " << endl;
    //cout << "size is : " << (*anit).second.size() << endl; 
    qiApp((*anit).second.size());
    for(set<int>::iterator sanit = (*anit).second.begin(); sanit != (*anit).second.end(); sanit++){
      infit = pSet->userAnnotation.find(*sanit);
      qiApp((*infit).second.index);
      qiApp((*infit).second.owner);
      sApp((*infit).second.annotation);
      qiApp((*infit).second.members.size());
      for(set<int>::iterator mit = (*infit).second.members.begin(); mit != (*infit).second.members.end(); mit++){
	qiApp(*mit);
      }
    }
  }else{
    qiApp(0);
  }
  pSet->annotationMutex->unlock();

  ///// which is pretty much everything for that.. 
  sApp(pdata->afdes);
  sApp(pdata->tigrDescription);
  qiApp(pdata->go.size());
  for(int i=0; i < pdata->go.size(); i++){
    qiApp(pdata->go[i].size());
    for(int j=0; j < pdata->go[i].size(); j++){
      sApp(pdata->go[i][j]);
    }
  }
  sApp(string("<ProbeDataEnd>"));
  writeArray();
  //cout << "\tFinished writing probe data to the socket" << endl; 
}

// bool ConnectionObject::connected(){
//   if(socket->error() == QSocketDevice::NoError){      // I need to change this somewhere.. 
//     return(true);
//   }
//   return(false);
// }

// Some statistical stuff::

void ConnectionObject::zScore(float* v, int s){  //s for the size.. 
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

void ConnectionObject::stdDeviation(float* v, uint s, float& mean, float& std, bool sampling){
  if(s < 2 && sampling || s < 1){
    return;                       // can't do this..
  }
  std = 0;
  float xsum = 0;
  float xsquaredSum = 0;
  for(int i=0; i < s; i++){
    xsum += v[i];
    xsquaredSum += (v[i] * v[i]);
  }
  mean = xsum / s;
  float SS = xsquaredSum - (xsum*xsum)/(float)s;
  if(sampling){
    std = sqrt(SS/(float)(s-1));
  }else{
    std = sqrt(SS/(float)(s));
  }
}

void ConnectionObject::squaredSum(float* v, uint s, float& mean, float& sqsum){
  if(s < 1){
    return;
  }
  mean = 0;
  sqsum = 0;
  for(int i=0; i < s; i++){
    mean += v[i];
    sqsum += (v[i] * v[i]);
  }
  //  cout << "\t\t\tsqsum : " << sqsum << "\t\tmean : " << mean << endl; 
  sqsum = (sqsum - (mean * mean)/(float)s);
  mean = mean / (float)s;
  //cout << "\t\t\tsqsum : " << sqsum << "\t\tmean : " << mean << endl; 
  // and that's it..
}
	

void ConnectionObject::doDiffSort(){
  // parse the lastMessage, and if it's ok, then run diffSort with the appropriate
  // parameters..
  //  cout << "doDiffSort: lastMessage: " << endl;
  int first = lastMessage.find(":", 0);
  int second = lastMessage.find(":", first+1);
  if(first == -1 || second == -1){ 
    return; }
  QString upString = lastMessage.left(first);
  QString downString = lastMessage.mid(first+1, second-(first+1));
  QString normString = lastMessage.right(1);   // just the last character.. 
  // convert to ints..
  bool u_ok, d_ok, b_ok;
  uint up = upString.toInt(&u_ok);
  uint down = downString.toInt(&d_ok);
  int normalised = normString.toInt(&b_ok);
  //  cout << "doDiffSort up: " << up << "\tdown: " << down << "\tnormed: " << normalised << endl;
  if(!u_ok || !d_ok || !b_ok){ 
    return; }
  if(normalised > 1 || normalised < 0) { 
    return; }          // the boolean identifier should be 1 or 0..
  bool normed;
  if(normalised == 0){ normed = false; }else{ normed = true; }
  //  cout << "calling diffSort up: " << up << endl;
  diffSort(up, down, normed);
}

void ConnectionObject::diffSort(uint up, uint down, bool normalised){
  // go trough and sort by the diff score... use values normalised across the whole 
  // experimental set if normalised is true..
  // (but for the implementation just pass it on to the diff function..)..
  //  cout << "beginning of diffSort clientIndex.size: " << clientIndex.size() << endl;
  if(clientIndex.size() == 0){
    return;
  }
  vector<dist_set> scores(clientIndex.size());
  for(int i=0; i < clientIndex.size(); i++){
    scores[i].index = pSet->data[clientIndex[i]]->index;
    scores[i].value = diff(pSet->data[clientIndex[i]], up, down, normalised);
  }
  //  cout << "got the scores.. " << endl;
  sort(scores.begin(), scores.end(), r_comp_set());     // sort biggest first..
  vector<int> tIndex(scores.size());
  for(int i=0; i < scores.size(); i++){
    tIndex[i] = scores[i].index;
  }
  //  cout << "just before calling writeIndex" << endl; 
  writeIndex(tIndex, "diffSort");
}

float** ConnectionObject::copyProbes(float** p, int ps, int es){
  float** cp = new float*[ps];
  for(int i =0; i < ps; i++){
    cp[i] = new float[es];
    for(int j=0; j < es; j++){
      cp[i][j] = p[i][j];
    }
  }
  return(cp);
}

void ConnectionObject::delProbes(float** p, int ps){
  for(int i=0; i < ps; i++){
    delete []p[i];
  }
  delete []p;
}

float ConnectionObject::diff(probe_set* p, uint up, uint down, bool normalised){
  // do a simple comparison,, of the two points.. 
  // returns -1 on error.
  float score = 0;
  float error = -1;
  if(p->probeNo < 1){ return(error); }
  if(up > p->allExptNo || down > p->allExptNo) { return(error); }
  // determine if the numbers are actually defined,, -through the exptLookup,,
  if(p->exptLookup[up] != -1 && p->exptLookup[down] != -1){
    up = p->exptLookup[up];
    down = p->exptLookup[down];
  }else{
    return(error);
  }
  // SS = (sum(X^2)) - ((sum(X)^2)/N
  // and variance = SS/N-1
  float mean = 0; 
  float sq_sum =0;
  float sum_sq = 0;
  float delta;
  float** probes = copyProbes(p->probes, p->probeNo, p->exptSize);
  //  float** probes = new float*[p->probeNo];
  //for(int i=0; i < p->probeNo; i++){
  // probes[i] = new float[p->exptSize];
  // for(int j=0; j < p->exptSize; j++){
  //   probes[i][j] = p->probes[i][j];
  // }
  //}   /// hmmmm, there should be a better way of doing that. -- can't I just do *probes ?? 

  //vector< vector<float> > probes = p->probes;
  if(normalised){
    for(int i=0; i < p->probeNo; i++){
      zScore(probes[i], p->exptSize);
    }
  }
  for(int i=0; i < p->probeNo; i++){
    delta = probes[i][up] - probes[i][down];
    mean += delta;
    sq_sum += (delta * delta);
  }
  sum_sq = (mean*mean);
  mean = mean/p->probeNo;
  //  mean = mean/probes.size();
  float SS = sq_sum - (sum_sq/p->probeNo);
  float std = sqrt(SS/(p->probeNo-1));
  // and delete probes..
  delProbes(probes, p->probeNo);
  //  for(int i=0; i < p->probeNo; i++){
  //  delete []probes[i];
  //}
  //delete []probes;     // mmm, maybe.. 
  return(mean/std);
}


float ConnectionObject::euclidean(float* v1, float* v2, uint s){
  //cout << "beginning of euclidean s is " << s << endl;
  float distance = 0;
  //  if(v1.size() != v2.size() || v1.size() == 0) { 
  //  cout << "bugger, BUGGER,, euclidean functions, the vectors are of different size, returning -1" << endl;
  //  return(-1); 
  //}
  if(s < 1){ return(-1);}
  for(int i=0; i < s; i++){
    distance += ((v1[i]-v2[i]) * (v1[i]-v2[i]));
    //    distance += pow((v1[i]-v2[i]), 2);
  }
  // NOTE PLEASE,, note that distance is now the square of the euclidean distance, and 
  // for this to be real, I should just return the square root,, this however is a little troublesome,
  // as the larger the number of dimensions, the larger the distance, and this would then screw things up
  // if I compare distances with many measurements and distances with few. So I just divide by the number of 
  // dimensions. Maybe this should be the number of degrees of freedome, or something, but I don't know!!
  //cout << "euclidean function returning square root of " << distance << " / " << s << endl;
  return((sqrt(distance))/s);
}


float ConnectionObject::anova(probe_set* pset, uint* expts, uint es){
  // calculates the anova score across th  experimental values given in the expts variable..
  // this is essentiall equal to the (variance between groups/degrees of freedom) / (variance within groups/degrees of freedom)..
  // equations and stuff obtained from http://davidmlabe.com/hyperstat/B918101
  
  
  float t_mean=0;       // the mean of all of the values..
  // groups are essentially experimental time points..
  // so get the mean.. for all ..
  int counter = 0;

  float** probes = copyProbes(pset->probes, pset->probeNo, pset->exptSize);
  //  vector< vector<float> > probes = pset->probes; // then normalise it..
  for(int i=0; i < pset->probeNo; i++){
    zScore(probes[i], pset->exptSize);           
  }

  // OK slight complication.. the int* expts contain the experiment indices that we want to 
  // look at. These do not directly map to the indices of the probe set vectors, and indeed may 
  // not be represented at all. This means that we have to create a second index that is specifice
  // to this one by using the exptLookup map

  //   OK, this isn't entirely true. the other option is just to check for presence, i.e. doesn't equal -1
  //   in the int* exptLookup, but we would have to do this for each probe set, i.e. 16 times most of the time
  //   so it seems that it might be faster to make a temporary thingy, and check properly in that one.. 

  uint* eindex = new uint[es];
  int selEx = 0;      // the selected experiments.. 
  for(int i=0; i < es; i++){
    //if(expts[i] < pset->exptSize){   /// but this should surely be allExptNo, not exptSs
    if(expts[i] < pset->allExptNo){
      if(pset->exptLookup[expts[i]] != -1){
	eindex[selEx] = pset->exptLookup[expts[i]];
	selEx++;
      }
    }
  }
  //  vector<int> eindex;
  //eindex.reserve(expts.size());       // may be smaller but cannot be larger..
  //map<int, int>::iterator it;
  //for(int i=0; i < expts.size(); i++){
  //   it = pset->exptLookup.find(expts[i]);
  //  if(it != pset->exptLookup.end()){
  //    eindex.push_back((*it).second);
  //  }
  //}
    
  for(int i=0; i < pset->probeNo; i++){
    for(int j=0; j < selEx; j++){
      //cout << "expts : " << expts[j] << "   probes: " << i << endl;
      //      if(expts[j] < probes[i].size()){             // i.e. we are not overstepping ourselves..
      t_mean += probes[i][eindex[j]];
      counter++;
	// }
    }
  }
  t_mean = t_mean/counter;
  counter=0;          // then we can reuse it.. 
  // calculate the mean for each experimental time point... // SO THERE IS A MUCH BETTER WAY of doing this, but I'm not going to have time to check it in one hour..
  //  int probesUsed = probes.size();  // we don't have any excluded at the moment.
  uint probesUsed = pset->probeNo;  // we don't have any excluded at the moment.

  //vector<float> exp_means(expts.size());
  float* exp_means = new float[selEx];
  for(int i=0; i < selEx; i++){
    exp_means[i] = 0;
    for(int j=0; j < pset->probeNo; j++){
      exp_means[i] += probes[j][eindex[i]];
    }
    //    exp_means[i] = exp_means[i]/probes.size();
    exp_means[i] = exp_means[i]/probesUsed;
  }
  // and then.. for each experimental time point calculate the variance within the group... 
  int er_df =0;
  float er_variance = 0;
  float var = 0;
  for(int i=0; i < selEx; i++){
    er_variance = 0;
    for(int j=0; j < pset->probeNo; j++){
      var += ((probes[j][eindex[i]]-exp_means[i]) * (probes[j][eindex[i]]-exp_means[i]));
      //      var += pow((probes[j][eindex[i]]-exp_means[i]), 2);
    }
    //    er_variance += (var/probes.size());
    er_variance += (var/probesUsed);
    er_df += (probesUsed-1);
    //    er_df += (probes.size()-1);
  }
  // and then find the variance between the different experimental points..
  int b_df= es-1;           // -- between degrees of freedom//
  //  int b_df= expts.size()-1;           // -- between degrees of freedom//
  float b_variance = 0;
  for(int i=0; i < selEx; i++){
    b_variance += ((exp_means[i]-t_mean) * (exp_means[i]-t_mean));
    //    b_variance += pow((exp_means[i]-t_mean), 2);
  }
  b_variance = b_variance/es;
  //  b_variance = b_variance/expts.size();
  // so I should be able to just calculate the anova scrore at this point....
  float a_score = (b_variance/(float)b_df)/(er_variance/(float)er_df);
  delProbes(probes, pset->probeNo);
  delete []eindex;
  delete []exp_means;
  return(a_score);
}

void ConnectionObject::anovaSelectSort(){
  if(clientIndex.size() == 0){ return; }
  vector<uint> expts; // simpler to create.. -then convert to an int*
  bool ok;  // for conversions.
  int start = 0;
  int end = lastMessage.find(":", start);
  int ei;
  QString temp;
  while(end != -1 && end != 0){
    temp = lastMessage.mid(start, (end-start));
    start = end+1;
    end = lastMessage.find(":", start);
    // ok and lets see if we can get a reasonable int out of this on
    ei = temp.toInt(&ok);
    if(ok){
      expts.push_back(ei);
    }
  }
  if(expts.size() == 0){ return; }
  uint* exptArray = new uint[expts.size()]; // ugly
  for(int i=0; i < expts.size(); i++) { exptArray[i] = expts[i]; } // should be a better way, but I don't know it.
  // lets make a couple of AnovaProcessors..
  vector<dist_set> scores(clientIndex.size());
  //////////////////// REMEMBER .. clientIndex has already had 1 subtracted to make sure that
  //////////////////// that it fits the pset-Data.. something..
  
  // split clientIndex into 2 parts, and work out the starts and stops. 
  uint mid = clientIndex.size()/2;      // note this is used as the size in the loop, and as such, will not be duplicated. 
  //  cout << "the processors started .. " << endl;
  AnovaProcessor* process1 = new AnovaProcessor(pSet, 0, mid, &clientIndex, &scores, exptArray, expts.size());
  AnovaProcessor* process2 = new AnovaProcessor(pSet, mid, clientIndex.size(), &clientIndex, &scores, exptArray, expts.size());

  process1->start();
  process2->start();     // hmm, I wonder.. 
  process1->wait();
  process2->wait();

  //cout << "the processors finished .. " << endl;
  // and now everything is just as before..
  sort(scores.begin(), scores.end(), r_comp_set());
  vector<int> tIndex;
  //  tIndex.reserve(scores.size());
  for(int i=0; i < scores.size(); i++){
    if(scores[i].value > 0){
      tIndex.push_back(scores[i].index);
    }
  }
  //  for(rit = scores.rbegin(); rit != scores.rend(); rit++){
  //  tIndex.push_back((*rit).second);
  //}
  delete process1;
  delete process2;
  delete exptArray;  // don't forget !!. 
  writeIndex(tIndex, "Anova (Select)");
}

void ConnectionObject::doKCluster(){
  // create a KClusterObject, put this into the set<kclusterProcesses*> and then go back.
  // should be OK..
  cout << "at the beginning of kCluster " << endl;
  vector<QString> words = splitString(lastMessage, '|');
  // ok, what information do I need.. 
  // 1. k -- the number of clusters to create
  // 2. meanNorm..          -- normalise mean values.. 
  // 3. localNorm   boolean, -base normalisation on local normalisation or global normalisation. 
  // 4. individual Norm..   -- normalise individual members before taking the mean..
  // 5. experimental selection  .. -which experiments to base the clustering on..
  //
  // 4 and 5 are kind of strange options, and for now the kClusterProcess assumes they are both true anyway. 
  // 3 is also a little strange, but I think it can be argued that it might be useful to do global normalisation 
  // from time to time.. 
  
  // use the current Client Index to sort this out.... 
  if(words.size() < 6){           // doesn't make any sense to cluster based on only one experiment.. 
    return;
  }
  bool ok;
  uint k = words[0].toInt(&ok);
  if(!ok){
    return;
  }
  int b = words[1].toInt(&ok);
  if(!ok){
    return;
  }
  bool localNorm = (bool)b;
  b = words[2].toInt(&ok);
  if(!ok){
    return;
  }
  bool individualNorm = (bool)b;
  b = words[3].toInt(&ok);
  if(!ok){
    return;
  }
  uint e;
  bool meanNorm = (bool)b;
  uint N = words.size()-4;
  uint* expts = new uint[words.size()-4];
  for(int i=4; i < words.size(); i++){
    e = words[i].toInt(&ok);
    if(!ok){
      delete []expts;
      return;
    }
    expts[i-4] = e;
  }
  cout << "just before creating the KClusterProcess " << endl;
  cout << " N is " << N << endl;
  // send an int* for the client index instead..
  uint* clientI = new uint[clientIndex.size()];
  for(int i=0; i < clientIndex.size(); i++){ clientI[i] = clientIndex[i]; }
  
  // and now just make the cluster with the stuff in it.. 
  KClusterProcess* cluster = new KClusterProcess(k, clientI, clientIndex.size(), &pSet->data, expts, N, localNorm, &clusterProcesses, &clusterMutex);
  //  KClusterProcess* cluster = new KClusterProcess(k, clientI, clientIndex.size(), pSet, expts, N, localNorm, &clusterSets, &clusterMutex);
  //clusterProcesses.insert(cluster);
  cluster->start();        // and then we can leave.. hooo yeahh.. 
  //cluster->wait();
}

void ConnectionObject::sendKClusters(){
  // go through the clusterSets and send the information in some sort of reasonable manner.. remove all of the ones..
  // then delete the clusterSets.. too much information to keep around.. really..
  cout << "Finished Kclustering will now delete all the information in the cluster Set " << endl
       << "if it's appropriate " << endl;
  set<void*>::iterator cvit;
  //  set<KClusterProcess*>::iterator psit;
  clusterMutex.lock();
  while(clusterProcesses.size()){
    cvit = clusterProcesses.begin();
    //  for(cvit = clusterProcesses.begin(); cvit != clusterProcesses.end(); cvit++){
    KClusterProcess* tempP = (KClusterProcess*)(*cvit);
    tempP->wait();
    //cout << "We have now waited for the cluster process and can safely delete it" << endl;
    sApp("<KClusters>");
    // first whether localNorm or not..
    qiApp((int)tempP->localNorm);
    // first send the number of experiments...
    qiApp(tempP->N);
    for(int i=0; i < tempP->N; i++){
      qiApp(tempP->selectedExperiments[i]);
    }
    // then the number of clusters
    qiApp(tempP->k);
    // and the cluster centers, and the members. 
    for(int i=0; i < tempP->k; i++){
      // first the center
      for(int j=0; j < tempP->N; j++){
	fApp(tempP->centers[i][j]);
      }
      // then the number of cluster members..
      qiApp(tempP->clusterSizes[i]);
      for(int j=0; j < tempP->clusterSizes[i]; j++){
	qiApp(tempP->probeIndices[tempP->clusters[i][j]]);      // the index of the member
	for(int k=0; k < tempP->N; k++){
	  fApp(tempP->points[tempP->clusters[i][j]][k]);
	}
      }
    }
    sApp("<KClusterEnd>");
    writeArray();
    doFlush();
    clusterProcesses.erase(cvit);
    delete tempP;
  }
  clusterMutex.unlock();
      
  //  clusterMutex.lock();
  //while(clusterSets.size()){
  //  clusterSets.erase(clusterSets.begin());
  //}
  //clusterMutex.unlock();
//   set<KClusterProcess*>::iterator cpit;
//   // and 
//   set<clusterSet*> tempClusters;
//   set<KClusterProcess*> tempProcs;
//   clusterMutex.lock();    // 
//   for(csit = clusterSets.begin(); csit != clusterSets.end(); csit++){
//     tempClusters.insert(*csit);
//     clusterSets.erase(csit);
//   }
//   clusterMutex.unlock();
//   // now. let's see if we can find the cluster processes and delete those..
//   for(cpit = clusterProcesses.begin(); cpit != clusterProcesses.end(); cpit++){
//     if(tempClusters.count((*cpit)->clusters)){
//       tempProcs.insert(*cpit);
//       clusterProcesses.erase(cpit);
//     }
//   }
//   ///  now we have a potential problem.. but I think it should be ok. It's messy though.
//   ///  send the data from tempClusters.. and then delete everything.. OK.. 
//   // if I delete the process does that delete the cluster ?? shouldn't do, but I'm not entirely sure, but we will find out I guess.
//   // delete the processes..
//   for(cpit = tempProcs.begin(); cpit != tempProcs.end(); cpit++){
//     cout << "deleting a process how is it going, crash, crash and burn" << endl;
//     (*cpit)->wait();
//     //delete *cpit;
//     tempProcs.erase(cpit);
//   }
//   for(csit = tempClusters.begin(); csit != tempClusters.end(); csit++){
//     cout << "Sending one cluster how is it looking? " << endl;
//   }
//   for(csit = tempClusters.begin(); csit != tempClusters.end(); csit++){
//     //    delete *csit;
//     tempClusters.erase(csit);
//   }
}

void ConnectionObject::traceExperiments(){
  cout << "Tracing the experiments.. " << endl;
  // although it is terribly inefficient to keep rewriting the code for getting the experiment indices from the 
  // message, I keep doing this, such that if I want to change the protocol to get other things, it will be easy. Not good, 
  // but,,, there is bound to be a better structural way for this somewhere, but..
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() < 4){
    cerr << "connection object, trace experiments, word size is only " << words.size() << " this is not very useful is it? " << endl;
    return;
  }
  QString sigmaString = words.back();
  words.pop_back();
  bool ok;
  float sigmaMultiplier = sigmaString.toFloat(&ok);
  if(!ok){
    cerr << "Couldn't get a sigma multiplier from : " << sigmaString << "\treturning from traceExperiments.. " << endl;
    return;
  }
  vector<uint> expts;
  for(uint i=0; i < words.size(); i++){
    int n = words[i].toInt(&ok);
    if(ok){
      expts.push_back(n);
    }else{
      cerr << "connection object trace Experiments : one of the experiments doesn't look like an int returning without action" << endl;
      return;
    }
  }
  // and I need to collect a vector of probe set thingies as well..
  vector<probe_set*> ps;
  ps.reserve(clientIndex.size());
  for(uint i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size()){
      ps.push_back(pSet->data[clientIndex[i]]);
    }
  }
  // and we can now make the thingy without too much problem.. 
  cout << "making the experiment tracer " << endl;
  ExperimentTracer tracer(ps, expts, &exptTracers, &exptTracerMutex, sigmaMultiplier);
  cout << "experimnent tracer made, and we're off somewhere else.. " << endl;
  // which will take care of itself as far as I know.. 
}

void ConnectionObject::doExperimentCompare(){
  vector<QString> words = splitString(lastMessage, '|');
  // all we are looking for is a series of experiments, at least two..
  if(words.size() < 2){
    cerr << "doExperiemntCompare, no point comparing " << words.size() << "  experiments is there ? " << endl;
    return;
  }
  vector<uint> expts;   
  for(uint i=0; i < words.size(); i++){
    bool ok;
    int n = words[i].toInt(&ok);
    if(ok){
      expts.push_back(n);
    }else{
      cerr << "doExperimentCompare: one of the experiments doesn't look like an int returning without action" << endl;
      return;
    }
  }
  // now make a couple of data structs more useful threaded stuff,, 
  //cout << "just before making data structs for experiments and other things.. " << endl;
  uint* experiments = new uint[expts.size()];
  for(int i=0; i < expts.size(); i++){
    experiments[i] = expts[i];
  }
  uint* genes = new uint[clientIndex.size()];
  for(int i=0; i < clientIndex.size(); i++){
    genes[i] = clientIndex[i];
  }
  cout << "just before making the processor" << endl;
  /// make the processor object..
  ExperimentCompareProcess* processor = new ExperimentCompareProcess(genes, clientIndex.size(), &pSet->data, experiments, expts.size(), &compareExperimentProcesses, &experimentComparisonMutex);
  //cout << "just before calling start " << endl;
  processor->start();
  //cout << "started the ExperimentCompareProcess : should do something with the pointer and things, I wonder if there is an auto delete function somewhere" << endl;
}

void ConnectionObject::sendExperimentDistances(){
  //cout << "Sending Experimental Distances " << endl;
  set<void*>::iterator it;
  experimentComparisonMutex.lock();     // the processes must finish first..
  while(compareExperimentProcesses.size()){
    it = compareExperimentProcesses.begin();
    //  for(it = compareExperimentProcesses.begin(); it != compareExperimentProcesses.end(); it++){
    ExperimentCompareProcess* temp = (ExperimentCompareProcess*)(*it);
    temp->wait();  
    sApp("<ExperimentDistances>");
    qiApp(temp->probeCounter);    // number of genes included.. 
    qiApp(temp->exptNo);
    for(int i=0; i < temp->exptNo; i++){
      qiApp(temp->experiments[i]);
      for(int j=i+1; j < temp->exptNo; j++){ 
	//qiApp(temp->experiments[j]);   // this is sort of redundant we can work out things.. 
	fApp(temp->distances[i][j]);
      }
    }
    sApp("<ExperimentDistancesEnd>");
    writeArray();
    compareExperimentProcesses.erase(it);
    delete temp;
  }
  experimentComparisonMutex.unlock();
}

void ConnectionObject::sendExperimentTrace(PathTracer* tracer){
  pointLink* pl = tracer->chain();
  //cout << "Sending a trace.. " << endl;
  sApp("<ExperimentTrace>");
  while(pl){
    qiApp(1);
    npoint* p = pl->point;
    cout << "id: " << p->id << "\td: " << p->d << "\tx: " << p->x << "\ty: " << p->y << endl;
    qiApp(p->id);
    fApp(p->d);   // the distance to the next point..
    fApp(p->x);
    fApp(p->y);
    pl = pl->next;
  }
  qiApp(0);
  sApp("<ExperimentTraceEnd>");
  writeArray();
}
    

void ConnectionObject::doFlatExptCompare(){
  vector<QString> words = splitString(lastMessage, '|');
  // all we are looking for is a series of experiments, at least two..
  if(words.size() < 4){
    cerr << "doFlatExperiemntCompare, no point comparing " << words.size()-2 << "  experiments is there ? " << endl;
    return;
  }
  float sigma, order;
  bool fok1, fok2;
  sigma = words[0].toFloat(&fok1);
  order = words[1].toFloat(&fok2);
  if(!fok1 || !fok2){
    cerr << "Couldn't get order or sigma from the message,, " << endl;
    return;
  }
  vector<uint> expts;   
  for(uint i=2; i < words.size(); i++){
    bool ok;
    int n = words[i].toInt(&ok);
    if(ok){
      expts.push_back(n);
    }else{
      cerr << "doFloatCompare: one of the experiments doesn't look like an int returning without action" << endl;
      return;
    }
  }
  // now make a couple of data structs more useful threaded stuff,, 
  //cout << "just before making data structs for experiments and other things.. " << endl;
  uint* experiments = new uint[expts.size()];
  for(int i=0; i < expts.size(); i++){
    experiments[i] = expts[i];
  }
  uint* genes = new uint[clientIndex.size()];
  for(int i=0; i < clientIndex.size(); i++){
    genes[i] = clientIndex[i];
  }
  //cout << "just before making the processor" << endl;
  FlatExptCompare* process = new FlatExptCompare(genes, clientIndex.size(), &pSet->data, experiments, expts.size(), sigma, order, 
						 &flatExptComparers, &flatExptCompareMutex);
  process->start();
}

float ConnectionObject::doFlatExptCompare(int* a, uint as, int* b, uint bs, float* v, uint vs, float order, float sigma){
  // a experiments a, as a size
  // b experiments b, bs b size
  // v the values, already normalised. -we'll need to find the minimum, but we shouldn't modify thise
  // vs the size of the values..
  // order -the order of the sigmoid curve used
  // sigma - size of sigma.. (i.e. spread of the curve..)
  if(sigma <= 0 || !as || !bs){
    cerr << "doFlatExptCompare thingy, sigma is less than or equal to 0, doesn't make sense " << endl;
    return(0);
  }
  // get the minimum mean value..
  float minMean = 0;    // we know the minimum will be less than 0, so this is OK.
  for(int i=0; i < vs; i++){
    if(minMean > v[i]){ minMean = v[i]; }
  }
  //cout << "\t\tminMean is : " << minMean << endl;
  // then go through the values and calculate the delta value;.
  float aMean = 0;
  float bMean = 0;
  for(int i=0; i < as; i++){
    //cout << "\t\t\ta: " << a[i] << "\tvalue :" << v[a[i]] << endl;
    float tv =  (v[a[i]] - minMean);    // a temporary variable to make the next line more readable..
    aMean += (pow(tv/sigma, order))/ (1 + pow(tv/sigma, order));
  }
  //cout << "\t\tas : " << as << "\taMean : " << aMean << endl;
  aMean = aMean/(float)as;
  for(int i=0; i < bs; i++){
    float tv = v[b[i]] - minMean;
    //cout << "\t\t\tb: " << b[i] << "\tvalue :" << v[b[i]] << endl;
    //cout << "\t\t\tsigma : " << sigma << "\torder : " << order << "\ttv : " << tv << endl;
    bMean += (pow((tv/sigma), order)/ (1 + pow(tv/sigma, order)));
  }
  //cout << "\t\tbs : " << bs << "\tbMean : " << bMean << endl;
  bMean = bMean/(float)bs;
  //cout << "\t\taMean : " << aMean << "\tbMean: " << bMean << endl;
  // and return the difference between bMean and aMean.. -very simple indeed.    --making it absolute.. well, it loses information, but means we don't have to worry about 
  // other stuff. 
  //return(fabs(aMean - bMean));
  return(aMean-bMean);
}
  

void ConnectionObject::sendFlatExptDistances(){
  cout << "Sending Norm Experimental Distances " << endl;
  set<void*>::iterator it;
  flatExptCompareMutex.lock();     // the processes must finish first..
  while(flatExptComparers.size()){
    //cout << "The first comparer being sent.. " << endl;
    it = flatExptComparers.begin();
    //  for(it = flatExptComparers.begin(); it != flatExptComparers.end(); it++){
    FlatExptCompare* temp = (FlatExptCompare*)(*it);
    temp->wait();  
    //cout << "wait returned should be OK " << endl;
    sApp("<FlatExperimentDistances>");  
    qiApp(temp->probeCounter);
    fApp(temp->sigma);
    fApp(temp->order);
    qiApp(temp->exptNo);
    for(int i=0; i < temp->exptNo; i++){
      qiApp(temp->experiments[i]);
      for(int j=i+1; j < temp->exptNo; j++){ 
	//qiApp(temp->experiments[j]);   // this is sort of redundant we can work out things.. 
	fApp(temp->distances[i][j]);
      }
    }
    sApp("<FlatExperimentDistancesEnd>");
    writeArray();
    flatExptComparers.erase(it);
    //    cout << "calling delete on the comparer" << endl;
    delete temp;
    //cout << "delete called on the comparer " << endl;
  }
  flatExptCompareMutex.unlock();
}
  
void ConnectionObject::devsFromMean(){
  //cout << "This is devs From Mean, and we are getting serious.. " << endl;
  vector<QString> words = splitString(lastMessage, '|');
  // hmm, we just need a set of experiments to do this with, right..
  vector<unsigned int> expts;
  bool ok;
  //cout << "And words size is " << words.size() << endl;
  for(int i=0; i < words.size(); i++){
    expts.push_back(words[i].toUInt(&ok));
    if(!ok){
      cerr << "We couldn't get what we were looking for" << endl;
      return;
    }
    //cout << "Added " << expts.back() << "  to expts " << endl;
  }
  //  cout << "devs From Mean, and expts. size is : " << expts.size() << endl;
  // because threads don't like vectors..
  unsigned int* e = new unsigned int[expts.size()];
  for(int i=0; i < expts.size(); i++){ e[i] = expts[i]; }
  // and let's make a couple of things
  DataExtractor de;
  ProbeStats ps;   // for calculating the things..
  // and go through evertying .. 
  vector<float> scores;            // this will be huge.. I'm sure..
  for(int i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size()){
      ExData* exd = de.extract(pSet->data[clientIndex[i]], e, expts.size());
      if(!exd->exptNo){
	cerr << "Didn't get an exd data set from data extractor.. " << endl;
	delete exd;
	continue;
      }
      float* v = ps.devFromMean(exd);
      if(!v){
	cerr << "Didn't get a set of devs from mean.. " << endl;
	delete exd;
	continue;
      }
      for(int j=0; j < exd->probeNo; j++){
	scores.push_back(v[j]);
      }
      delete []v;
      delete exd;
    }
  }
  // at which point I need to send the stats.. in some way or another.. 
  // cout << "devs from mean scores size is : " << scores.size() << endl;
  if(scores.size()){
    sApp("<statCollection>");
    qiApp(1);
    sApp("P Dev From Mean");
    qiApp(scores.size());
    for(int i=0; i < scores.size(); i++){
      qiApp(-1);
      fApp(scores[i]);
    }
    sApp("<statCollectionEnd>");
    writeArray();
  }else{
    cerr << "Didn't get any scores for some strange reason " << endl;
  }
  delete []e;
}

void ConnectionObject::exportMeans(){
  // there are lots of options by which we can do this, but first let's just do this for the 
  // the selected experiment sets and the selected index, as this should be simpler..
  vector<QString> words = splitString(lastMessage, '|');
  // hmm, we just need a set of experiments to do this with, right..
  vector<unsigned int> expts;
  bool ok;
  if(words.size() < 5){     // we don't allow export of less as this gives us trouble with normalisation.. 
    cerr << "Export Means couldn't find any experiments to include.. " << endl;
    return;
  }
  string fileName = words[0];
  int normMethod = words[1].toInt(&ok);
  if(!ok){
    cerr << "Export means coulnd't get a normalisation method from word " << words[1] << endl;
    return;
  }
  // we know about methods 1, 2, 3.. so 
  if(!(normMethod < 4 && normMethod > 0)){
    cerr << "Export Means : Unknown normalisation method .. " << endl;
    return;
  }
  for(int i=2; i < words.size(); i++){
    expts.push_back(words[i].toUInt(&ok));
    if(!ok){
      cerr << "Export Means couldn't parse an experiment index. " << endl;
      return;
    }
  }
  //  cout << "devs From Mean, and expts. size is : " << expts.size() << endl;
  // because threads don't like vectors..
  unsigned int* e = new unsigned int[expts.size()];
  for(int i=0; i < expts.size(); i++){ e[i] = expts[i]; }
  // and let's make a couple of things
  DataExtractor de;
  Normaliser norm;     // for normalisation.. 
  // and then let's go through the values and see how it goes.. 
  sApp("<ExportMeans>");
  sApp(fileName);     // so the client knows what to do. 
  // doing it this way is potentially a security problem as a malicious server could try to write files to 
  // somewhere other than the client intended. -- Probably good for the client to remember what it requested
  // and refuse to write to anywhere else.
  ostringstream header;
  header << "probe_index";
  for(int i=0; i < expts.size(); i++){
    header << "\t" << expts[i];
  }
  sApp(header.str());
  for(int i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size()){
      ExData* exd = de.extract(pSet->data[clientIndex[i]], e, expts.size(), true);
      if(!exd->exptNo){
	delete exd;
	continue;
      }
      // since we are exporting to text we are just going to send one string per line..
      // which all the client has to do is stick into some file..
      // just end by having the end marker <exportDataEnd> or something.. should be quite simple
      
      // First normalise all the values in my exData..
      // we should allow the option of normalising using 

      switch(normMethod){
      case 1:
	norm.zScore(exd->values, exd->probeNo, exd->exptNo);
	break;
      case 2:
	norm.mScore(exd->values, exd->probeNo, exd->exptNo);
	break;
      default :
	cerr << "no normalisation performed in export: " << endl;
      }
      // get a mean if possible..
      float* meanValues = exd->mean();
      if(!meanValues){
	cerr << "ExData mean() returned 0 in export Mean function, bugger " << endl;
	delete exd;
	continue;
      }
      // make a line.. 
      ostringstream line;
      line << pSet->data[clientIndex[i]]->index;    // make sure that this is indeed the index.
      for(int j=0; j < exd->exptNo; j++){
	line << "\t" << meanValues[j];
      }
      // and just write the string..
      sApp(line.str());   // which should be fine..
      delete exd;
      delete []meanValues;
    }
  }
  sApp("<ExportMeansEnd>");
  writeArray();
}
  

void ConnectionObject::parseGenomicRegionRequest(){
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() < 4){
    return;
  }
  // and then split the first word
  bool probeSet = false;
  bool ensGenes = false;
  // and we can increase the number..
  vector<QString> types = splitString(words[0], ',');    // get the types..
  for(int i=0; i < types.size(); i++){
    if(types[i] == "ps"){ probeSet = true; }
    if(types[i] == "egene"){ ensGenes = true; }
  }
  // and then we can do some stuff depending on the above.
  string chrom = words[1].latin1();
  bool ok1, ok2;
  int start = words[2].toInt(&ok1);
  int end = words[3].toInt(&ok2);
  if(!ok1 || !ok2 || start > end){
    cerr << "couldn't get good numbers, or inverted range" << endl;
    return;
  }
  // and now lets just send the bugger..
  if(probeSet){
    sendGenomicRegionProbeSetMatches(chrom, start, end);
  }
  if(ensGenes){
    int requestId = 0;
    if(words.size() > 4){
      bool ok3;
      requestId = words[4].toInt(&ok3);
      if(!ok3){
	//cout << "couldn't get a request id.. from " << words[4] << endl;
	requestId = 0;
      }
    }
    sendGenomicRegionEnsemblGenes(requestId, chrom, start, end);
  }
  // and as I don't have any methods for anything else forget..
  // -- how to handle combinataions.. I'm not sure of yet anyway..
}
  
void ConnectionObject::sendIshThumbnails(){
  // first get the ish probe index from the lastmessage
  bool ok;
  int regionSize = 200000;    /// static default for now. fix later.. 
  int index = lastMessage.toInt(&ok);
  if(!ok){
    return;
  }
  ///////////// seems that we cannot actually do the whole thing with just a large object but that we must make
  ///////////// another connection to the database as well..
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "sendIshThumbNails connection bad : " << conn.ErrorMessage() << endl;
    return;
  }
  ostringstream query;    // probably not threadsafe, but for god's sake do I have to write my own class for everything ??
  query << "select index, thumbnail from ish_images where probe = " << index; 
  if(!conn.Exec(query.str().c_str())){
    cerr << "Couldn't select thumbnail from ish_images, not sure why.. " << endl 
	 << query.str() << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  // send the probe data first so that the client can put some appropriate labels on the things.. 
  sendIshProbeData(index);     
  for(int i=0; i < conn.Tuples(); i++){
    //cout << "Sending an Image using the sendIshThumbnail function " << endl;
    sendIshThumbnail(index, atoi(conn.GetValue(i, 0)), atoi(conn.GetValue(i, 1))  );
    //// note by doing this in a seperate function we are forced to create a new large object connection to the database for each 
    //// thumbnail image. This is probably slower than creating one, and then using this one large object by opening and closing..
    ///  the objects. However, should we have some sort of corrupted object and we don't recover correctly, then we may have an easier time
    ///  with sepearte objects... -- There is a problem at the moment though, which is that the client programme doens't know how many objects 
    /// it is likely to receive.. -- and it might be an idea to send a little message to the client telling it the number of images which 
    /// will be sent..
    // perhaps later..
  }


  // now let's look up the data for this probe and see what genomic regions we can send to the client.. well, you never know, it might look good.. eh..
  //cout << "Sending the IshImageRegions, but I think there may be lock somewhere that's locked at the moment " << endl;
  sendIshImageRegions(index, regionSize);     // well acutally that's regions, but never mind.. 
 //  map<int, ishProbeData>::iterator it = pSet->ishProbes.find(index);
//   if(it != pSet->ishProbes.end()){
//     sApp("<sendingRegions>");
//     qiApp(2);      // the target
//     qiApp(index);    // the index
//     sApp("<sendingRegionsEnd>");
//     writeArray();
//     vector<ishProbeMatchSet>::iterator vit;
//     for(vit =(*it).second.probeMatches.begin(); vit != (*it).second.probeMatches.end(); vit++){
//       int range = (*vit).maxPos - (*vit).minPos;
//       int margin = 100;
//       if(range < regionSize){
// 	margin = (regionSize - range)/2;
//       }
//       sendGenomicRegion((*vit).chromosome, (*vit).minPos-margin, (*vit).maxPos+margin, index, 2);
//     }
//     sApp("<finishedSendingRegions>");
//     qiApp(2);      // the target
//     qiApp(index);    // the index
//     sApp("<finishedSendingRegionsEnd>");
//     writeArray();
//   }
  
}

void ConnectionObject::parseIshProbeRegionRequest(){
  bool ok;
  int regionSize = 200000;
  int index = lastMessage.toInt(&ok);
  if(!ok){
    return;
  }
  sendIshProbeData(index);      // it's a little redundant.. but it's the easiest way at the moment.. -- client should make an explicit request
  sendIshImageRegions(index, regionSize);
}

void ConnectionObject::sendIshImageRegions(int index, int regionSize){
  //cout << "beginning of sendIshImageRegions " << endl;
  pSet->ishProbeDataMutex->lock();
  //cout << "managed to lock the mutex " << endl;
  map<int, ishProbeData>::iterator it = pSet->ishProbes.find(index);
  if(it != pSet->ishProbes.end()){
    sApp("<sendingRegions>");
    qiApp(2);      // the target
    qiApp(index);    // the index
    sApp("<sendingRegionsEnd>");
    writeArray();
    vector<ishProbeMatchSet*>::iterator vit;
    for(vit =(*it).second.probeMatches.begin(); vit != (*it).second.probeMatches.end(); vit++){
      int range = (*vit)->maxPos - (*vit)->minPos;
      int margin = 100;
      if(range < regionSize){
	margin = (regionSize - range)/2;
      }
      sendGenomicRegion((*vit)->chromosome, (*vit)->minPos-margin, (*vit)->maxPos+margin, index, 2);
    }
    pSet->ishProbeDataMutex->unlock();
    sApp("<finishedSendingRegions>");
    qiApp(2);      // the target
    qiApp(index);    // the index
    sApp("<finishedSendingRegionsEnd>");
    writeArray();
  }else{
    pSet->ishProbeDataMutex->unlock();
  }
}

void ConnectionObject::sendIshProbeData(int index){
  cout << "beginning of sendIshProbeData " << endl;
  pSet->ishProbeDataMutex->lock();
  map<int, ishProbeData>::iterator it = pSet->ishProbes.find(index);      // thread safe or not.. probably I should lock the mutex before ..
  if(it == pSet->ishProbes.end()){
    cerr << "No such ishProbe data for index : " << index << endl;
    pSet->ishProbeDataMutex->unlock();
    return;
  }
  ishProbeData* d = &(*it).second;
  sApp("<ishProbeData>");
  qiApp(index);
  sApp(d->consensusSequence);
  sApp(d->antisensePromoter);
  qiApp(d->afIdAssociation);
  sApp(d->vectorName);
  qiApp(d->designLength);
  sApp(d->probeName);
  sApp(d->probeIdentifier);
  qiApp(d->ensemblGuess);
  //// and the annotation
  qiApp(d->textAnnotation.size());
  map<int, ish_annotation>::iterator mit;
  for(mit = d->textAnnotation.begin(); mit != d->textAnnotation.end(); mit++){
    sApp((*mit).second.fieldName);
    appendIshAnnotation((*mit).second);
  }
  qiApp(d->numberAnnotation.size());
  for(mit = d->numberAnnotation.begin(); mit != d->numberAnnotation.end(); mit++){
    sApp((*mit).second.fieldName);
    appendIshAnnotation((*mit).second);
  }
  qiApp(d->classification.size());
  for(mit = d->classification.begin(); mit != d->classification.end(); mit++){
    sApp((*mit).second.fieldName);
    appendIshAnnotation((*mit).second);
  }
  /// and unlock the mutex..
  pSet->ishProbeDataMutex->unlock();          // necessary since I am reading from what is potentially a dynamic structure.. -- very unlikely to cause 
                                              // a problem, but nevertheless possible.. -- I need to something like this for the other ones as well. 
  sApp("<ishProbeDataEnd>");
  writeArray();
  //cout << "end of sendIshProbeData " << endl;
}

void ConnectionObject::sendIshTextFields(){
  set<string>::iterator it;
  pSet->ishProbeDataMutex->lock();
  sApp("<ishTextFields>");
  qiApp(pSet->ishTextFields.size());
  for(it = pSet->ishTextFields.begin(); it != pSet->ishTextFields.end(); it++){
    sApp((*it));
  }
  pSet->ishProbeDataMutex->unlock();
  sApp("<ishTextFieldsEnd>");
  writeArray();
}

void ConnectionObject::sendIshFloatFields(){
  set<string>::iterator it;
  pSet->ishProbeDataMutex->lock();
  sApp("<ishFloatFields>");
  qiApp(pSet->ishFloatFields.size());
  for(it = pSet->ishFloatFields.begin(); it != pSet->ishFloatFields.end(); it++){
    sApp((*it));
  }
  pSet->ishProbeDataMutex->unlock();
  sApp("<ishFloatFieldsEnd>");
  writeArray();
}

void ConnectionObject::sendIshClasses(){
  set<string>::iterator it;
  pSet->ishProbeDataMutex->lock();
  sApp("<ishClasses>");
  qiApp(pSet->ishClasses.size());
  for(it = pSet->ishClasses.begin(); it != pSet->ishClasses.end(); it++){
    sApp(*it);
  }
  pSet->ishProbeDataMutex->unlock();
  sApp("<ishClassesEnd>");
  writeArray();
}

void ConnectionObject::sendProtocols(){
  /// create a protocol collection, and send this.. 
  ProtocolCollection collection(conninfo);   // this does all the database lookups and things..
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 2){
    cerr << "sendProtocols words size should be 2 but is : " << words.size() << endl;
    return;
  }
  sApp("<ProtocolCollection>");
  bool ok;
  int requesterId = words[0].toInt(&ok);
  if(ok){
    qiApp(requesterId);
  }else{
    qiApp(0);
  }
  int idrequest = words[1].toInt(&ok);
  if(ok){
    qiApp(idrequest);
  }else{
    qiApp(0);
  }
  collection.serialise(sendData);
  sApp("<ProtocolCollectionEnd>");
  writeArray();
}

void ConnectionObject::sendProtocol(){
  bool ok;
  int id = lastMessage.toInt(&ok);
  if(!ok){
    return;
  }
  Protocol protocol(conninfo, id, true);
  sApp("<Protocol>");
  //cout << "calling serialise " << endl;
  protocol.serialise(sendData, true);
  //cout << "serialise returned : " << endl;
  sApp("<ProtocolEnd>");
  writeArray();
  //cout << "wrote Array " << endl;
}

void ConnectionObject::receiveNewProtocol(){
  cout << "receiveNewprotococol" << endl;
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() < 5){
    cerr << "receiveNewProtocol, not enough words need at least 4 have only : " << words.size() << endl;
    return;
  }
  // first get the various bits and pieces...
  bool ok;
  int requestId = words[0].toInt(&ok);
  if(!ok){ 
    cerr << "coulnd't get a request Id : " << endl;
    return; }
  int parentId = words[1].toInt(&ok);
  if(!ok){ 
    cerr << "couldn't get a parent id : " << endl;
    return; }
  string protocolName = words[2].latin1();
  string protocolDescription = words[3].latin1();
  int stepNo = words[4].toInt(&ok);
  if(!ok){
    cerr << "receiveNewProtocol couldn't get a stepNo from the lastMessage .. " << endl;
    return;
  }
  int wc = 5;
  if(words.size() < (5 + (4 * stepNo)) ){
    cerr << "receiveNewProtocol words size is too small need at least : " << (4 + (4 * stepNo)) << "  for " << stepNo <<  "  steps. Have only : " << words.size() << endl;
    return;
  }
  vector<ProtocolStep> steps;
  steps.reserve(stepNo);
  for(int i=0; i < stepNo; i++){
    bool ok1, ok2, ok3;
    int stepId = words[wc++].toInt(&ok1);      // if this is 0, then it is a new step.. 
    int parentId = words[wc++].toInt(&ok2);
    string description = words[wc++].latin1();
    int creatorId = words[wc++].toInt(&ok3);   // if this is 0, then set to userId  
    if(!creatorId){ creatorId = userId; }
    if(ok1 && ok2 && ok3){
      steps.push_back(ProtocolStep(stepId, parentId, creatorId, description, conninfo));
    }
  }
  //cout << "about to make the protocol.. " << endl;
  /// at this point I should be able to create a protocol,, and ask for it to go into the database..
  Protocol protocol(steps, parentId, userId, currentUserName, protocolName, protocolDescription, conninfo, true);
  // true at the end means it gets automatically entered, as long as it has the appropriate parent id and all this.. OK la.. 
  //cout << "after the protocol made " << endl;
  // At this point we should send the protocol.. back for the guy to load with id, and time of creation and stuff like that..
  sApp("<Protocol>");
  qiApp(-1);          // indicate this is a new protocol..
  qiApp(requestId);
  protocol.serialise(sendData, true);
  sApp("<ProtocolEnd>");
  writeArray();
  //cout << "New protocol made and copy sent to the client.. OK la.. " << endl;
}

void ConnectionObject::commitExperimentToDB(){
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 5){
    cerr << "commit Experiment to DB, words size is not so great.. should be 5 but is : " << words.size() << endl;
    return;
  }
  bool rok, uok, tok, pok;
  int requestId = words[0].toInt(&rok);
  int userid = words[1].toInt(&uok);    // member value is userId... 
  uint time = words[2].toUInt(&tok);
  int protocol = words[3].toInt(&pok);
  if(!(rok && uok && tok && pok)){
    cerr << "commit experiment to db, couldn't get one of the integer variables, giving up " << endl;
    return;
  }
  string comment = words[4].latin1();    // and let's make an object,, with the various things.. 
                                         // and let that object enter itself if it can, if it can't 
                                         // then it just doesn't get the thingy.. OK.. 
  if(userId != userid){
    cerr << "Conflicting user id's better not do anything.. " << endl;
    return;
  }
  Experiment expt(userId, time, protocol, comment, conninfo);
  // which should do everything,, -- if it is set appropriately.. the user id is not -1
  if(expt.isEntered()){
    //cout << "sending the experiment details to the client" << endl;
    sApp("<Experiment>");
    qiApp(requestId);
    expt.serialise(sendData);
    sApp("<ExperimentEnd>");
  }else{
    cerr << "Ehh, what happendned mate " << endl;
  }
  writeArray();
  /// voila .. 
}

void ConnectionObject::sendExperiments(){
  ExperimentCollection collection(conninfo);
  sApp("<ExperimentCollection>");
  collection.serialise(sendData);
  sApp("<ExperimentCollectionEnd>");
  writeArray();
}

void ConnectionObject::sendTissues(){
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "sendTissues connection is bad " << endl;
    return;
  }
  ostringstream query;
  query << "select * from ish_tissue";
  if(!conn.Exec(query.str().c_str())){
    cerr << "sendTissues couldn't do query " << endl
	 << conn.ErrorMessage();
    return;
  }
  sApp("<Tissues>");
  qiApp(conn.Tuples());
  for(int i=0; i < conn.Tuples(); i++){
    qiApp(atoi(conn.GetValue(i, 0)));
    sApp(conn.GetValue(i, 1));
    fApp(atof(conn.GetValue(i, 2)));
  }
  sApp("<TissuesEnd>");
  writeArray();
}
	 
void ConnectionObject::sendIshAnnotationFields(){
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "sendTissues connection is bad " << endl;
    return;
  }
  const char* query = "select * from ish_annotation_fields";
  if(!conn.Exec(query)){
    cerr << "coudn't execute ish annotation field query " << endl;
    return;
  }
  sApp("<IshAnnotationFields>");
  qiApp(conn.Tuples());
  for(int i=0; i < conn.Tuples(); i++){
    qiApp(atoi(conn.GetValue(i, 0)));
    sApp(conn.GetValue(i, 1));
  }
  sApp("<IshAnnotationFieldsEnd>");
writeArray();
}

void ConnectionObject::makeIshTissue(){
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() != 2){
    cerr << "MakeIshTissue : Words size is wrong should be 2 but is : " << words.size() << endl;
    return;
  }
  bool ok;
  float age = words[1].toFloat(&ok);  // it's a little redundat as we will convert to text again.. oh well.. 
  if(!ok){
    cerr << "MakeTissue Coulnd't make a float from " << words[1] << endl;
    return;
  }
  char* cleanTissue = new char[words[0].length() * 2 + 1];
  PQescapeString(cleanTissue, words[0].latin1(), words[0].length());
  // make a connection to the database backend..
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "makeIshTissue, db connection bad" << endl;
    return;
  }
  ostringstream entry;
  entry << "insert into ish_tissue (name, age) values ('" << cleanTissue << "', " << age << ")";
  cout << "MakeTissue : " << entry.str() << endl;
  delete cleanTissue;
  if(!conn.Exec(entry.str().c_str())){
    cerr << "makeTissue couldn't enter into database, bugger : " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  // rather than checking let's just force an update by calling sendTissues..
  sendTissues();
}

void ConnectionObject::makeIshAnnotationField(){
  // firs make a string for the field name..
  char* field = new char[lastMessage.length()*2 + 1];
  PQescapeString(field, lastMessage.latin1(), lastMessage.length());
  ostringstream entry;
  entry << "insert into ish_annotation_fields (field_name) values ('" << field << "')";
  cout << "makeIshAnnotatinField : " << entry.str() << endl;
  delete field;
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "makeIshAnnotationField couldn't connect to database : " << endl;
    return;
  }
  if(!conn.Exec(entry.str().c_str())){
    cerr << "MakeIshAnnotationField couldn't enter field " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  sendIshAnnotationFields();
}
    
void ConnectionObject::commitIshImageToDB(){
  /// arggh, this is going to be painful..
  vector<QString> words = splitString(lastMessage, '|');  // as usual...

  statusMessage status; // set identifier below..
  // the size we want is ..
  if(words.size() < 8){    // size depends on the number of comments.. 
    cerr << "commitIshImageToDB words size should be larger than 8 but is : " << words.size() << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, wrong word number");
    writeStatus(status);
    return;
  }
  bool ok;
  int probeId = words[0].toInt(&ok);
  if(!ok){
    cerr << "commitIshImageToDB couldn't get probe id" << endl; 
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get probe id");
    writeStatus(status);
    return;
  }
  QString promoter = words[1];
  int experimentId = words[2].toInt(&ok);
  if(!ok){
    cerr << "commitIshImageToDB coulnd't get experimentId" << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get experiment id");
    writeStatus(status);
    return;
  }
  int tissueId = words[3].toInt(&ok);
  if(!ok){
    cerr << "commitIshImageToDB coulnd't get tissueId" << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get tissue id");
    writeStatus(status);
    return;
  }
  QString fileName = words[4];
  int commentNo = words[5].toInt(&ok);
  if(!ok){
    cerr << "couldn't parse the number of comments for ish image data " << endl;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get comment no");
    status.ok = false;
    writeStatus(status);
    return;
  }
  //// words size should be, 7 + commentNo * 2, check that this is correct.. 
  if(words.size() != (8 + commentNo*2)){
    cerr << "ish image data comments, comment No is " << commentNo << "  so words size should be " << (8 + commentNo*2) << "  but in fact is " << words.size() << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, wrong word number after getting commentNo");    
    writeStatus(status);
    return;
  }
  int ci = 6;   // the current index that we need to use..
  vector<Comment> comments(commentNo);
  for(int i=0; i < commentNo; i++){
    int commentField = words[ci].toInt(&ok);
    if(!ok){
      cerr << "Couldn't get commentfield " << endl;
      status.ok = false;
      status.errorMessages.push_back("commitIshImageToDB function, couldn't get experiment id");
      writeStatus(status);
      return;
    }
    ci++;
    QString comment = words[ci];
    ci++;
    comments[i] = Comment(commentField, comment);
  }
  int imageSize = words[ci].toInt(&ok);
  if(!ok){
    cerr << "commitIshImageToDB couldn't get the size of the image.. bugger " << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get image size");
    writeStatus(status);
    return;
  }
  // and last should be the identifier of the request..
  // better to put this first, but as I'm adding on ere..
  ci++;
  unsigned int identifier = words[ci].toUInt(&ok);
  if(!ok){
    cerr << "Could not get identifier for the thingy. This is a bit trouble some as I really can't do anything at all " << endl;
    status.ok = false;
    status.errorMessages.push_back("commitIshImageToDB function, couldn't get identifier");
    writeStatus(status);
    return;
  }
  // make a status message with this
  status.id = identifier;
  //  statusMessage status(identifier);        // later we can change this if I need to.. set bool ok to false if some problem.. 
  // and now we'll just read the data from the netArray.. and see what we get..
  char* imageData = inData->readBinaryFromSocket(socketNumber, imageSize);
  if(!imageData){
    cerr << "Couldn't get image Data : " << endl;
    status.ok = false;
    status.errorMessages.push_back("Couldn't read the image data from the stream");
    writeStatus(status);
    return;
  }
  // then we need some data as to how big the next file is...
  char* sizeInfo = inData->readBinaryFromSocket(socketNumber, 12);    // read 3 ints, each one 4 bytes.. 
  if(!sizeInfo){
    cerr << "couldn't get size info .. " << endl;
    delete []imageData;
    status.ok = false;
    status.errorMessages.push_back("Couldn't get the appropriate size information for the thumbnail image");
    writeStatus(status);
    return;
  }
  int first = ntohl(*(uint*)sizeInfo);
  int second = ntohl(*(uint*)(sizeInfo + 4));
  int third = ntohl(*(uint*)(sizeInfo + 8));         // so this isn't very good as it assumes the size of things.. 
  if(second != (third/2) || first != 0 || third <= 0){
    cerr << "Thumbnail data size specification appears wrong : " << endl
	 << "first : " << first << endl << "second : " << second << endl
	 << "third : " << third << endl;
    status.ok = false;
    status.errorMessages.push_back("Thumbnail data size specification appears inconsistent");
    writeStatus(status);
    return;
  }
  uint thumbSize = third;
  char* thumbData = inData->readBinaryFromSocket(socketNumber, thumbSize);
  if(!thumbData){
    delete []imageData;
    delete []sizeInfo;
    status.ok = false;
    status.errorMessages.push_back("Couldn't read thumbnail image from stream");
    writeStatus(status);
    cerr << "coulnd't get thumbdata .. returning .. " << endl;
    return;
  }
  // which should be ok,, let's see if that works before we start trying to stick things into the database..
  //cout << "Got the data for the thumb " << thumbSize << "  bytes " << endl
  //    << "And got the main image data " << imageSize << "  bytes " << endl;

  /// and now we have to get some pqlargeobjects and work out how to use them.. 
  PgLargeObject imageObject(conninfo);
  PgLargeObject thumbObject(conninfo);
  if(imageObject.ConnectionBad()){
    cerr << "couldn't make large object connection " << endl;
    delete []imageData;
    delete []thumbData;
    delete []sizeInfo;
    status.ok = false;
    status.errorMessages.push_back("Couldn't create a large object connection");
    writeStatus(status);
    return;
  }
  imageObject.Exec("begin");
  thumbObject.Exec("begin");
  cerr << "After exec begin error Message : " << imageObject.ErrorMessage() << endl;
  
  // cout << "now the image object has an oid of : " << imageObject.LOid() << endl;
  // create an object ..
  imageObject.Create();
  Oid imageOid = imageObject.LOid();
  //cout << "created and now the imageOid  is " << imageOid << endl;
  //cerr << "after checking the oid error Message : " << imageObject.ErrorMessage() << endl;

  imageObject.Open();        // does what .. 
  int imageWritten = imageObject.Write(imageData, imageSize);
  cerr << "after writing the first time .. error Message : " << imageObject.ErrorMessage() << endl;

  // then create and write the data to the other thingy.. 
  thumbObject.Create();
  thumbObject.Open();      // ?? 
  Oid thumbOid = thumbObject.LOid();
  //cout << "created and now the imageOid  is " << imageOid << endl;
  int thumbWritten = thumbObject.Write(thumbData, thumbSize);

  delete []imageData;
  delete []thumbData;
  delete []sizeInfo;
  if(thumbWritten != thumbSize || imageWritten != imageSize){
    cerr << "It seems that we didn't write enough for the large objects bugger " << endl;
    cerr << "imageWritten is : " << imageWritten << endl
	 << "thumbWritten is : " << thumbWritten << endl;
    cerr << "error Message : " << imageObject.ErrorMessage() << endl;
    status.ok = false;
    status.errorMessages.push_back("Problems writing the image data into the database");
    writeStatus(status);
    return;
  }
  ostringstream entry;    // now we need to make an entry for the table..
  char* cleanFile = new char[fileName.length()*2 + 1];
  PQescapeString(cleanFile, fileName.latin1(), fileName.length());
  char* cleanProm = new char[promoter.length()*2 +1];
  PQescapeString(cleanProm, promoter.latin1(), promoter.length());
  entry << "insert into ish_images (tissue, file_name, promoter, image, thumbnail, probe, experiment, user_id) "
	<< "values (" << tissueId << ", '" << cleanFile << "', '" << cleanProm << "', "
	<< imageOid << ", " << thumbOid << ", " << probeId << ", " << experimentId << ", " << userId << ")";
  delete cleanFile;
  delete cleanProm;
  if(imageObject.ExecCommandOk(entry.str().c_str()) != 1){
    cerr << "couln't enter the stuff into the ish images table " << endl; 
    status.ok = false;
    status.errorMessages.push_back("Couldn't enter data into the ish images table");
    writeStatus(status);
    return;
  }
  //if(imageObject.Tuples() != 1){
  // cerr << "tuples number is wrong... should be 1 but is " << imageObject.Tuples() << endl;
  // return;;
  //}
  // if here let's commit..
  if(!imageObject.Exec("commit")){
    cerr << "Couldn't commit the images to the database.. " << endl;
    status.ok = false;
    status.errorMessages.push_back("Couldn't commit the image to the database");
  }
  if(!thumbObject.Exec("commit")){
    cerr << "couldn't exec for the thumb object, bugger " << endl;
    status.ok = false;
    status.errorMessages.push_back("Couldn't commit the thumbnail image to the database");
  }
  /// take care of the comments as well,,, but maybe that is most of the stuff done...
  // this is a real pain, but I need to find out what the index of the image we just added is. Fortunately, I can use the imageOid as this 
  // should be unique. But this is a real pain, to do, and I'm sure that there must be a better way of doing this. Probably this can all be
  // written as an SQL function,  -- could do with some better SQL programming skills.
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "commit Ish Image to DB, -- I'm giving up on adding comments as I can't get a connection to check the image index " << endl;
    status.ok = false;
    status.errorMessages.push_back("Failed to write the annotation to the database due to connection failure");
    writeStatus(status);
    return;
  }
  ostringstream check;
  check << "select index from ish_images where image = " << imageOid;
  if(!conn.Exec(check.str().c_str())){
    cerr << "Couldn't do a select on the images data base ,, so can't add the comments, very sorry.. " << endl;
    status.ok = false;
    status.errorMessages.push_back("Failed to write the annotation to the database due to failure to obtain the image oid");
    writeStatus(status);
    return;
  }
  if(!conn.Tuples() == 1){
    cerr << "Commit ish image to database, -- I coulnd't find the index for the image, so can't add comments.." << endl;
    status.ok = false;
    status.errorMessages.push_back("Failed to write the annotation to the database due to failure to obtain the image oid");
    writeStatus(status);
    return;
  }
  int imageIndex = atoi(conn.GetValue(0, 0));
  if(commitIshImageCommentsToDB(comments, imageIndex)){
    cout << "Ish Image Annotation succesfully added to the database " << endl;
  }else{
    status.ok = false;
    status.errorMessages.push_back("Failed to write the annotation to the database for unknown reasons");
  }
  writeStatus(status);
}
  
bool ConnectionObject::commitIshImageCommentsToDB(vector<Comment> comments, int imId){
  // 
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "CommitIshImageComments to DB couldn't get a datbase connection, nothing added, very sorry.. " << endl;
    return(false);
  }
  int counter = 0;
  for(int i=0; i < comments.size(); i++){
    char* cleanComment = new char[comments[i].comment.length() * 2 + 1];
    PQescapeString(cleanComment, comments[i].comment.latin1(), comments[i].comment.length());
    ostringstream entry;
    entry << "insert into ish_annotation (user_index, image, field, annotation) values (" << userId << ", " << imId << ", " << comments[i].index << ", '"
	  << cleanComment << "')";
    //cout << "entry string : " << entry.str() << endl;
    if(!conn.Exec(entry.str().c_str())){
      cerr << "Couldn't insert the comment " << comments[i].comment.latin1() << endl;
      cerr << conn.ErrorMessage() << endl;
    }
    if(conn.CmdTuples() != 1){
      cerr << "Didn't insert, not sure why, but error message follows " << endl
	   << conn.ErrorMessage() << endl;
    }else{
      counter += conn.CmdTuples();
    }
    delete []cleanComment;
  }
  cout << "commitIshImageCommentsToDB inserted a total of  " << counter << "  lines into table ish_annotation " << endl;
  return(true);
}

void ConnectionObject::doBlast(){
  QString sequence(">testSequence initially from hes7\nATGGTCACCCGGGATCGAGCTGAGAATAGGGACGGCCCCAAGATGCTCAAGCCGCTTGTGGAGAAGCGGCGCCGGGACCGCATCAACCGCAGCCTGGAAGAGCTGAGGCTGCTGCTGCTGGAGCGGACCCGGGACCAGAACCTCCGGAACCCGAAGCTGGAGAAAGCGGAGATATTGGAGTTCGCCGTGGGCTACTTGAGGGAGCGAAGCCGGGTGGAGCCCCCGGGGGTTCCCCGGTCCCCAGTCCAGGACGCCAAGGCGCTCGCCAGCTGCTACTTGTCCGGTTTCCGCGAGTGCCTGCTTCGCTTGGCGGCCATCGCGCACGACGCCAGCCCGACCGCCCGCGCCCAGCTCTTCTCCGCGCTGCACGGCTATCTGCGCCCCAAACCGCCCCGGCCCAAGCCGGTAGATCCGAGGCCTCCAGCGCCGCGCCCATCCCTGGACCCCGCCGCACCGGCCCTTGGCCCTGCGCTGCACCAGCGCCCCCCAGTGCACCAGGGCCACCCTAGCCCGCGCTGCGCATGGTCCCCATCCCTCTGCTCCCCGCGCGCCGGGGATTCTGGCGCGCCGGCGCCCCTCACCGGACTGCTGCCGCCGCCACCGCCGCCTCACAGACAAGACGGGGCGCCCAAGGCCCCGCTGCCCCCGCCGCCCGCTTTCTGGAGACCTTGGCCCTGA");
  BlastClient* blastClient = new BlastClient(sequence);
  blastClients.insert(blastClient);
  blastClient->start();
  ////////  which gives us a memory leak, as we have no way of knowing when the blast client is done.. and hence when to delete it
  ///////   later fix this in some way or other.. not too difficult..
}

void ConnectionObject::appendIshAnnotation(ish_annotation& a){    // no need to copy it..
  qiApp(a.annotation_id);
  qiApp(a.userId);
  sApp(a.userName);
  if(a.numerical){
    fApp(a.value);
  }else{
    sApp(a.annotation);
  }
}
  

void ConnectionObject::sendIshThumbnail(int index, int imageIndex, int oid){
  // index is the ish_probe_index,, -- we need to tell the client what this is.. after all.. 
  PgLargeObject lo(oid, conninfo);
  if(lo.ConnectionBad()){
    cerr << "Couldn't connect to large object oid : " << index << endl
	 << lo.ErrorMessage() << endl;
  }
  
  lo.Exec("BEGIN");    // do I need this ?? I'm not actually sure,, but it may be useful.. 
  lo.Open();
  int seekValue = lo.LSeek(0, SEEK_END);
  int end = lo.Tell();
  //cout << "After seeking to the end. The seekValue returned is : " << seekValue << "  and the tell value returned is " << end << endl;
  if(end < 0){
    cerr << "lo.LSeek returns negative value, return without reading " << endl;
    lo.Exec("COMMIT");      // maybe unnecessary but seems 
    lo.Close();              // presumably closes on destruct, but what the hell. 
    return;
  }
  int returnValue = lo.LSeek(0, 0);
  int begin = lo.Tell();
  //cout << "After seeking back to the beginning. The value returned by seek is : " << returnValue << "  and tell says " << begin << endl;
  // length of the buffer should be equal to end.. so let's make a buffer for reading into..
  char* buf = new char[end];
  int bytesRead = lo.Read(buf, end);
  lo.Exec("COMMIT");  // ?? 
  //cout << "read " << bytesRead << "  bytes into buffer" << endl;
  // and now let's make a message containing buffer,, and send it.. then delete the buffer or something.
  sApp("<ishThumbNail>");
  qiApp(index);   // the index of the probe..
  qiApp(imageIndex);  // the index of the image,, so we know which one it is.. 
  //qiApp(end);     // the length of the data.. // not needed the sApp function anyway provides a function for this.. 
  sApp(buf, end);  // just appends the buffer ,, -don't forget to delete it... 
  sApp("<ishThumbNailEnd>");
  delete buf;
  lo.Close();    // hopefully that is OK.. 
  writeArray();
  ///
}

void ConnectionObject::sendIshImage(){
  // first get the index.. from the last message..
  // this index should be the database image line index.. 
  bool ok;
  int index = lastMessage.toInt(&ok);
  if(!ok){
    cerr << "coulnd't get an index from the last message (sendIshImage), message : " << endl; //lastMessage << endl;
    return;
  }
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "sendIshImage : connection bad" << endl << conn.ErrorMessage() << endl;
    return;
  }
  ostringstream query;
  query << "select probe, image from ish_images where index=" << index;
  if(!conn.Exec(query.str().c_str())){
    cerr << "couldn't query for the image " << endl << conn.ErrorMessage() << endl;
    return;
  }
  if(!conn.Tuples()){
    cerr << "query didnt' return any images getting out of here" << endl;
    return;
  }
  int probeIndex = atoi(conn.GetValue(0, 0));
  int oid = atoi(conn.GetValue(0, 1));
  // ok.. now we have to make a large object and fill it with stuff..
  PgLargeObject lo(oid, conninfo);
  if(lo.ConnectionBad()){
    cerr << "couldn't connect to large object oid : " << oid << endl
	 << lo.ErrorMessage() << endl;
    return;
  }
  lo.Exec("BEGIN");
  lo.Open();
  int end = lo.LSeek(0, SEEK_END);
  //cout << "sendIshImage after seeking to end end is : " << end << endl;
  // and seek back to the beginning..
  int begin = lo.LSeek(0, 0);
  if(end < 0){
    cerr << "lo.LSeek returning negative number  (sendIshImage)" << endl;
    lo.Exec("COMMIT");
    lo.Close();
  }
  char* buffer = new char[end];
  int bytesRead = lo.Read(buffer, end);
  lo.Exec("COMMIT");
  if(bytesRead != end){
    cerr << "bytes read is not as requested, perhaps we should get out of here.. " << endl;
    delete buffer;
    lo.Close();
    return;
  }
  sApp("<ishFullImage>");
  qiApp(probeIndex);
  qiApp(index);
  sApp(buffer, end);
  sApp("<ishFullImageEnd>");
  delete buffer;
  lo.Close();
  writeArray();
}

void ConnectionObject::sendGenomicRegionProbeSetMatches(string chr, int start, int end){
  map<string, chromAnnotation*>::iterator cait = pSet->chromosomeAnnotation.find(chr);
  if(cait == pSet->chromosomeAnnotation.end()){
    cerr << "No such chromsome known" << endl;
    return;
  }
  chromAnnotation* annot = (*cait).second;   // so much easier to deal with.. 
  // now lets do some work on the start and end
  if( end < start){
    cerr << "end smaller than start, can't be bothered, sort it out man" << endl;
    return;
  }
  if(start < 1) { start = 1; }
  if(end < 1 || end > annot->chromosomeSize){ end = annot->chromosomeSize; }
  //if(end > annot->chromosomeSize || end < 1){ end = annot->chromosomeSize; }
  // ok,, it is still possible that start is larger than chromsomeSize,,
  if(start > end){ // annot->chromosomeSize){
    cerr << "stupid bloody attmept, get over it, what do you think you are doing ?? " << endl;
    return;
  }
  // but now, start is bigger than 0, smaller than end, and end is smaller than or equal to the size.. 
  // work out where to start..
  int sRegion = (start-1)/annot->regionSize;
  int eRegion = (end-1)/annot->regionSize;   // which is all very pretty.. 
  if(eRegion >= annot->regionNo){     // something is not so good..
    cerr << "eRegion is larger than the annot->regionNo, this is not good" << endl;
    return;
  }
  // ok, if we got here we better start writing something to the poor client..
  sApp("<chromRegionProbeSetMatches>");       // catchy header that!! 
  // include some description of the matches here, i.e. the chromsome and start - stop
  sApp(chr);
  qiApp(start);
  qiApp(end);
  // and lets include the chromsome size..
  qiApp(annot->chromosomeSize);

  for(int i=sRegion; i <= eRegion; i++){
    // unfortunately, we now have to go through everything..
    //cout << "going through regions from sRegion " << sRegion << "  to eRegion " << eRegion << "  currently at " << i << endl;
    for(int j=0; j < annot->regions[i]->pMatchNo; j++){
      //cout << "and j is " << j << endl; 
      // do we have an overlap between our end and start or not..
      probeSetMatch* match = annot->regions[i]->pMatches[j];    // otherwise the semantics will be too bad..
      // check to see if we've already covered this one..
      if(i > sRegion){
	if(match->cStart < annot->regions[i]->start || match->cEnd < annot->regions[i]->start){
	  // i.e. it overlaps with the preceding one.. 
	  continue;
	}
      }
      // determine if overlap ..
      if(  ((start <= match->cStart) != (end <= match->cEnd)) || ((start <= match->cEnd) != (end <= match->cStart)) ){
	// we have ovelap.. Ho yea..
	//cout << "we have an overlap between a feature and a request, we should be sending something... hooo yeahh." << endl;
	// so if we have an overlap, write 1..
	qiApp(1);
	// first the coordinates of the match..
	qiApp(match->cStart);
	qiApp(match->cEnd);
	// and then the details of the probe set match.. essentially a whole load of numbers..
	qiApp(match->dbIndex);
	qiApp(match->afStart);
	qiApp(match->afEnd);
	qiApp(match->afLength);
	qiApp(match->alignLength);

	//qiApp(match->af_n_count);

	qiApp(match->match);
	dApp(match->expectation);
      }
    }
  }
  // and if we get here, then we have no more probeSetMatches to send so.. let's write a -1;
  qiApp(-1);
  // and then the tail..
  sApp("<chromRegionProbeSetMatchesEnd>");
  // call writeArray..
  writeArray();
  // and we are done... ho yeahh... 
}

/// ok how about sending a genomicRegion with things.. in it..
void ConnectionObject::sendGenomicRegionEnsemblGenes(int requestId, string chr, int start, int end){
  map<string, chromAnnotation*>::iterator cait = pSet->chromosomeAnnotation.find(chr);
  if(cait == pSet->chromosomeAnnotation.end()){
    cerr << "No such chromsome known" << endl;
    return;
  }
  chromAnnotation* annot = (*cait).second;   // so much easier to deal with.. 
  // now lets do some work on the start and end
  if( end < start){
    cerr << "end smaller than start, can't be bothered, sort it out man" << endl;
    return;
  }
  if(start < 1) { start = 1; }
  if(end > annot->chromosomeSize || end < 1){ end = annot->chromosomeSize; }
  // ok,, it is still possible that start is larger than chromsomeSize,,
  if(start > annot->chromosomeSize){
    cerr << "stupid bloody attmept, get over it, what do you think you are doing ?? " << endl;
    return;
  }
  // but now, start is bigger than 0, smaller than end, and end is smaller than or equal to the size.. 
  // work out where to start..
  int sRegion = (start-1)/annot->regionSize;
  int eRegion = (end-1)/annot->regionSize;   // which is all very pretty.. 
  if(eRegion >= annot->regionNo){     // something is not so good..
    cerr << "eRegion is larger than the annot->regionNo, this is not good" << endl;
    return;
  }
  /// tiem to start writing things to the client..
  sApp("<chromRegionEnsemblGenes>");
  qiApp(requestId);
  sApp(chr);
  qiApp(start);
  qiApp(end);
  qiApp(annot->chromosomeSize);
  for(int i=sRegion; i <= eRegion; i++){
    // unfortunately, we now have to go through everything.. 
    for(int j=0; j < annot->regions[i]->ensGeneNo; j++){
      // do we have an overlap between our end and start or not..
      ensemblGene* gene = annot->regions[i]->ensGenes[j];    // otherwise the semantics will be too bad..
      // check to see if we've already covered this one..
      if(i > sRegion){
	if(gene->start < annot->regions[i]->start){    // for the ens genes start is always smaller than end anyway.. 
	  // i.e. it overlaps with the preceding one.. 
	  // and we should have sent it the last time.. 
	  continue;
	}
      }
      /// determine if we have overlap. This is incredibly wasteful as we only need to do this for the first and the last 
      // genomic regions, but that would be more code, so what the helll....
      // something to improve in the future.. 
      if(  ((start <= gene->start) != (end <= gene->stop)) || ((start <= gene->stop) != (end <= gene->start)) ){
	/// just write the thingy like before..
	qiApp(1);   // so I know to keep reading..
	qiApp(gene->dbIndex);
	sApp(gene->ensemblId);

	//sApp(gene->externalId);
	//sApp(gene->description);

	sApp(gene->chromosome);
	qiApp(gene->start);
	qiApp(gene->stop);
	qiApp(gene->strand);
      }
    }
  }
  qiApp(-1);
  sApp("<chromRegionEnsemblGenesEnd>");
  writeArray();
}


void ConnectionObject::sendGenomicRegion(string chrom, int start, int stop, int association, int target){
  vector<int> temp;     // we won't use it anyway..
  sendGenomicRegion(chrom, start, stop, association, target, temp, false);
}

void ConnectionObject::sendGenomicRegion(string chrom, int start, int stop, int association, int target, vector<int>& psets, bool useT){
  // the target is related to who wants to see the genomic region down at the client side,, and generally is 
  // just one or 2 . at the moment.. 

  // send everything we've got associated with a given region. Not particularly flexible I know, but primarily 
  // aimed at sending small regions with each probe set, -and then perhaps several different loci at the same time
  // as is quite commonly the case. Hope this isn't going to get to slow .. 
  
  // first some checks to make sure that chrom, stop, and start are kind of reasonable..
  int maxRange = 100000000;   // no more than 500 regions sent.. Not good.. but maybe it will stop the crashes..  
  cout << "beginning of sendgenomic region " << endl;
  map<string, chromAnnotation*>::iterator cait = pSet->chromosomeAnnotation.find(chrom);
  if(cait == pSet->chromosomeAnnotation.end()){
    cerr << "No such chromosome : " << endl; //chrom << endl;
    return;
  }
  chromAnnotation* annot = (*cait).second;
  // then some sanity checks on the numbers..
  if(start < 1){ start = 1; }     // for whole chromosome give two negative numbers 
  if(stop < 1 || stop > annot->chromosomeSize){ stop = annot->chromosomeSize; }

  if(stop < start || stop == start){
    cerr << "give me a break, stop is smaller than or equal to start" << endl;
    return;
  }
  if(stop - start > maxRange){
    start = (start + stop + maxRange)/2;
    stop = start + maxRange;
  }
  // ok start and stop should both be within range now.. unless I'm forgetting something stupid..
  uint sRegion = (start-1)/annot->regionSize;
  uint eRegion = (stop-1)/annot->regionSize;
  //// eRegion should not be able to be larger than the regionNo, nevertheless..
  if(eRegion >= annot->regionNo || sRegion >= annot->regionNo || sRegion > eRegion){
    cerr << "socketNumber : " << socketNumber <<  "eRegion is larger than regionNo, this is bad, fix it!! " << endl;
    return;
  }
  /// and now we should be able to just go through everything. As our new data structure is more in line with 
  /// ensembl, -and gene begin is always the lower number, this simplifies the overlap determinations a lot, nice..
  sApp("<chromRegion>");
  qiApp(target);         // who want to see this map.. -- make sure to update the client on this one.. 
  qiApp(association);    // this should be the dbIndex of the afid this is associated with.. this is important, so that I have some idea of what to do with this at the other end. 
  sApp(chrom);
  qiApp(start);
  qiApp(stop);
  qiApp(annot->chromosomeSize);
  // and then just go through everything and see if it overlaps. This could be speeded up by only running the overlap functions
  // on the last go, but for now, keep the code simple..
  for(int i=sRegion; i <= eRegion; i++){
    for(int j=0; j < annot->regions[i]->ensGeneNo; j++){
      ensemblGene* gene = annot->regions[i]->ensGenes[j];
      // check to see if we've covered this one already..
      if(i > sRegion && gene->start < annot->regions[i]->start){  // should already have been sent. 
	continue;
      }
      //cout << "Sending Ensembl Gene ... " << endl;
      // Determine overlap.. 
      if(gene->start < stop && gene->stop > start){
	qiApp(1);
	qiApp(gene->dbIndex);
	sApp(gene->ensemblId);
	qiApp(gene->start);
	qiApp(gene->stop);
	qiApp(gene->strand);  // should be 1 or -1..
	// and then for each transcript we need to send individual information..
	qiApp(gene->transcriptNo);
	for(int k=0; k < gene->transcriptNo; k++){
	  qiApp(gene->transcripts[k]->index);
	  sApp(gene->transcripts[k]->id);
	  qiApp(gene->transcripts[k]->start);
	  qiApp(gene->transcripts[k]->stop);
	  qiApp(gene->transcripts[k]->exonNo);
	  // and then the exons..
	  for(int m=0; m < gene->transcripts[k]->exonNo; m++){
	    sApp(gene->transcripts[k]->exons[m]->id);
	    qiApp(gene->transcripts[k]->exons[m]->start);
	    qiApp(gene->transcripts[k]->exons[m]->stop);
	    qiApp(gene->transcripts[k]->exons[m]->codeStart);
	    qiApp(gene->transcripts[k]->exons[m]->codeStop);
	  }
	}
	//cout << "socketNumber " << socketNumber << "\t\t\t: i " << i << "  j: " << j << "  after " << endl;
	/////////// let's try sending all of the information we have on this gene as well. -which is a fair bit of data.. 
	//////////  just use the ensemblAnnotation and the index..
	map<int, multimap<int, string> >::iterator ait = pSet->ensemblAnnotation.find(gene->dbIndex);
	map<int, string>::iterator fit;
	//cout << "socketNumber " << socketNumber << "\t\t\t: got the ait iterator " << endl;
	for(multimap<int, string>::iterator mit=(*ait).second.begin(); mit != (*ait).second.end(); mit++){
	  fit = pSet->ensemblFields.find((*mit).first);
	  if(fit != pSet->ensemblFields.end()){
	    // everything is OK, we can now go through and write out the annotation as it comes out.
	    //    cout << "\t\tSocket Number : " << socketNumber << "\t\t before appending to data" << endl;
	    qiApp(2);  // it's not 1, or we can't distinguish between here and above.. 
	    sApp((*fit).second);
	    sApp((*mit).second);
	    //cout << "\t\tSocket Number : " << socketNumber << "\t\t after appending to data" << endl;
	  }
	}
	qiApp(-1);  // don't do this anymore.. 
      }
    }
    //cout << "socketNumber : " << socketNumber << "\t\t\t i : " << i << endl;
  }
  qiApp(-1);   // indicate that we're done. with ensemblGenes.

  // Send other types of transcripts. These will be identified by some sort of string -- but basically these are not definite transcripts..
  // but blast matches to some sort of transcript. We will work out some way of looking at them later.. Ok..
  for(int i=sRegion; i <= eRegion; i++){
    for(int j=0; j < annot->regions[i]->transcriptNo; j++){
      Transcript* trans = annot->regions[i]->transcripts[j];
      if(i > sRegion && trans->start < annot->regions[i]->start){   // we should have sent this already..
	continue;
      }
      // determine if it overlaps with the requested range..
      if(trans->start < stop && trans->stop > start){
	qiApp(3);          // 3 indicates transcript... 
	sApp(trans->source);   // what kind is it ..
	sApp(trans->id);
	//cout << "\t\tTRANS CHROMOSOME " << trans->chromosome << endl;
	sApp(trans->chromosome);
	qiApp(trans->start);
	qiApp(trans->stop);
	qiApp(trans->strand);
	qiApp(trans->length);   // the length of the transcript.. -which may not be fully covered by the alignments..
	// and then the number of exons..
	qiApp(trans->exonNo);
	for(int k=0; k < trans->exonNo; k++){
	  qiApp(trans->exons[k]->start);
	  qiApp(trans->exons[k]->stop);
	  qiApp(trans->exons[k]->tbegin);
	  qiApp(trans->exons[k]->tend);
	  
	}
      }
    }
  }
  qiApp(-1);
		  
  /// set a maxExpectation for stuff we're sending. Again this should be user controlled, but we can set quite high, 
  /// as we would probably want to include things from the same probeSetMatch.. 
  double maxExpectation = 1;   // to see if this is causing strange behaviour... 
  //  double maxExpectation = 1e-6;
  //cout << "\t\t\tprobe set matches " << endl;
  set<int> alreadyDone;    // so we don't include the same probe set more than once in the vector psets.. 
  for(int i=sRegion; i <= eRegion; i++){
    for(int j=0; j < annot->regions[i]->pMatchNo; j++){
      probeSetMatch* match = annot->regions[i]->pMatches[j];
      //cout << "trying to send a probe set match, cStart is " << match->cStart << endl
      //	   << "                                  cEnd   is " << match->cEnd << endl;
      if(i > sRegion && match->cStart < annot->regions[i]->start){
	continue;
      }
      if(match->cStart < stop && match->cEnd > start && match->expectation < maxExpectation){
	int logExp = -500;
	if(match->expectation){
	  logExp = (int)log10(match->expectation);
	}
	if(logExp <= pSetThresholds.maxExpect && match->alignLength > pSetThresholds.minLength && (float)match->match/(float)match->alignLength > pSetThresholds.minMatch){
	  if(!alreadyDone.count(match->dbIndex)){
	    psets.push_back(match->dbIndex);
	    alreadyDone.insert(match->dbIndex);
	  }
	}
	qiApp(1);
	// and then the information describing the match...
	//cout << "\tsending probe Set match id : " << match->dbIndex << "\t start: " << match->cStart << "\t end: " << match->cEnd << "\t expect: " << match->expectation << endl;
	qiApp(match->dbIndex);
	qiApp(match->cStart);
	qiApp(match->cEnd);
	qiApp(match->afStart);
	qiApp(match->afEnd);
	qiApp(match->afLength);
	qiApp(match->alignLength);
	//qiApp(match->alignLength);
	qiApp(match->match);
	dApp(match->expectation);
	qiApp(match->strand);
      }
    }
  }  
  qiApp(-1);
  /// and finallly let's send any in_situ_probe_blast matches if there are any.. pretty much the same idea...
  /// really should sort this out sometime..
  //cout << "\t\t\tishMatches : " << endl;
  for(int i=sRegion; i <= eRegion; i++){
    for(int j=0; j < annot->regions[i]->ishMatchNo; j++){
      qiApp(1);
      ishProbeMatchSet* match = annot->regions[i]->ishMatches[j];
      qiApp(match->dbIndex);
      qiApp(match->minPos);
      qiApp(match->maxPos);
      qiApp(match->pLength);
      //qiApp(match->alignLength);
      //qiApp(match->match);
      //qiApp(match->cStart);
      //qiApp(match->cEnd);
      fApp(match->score);
      //      dApp(match->expectation);
      qiApp(match->strand);
      sApp(match->chromosome);   // don't send the sequences, not really needed. eh... !!..
      qiApp(match->assemblyId);
      qiApp(match->matches.size());
      set<ishProbeMatch*>::iterator imit;
      for(imit = match->matches.begin(); imit != match->matches.end(); imit++){
	qiApp((*imit)->pStart);
	qiApp((*imit)->pEnd);
	qiApp((*imit)->cStart);
	qiApp((*imit)->cEnd);
	fApp((*imit)->percent);
      }
    }
  }
  //cout << "\t\t\twritten all of the stuff" << endl;
  qiApp(-1);    // end of thingy.. hooo yeahh. 
  sApp("<chromRegionEnd>");
  writeArray();    // and that is all we really have to do..
  /// and this is the time when one realises that a little bit of object orientated inheritance and stuff would really server me well...
  /// but there you go, it would be a little more thinking, and I have to many other things to think about..
  // ---- and now let's send the probe Set matches on their way in the same manner. These thoug
}
/// hoohooo yeaah. 


void ConnectionObject::anovaSort(){
  //  cout << "beginning of anova sort" << endl;
  if(clientIndex.size() < 2) { return; }
  //map<float, int> scores;              // will automatically sort by the float,, -but 0 is 0 is 0.. not good.
  vector<dist_set> scores(clientIndex.size());
  //////////////////// REMEMBER .. clientIndex has already had 1 subtracted to make sure that
  //////////////////// that it fits the pset-Data.. something..

  // split clientIndex into 2 parts, and work out the starts and stops. 
  uint mid = clientIndex.size()/2;      // note this is used as the size in the loop, and as such, will not be duplicated. 
  //cout << "\tanovaSort mid set to : " << mid << endl
  //   << "\tand sizeo f clientIndex is : " << clientIndex.size() << endl;
  AnovaProcessor* process1 = new AnovaProcessor(pSet, 0, mid, &clientIndex, &scores);
  AnovaProcessor* process2 = new AnovaProcessor(pSet, mid, clientIndex.size(), &clientIndex, &scores);
  //cout << "\tanovaprocessors created and ready to roll" << endl;
  process1->start();
  process2->start();     // hmm, I wonder.. 
  //cout << "\tboth processes started" << endl;
  process1->wait();
  process2->wait();

  // create a vector and fill it using a reverse iterator..
  sort(scores.begin(), scores.end(), r_comp_set());
  //cout << "\tscores sorted " << endl;
  vector<int> tIndex;
  //  tIndex.reserve(scores.size());
  for(int i=0; i < scores.size(); i++){
    if(scores[i].value > 0){
      tIndex.push_back(scores[i].index);
    }
  }
  //cout << "\ttempIndex created " << endl;
  //  for(rit = scores.rbegin(); rit != scores.rend(); rit++){
  //  tIndex.push_back((*rit).second);
  //}
  delete process1;
  delete process2;
  //cout << "processes deleted " << endl;
  writeIndex(tIndex, "Anova");
  //cout << "end of process " << endl;
}

//QMutex tLock;
void ConnectionObject::euclidSort(probe_set* pset1, vector<uint> expts){
  // compare using all of the available experiments for pset1, against all others..
  //  cout << "beginning of euclid compare" << endl;
  vector<dist_set> distances(pSet->data.size());
  // find the midpoint and make a couple of euclidSortProcessors!!
  uint mid = pSet->data.size()/2;

//   cout << "Creating euclid sort processors for experiments : ";
//   for(uint i=0; i < expts.size(); i++){
//     cout << "\t" << expts[i] << endl;
//   }
//   cout << endl;

  EuclidSortProcessor* proc1 = new EuclidSortProcessor(pSet, pset1, 0, mid, &distances, clientChips, expts);
  EuclidSortProcessor* proc2 = new EuclidSortProcessor(pSet, pset1, mid, pSet->data.size(), &distances, clientChips, expts);
  
  proc1->start();
  proc2->start();
  proc1->wait();
  proc2->wait();
  //cout << "initialised distances " << endl;
  //probe_set probe1 = *pset1;
  //probe_set probe2;
  //  int size = pSet->data.size();
  //QMutex tLock;
  //qApp->lock();
  //tLock.lock();
  //for(int i=0; i < size; i++){
    //    probe2 = pSet->data[i];
    //tLock.unlock();
    //distances[i].index = pSet->data[i]->index;
    //distances[i].index = probe2.index;
    //distances[i].value = heavyCompare(&probe1, &probe2, pset1->exptIndex);
    //distances[i].value = heavyCompare(pset1, pSet->data[i], pset1->exptIndex, pset1->exptSize);
  //}
  //tLock.unlock();
  //qApp->unlock();
  //cout << "just before sorting the distances" << endl;
  sort(distances.begin(), distances.end(), comp_set());
  vector<int> tIndex;
  for(int i=0; i < distances.size(); i++){
    if(distances[i].value > 0){
      tIndex.push_back(distances[i].index);
    }
  }
  delete proc1;
  delete proc2;
  //cout << "just before writing the index " << endl;
  writeIndex(tIndex, "Compare to probe set");
}

//QMutex hLock;

float ConnectionObject::heavyCompare(probe_set* pset1, probe_set* pset2, uint* expts, uint es){
  //
  //   Calculate the euclidean distance between the probe set 1 and the probe set 2 using a heavy
  //   algorithm which calculates the mean euclidean distance of all against all probe sets.
  //   
  //   Hence passing the same probe set in one and 2 will not give a 0 distance, but will give a 
  //   value which is related to how closely the different probe sets correlate with each other.
  //   
  //   Obviously I can only make a comparison for the same experimental points, so I have 
  //   to take that into account. Use the slow approach of actually creating a new set of vectors
  //   and manually filling those where both have got values. Then z-score normalise these..
  //   
  // OK, here goes..

  //  vector< vector<float> > v1(pset1->probes.size());
  //vector< vector<float> > v2(pset2->probes.size());
  
  float** v1 = new float*[pset1->probeNo];  
  float** v2 = new float*[pset2->probeNo];
  // and reserve enough memory for each one..
  for(int i=0; i < pset1->probeNo; i++){ v1[i] = new float[es]; }
  for(int i=0; i < pset2->probeNo; i++){ v2[i] = new float[es]; }
  int selEx = 0; 


  //map<int, int>::iterator it1, it2;
  for(int i=0; i < es; i++){
    if(expts[i] < pset1->allExptNo){           // allExptNo should be the same for both !!
      if(pset1->exptLookup[expts[i]] != -1 && pset2->exptLookup[expts[i]] != -1){
	for(int j=0; j < pset1->probeNo; j++){
	  v1[j][selEx] = pset1->probes[j][pset1->exptLookup[expts[i]]];
	}
	for(int j=0; j < pset2->probeNo; j++){
	  v2[j][selEx] = pset2->probes[j][pset2->exptLookup[expts[i]]];
	}
	selEx++;
      }
    }
  }
//     it1 = pset1->exptLookup.find(expts[i]);
//     it2 = pset2->exptLookup.find(expts[i]);
//     if(it1 != pset1->exptLookup.end() && it2 != pset2->exptLookup.end()){
//       //  index1.push_back((*it1).second);
//       //index2.push_back((*it2).first);
//       for(int j=0; j < pset1->probes.size(); j++){
// 	v1[j].push_back(pset1->probes[j][(*it1).second]);
//       }
//       for(int j=0; j < pset2->probes.size(); j++){
// 	v2[j].push_back(pset2->probes[j][(*it2).second]);
//       }
//     }
//   }

  // AND now just do an all against all comparison, and take the mean value of the comparisons.. 
  float distance = 0;
  //cout << "beginning the distance calculations" << endl;
  int v1size = pset1->probeNo;
  int v2size = pset2->probeNo;

  for(int i=0; i < v1size; i++){
    for(int j=0; j < v2size; j++){
      distance += euclidean(v1[i], v2[j], selEx);  
    }
  }
  distance = distance / (float)(v1size * v2size);
  // and lets delete the two float**,, we can use the function..
  delProbes(v1, v1size);
  delProbes(v2, v2size);
  //  distance = distance / (float)(pset1->probes.size() * pset2->probes.size());
  return(distance);
}

float ConnectionObject::meanEuclidCompare(probe_set* p, float* target, uint* expts, uint es, bool normed){
  // like the function below, except that it compares the target to either the normalised mean of the 
  // raw data, or the normalised mean of locally normalised data... depending on the boolean normed.
  // quite a mouthful that... eh..
  
  if(!p->index){
    return(-1.0);
  }
  //if(target.size() != expts.size()) { return(-1); }
  float* meanValues = new float[es];
  //  vector<float> meanValues(expts.size(), 0);     // will puke if it doesn't have the right sizes..

  float** tempValues = copyProbes(p->probes, p->probeNo, p->exptSize);
  //  vector< vector<float> > tempValues;            // actually need this as we may want to normalise.. 
  //map<int, int>::iterator it;
  //for(int i=0; i < meanValues.size(); i++){
  //  cout << "\t\tInitial value of meanValues : " << i << "\t" << meanValues[i] << endl;
  //}
  //// Just Checking, I'm not sure if the constructor is correct!! 
  //tempValues = p->probes;
  
  if(normed){
    for(int i=0; i < p->probeNo; i++){
      zScore(tempValues[i], p->exptSize);
      //      zScore(tempValues[i], p->probeNo);
    }
  }
  
  for(int i=0; i < es; i++){
    if(expts[i] < p->allExptNo){
      if(p->exptLookup[expts[i]] != -1){
	for(int j=0; j < p->probeNo; j++){
	  meanValues[i] += tempValues[j][ p->exptLookup[expts[i]] ];
	  //	  meanValues[i] += p->probes[j][ p->exptLookup[expts[i]] ];
	}
      }else{
	delete []meanValues;
	delProbes(tempValues, p->probeNo);
	return(-1);
      }
    }else{
      delete []meanValues;
      delProbes(tempValues, p->probeNo);
      return(-1);
    }
  }
  //  it = p->exptLookup.find(expts[i]);
  //  if(it == p->exptLookup.end()){
  //    return(-1);
  //  }
  //  for(int j=0; j < p->probes.size(); j++){
  //    meanValues[i] += tempValues[j][(*it).second];
  //  }
  //}
  //// hmm, really I should divide by the meanValues, --but because in the long run we are going to normalise 
  //// them anyway, there isn't actually any point.. 
  //// Assume that the target is already normalised, and that just run the euclidean function...
  zScore(meanValues, es);
  float score = euclidean(meanValues, target, es);
  delete []meanValues;
  delProbes(tempValues, p->probeNo);
  return(score);
  //return(euclidean(meanValues, target));
}
	 

float ConnectionObject::singleEuclidCompare(probe_set* pset1, float* target, uint* expts, uint es){
  // compare a each probe in a probe set against a single vector of already normalised values
  // for the given experiments.. -- first look up the experiments..
  if(pset1->probeNo < 1) {
    cout << "single euclid compare. probe set has no probes " << endl;
    return(-1); 
  }
  float** v1 = new float*[pset1->probeNo];
  for(int i=0; i < pset1->probeNo; i++){ v1[i] = new float[es]; }

  //  vector< vector<float> > v1(pset1->probes.size());   // for keeping normalised values in.
  //vector<float> v2;                                  // for keeping the target values that correspond to the expt. 
  //// Here is the QUESTION...........
  ////         If there is little overlap between the expts specified in the target, and those present for 
  ///          the gene we are currently looking at.. then What should we do????????
  ///          Sounds like arbritrary threshold time to me.. -- but as the user can check this,, and select the 
  ///          points, I will for the moment.. be rather brutal, and say if there is not a perfect overlap, then
  ///          I will return a -1, which I can check for in the thingy. This means that the index returned will 
  ///          only contain those that have a perfect match for the things we are looking for. Which in our case always
  ///          will include all the genes from chip A... -- at least its a definite answer as opposed to ugy 
  ///          prevarication...

  //map<int, int>::iterator it;
  //cout << "probe index : " << pset1->index << endl;
  for(int i=0; i < es; i++){
    //it = pset1->exptLookup.find(expts[i]);
    if(expts[i] >= pset1->allExptNo){
      //cout << "experiment index is larger than the exptMap " << endl
      //	   << "index : " << expts[i] << endl;
      delProbes(v1, pset1->probeNo);
      return(-1);
    }
    if(pset1->exptLookup[expts[i]] == -1){
      //cout << "exptLookup[" << expts[i] << "] : " << pset1->exptLookup[expts[i]] << endl;
      //cout << "expts[" << i << "] : " << expts[i] << endl;
      delProbes(v1, pset1->probeNo);
      return(-1);
    }
    //
    //if(it == pset1->exptLookup.end()){
    //  cout << "BY BY baby, I'm ignoring you because you don't have the answers I'm looking for" << endl;
    //  return(-1);
    //}
    for(int j=0; j < pset1->probeNo; j++){
      v1[j][i] = pset1->probes[j][pset1->exptLookup[expts[i]]];
      //      v1[j].push_back(pset1->probes[j][(*it).second]);
    }
  }
  // NORMALISE and get the distances..
  float distance = 0;
  //if(v1.size() < 1) { return(-1); }
  for(int i=0; i < pset1->probeNo; i++) { 
    //cout << "doint the z-score " << endl;
    zScore(v1[i], es); 
    //cout << "calling for some distance " << endl;
    distance += euclidean(target, v1[i], es);
  }
  // and go through and get the actual distance..
  delProbes(v1, pset1->probeNo);
  //cout << "Total distance is " << distance << endl;
  return(distance/pset1->probeNo);
}

void ConnectionObject::doCreateNewUser(){
  vector<QString> words = splitString(lastMessage, ':'); // should be OK.. if I remember correctly.
  if(words.size() != 8){
    cout << "createNewUser, words size is wrong. Forget this.. " << endl;
    return;
  }
  if(userId != 1){
    cout << "You don't have the privilege byby.. " << endl;
  }
  bool ok;
  vector<int> oldh(3);
  vector<int> newh(3);
  for(int i=0; i < 3; i++){
    oldh[i] = words[i+1].toInt(&ok);
    newh[i] = words[i+5].toInt(&ok);  // should check the Ok,, but it doesn't matter..
  }
  // first we need to check whether we have the appropriate privileges.. 
  ostringstream query;
  query << "select index from users where index = 1 and key1 = " << oldh[0] << " and key2=" << oldh[1]
									      << " and key3=" << oldh[2];
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "connection not so good.. " << endl;
    return;
  }
  if(!conn.Exec(query.str().c_str())){
    cerr << "Command didn't work.. " << endl;
    cerr << conn.ErrorMessage() << endl;
    return;
  }
  if(conn.Tuples() != 1){
    cerr << "Wrong ID, or wrong passwords, or something like that " << endl;
    return;
  }
  // ok. so we have the right user id. and we also have the appropriate password.. well.. we can now make the insert string
  ostringstream insertstring;
  insertstring << "insert into users (index, user_name, key1, key2, key3) values ("
	       << "nextval('user_ids'), '" << words[4].latin1() << "', " << newh[0] << ", " << newh[1] << ", " << newh[2] << ")";
  if(!conn.Exec(insertstring.str().c_str())){
    cerr << "couldn't insert new user, not sure why.. but something or other.. " << endl;
  }
  if(conn.CmdTuples() != 1){
    cerr << "same again, couldn't insert.. what's going on here.. " << endl;
  }else{
    // it worked and we need to update the the map<int, userInformation> userTable..
    //    pSet->userTable.insert(make_pair(
    ostringstream idQueryString;
    idQueryString << "select index from users where user_name = '" << words[4].latin1() << "'";
    if(!conn.Exec(idQueryString.str().c_str())){
      cerr << "Coulnd't execute id lookup. " << endl;
      cerr << conn.ErrorMessage() << endl;
      return;
    }
    if(conn.Tuples() != 1){
      cerr << "conn.Tuples isn't 1. This doesn't make much sense " << endl
	   << "Tuples : " << conn.Tuples() << endl;
      return;
    }
    int index = atoi(conn.GetValue(0, 0));
    // and now we do the insert.. first make sure it's not taken.. hmm 
    // hmm oh bugger that.
    string uName(words[4].latin1());
    string fName;
    string lName;
    pSet->userMutex->lock();
    pSet->userTable.insert(make_pair(index, userInformation(index, uName, fName, lName)));
    //map<int, userInformation>::iterator it = pSet->userTable.find(index); // ok.. --- I don't user this after all.. 
    pSet->userMutex->unlock();
    QCustomEvent* updateUsers = new QCustomEvent(1001);     // user is 1000.. 
    QThread::postEvent(server, updateUsers);
    sApp("<newUser>");
    qiApp(index);
    sApp(uName);
    sApp(fName);
    sApp(lName);
    sApp("<newUserEnd>");
    writeArray();  // it's quite easy after all.. 
  }
  
  // write some message to indicate if it worked or not. ADD later when we have proper mechanism for doing this..
}


void ConnectionObject::doChangePassword(){
  bool ok;
  int start = 0;
  int end = lastMessage.find(":", start);
  QString userName = lastMessage.left(end);
  cout << "user Name is " << userName << endl;
  start = end + 1;
  vector<int> numbers; 
  end = lastMessage.find(":", start);
  while(end > 0){
    QString numstring = lastMessage.mid(start, end-start);
    //cout << "numstring is " << numstring << endl;
    numbers.push_back(numstring.toInt(&ok));  // hmm 
    start = end+1;
    end = lastMessage.find(":", start);
  }
  if(numbers.size() != 6){
    // I need some sort of error message here, but am not sure how to implement, could just 
    // roll a message.. but there you go..
    cerr << "numbers.size is not 6, but : " << numbers.size() << endl;
    return;
  }
  // first three numbers are the old key1, key2, and key3.. -- we can do this update in one go.. if I remember
  // how to.. 
  ostringstream query;
  query << "update users set key1 = " << numbers[3] << ", key2 = " << numbers[4] << ", key3 = " << numbers[5] << " where "
	<< "user_name = '" << userName << "' and key1 = " << numbers[0] << " and key2 = " << numbers[1] 
	<< " and key3 = " << numbers[2];
  //cout << "query is: \n" << query.str() << endl;
  //const char* conninfo = "dbname=expression";
  //PgConnection conn(conninfo);
  PgDatabase conn(conninfo);
  //  PgCursor conn(conninfo, "portal");
  if(conn.ConnectionBad()){
    cerr << "connection not good" << endl;
    return;
  }
  if(!conn.Exec(query.str().c_str())){
    sApp(string("<ChangePassWordError>"));
    sApp(string("<ChangePassWordErrorEnd>"));
    writeArray();
    cerr << "command did not work: " << query << endl;
    cerr << "Error: " << conn.ErrorMessage() << endl;
    return;
  }
  int changeOK;
  if(conn.CmdTuples() != 1){
    cerr << "Password Change failed or went multiplied.. cmdTuples: " << conn.CmdTuples() << endl;
    changeOK = 0;
  }else{
    changeOK = 1;
  }
  sApp(string("<passWordChanged>"));
  qiApp(changeOK);
  sApp(string("<passWordChangedEnd>"));
  writeArray();
}

void ConnectionObject::meanComparisonSort(){
  /// take a vector of floats (already normalised).. -
  /// then compare these against everything that is in client Index..
  /// and  sort.. do a comparison against either the mean of normalised or 
  /// the mean of the raw (in both cases the mean is normalised across the series!!). 
  /// depending on the normed boolean.. (0 or 1).. 
  vector<uint> eIndex;
  vector<float> values;
  
  bool v_ok;         // value -float
  bool i_ok;         // int, experimental Index,, 
  float v;
  uint ei;

  vector<QString> words = splitString(lastMessage, ':');
  if(words.size() % 2 == 1 || words.size() == 0){    // need an odd number.. 
    cerr << "meanComparisonSort , words size doesn't make any sense: " << endl;
    return;
  }
  for(int i=0; i < (words.size()-2); i += 2){
    v = words[i].toFloat(&v_ok);
    ei = words[i+1].toInt(&i_ok);
    if(v_ok && i_ok){
      values.push_back(v);
      eIndex.push_back(ei);
    }
  }
  QString normString = words[words.size()-2];
  int normInt = normString.toInt(&i_ok);
  bool distribution = (words[words.size()-1] != "0");
  //cout << "Do mean comparison sort, and the user is looking for a distribution " << distribution << endl;
//   int start = 0; 
//   int end = lastMessage.find(":", start);
//   while(end != -1 && end != 0){
//     QString vs = lastMessage.mid(start, (end-start));
//     start = end+1;
//     end = lastMessage.find(":", start);
//     QString is = lastMessage.mid(start, (end-start));
//     v = vs.toFloat(&v_ok);
//     i = is.toInt(&i_ok);
//     if(v_ok && i_ok){
//       values.push_back(v);
//       eIndex.push_back(i);
//     }
//     start = end +1;
//     end = lastMessage.find(":", start);
//   }
//   QString normString = lastMessage.right(1); 
  //cout << endl << "normString:\t" << normString << "\tnormInt: " << normInt << endl;
  //cout << "In the mean comparison sort function, last Message : " << lastMessage << "\ti_ok\t" << i_ok << endl;
  //for(int i=0; i < values.size(); i++){
    //cout << eIndex[i] << "\t: " << values[i] << endl;
  //}
  if(i_ok){
    if(normInt < 0 || normInt > 1){
      return;
    }
  }else{
    return;
  }
  // i is either 0 or 1. just cast it to a boolean. 
  if(values.size() < 1) { return; }
  // and now we can just use the mean Euclid Compare thingy, and do the sort.. ho yeah.. 
  
  vector<dist_set> distances;
  distances.reserve(clientIndex.size());
  float d;
  float* target = new float[values.size()];
  uint* expts = new uint[eIndex.size()];
  for(int i=0; i < values.size(); i++){ 
    expts[i] = eIndex[i];
    target[i] = values[i]; 
  }
  //cout << "normalisation is " << (bool)normInt << endl;
  for(int i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size()){
      d = meanEuclidCompare(pSet->data[clientIndex[i]],  target, expts, eIndex.size(), (bool)normInt);
      //      d = meanEuclidCompare(&pSet->data[clientIndex[i]],  values, eIndex, (bool)normInt);
      if(d >= 0){
	distances.push_back(dist_set(pSet->data[clientIndex[i]]->index, d));
      }
    }
  }
  // and if we actually get something after all of that,, we have to check it out..
  sort(distances.begin(), distances.end(), comp_set());
  vector<int> tIndex(distances.size());
  for(int i=0; i < tIndex.size(); i++){
    tIndex[i] = distances[i].index;
  }
  writeIndex(tIndex, "Mean Profile Sort");
  // if distribution, then send the actual values...
  if(distribution){
    sApp("<statCollection>");
    qiApp(1);
    sApp("Mean E Distance");
    qiApp(distances.size());
    for(int i=0; i < distances.size(); i++){
      qiApp(distances[i].index-1);
      fApp(distances[i].value);
    }
    sApp("<statCollectionEnd>");
    writeArray();
  }
  delete []target;
  delete []expts;
}

float ConnectionObject::variationCoefficient(probe_set* pset, uint* expts, uint es){
  // calculate the coefficient of variation for the given probe set and 
  // the given experimental indices..
  // There are actually several ways in which we could do this, what I will calculate here
  // is the mean value of the individual probe pair profile correlations..
  //
  // the coefficient of variation is simply the standard deviation / mean value
  // which obviously will be calculated separately for each probe pair... 
  // however,,, we could end up with situations where we have mean values which are close to 
  // 0. This is problematic,, and for this kind of data there really isn't much we can do 
  // about that.. hmmmm... 
  //
  // we could do something rather silly, like add a pseudocount when dealing with weight matrices,
  // -- we could find the minimum value,, and calculate something to add to all values which makes 
  // sure that they are all above 0, but that reduces the measure for low values.. 
  // -- the real problem is that although it is simple to say we don't want negative or 0 means
  // 
  // in the long run I may include the pm scores in the probe sets, -it's quite easy to incorporate this
  // change, but for now I think I will do the simple thing and just stick in a couple of thresholds..

  // ignore the data from the probe pair if the mean is not above 0
  // ignore value if the coefficient is above 10. 10 is pretty high anyway, but should exclude
  // the influence of stupid numbers..
  //
  // ok first work out what the appropriate experimental indices are going to be..
  //cout << "making a data Extractor " << endl;
  DataExtractor de;
  //cout << "getting the data " << endl;
  ExData* ed = de.extract(pset, expts, es);
  if(ed->exptNo < 2){
    cerr << "Connection Object variation Coefficient exptNo is less than 2" << endl;
    return(0);
  }
  int exptNo = ed->exptNo;
  // do an mScore normalisation of the data..
  //cout << "making a normaliser, exptNo is now : " << exptNo << endl;
  Normaliser norm;
  norm.mScore(ed->values, ed->probeNo, exptNo);
  //cout << "Normalised.. " << endl;

  // find the mean of the above and calcuate the appropriate values..
  float* meanV = new float[exptNo];
  for(int i=0; i < exptNo; i++){ meanV[i] = 0; };
  for(int i=0; i < ed->probeNo; i++){
    for(int j=0; j < exptNo; j++){
      meanV[j] += (ed->values[i][j]/((float)ed->probeNo));
    }
  }
  //cout << "Made the mean : " << endl;
  // ok, now find the mean and standard deviation of the mean, and return the coefficient.. 
  float sum = 0;
  float sq_sum = 0;
  for(int i=0; i < exptNo; i++){
    sum += meanV[i];
    sq_sum += (meanV[i] * meanV[i]);
  }
  float SS = sq_sum - (sum * sum)/((float)(exptNo));
  float std = sqrt(SS/((float)(exptNo - 1)));
  float mean_mean = sum / (float)exptNo;
  //cout << "made lots of stuff " << endl;
  delete meanV;
  //cout << "deleted mean,, or not as the case may be.. " << endl;
  delete ed;
  //cout << "deleted ed or not ..  " << endl;
  if(mean_mean > 0){
    return(std/mean_mean);
  }else{
    if(mean_mean == 0){
      return(0);
    }
    return(-std/mean_mean);
  }

//   uint* eI = new uint[es];
//   int selEx = 0;
  
//   //vector<int> eI;       // for the local indices..
//   float meanCoefficient = 0;        // so we can just return this if something goes wrong.. 
//   //map<int, int>::iterator it;
//   for(int i=0; i < es; i++){
//     if(expts[i] < pset->allExptNo){
//       if(pset->exptLookup[expts[i]] != -1){
// 	eI[selEx] = pset->exptLookup[expts[i]];
// 	selEx++;
//       }
//     }
//   }
//   // it = pset->exptLookup.find(expts[i]);
//   //  if(it != pset->exptLookup.end()){
//   //    eI.push_back((*it).second);
//   //  }
//   //}
//   if(selEx < 2){         // then it doesn't make any sense,, and we can't calculate any of these..
//     return(meanCoefficient);
//   }
//   // and then just go throught the different things and calculate
//   int counter = 0;
//   float mean_mean = 0;
//   float mean_std = 0;
//   for(int i=0; i < pset->probeNo; i++){
//     float mean = 0; 
//     float std = 0; 
//     float sum = 0;
//     float squaredSum = 0;     // a bit much, but what the hell
//     for(int j=0; j < selEx; j++){
//       sum += pset->probes[i][eI[j]];
//       squaredSum += (pset->probes[i][eI[j]] * pset->probes[i][eI[j]]);
//     }
//     mean = sum / (float)selEx;
//     float SS = squaredSum - (sum * sum)/(float)selEx;
//     std = sqrt(SS/(float)(selEx -1));
//     mean_mean += mean;
//     mean_std += std;
//     // now check if the mean is above 0,,
//     //if(mean > 0){
//       //      if(std / mean < 10){
//       //meanCoefficient += (std/mean);
//       //counter++;
// 	//}
//     //}
//   }
//   delete []eI;
//   //  return(mean_std / (float)pset->probeNo);
//   if(mean_mean > 0){
//     return(mean_std/mean_mean);
//   }else{
//     if(mean_mean != 0){
//       return(-mean_std/mean_mean);   // it makes more sense,, but hopefully 
//     }else{
//       return(0);
//     }
//   }
//   //   if(counter > 0){
//   //     delete []eI;
// //     return(meanCoefficient/counter);
// //   }
// //   delete []eI;
// //   return(0);    // just in case.. !
}
	
	

void ConnectionObject::collectStats(){
  // first get the experimental indices from the lastMessage and then 
  // go through and obtain a range of stats for everything... -hmmmmmm
  vector<uint> eIndex;
  bool ok;
  uint e;
  int start = 0;
  int end = lastMessage.find(":", start);
  while(end != -1 && end != 0){
    QString is = lastMessage.mid(start, (end-start));
    start = end+1;
    e = is.toInt(&ok);
    if(ok){
      eIndex.push_back(e);
    }
    start = end +1;
    end = lastMessage.find(":", start);
  }
  if(eIndex.size() == 0){
    return;
  }
  uint* expts = new uint[eIndex.size()];
  for(int i=0; i < eIndex.size(); i++){ expts[i] = eIndex[i]; }
  ///////////// now simply go through the currently selected set of probe sets
  ///////////// and get the anova scores and all the other stuff... who yeahh.
  vector<float> anovaScores;             // the anova scores
  vector<float> varCoefficients;          // coefficients of variation
  vector<float> euclidQual;               // euclidean qualitites (self comparison)...
  vector<int> doneIndices;
  vector<string> statNames(3);      // rationalise at some point to make variable
  statNames[0] = "anova";
  statNames[1] = "euclidean";
  statNames[2] = "variation";     // to send to thingy.. 
  doneIndices.reserve(clientIndex.size());   // should really do the same with the others, but what the hell..
  anovaScores.reserve(clientIndex.size());
  varCoefficients.reserve(clientIndex.size());
  euclidQual.reserve(clientIndex.size());
  // don't like 0's or negative values.. 
  float minAnova = 99999;
  float minEuclid = 99999;
  float minVariation = 9999;

  for(int i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size() && pSet->data[clientIndex[i]]->index){
      //cout << "Getting stats for index : " << i << endl;
      doneIndices.push_back(clientIndex[i]);      // and then use this for thingy.. 
      anovaScores.push_back(anova(pSet->data[clientIndex[i]], expts, eIndex.size()));
      euclidQual.push_back(heavyCompare(pSet->data[clientIndex[i]], pSet->data[clientIndex[i]], expts, eIndex.size()));
      varCoefficients.push_back(variationCoefficient(pSet->data[clientIndex[i]], expts, eIndex.size()));

      //cout << anovaScores[i] << endl;
      // check if we need to update minimums..
      if(anovaScores.back() > 0 && anovaScores.back() < minAnova){ minAnova = anovaScores[i]; }
      if(euclidQual.back() > 0 && euclidQual.back() < minEuclid){ minEuclid = euclidQual[i]; }
      if(varCoefficients.back() > 0 && varCoefficients.back() < minVariation){ minVariation = varCoefficients[i]; }
    }
  }
  // then write the stuff to the netArray and send it off.. have to check how to do that, but lets see if it 
  // compiles and then check that I'm using the right functions!!
  sApp(string("<statCollection>"));
  qiApp(statNames.size());
  for(int i=0; i < statNames.size(); i++){
    sApp(statNames[i]);
  }
  qiApp(doneIndices.size());       // so we know how many to read..
  for(int i=0; i < doneIndices.size(); i++){
    qiApp(doneIndices[i]);
    if(anovaScores[i] > 0){
      fApp(anovaScores[i]);
    }else{
      fApp(minAnova);
    }
    if(euclidQual[i] > 0){
      fApp(euclidQual[i]);
    }else{
      fApp(minEuclid);
    }
    if(varCoefficients[i] > 0){
      fApp(varCoefficients[i]);
    }else{
      fApp(minVariation);
    }
  }
  sApp(string("<statCollectionEnd>"));
  // and send it off.. 
  writeArray();  
  delete []expts;
}

void ConnectionObject::compareCellGroups(){
  // First work out how many members in group a..
  vector<QString> words = splitString(lastMessage, '|');
  if(words.size() < 2){
    cerr << "compareCellGroups words is smaller than two,, eh.. no sense eh " << endl;
    return;
  }
  bool ok;
  int aSize = 0;
  aSize = words[0].toInt(&ok);
  if(!ok || aSize > pSet->experiments.size() || words.size() < 2+aSize || aSize < 1){
    cerr << "compareCellGroups couldn't get the size for the first image,, should make sure that it's not to big as well" << endl;
    cerr << "words size : " << words.size() << endl
	 << "aSize      : " << aSize << endl
	 << "words are : " << endl;
    for(int i=0; i < words.size(); i++){
      cerr << words[i] << "  ";
    }
    cerr << endl;
    return;
  }
  vector<int> a(aSize);
  for(int i=0; i < aSize; i++){
    int exp = words[i+1].toInt(&ok);
    if(!ok){
      cerr << "Couldn't get one of the experiments for the a set returning  " << endl;
      return;
    }
    a[i] = exp;
  }
  int bSize = 0;
  bSize = words[aSize+1].toInt(&ok);
  if(!ok || bSize > pSet->experiments.size() || words.size() < 2+aSize+bSize || bSize < 1){
    cerr << "compareCellGroups, something not so good with bSize " << endl;
    cerr << "words size : " << words.size() << endl
	 << "aSize      : " << aSize << endl
	 << "bSize      : " << bSize << endl
	 << "words are : " << endl;
    for(int i=0; i < words.size(); i++){
      cerr << words[i] << "  ";
    }
    cerr << endl;
    return;
  }
  vector<int> b(bSize);
  for(int i=0; i < bSize; i++){
    int exp = words[i+2+aSize].toInt(&ok);
    if(!ok){
      cerr << "Couldn't get one of the experiments for the a set returning  " << endl;
      return;
    }
    b[i] = exp;
  }
  // Ok, I think that we can now go through and actually do the comparison,, ofcourse using all of the careful things that
  // we always use...
  int* aIndex = new int[a.size()];
  int* bIndex = new int[b.size()];   // could probably use vectors without a problem.. 
  float* avalues = new float[a.size()];
  float* bvalues = new float[b.size()];
  vector<float> correlates;          // the factors.. 
  vector<float> deltas;
  vector<float> compDeltas;     // mean square difference divided by the sum of the values.. or so
  vector<float> aCorrelates;    // a funky correlation score..
  vector<float> bCorrelates;    // a funky internal correlation score.. 
  vector<float> flatDistances;  // more or less the one used by the flat expt compare process.. 
 
  vector<int> doneIndices;      // as if we are dealing with distances we can't just assign 0 for nonfunctional ones.. 
  deltas.reserve(clientIndex.size());
  compDeltas.reserve(clientIndex.size());
  doneIndices.reserve(clientIndex.size());
  correlates.reserve(clientIndex.size());
  aCorrelates.reserve(clientIndex.size());
  bCorrelates.reserve(clientIndex.size());
  flatDistances.reserve(clientIndex.size());
  cout << "Created the correlates vector " << endl;
  for(int i=0; i < clientIndex.size(); i++){
    //cout << "client Index " << i << "   " << clientIndex[i] << endl;
    //correlates[i] = 0;        // the default, then increment from here..
    // first we have to work out which indices we are going to use.. 
    int aUsed = 0;
    int bUsed = 0;
    probe_set* probe = pSet->data[clientIndex[i]];     // this is much easier to read.. 
    //cout << " a : ";
    for(int j = 0; j < a.size(); j++){
      //cout << a[j] << " -> ";
      if(a[j] < probe->allExptNo){
	if(probe->exptLookup[a[j]] != -1){
	  aIndex[aUsed] = probe->exptLookup[a[j]];
	  //cout << aIndex[aUsed] << "   ";
	  aUsed++;
	}
      }
    }
    //    cout << endl << " b : ";
    for(int j = 0; j < b.size(); j++){
      //cout << b[j] << " -> ";
      if(b[j] < probe->allExptNo){
	if(probe->exptLookup[b[j]] != -1){
	  bIndex[bUsed] = probe->exptLookup[b[j]];
	  //cout << bIndex[bUsed] << "   ";
	  bUsed++;
	}
      }
    }
    //    cout << endl;
    //    cout << "got the indices... aUsed : " << aUsed << "  bUsed : " << bUsed << endl;
    ///// if everything ok, we have to copy the probes and normalise them, and then take the mean values.. 
    ////  there is a more effective way of doing this, but using the function will do for now..
    if(!aUsed || !bUsed){
      cerr << "Didnd't get any values for this one, no point in doint the rest for god's sake... " << endl;
      continue;
    }
    float** probes = copyProbes(probe->probes, probe->probeNo, probe->exptSize);
    // remember to delete when we are done.....
    // normalise the values in the probes and make ourselves a mean vector...
    float* mean = new float[probe->exptSize];
    // set to 0..
    //cout << "made probes and mean : " << endl;
    for(int j=0; j < probe->exptSize; j++){
      mean[j] = 0;
    }
    //cout << "set the means to 0:L " << endl;
    for(int j=0; j < probe->probeNo; j++){
      zScore(probes[j], probe->exptSize);
      for(int k=0; k < probe->exptSize; k++){
	mean[k] += probes[j][k];
      }
    }
    for(int j=0; j < probe->exptSize; j++){
      mean[j] = mean[j]/(float)probe->probeNo;    // rather important to make sure I cast it to a float..
    }
    //cout << "Got the mean values will now go and do something else : " << endl; 
    //// OK,, now I have the data that I want to look at, I can now look at doing the actual math...
    if(aUsed && bUsed){
      int pos = correlates.size();
      correlates.push_back(0);
      for(int j=0; j < aUsed; j++){
	for(int k=0; k < bUsed; k++){
	  //cout << aIndex[j] << "x" << bIndex[k] << "   ";
	  correlates[pos] += (mean[aIndex[j]] * mean[bIndex[k]]);
	}
      }
      correlates[pos] = correlates[pos]/(float)(aUsed*bUsed);
      // the above calculates some sort of correlation coefficient, but it isn't exactly what we want...
      // try calculating the T -statistic for the two groups. The T statistic for uncoupled difference or two way thingy
      // is simply the difference in the mean of the two groups divided by the pooled variance of the two groups.. 
      // the pooled variance in the group is the sum of the variances adjusted by the degree of fredoom each has..
      // Note that I am going to use a bastardised equation in this case,, as we want to be able to cope with a group of genes containing
      // only one sample.. -- so rather than dividing by N-1, I will divide by N.
      // the equatione for pooling variance for a and b (as and bs) is ..
      float amean, bmean, asqsum, bsqsum;    // get these assigned by the wonderful thingy..
      //float* avalues = new float[bUsed];
      //float* bvalues = new float[aUsed];
      for(int j=0; j < aUsed; j++){ avalues[j] = mean[aIndex[j]]; }
      for(int j=0; j < bUsed; j++){ bvalues[j] = mean[bIndex[j]]; }      // ok, so that's ugly, and stuff, but there you go.. never mind..
      squaredSum(avalues, aUsed, amean, asqsum);
      squaredSum(bvalues, bUsed, bmean, bsqsum);
      // and lets work out the sum of the squares for each one.. already done this, but..
      float aSquared =0;
      float bSquared= 0;
      float adeltaSquared = 0;
      float bdeltaSquared = 0;
      for(int j=0; j < aUsed; j++){
	aSquared += (avalues[j] * avalues[j] * avalues[j] * avalues[j]);
	//	adeltaSquared += (avalues[j] - amean) * (avalues[j] - amean);
	adeltaSquared += fabs(avalues[j] - amean);
      }
      for(int j=0; j < bUsed; j++){
	bSquared += (bvalues[j] * bvalues[j] * bvalues[j] * bvalues[j]);
	//	bdeltaSquared += (bvalues[j] - bmean) * (bvalues[j] - bmean);
	bdeltaSquared += fabs(bvalues[j] - bmean);
      }
      
      aCorrelates.push_back(adeltaSquared/aSquared);
      bCorrelates.push_back(bdeltaSquared/bSquared);  // hmmm.. maybe interesting.. 
      // and then divide the difference in means by thingies.........
      float tscore = amean-bmean/ ( sqrt( ((asqsum + bsqsum)/(aUsed + bUsed)) * ( (float)1/(float)aUsed + (float)1/(float)bUsed) ));  // which looks pretty ugly, but.. 
      //      cout << "amean : " << amean << "\t bmean :" << bmean << "\tasqsum : " << asqsum << "\tbsqsum : " << bsqsum << "\t\t T score : " << tscore << endl;
      // I just want to have a look at some values ..
      float pooledStd = sqrt( (((asqsum * (float)aUsed) + (bsqsum * (float)bUsed))/(aUsed + bUsed)) * ((float)1/(float)aUsed + (float)1/(float)bUsed) );
//       cout << "amean                         : " << amean << endl;
//       cout << "bmean                         : " << bmean << endl;
//       cout << "Pooled Std Deviation          : " << pooledStd << endl;
//       cout << "std deviation over difference : " << pooledStd / (amean-bmean) << endl;
      float ratio = pooledStd / fabs(amean-bmean);
      //// now the tscore is not the ideal score for us. We are rather looking for something that indicates information content, which doesn't do.
      /// hence we will modify it in the future, but for now, let's just return that and see how it behaves. but Since I don't really want 
      deltas.push_back(fabs( (amean-bmean)/pooledStd) );
      // comp deltas is the mean difference multiplied by the sum of the squares..
      float aSumSquared = 0;
      float bSumSquared = 0;
      //     cout << "Adding up stuff for things " << endl;
      for(int j=0; j < aUsed; j++){ aSumSquared += (avalues[j]  * avalues[j]); }
      // cout << "Added up for a, and about to do for b " << endl;
      for(int j=0; j < bUsed; j++){ bSumSquared += (bvalues[j] * bvalues[j]); }
      //cout << "got the thingy.. " << endl;
      //compDeltas.push_back((amean-bmean)*(amean-bmean)*(aSumSquared + bSumSquared));    // hmm, should be ok.. 
      compDeltas.push_back((amean-bmean)*(amean-bmean) * (1 - ratio/(1+ratio)) );    // hmm, should be ok.. maybe ... hmm..  
      //cout << "pushed back compdeltas value : " << (amean-bmean)*(amean-bmean) * (1 - ratio/(1+ratio)) << endl;
      float flatD = doFlatExptCompare(aIndex, aUsed, bIndex, bUsed, mean, probe->exptSize, 3.667, 0.75);
      cout << "Flat distance i : " << i << "   value : " << flatD << endl;
      flatDistances.push_back(flatD);
      //   flatDistances.push_back(doFlatExptCompare(aIndex, aUsed, bIndex, bUsed, mean, probe->exptSize, 3.667, 0.75));
      doneIndices.push_back(clientIndex[i]);
      //cout << "pushed back done indices with client index for i : " << i << endl;
    }
    //cout << "correlates  is : " << correlates[i] << endl;
    // and make sure to delete probes and mean..
    // this is a bit wasteful, but probably the safest option..
    delete []mean;
    //cout << "deleted mean" << endl;
    for(int j=0; j < probe->probeNo; j++){
      delete []probes[j];
    }
    //cout << "delete probe things" << endl;
    delete []probes;
    //cout << "deleted probes : " << endl;
  }
  delete []aIndex;
  delete []bIndex;
  delete []avalues;
  delete []bvalues;
  // and then we just have to write the message to the client, and let it know what we found out... !! ho ho  ho yeaah..
  // have a feeling that there is something that I need to sort out, but can't think what...
  sApp("<statCollection>");
  qiApp(6);
  sApp("Cell gene correlation stuff");
  sApp("T score ");
  sApp("Compensated Difference");
  sApp("A internal");
  sApp("B internal");
  sApp("Flat delta");
  qiApp(doneIndices.size());
  for(int i=0; i < doneIndices.size(); i++){
    qiApp(doneIndices[i]);
    fApp(correlates[i]);
    fApp(deltas[i]);
    fApp(compDeltas[i]);
    fApp(aCorrelates[i]);
    fApp(bCorrelates[i]);
    fApp(flatDistances[i]);
  }
  sApp("<statCollectionEnd>");
  writeArray();
}


void ConnectionObject::rawComparisonSort(){
  /// take a vector of floats (already normalised).. -
  /// then compare these against everything that is in client Index..
  /// and  sort.. Do a comparison against individual probes,, (normalised).

  // first parse the lastMessage..
  vector<uint> eIndex;
  vector<float> values;
  bool v_ok;
  bool i_ok;
  float v;
  uint ei;
  vector<QString> words = splitString(lastMessage, ':');
  //  cout << "rawComparisonSort words size is " << 
  if(words.size() % 2 == 0 || words.size() == 0){
    cerr << "rawComparisonSort , words size doesn't make any sense: " << endl;
    return;
  }
  for(int i=0; i < words.size()-1; i += 2){
    v = words[i].toFloat(&v_ok);
    ei = words[i+1].toInt(&i_ok);
    if(v_ok && i_ok){
      values.push_back(v);
      eIndex.push_back(ei);
    }
  }
  // The last word should be a bool determining whether or not we return the distribution.. hmm
  // hmmmm... 
  bool distribution = (words[words.size()-1] != "0");
  cout << "values.size is " << values.size() << endl
       << "eIndex size is " << eIndex.size() << endl
       << "and the guy wants the distribution ? " << distribution << endl;
//   int start = 0; 
//   int end = lastMessage.find(":", start);
//   while(end != -1 && end != 0){
//     QString vs = lastMessage.mid(start, (end-start));
//     start = end+1;
//     end = lastMessage.find(":", start);
//     QString is = lastMessage.mid(start, (end-start));
//     v = vs.toFloat(&v_ok);
//     i = is.toInt(&i_ok);
//     if(v_ok && i_ok){
//       values.push_back(v);
//       eIndex.push_back(i);
//     }
//     start = end +1;
//     end = lastMessage.find(":", start);
//     cout << "rawComparisonSort values size : " << values.size() << endl;
//   }
  ////////// which I suppose could fall over and do bad stuff, but hopefully it will be OK
  /////////  Now to actually work out what is what.. Remembering that we have to look up the experimental
  /////////  indices as we go along..
  if(values.size() < 1) { return; }
  
  float* va = new float[values.size()];
  uint* e = new uint[eIndex.size()];   // which should be the same..
  for(int i=0; i < values.size(); i++){
    va[i] = values[i];
    e[i] = eIndex[i];
  }
  //////////// ugly, I have to fix it later on.. hmmm... 

  vector<dist_set> distances;
  distances.reserve(clientIndex.size());
  float d;
  for(int i=0; i < clientIndex.size(); i++){
    if(clientIndex[i] < pSet->data.size()){
      d = singleEuclidCompare(pSet->data[clientIndex[i]],  va, e, values.size());
      //      d = singleEuclidCompare(&pSet->data[clientIndex[i]],  values, eIndex);
      if(d >= 0){
	distances.push_back(dist_set(pSet->data[clientIndex[i]]->index, d));
      }
    }
  }
  // and if we actually get something after all of that,, we have to check it out..
  sort(distances.begin(), distances.end(), comp_set());
  vector<int> tIndex(distances.size());
  for(int i=0; i < tIndex.size(); i++){
    tIndex[i] = distances[i].index;
  }
  writeIndex(tIndex, "Raw Profile Search");
  if(distribution){
    cout << "Sending the Stat collection " << endl;
    sApp("<statCollection>");
    qiApp(1);
    sApp("E. Distance");
    qiApp(distances.size());
    for(int i=0; i < distances.size(); i++){
      qiApp(distances[i].index-1);   // now that is a bit ugly.. hmm 
      fApp(distances[i].value);
    }
    sApp("<statCollectionEnd>");
    writeArray();
  }
  
  delete []va;
  delete []e;
}

void ConnectionObject::qiApp(int value){
  //  Q_INT32 q = (Q_INT32)value;
  int q = htonl(value);
  char* c = (char*)&q;
  //cout << "value is : " << value << endl;
  for(int i=0; i < 4; i++){
    //cout << "\t" << i << " : " << (int)*(c+i) << endl;
    sendData->app(*(c + i));
  }
}

void ConnectionObject::iApp(int value){
  char* c = (char*)&value;
  for(int i=0; i < 4; i++){
    sendData->app(*(c+i));
  }
}

void ConnectionObject::sApp(char* c, int l){
  //// please check before using this function. easy to segment fault if incorrectly used. 
  qiApp(l);
  for(int i=0; i < l; i++){
    sendData->app(c[i]);
  }
}

void ConnectionObject::sApp(string s){
  // add the size including the thingy.. 0 terminator
  // and put the size at the beginning.. (as a Q_INT32)
  //Q_UINT32 wl = (Q_UINT32)(s.size());
  //char* c = (char*)&wl;
  //for(int i=0; i < 4; i++){
  //  sendData->app(*(c+i));
  //}
  qiApp(s.size());
  //  qiApp(s.size()+1);  // may need to make this a Q_UINT32 ,, but lets see for now
  //cout << "socketNumber << " << socketNumber << "  appending string : " << s << endl;
  for(int i=0; i < s.size(); i++){
    //cout << " " << s[i]; // << endl;
    sendData->app(s[i]);
  }
  //cout << endl << "socketNumber : " << socketNumber << "  appended string done " << endl;
  //  sendData.push_back('\0');
}

void ConnectionObject::fApp(float f){
  char* c = (char*)&f;
  for(int i=0; i < 4; i++){
    sendData->app(*(c + i));
  }
}

void ConnectionObject::dApp(double d){
  char* c = (char*)&d;
  for(int i=0; i < 8; i++){
    sendData->app(*(c + i));
  }
}

void ConnectionObject::catchAlarm(int value){
  //cout << "Alarm caught value is : " << value << endl;
  //cout << "my socket number is " << socketNumber << endl;
}
