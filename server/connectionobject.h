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

#ifndef CONNECTIONOBJECT_H
#define CONNECTIONOBJECT_H

//#include <qsocket.h>
#include <qsocketdevice.h>
#include <qstring.h>
#include <qthread.h>
#include <qmutex.h>
#include <qcstring.h>
#include <set>
#include <string>
#include <vector>
#include <signal.h>            // need to be able to handle or ignore signals of different sorts. ?? 
#include <time.h>              // since we store the sessionId as a value of time_t
//#include "processor.h"
#include "../netArray/netArray.h"
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include "../util/pathTracer.h"
#include "kClusterProcess.h"
#include "blastClient.h"
#include <qwidget.h>            // so I can send events to the server process. 

struct ProbeSetThresholds {
  int maxExpect;
  int minLength;
  float minMatch;
  ProbeSetThresholds(){
    maxExpect = -10;
    minLength = 35;
    minMatch = 0.9;
  }
  ProbeSetThresholds(int me, int ml, float mm){
    maxExpect = me;
    minLength = ml;
    minMatch = mm;
  }
  void setThresholds(int me, int ml, float mm){
    maxExpect = me;
    minLength = ml;
    minMatch = mm;
  }
};

struct dist_set {
  int index;
  float value;
  dist_set(int i, float v){
    index = i;
    value = v;
  }
  dist_set(){
  }
};

class ConnectionObject : public QObject, public QThread
{
  //  Q_OBJECT
  // a class that holds and processes the communication with 
    public:
  ConnectionObject(QSocketDevice* s, ProbeSetSet2* p, QWidget* sp);
  ~ConnectionObject();        // not sure what to put in there, but probably we will need to put something..

  // some boolean flags that can be changed by the parent server process..
  bool updateUsers;    // send updated usertable to the thingy.. 
  int commandCounter;  // number of commands served.. -public so that the server process can check when it closes the client. 

 private:

  typedef struct sessionData {
    sessionData(){      // do nothing..
    }
    int sessionIndex;
    int userId;     // important
    string sessionTitle;
    string description;
    set<string> keyWords;
    set<int> geneIds;
  };
  QSocketDevice* socket;    // not object based so may be better.. than the socket.. -- but anyway, I think that I can soon get rid of this
                            // and just use the socket number,, - and close it manually or something like that.
  int socketNumber;         // the number of the socket, just for convenience.. 
  struct sigaction sa;             // a signal handler,, set up the things in the constructor... 
                                   // this definition could just as well have been made in the constructor itself, but maybe
                                   // there is some advantage for it to be known within the whole object. 
  // one wonders if one could also set up a signal handler to handle segmentation faults gracefully.. hmm -- questionable if one wants to do this though. 
  struct sigaction alarmAction;      /// I think I don't use this,, as it seems there is no easy way of combining the use of signals and threads.
  ProbeSetSet2* pSet;        // this is shared between clients, I should probably provide accessor functions that check 
                             // relevant mutexes, rather than the direct access I use at the moment.
  //  vector<probe_data>* data; 
  //QString name;        // it seems that we don't use this !!.. 

  int userId;       // the user id from the users table..
  time_t logInTime; // the time when we logged in. Use this as a unique identifier for this person.
  string currentUserName;  // the current User Name.. long name as I use userName all over the place.. bugger.. that.. 
  QString addressString;  // the address as a string.. 
  QWidget* server;   // the parent server process for sending events to.. 
  
  const char* conninfo;
  int bufSize;      // size of sending things.. 

  QString currentMessage;
  QString lastMessage;
  QString currentCommand;
  QString lastCommand;
  vector<uint> clientIndex;     // read from the client and used for statistical operations carried out on a subset of genes.
  set<int> clientChips;         // active chips, send only data from chips which are active.. (when setting the index for the client).

  vector<RegionSpecification> clientRegions;    // a set of region Specification returned by some sort of database operation.. 
  ProbeSetThresholds pSetThresholds;            // determine whether or not to send probe sets contained in a region.. 
  //map<int, vector<QByteArray> > clientRequests; // if not 
  //set<void*> clusterSets;                 // can be created and things. should always be empty, but if not, then send one.. 
  set<void*> clusterProcesses;             // so I can delete the process when I take care of the data in the clusterSet. 

  set<PathTracer*> exptTracers;        
  QMutex exptTracerMutex;                  // hmm, why do I use void pointers below, not really sure.. 
  set<void*> compareExperimentProcesses;   // as above.. -both of these should probably be private.. but have some reason perhaps for keeping them here. 
  set<void*> flatExptComparers;        // hmmmm.. -- not sure if we need to separate these different things.. but who knows.. 
  set<BlastClient*> blastClients;      // just go through these when checking flags, and see if any are done, if done, delete and get rid of.. 
  QMutex clusterMutex;
  QMutex experimentComparisonMutex;
  QMutex flatExptCompareMutex;
  //vector<char> sendData;                       // pack the data to send in this.. !! instead.. 
  NetArray* sendData;         // NO use this instead.. 
  NetArray* inData;           // for incoming data.. 

  int commandState;            // command processing state. 
  bool authenticated;
  QWaitCondition dataWait;
  bool stayConnected;

  // Log in functions ..
  void logIn();
  void logOut();  // these change the tables user_log
  int dbCommand(string cmd);  // connect to the database and execute the command. Return CmdTuples or -1 if some problem

  // some commands..   -- 
  void checkFlags();              // check the boolean flags to see if we need to update the user with something.. 
  void parseCommand();                        // parse the contents and do one of the following things..
  void sendProbeSet();                     // send a probe Set to the client
  void sendRegion();
  void doDBLookup();       // look up query in database, write as separate function !!
  void doGenDBLookup();    // look up genomic regions as opposed to database positions.. 

  void writeStatus(statusMessage message);      // Tell the client if something worked or not, and give information as to why. 
  void writeDBChoices();   // tell the client what kinf of databaes lookups it can do.. 
  void writeProbeSet(probe_set* pset);
  void writeProbeData(probe_data* pdata);     // write the probe data to a socket.. 
  void writeFileInfo();
  void writeExptInfo();
  void writeChipInfo();    // send the data for the chips.. 
  //void writeString(string s);                 // convenience, write a string to s socket..
  //void writeInt(int i);
  //void writeFloat(float f);
  void writeIndex(vector<uint> v, string term, bool setClientIndex=1);              // used by functions that changes the index.. 
  void writeIndex(vector<int> v, string term, bool setClientIndex=1);
  void doFlush();                              // just send some bits down the pipe to flush things out.. 
  void doChangePassword();                     // try to change the password.. 
  void doCreateNewUser();                      // try to create a new user.. 
  bool createSession();                        // create a new session.. should be quite simple.. 
  bool newAnnotation();                        // insert a new user annotation.. 
  bool updateSessionInformation();             // update the description for a session, only for the owner of that session,, hmmm.. 
  bool updateAnnotation();
  void sendSessionInformation();               // get the information on the sessions present in the database and send it to the client.
  void sendUserInformation();
  void sendGenomeSequence();
  void sendProbeSetSequence();
  void sendEnsemblPeptide();
  void sendEnsemblTranscriptSequence();
  void doKCluster();                          // start a kClustering thing.. 
  void sendKClusters();
  void traceExperiments();                    // find a reasonable trace through the different experiments.. 
  void doExperimentCompare();                 // compare experiments.. -- maybe it's not that heavy.. !!  (although
  void sendExperimentDistances();             // send the distances.. 
  void sendExperimentTrace(PathTracer* tracer);           // se
  void doFlatExptCompare();                   // use a different algorithm,, I should really merge these but my head is no good today.. 
  float doFlatExptCompare(int* a, uint as, int* b, uint bs, float* v, uint vs, float order, float sigma);   // compare the indicated experiments as groups and return 
  void sendFlatExptDistances();               // ugly.. should use same function, but I'm a bit .. off colour today. 
  void devsFromMean();                       // calculate standard deviations from the mean profile for all probe pair profiles --> send to statviewer..
  void exportMeans();                        // export means of z-score normalised stuff. later we can redo this to read the lastMessage, but first let's make simple..
  

  void sendGenomicRegionProbeSetMatches(string chr, int start, int end);    // send the data for the things.. !
  void sendGenomicRegion(string chrom, int start, int stop, int association, int target);
  void sendGenomicRegion(string chrom, int start, int stop, int association, int target, vector<int>& psets, bool useT);   // fill in which probe sets have been included
  void sendIshThumbnails();    // parse the string and send in one go. it should be quite easy.. 
  void sendIshThumbnail(int index, int ImageIndex, int oid);
  void sendIshImage();        // parse the string, and send.. in the same function,, should be only one image.. 
  void parseIshProbeRegionRequest();  // just get the int..
  void sendIshImageRegions(int index, int regionSize);   // separate into two functions so that I can reuse this function and clean up the code (fat bloody chance of that eh??).
  void sendIshProbeData(int index);     // look up the data struct and send the ish image.. 
  void appendIshAnnotation(ish_annotation& a);  // just append,, -- convenience function.. 
  void associateIshProbeWithEnsemblGene();
  void setIshProbeName();         
  void insertIshProbeTextAnnotation();
  void updateIshProbeTextAnnotation();
  void insertIshProbeFloatAnnotation();
  void insertIshProbeClassification();
  void sendIshTextFields();
  void sendIshFloatFields();
  void sendIshClasses();

  //// some stuff for handling protocols..
  void sendProtocols();        // sends all the protocols in the database..
  void sendProtocol();         // parses the request and sends individual protocols.. 
  void receiveNewProtocol();   // parse a new protocol and do something with it.. 

  /// and for experiments  (mostly ish_experiments, but the name of the table isn't too important.. 
  void commitExperimentToDB();
  void sendExperiments();

  void sendTissues();
  void sendIshAnnotationFields();

  void makeIshTissue();
  void makeIshAnnotationField();

  void commitIshImageToDB();
  bool commitIshImageCommentsToDB(vector<Comment> comments, int imId);    // and add comments to the appropriate tables.. if I can work that out.

  /// some stuff for handling requests to blast things.. 
  void doBlast();       // initially just for testing.. 

  void parseGenomicRegionRequest();
  void sendGenomicRegionEnsemblGenes(int requestId, string chr, int start, int end);
  vector<QString> splitString(QString& word, char delimiter);

  void iApp(int value);       // append the value to the vector.. SEND as Q_INT32 ??
  void qiApp(int value);   // do the same with a 
  void sApp(string s);        // append a string..
  void sApp(char* c, int l);  // append a char array of length l, starting at c,, -kind of dangerous really.. 
  void fApp(float f);           // append a float... (do not resize.. but do suggest reserver..)
  void dApp(double d);

  void writeArray();            // writes the array to the socket. the only thing that does!!
  ssize_t writen(const char* dptr, size_t n);
  //void sendArray(vector<QByteArray> data);
  //  void sendArray(int n);

  void readClientIndex();                    // not really read as its more a case of parsing the lastCommand variable.. and setting the clientIndex vector
  void setClientChips();                      // set all the client chips. erase everything.. then insert the stuff if appropriate. 
  bool authenticate();
  bool connected();         // just a convenience function,, checks the state of the socket..
  
  float** copyProbes(float** p, int ps, int es);   // ps- probe no,  es -experiment size, 
  void delProbes(float** p, int ps);

  void stdDeviation(float* v, uint s, float& mean, float& std, bool sampling=true);   // assign mean and std deviation to the references..
                                                                               // if sampling is true, then do proper thing and
                                                                               // divide by n-1, otherwise divide by n (in case we may have groups of 1).
  void squaredSum(float* v, uint s, float& mean, float& sqsum); 
  void zScore(float* v, int s);
  float euclidean(float* v1, float* v2, uint s);
  float heavyCompare(probe_set* pset1, probe_set* pset2, uint* expts, uint es);
  float diff(probe_set* p,  uint up, uint down, bool normalised);     // do a binary comparison between two experimental points.. 
  void diffSort(uint up, uint down, bool normalised);
  void doDiffSort();   // for parsing the last message.. 
  float singleEuclidCompare(probe_set* pset1, float* target, uint* expts, uint es);
  float meanEuclidCompare(probe_set* p, float* target, uint* expts, uint es, bool normed);
  void meanComparisonSort();    
  void euclidSort(probe_set* pset1, vector<uint> expts);
  void rawComparisonSort();
  float anova(probe_set* pSet, uint* expts, uint es);
  float variationCoefficient(probe_set* pset, uint* expts, uint es);
  void anovaSort();        // uses clientIndex, -and uses all of the experiments,, -so just passes the thing.. to the thing.. 
  void anovaSelectSort();  // looks at the lastMessage to determine which experiments to utilize.. 
  void collectStats();     // reads the last message to determine which experiments to take into consideration. Then obtains a number of stats for each one which it writes to the thingy.. 
  void doEuclidCompare();
  void compareCellGroups();              // calculate the correlation coefficient for each gene for two groups of cells.. ahh, just read the source to work out what it means.. 
         
  void readClient();                        // read from the client, and try to work out what's going on..

  static void catchAlarm(int value);               // catches the alarm signal, but doesn't do very much with it.. 
  void reportError(string errorMessage, int errorNumber);

 protected:
 void run();

};

#endif
 
 
 
