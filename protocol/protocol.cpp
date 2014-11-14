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

#include "protocol.h"
#include "../netArray/netArray.h"
#include <qdatetime.h>
#include <vector>
#include <string>
//#include <libpq++.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <libpq-fe.h>

using namespace std;

ProtocolStep::ProtocolStep(){
  stepId = 0;
  stepDescription = "empty Description";
  parentId = 0;
  creatorId = 0;
  creatorName = "NULL";
  filled = false;
  newStep = true;   // no sure what I want to do with this, but what the hell.. 
}

ProtocolStep::ProtocolStep(int sid, const char* conninfo){
  stepId = sid;
  newStep = false;
  filled = false;
 //// This step should already exist. if it doesn't then we have some sort of a problem as we will not be able to 
  //// to insert it into the database. Hence we can argue that regardless of whether or not we want to fill it in, we might
  //// want to check whether or not it is in the database... -- and if it is not there.. then do something.. ---
  //// however, that is likely to take some time, and as the mistake is likely to be discovered at some later stage during an 
  //// an insert maybe we can just leave it.. and see if we catch the error later on. The trouble here is likely to be with 
  //// error reporting...  .. 
  if(conninfo){
    if(!fillFromDB(conninfo)){
      cerr << "Some trouble filling in the data base details " << endl;
    }
  }  
}

ProtocolStep::ProtocolStep(string des, int pid, int cid, string cname, QDateTime ctime, int sid, const char* conninfo){
  stepId = sid;
  stepDescription = des;
  parentId = pid;
  creatorId = cid;
  creatorName = cname;    // this is likely to be redundant.. as we may not need this anyway.. 
  creationTime = ctime;
  filled = true;
  if(!stepId){
    newStep = true;
  }
  // which is pretty much all for this.. -- generally this will be followed by an entry into the database.. 
  if(conninfo){
    //    PgDatabase conn(conninfo);
    if(!dbEnter(conninfo)){
      cerr << "some trouble entering into the database " << endl;
    }
  }
}


ProtocolStep::ProtocolStep(int sid, int pid, int cid, string des, const char* conninfo){
  // this is used for entering into the database..
  stepId = sid;
  parentId = pid;
  if(!parentId){ parentId = 1; }    // I am being redundant.. -at the moment I'm checking for this in the 
                                    // the dbEnter function as well.. -- as I'm not sure what the best place
                                    // to check is. One could argue this is up to the user.. but .. might as 
                                    // well be a bit more flexible.. 
  creatorId = cid;
  stepDescription = des;
  filled = false;              // as it doesn't actually have all the values.... wouldn't be any good to send back.. 
  if(!sid){        // i.e. 0 value..
    //    PgDatabase conn(conninfo);
    if(dbEnter(conninfo)){
      newStep = false;
    }else{
      cerr << "couldn't enter step into the database what's up? " << endl;
      newStep = true;
    }
  }else{
    newStep = false;
  }      ////////// HORRIBLY BADLY WRITTEN.. 
}

bool ProtocolStep::fillFromDB(const char* conninfo){
  ostringstream query;
  query << "select a.*, b.user_name from protocol_step where step=" << stepId << " and a.user_id=b.index";
  PGconn* conn = PQconnectdb(conninfo);
  if(!conn || PQstatus(conn) != CONNECTION_OK){
    PQfinish(conn);
    cerr << "fillFromDb unable to connect to database" << endl;
    return(false);
  }
  PGresult* res = PQexec(conn, query.str().c_str());
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    PQclear(res);
    PQfinish(conn);
    return(false);
  }

  /*
  PgDatabase conn(conninfo);
  if(!conn.Exec(query.str().c_str())){
    cerr << "Couldn't execute query for step information.. don't know why.. " << endl;
    return(false);
  }
  if(conn.Tuples() != 1){
    cerr << "Conn tuples is an inappropriate number .. " << endl;
    return(false);
  }
  */
  // and now try to get all the appropriate values..
  parentId = atoi(PQgetvalue(res, 0, 1));
  stepDescription = PQgetvalue(res, 0, 2);
  creatorId = atoi(PQgetvalue(res, 0, 3));
  creationTime = QDateTime::fromString( PQgetvalue(res, 0, 4) );
  creatorName = PQgetvalue(res, 0, 5);
  filled = true;
  return(true);

  /*
  parentId = atoi(conn.GetValue(0, 1));
  stepDescription = conn.GetValue(0, 2);
  creatorId = atoi(conn.GetValue(0, 3));
  creationTime = QDateTime::fromString(conn.GetValue(0, 4), Qt::ISODate);    // I hope.. -I'm not sure though..
  creatorName = conn.GetValue(0, 5);
  // and if we get here..
  filled = true;
  return(true);
  */
}

bool ProtocolStep::dbEnter(const char* conninfo){
  if(!conninfo){
    cerr << "ProtocolStep, no connection information" << endl;
    return(false);
  }
  PGconn* conn = PQconnectdb(conninfo);
  if(PQstatus(conn) != CONNECTION_OK){
    cerr << "ProtocolStep::dbEnter connection to DB failed" << endl;
    return(false);
  }

  ostringstream entry;
  char* escDescription = new char[stepDescription.size()*2 +1];
  PQescapeString(escDescription, stepDescription.c_str(), stepDescription.size());
  if(parentId == 0){ parentId = 1; }    // just so that the toplevel steps have their things.. going ok.. 
  entry << "insert into protocol_step (parent_step, description, user_id) values (" << parentId << ", '" << escDescription << "', " << creatorId << ")";
  delete escDescription;
  // the above could be written a bit clearer.
  PGresult* res = PQexec(conn, entry.str().c_str());
  int cmdTuples = atoi(PQcmdTuples(res));
  
  if(cmdTuples != 1){
    cerr << "CmdTuples != 1 but: " << cmdTuples;
    PQclear(res);
    PQfinish(conn);
    return(false);
  }

  /*
  if(!conn.Exec(entry.str().c_str())){
    cerr << "Couldn't exec the step dbentry " << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  if(conn.CmdTuples() != 1){
    cerr << "CmdTuples is not 1 but something else,, hmmm " << conn.ErrorMessage() << endl;
    return(false);
  }
  */
  Oid oid = PQoidValue(res);
  //  const char* oid = conn.OidStatus();     // which is bloody stupid, but I can't find a method that lets me get at this information.. any other way.. 
  // then use this to select..
  ostringstream query;
  query << "select step from protocol_step where oid = " << oid ;  // a number presumably
  
  res = PQexec(conn, query.str().c_str());
  if(PQresultStatus(res) == PGRES_TUPLES_OK){
    stepId = atoi( PQgetvalue(res, 0, 0) );
  }
  PQclear(res);
  PQfinish(conn);
  newStep=false;
  return(true);
  /*
  if(!conn.Exec(query.str().c_str())){
    cerr << "Couldn't exec the protocol step,, hmm, trouble man " << conn.ErrorMessage() << endl;
    return(false);
  }
  if(conn.Tuples() != 1){
    cerr << "Conn tuples is not 1 bugger : " << conn.ErrorMessage() << endl;
    return(false);
  }
  stepId = atoi(conn.GetValue(0, 0));
  newStep = false;
  return(true);
  */
}

void ProtocolStep::serialise(NetArray* na){
  // write the useful information to the thingy..
  na->iapp(stepId);
  na->iapp(parentId);
  na->iapp(creatorId);
  na->sapp(creatorName);
  na->sapp(stepDescription);
  na->iapp(creationTime.toTime_t());   
  // which is actually everything.. as the client anyway doesn't care about anything..
}

/////////////// and the Protocol functions... god, this is boring..

Protocol::Protocol(){
  protocolId = 0;
  creatorId = 0;
  conninfo = 0;
}

Protocol::Protocol(const char* cinfo, int pid, bool fillsteps){    // look up all the data in the database,
  conninfo = cinfo;     // why ?? 
  cout << "trying to make a Protocol for id : " << pid << endl;
  protocolId = 0;         // on success give it the right number .. 

  PGconn* conn = PQconnectdb(conninfo);
  if(!conn || PQstatus(conn) != CONNECTION_OK){
    PQfinish(conn);
    cerr << "Protocol constructor: unable to connect to database" << endl;
    return;
  }
  
  /*
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "Protocol creator, can't get database connection " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  */

  ostringstream query;
  query << "select a.*, b.user_name from protocols a, users b where a.user_id=b.index and a.index=" << pid;
  PGresult* res = PQexec(conn, query.str().c_str());
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    PQclear(res);
    PQfinish(conn);
    return;
  }
  
  /*
  if(!conn.Exec(query.str().c_str())){
    cerr << "coulnd't query for the protocol with id : " << pid << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  if(conn.Tuples() != 1){
    cerr << "Didn't find any tuples.. bugger that,, " << endl;
    return;
  }
  */

  protocolId = atoi(PQgetvalue(res, 0, 0));
  parentId = atoi(PQgetvalue(res, 0, 1));
  creatorId = atoi(PQgetvalue(res, 0, 2));
  creationTime = QDateTime::fromString(PQgetvalue(res, 0, 3), Qt::ISODate);
  protocolName = PQgetvalue(res, 0, 5);
  protocolDescription = PQgetvalue(res, 0, 6);
  creatorName = PQgetvalue(res, 0, 7);

  /*
  protocolId = atoi(conn.GetValue(0, 0));
  parentId = atoi(conn.GetValue(0, 1));
  creatorId = atoi(conn.GetValue(0, 2));
  creationTime = QDateTime::fromString(conn.GetValue(0, 3), Qt::ISODate);
  protocolName = conn.GetValue(0, 5);
  protocolDescription = conn.GetValue(0, 6);
  creatorName = conn.GetValue(0, 7);
  */
  cout << "it seems that we got a protocol with name : " << protocolName << "  from the database : " << endl;
  if(fillsteps){
    cout << "calling fetchSteps: " << endl;
    fetchSteps();
  }
}

Protocol::Protocol(int pid, int ppid, string pname, string pdes, QDateTime cTime, int cId, string cName, const char* cinfo, bool fillSteps){             // if fillSteps, then get the stepdescriptions and everything from the database .. 
  conninfo = cinfo;
  protocolId = pid;
  protocolName = pname;
  protocolDescription = pdes;
  creationTime = cTime;
  creatorId = cId;
  creatorName = cName;
  parentId = ppid;
  if(cinfo && fillSteps){
    if(!fetchSteps()){
      cerr << "Couldn't fill the steps.. bloody hell, bastards the lot of them " << endl;
    }
  }
}

Protocol::Protocol(vector<ProtocolStep> psteps, int ppid, int cid, string cName, string pName, string pdes, const char* cinfo, bool dbenter){
  cout << "Protocol making protocol " << endl;
  conninfo = cinfo;
  protocolSteps = psteps;
  parentId = ppid;
  if(!parentId){
    parentId = 1;
  }  // ugly.. but .. better here than elsewhere.. 
  creatorId = cid;
  creatorName = cName;
  protocolName = pName;
  protocolDescription = pdes;
  /// at this point we don't know the process id or the creation time (creation == db entry..)
  protocolId = 0;
  if(dbenter){
    cout << "Trying to enter into database " << endl;
    if(!dbEnter(conninfo)){
      cerr << "Couldn't enter protocol into database bugger .. " << endl;
    }
  }
}

bool Protocol::dbEnter(const char* conninfo){
  cout << "Protocol dbEnter" << endl;
  if(!conninfo){
    cerr << "Protocol::dbEnter, no connection information" << endl;
    return(false);
  }
  PGconn* conn = PQconnectdb(conninfo);
  if(PQstatus(conn) != CONNECTION_OK){
    cerr << "Protocol::dbEnter connection to DB failed" << endl;
    return(false);
  }


  /*
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "Coulnd't make the connection to the database " << endl;
    return(false);
  }
  /// start a work thingy..
  if(!conn.Exec("begin")){
    cerr << "Coulnd't begin the execution block for entering the database thingy.. " << endl;
    return(false);
  }
  */
  /// start by entering the descriptions and things, obtain an id for the thingy, then go through the 
  /// steps and enter these..
  /// at end if everything ok, then commit the changes..

  // WARNING:: this code doesn't use transaction control. It's unlikely to ever be used, and should be
  // removed anyway.
  ostringstream penter;
  char* escPName = new char[protocolName.size()*2+1];
  char* escPDes = new char[protocolDescription.size()*2+1];
  PQescapeString(escPName, protocolName.c_str(), protocolName.size());
  PQescapeString(escPDes, protocolDescription.c_str(), protocolDescription.size());
  penter << "insert into protocols (parent, user_id, name, description) values (" << parentId << ", " 
       << creatorId << ", '" << escPName << "', '" << escPDes << "')";
  delete escPName;
  delete escPDes;

  PGresult* res = PQexec(conn, penter.str().c_str());
  int cmdTuples = atoi(PQcmdTuples(res));
  if(cmdTuples != 1){
    cerr << "CmdTuples != 1 but: " << cmdTuples;
    PQclear(res);
    PQfinish(conn);
    return(false);
  }
  
  /*
  cout << "trying to call db to exec : " << penter.str() << endl;
  if(!conn.Exec(penter.str().c_str())){
    cerr << "Couldn't insert protocol into protcol table" << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  // we need to set the protocol id from the table now, this is painful but necessary unfortunately.. 
  if(conn.CmdTuples() != 1){
    cerr << "CmdTuples is not 1 but something else,, hmmm " << conn.ErrorMessage() << endl;
    return(false);
  }
  */
  
  //  const char* oid = conn.OidStatus();     // which is bloody stupid, but I can't find a method that lets me get at this information.. any other way.. 
  Oid oid = PQoidValue(res);
  // then use this to select..
  ostringstream query;
  query << "select index, entry_time from protocols where oid = " << oid ;

  res = PQexec(conn, query.str().c_str());
  if(PQresultStatus(res) == PGRES_TUPLES_OK){
    protocolId = atoi( PQgetvalue(res, 0, 0) );
    creationTime = QDateTime::fromString( PQgetvalue(res, 0, 1), Qt::ISODate);
    for(int i=0; i < protocolSteps.size(); ++i){
      if(protocolSteps[i].isNew() || protocolSteps[i].id() == 0){
	if(!protocolSteps[i].dbEnter(conninfo)){
	  cerr << "Couldn't enter a protocol step into the database.. may already exist.. but will quit trying. " << endl
	       << "Memory leak at this location" << endl;
	  return(false);
	}
      }
      ostringstream stepenter;
      stepenter << "insert into protocol_step_series values (" << protocolId << ", " << protocolSteps[i].id() << ", " << i+1 << ")";
      res = PQexec(conn, stepenter.str().c_str());
      if(atoi(PQcmdTuples(res)) != 1){
	cerr << "Trouble entering protocol step into protocol_step_series: " << i << endl;
      }
    }
  }
  PQclear(res);
  PQfinish(conn);
  return(true);   /// Should do a commit step here.. later though.
  
  /*
  if(!conn.Exec(query.str().c_str())){
    cerr << "Couldn't exec the query for the protocol id.. bugger.. ,, hmm, trouble man " << conn.ErrorMessage() << endl;
    return(false);
  }
  if(conn.Tuples() != 1){
    cerr << "Conn tuples is not 1 bugger : " << conn.ErrorMessage() << endl;
    return(false);
  }
  protocolId = atoi(conn.GetValue(0, 0));
  creationTime = QDateTime::fromString(conn.GetValue(0, 1), Qt::ISODate);
  // then just go through the steps in the protocol steps, and if not new, just insert into protocol_step_series.. 
  for(int i=0; i < protocolSteps.size(); i++){
    cout << "trying to enter the step.. " << endl;
    if(protocolSteps[i].isNew() || protocolSteps[i].id() == 0){
      if(!protocolSteps[i].dbEnter(conn)){
	cerr << "Couldn't enter a protocol step into the database.. may already exist.. but will quit trying. " << endl
	     << conn.ErrorMessage() << endl;
	return(false);
      }
    }
  */
  /*
    ostringstream stepenter;
    stepenter << "insert into protocol_step_series values (" << protocolId << ", " << protocolSteps[i].id() << ", " << i+1 << ")";
    cout << "Trying to call step series and to enter : " << stepenter.str() << endl;
    if(!conn.Exec(stepenter.str().c_str())){
      cerr << "Couldn't manage to enter the step into the series.. lets check if maybe it doesn't exist" << endl;
      cerr << conn.ErrorMessage() << endl;
      if(!protocolSteps[i].dbEnter(conn)){
	cerr << "Can't enter either bugger " << endl
	     << conn.ErrorMessage() << endl;
	return(false);
      }
      // if here, just try again,, in case it worked..
      if(!conn.Exec(stepenter.str().c_str())){
	cerr << "Really cant enter the step in the series.. " << endl
	     << conn.ErrorMessage() << endl;
	return(false);
      }
    }
  }
  if(!conn.Exec("commit")){
    return(false);
  }
  // if we get here, then we are done and we can just return true..
  return(true);
  */
}

bool Protocol::fetchSteps(){
  // get the steps from the database.. and actually create the steps..
  cout << "fetch steps : " << endl;
  PGconn* conn = PQconnectdb(conninfo);
  if(!conn || PQstatus(conn) != CONNECTION_OK){
    cerr << "fetchSteps: unable to connect to database" << endl;
    return(false);
  }
  /*
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "Coulnd't make connection to db in order to fetch the steps for the protocol" << endl;
    return(false);
  }
  */
  ostringstream query;
  query << "select b.*, c.user_name from protocol_step_series a, protocol_step b, users c where a.step=b.step and a.protocol=" << protocolId << " and b.user_id=c.index order by step_number";
  cout << "fetch steps query is : " << query.str() << endl;

  PGresult* res = PQexec(conn, query.str().c_str());
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    cerr << "fetchSteps: no tuples returned" << endl;
    return(false);
  }
  int nTuples = PQntuples(res);
  /*
  if(!conn.Exec(query.str().c_str())){
    cerr << "Couldn't execute the query for the protocol step series... " << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  */
  //  for(int i=0; i < conn.Tuples(); i++){
  for(int i=0; i < nTuples; i++){
    int stepid = atoi(PQgetvalue(res, i, 0));
    int pid = atoi(PQgetvalue(res, i, 1));
    string description = PQgetvalue(res, i, 2);
    int cid = atoi(PQgetvalue(res, i, 3));
    QDateTime ctime = QDateTime::fromString(PQgetvalue(res, i, 4), Qt::ISODate);
    string cname = PQgetvalue(res, i, 5);
    protocolSteps.push_back(ProtocolStep(description, pid, cid, cname, ctime, stepid));
  }
  /*
  for(int i=0; i < conn.Tuples(); i++){
    int stepid = atoi(conn.GetValue(i, 0));
    int pid = atoi(conn.GetValue(i, 1));
    string description = conn.GetValue(i, 2);
    int cid = atoi(conn.GetValue(i, 3));
    QDateTime ctime = QDateTime::fromString(conn.GetValue(i, 4), Qt::ISODate);
    string cname = conn.GetValue(i, 5);
    protocolSteps.push_back(ProtocolStep(description, pid, cid, cname, ctime, stepid));
  }
  */
  cout << "protocolSteps size is now : " << protocolSteps.size() << endl;
  return(true);
}

bool Protocol::fillSteps(){                            
  for(int i=0; i < protocolSteps.size(); i++){
    if(!protocolSteps[i].isFilled()){
      protocolSteps[i].fillFromDB(conninfo);
    }
  }
  return(true);   // make some check later.. 
}

void Protocol::serialise(NetArray* na, bool withSteps){
  na->iapp(protocolId);
  na->iapp(creationTime.toTime_t());
  na->iapp(parentId);
  na->iapp(creatorId);
  na->sapp(creatorName);
  na->sapp(protocolName);
  na->sapp(protocolDescription);
  if(!withSteps){
    na->iapp(0);
    return;
  }
  na->iapp(protocolSteps.size());
  for(int i=0; i < protocolSteps.size(); i++){
    protocolSteps[i].serialise(na);
  }
}

ProtocolCollection::ProtocolCollection(const char* conninfo){
  PGconn* conn = PQconnectdb(conninfo);
  if(!conn || PQstatus(conn) != CONNECTION_OK){
    PQfinish(conn);
    cerr << "ProtocolCollection constructor: unable to connect to database" << endl;
    return;
  }
  /*
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    return;
  }
  */
  const char* query = "select a.*, b.user_name from protocols a, users b where a.user_id=b.index";
  PGresult* res = PQexec(conn, query);
  if(PQresultStatus(res) != PGRES_TUPLES_OK){
    cerr << "ProtocolCollection constructor. Query did not return any tuples" << endl;
    PQclear(res);
    PQfinish(conn);
    return;
  }
  int nTuples = PQntuples(res);
  /*
  if(!conn.Exec(query)){
    cerr << "ProtocolCollection : couldn't make query returning" << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  */
  for(int i=0; i < nTuples; i++){
    //  for(int i=0; i < conn.Tuples(); i++){
    int pid = atoi(PQgetvalue(res, i, 0));
    int ppid = atoi(PQgetvalue(res, i, 1));
    int user_id = atoi(PQgetvalue(res, i, 2));
    QDateTime ctime = QDateTime::fromString(PQgetvalue(res, i, 3), Qt::ISODate);
    string protName = PQgetvalue(res, i, 5);
    string protDescription = PQgetvalue(res, i, 6);
    string user_name = PQgetvalue(res, i, 7);
    Protocol temp(pid, ppid, protName, protDescription, ctime, user_id, user_name, conninfo, false);

    // int pid = atoi(conn.GetValue(i, 0));
    // int ppid = atoi(conn.GetValue(i, 1));
    // int user_id = atoi(conn.GetValue(i, 2));
    // QDateTime ctime = QDateTime::fromString(conn.GetValue(i, 3), Qt::ISODate);
    // string protName = conn.GetValue(i, 5);
    // string protDescription = conn.GetValue(i, 6);
    // string user_name = conn.GetValue(i, 7);
    // Protocol temp(pid, ppid, protName, protDescription, ctime, user_id, user_name, conninfo, false);
    protocols.push_back(temp);
  }
}

void ProtocolCollection::serialise(NetArray* na){
  na->iapp(protocols.size());
  for(int i=0; i < protocols.size(); i++){
    protocols[i].serialise(na, false);
  }
}
