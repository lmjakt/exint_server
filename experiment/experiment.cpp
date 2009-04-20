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

#include "experiment.h"
#include <sstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

Experiment::Experiment(int uid, uint time, int prot, string com, const char* conninfo){
  id = -1;    // if sucessful, it will get a new value from the database..
  userId = uid;
  protocol = prot;   // this is the id,, -unless this is correct we'll not be able to enter into the database..
  comment = com;

  experimentTime.setTime_t(time);
  // now try to enter into the database, followed by an attempt to get fill the whole structure out..
  // which can then be serialised to the pope..
  if(!dbEnter(conninfo)){
    cerr << "Couldn't enter into database or some other problem.. " << endl;
  }else{
    cout << "Experiment Constructor : Managed to enter into datbase " << endl;
  }
}

bool Experiment::dbEnter(const char* conninfo){
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "Experiment dbEnter data base connection bad " << endl;
    return(false);
  }
  // ok let's escape the comment.
  char* cleanComment = new char[comment.size() * 2 +1];
  PQescapeString(cleanComment, comment.c_str(), comment.size());
  
  // start a transaction..
  if(!conn.Exec("begin")){
    cerr << "Experiment dbEnter couldn't start a transaction block, returning" << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  ostringstream enter;
  enter << "insert into ish_experiments (user_id, experiment_time, protocol, free_comment) values ("
	<< userId << ", '" << experimentTime.toString(Qt::ISODate) << "', " << protocol
	<< ", '" << cleanComment << "')";
  delete cleanComment;
  cout << enter.str() << endl;
  if(!conn.Exec(enter.str().c_str())){
    cerr << "Experiment dbEnter, couldn't call the entry thing : " << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  if(conn.CmdTuples() != 1){
    cerr << "Experiment dbEnter CmdTuples is not 1 but : " << conn.CmdTuples() << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }
  const char* oid = conn.OidStatus();     // hmm...
  /// and now we need a much stronger query across several tables..
  ostringstream query;
  query << "select a.experiment, a.user_id, a.experiment_time, a.entry_time, a.protocol, a.free_comment, "
	<< "b.name, b.description, c.user_name from ish_experiments a, protocols b, users c "
	<< "where a.oid=" << oid << " and a.user_id=c.index and a.protocol=b.index and a.user_id=b.user_id";
  cout << query.str() << endl;
  if(!conn.Exec(query.str().c_str())){
    cerr << "Experiment couldn't confirm the entry.. " << endl
	 << conn.ErrorMessage();
    return(false);
  }
  if(conn.Tuples() != 1){
    cerr << "conn.Tuples is not 1 but : " << conn.Tuples() << endl;
    return(false);
  }
  // and then I just have to make a thing for making a vector of thingies.. maybe a collection or something..
  id = atoi(conn.GetValue(0, 0));
  userId = atoi(conn.GetValue(0, 1));
  experimentTime = QDateTime::fromString(conn.GetValue(0, 2), Qt::ISODate);
  entryTime = QDateTime::fromString(conn.GetValue(0, 3), Qt::ISODate);
  protocol = atoi(conn.GetValue(0, 4));
  comment = conn.GetValue(0, 5);
  protocolName = conn.GetValue(0, 6);
  protocolDescription = conn.GetValue(0, 7);
  userName = conn.GetValue(0, 8);
  // and at this point we commit to the database..
  if(!conn.Exec("commit")){
    id = -1;
    cerr << "Experiment dbEnter : Couldn't commit to database " << endl
	 << conn.ErrorMessage() << endl;
    return(false);
  }  // and now finally we can just set everything to the correct value.. which should be pretty simple..
  return(true);
}

void Experiment::serialise(NetArray* na){
  na->iapp(id);
  na->iapp(userId);
  na->sapp(userName);
  na->iapp(experimentTime.toTime_t());
  na->iapp(entryTime.toTime_t());
  na->iapp(protocol);
  na->sapp(protocolName);
  na->sapp(protocolDescription);
  na->sapp(comment);
  // and that's pretty much it..
}

  
ExperimentCollection::ExperimentCollection(const char* conninfo){
  PgDatabase conn(conninfo);
  if(conn.ConnectionBad()){
    cerr << "Experiment dbEnter data base connection bad " << endl;
    return;
  }
  ostringstream query;
  query << "select a.experiment, a.user_id, a.experiment_time, a.entry_time, a.protocol, a.free_comment, "
	<< "b.name, b.description, c.user_name from ish_experiments a, protocols b, users c "
	<< "where a.user_id=c.index and a.protocol=b.index and a.user_id=b.user_id";
  if(!conn.Exec(query.str().c_str())){
    cerr << "ExperimentCollection couldn't make the query " << endl
	 << conn.ErrorMessage() << endl;
    return;
  }
  for(int i=0; i < conn.Tuples(); i++){
    int id = atoi(conn.GetValue(i, 0));
    int userId = atoi(conn.GetValue(i, 1));
    QDateTime experimentTime = QDateTime::fromString(conn.GetValue(i, 2), Qt::ISODate);
    QDateTime entryTime = QDateTime::fromString(conn.GetValue(i, 3), Qt::ISODate);
    int protocol = atoi(conn.GetValue(i, 4));
    string comment = conn.GetValue(i, 5);
    string protocolName = conn.GetValue(i, 6);
    string protocolDescription = conn.GetValue(i, 7);
    string userName = conn.GetValue(i, 8);
    experiments.push_back(Experiment(id, userId, userName, experimentTime, entryTime, protocol, protocolName, protocolDescription, comment));
  }
}

void ExperimentCollection::serialise(NetArray* na){
  na->iapp(experiments.size());
  for(int i=0; i < experiments.size(); i++){
    experiments[i].serialise(na);
  }
}
