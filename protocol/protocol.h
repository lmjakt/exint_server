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

#ifndef PROTOCOL_H
#define PROTOCOL_H

#include <string>
#include <vector>
#include <qdatetime.h>
#include "../netArray/netArray.h"
//#include <libpq++.h>

using namespace std;

class ProtocolStep {
  public :
    ProtocolStep();
  ProtocolStep(int sid, const char* conninfo=0);        // if conninfo is 0, just contains the id.. otherwise fills itself.. 
  ProtocolStep(string des, int pid, int cid, string cname, QDateTime ctime, int sid, const char* conninfo=0);
  ProtocolStep(int sid, int pid, int cid, string des, const char* conninfo);
  ~ProtocolStep(){}

  bool dbEnter(const char* conninfo);          // enter into the database,, return the id.. (also set the id given).. 
  //  bool dbEnter(PgDatabase& conn);          // enter into the database,, return the id.. (also set the id given).. 
  void serialise(NetArray* na);    
  bool fillFromDB(const char* conninfo);      // get the data from the database..

  int id(){
    return(stepId);
  }
  string description(){
    return(stepDescription);
  }
  bool isFilled(){
    return(filled);
  }
  bool isNew(){
    return(newStep);
  }
  private :
    int stepId;              // 0 if not described.. 
  string stepDescription;
  int parentId;
  int creatorId;              // who originiated the step.. (integer.. should be ok for most cases..)
  string creatorName;
  QDateTime creationTime;     // the time when the step was created.. 
  bool filled;             // true -is filled.. 
  bool newStep;
  // and some database methods...
};

class Protocol {
  public :
    Protocol();
  Protocol(int pid, int ppid, string pname, string pdes, QDateTime cTime, int cId, string cName, const char* cinfo,  bool fillSteps=false);
  Protocol(vector<ProtocolStep> psteps, int ppid, int cid, string cName, string pName, string pdes, const char* cinfo, bool dbenter);
  Protocol(const char* cinfo, int pid, bool fillsteps);
  ~Protocol(){}
  
  bool fillSteps();  // get or fill empty steps.. -- probably not very useful.. 
  bool fetchSteps();  // get and fill steps... 
  void serialise(NetArray* na, bool withSteps);  // serialise the data to a net array.. 
  bool dbEnter(const char* conninfo);  // stick it into the database if it works, then everyone happy..  

  vector<ProtocolStep> steps(){
    return(protocolSteps);
  }
  int id(){
    return(protocolId);
  }
  int pId(){
    return(parentId);
  }
  int cId(){
    return(creatorId);
  }
  
  private :
    vector<ProtocolStep> protocolSteps;
  int protocolId;
  string protocolName;
  string protocolDescription;
  QDateTime creationTime;
  int creatorId;
  string creatorName;
  const char* conninfo;

  // relationships..
  int parentId;          // actually we don't really need to build networks.. as anyway, we have to 
                         // serialise everything and send to the client.. not much point.. 
};

class ProtocolCollection {
  public :
    ProtocolCollection(const char* conninfo);       // just gets a series of protocols from the database,, and can write to thingy..
  ~ProtocolCollection(){}
  
  void serialise(NetArray* na);                     // serialise,, but just the id's.. right.. 
  
  private :
    vector<Protocol> protocols;
};

#endif 
  
    
