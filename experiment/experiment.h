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

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <string>
#include <qdatetime.h>
#include "../netArray/netArray.h"
//#include <libpq++.h>
#include <vector>

using namespace std;

class Experiment {
  public :
    Experiment(){ id = -1; }
  Experiment(int uid, uint time, int prot, string com, const char* conninfo);   // fill in details from the database.. 
  Experiment(int i, int uid, string uname, QDateTime expttime, QDateTime entime, int prot, string pName, string pDes, string com){
    id = i;
    userId = uid;
    userName = uname;
    experimentTime = expttime;
    entryTime = entime;
    protocol = prot;
    comment = com;
    protocolName = pName;
    protocolDescription = pDes;
  }
  void serialise(NetArray* na);     // for convenience.. 
  bool isEntered(){
    return(id != -1);
  }
  
  private :
    int id;   // -1 = not set yet. 
  int userId;
  string userName;
  QDateTime experimentTime;
  QDateTime entryTime;      // the time entered into the database.. 
  int protocol;             // to view the protocol, use a protocol viewer,, -it just needs to know the db index of the protocol..
  string protocolName;
  string protocolDescription;   // get these when loading experiments.. 
  string comment;           // extra information that the user might want to input into the system..
  
  bool dbEnter(const char* conninfo);   // for ease.. just return if some problem
};

class ExperimentCollection {
  public :
    ExperimentCollection(const char* conninfo);      // gets everything from the database.. 
  
  void serialise(NetArray* na);
  private :
    vector<Experiment> experiments;
};


#endif

