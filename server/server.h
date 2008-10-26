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

#ifndef SERVER_H
#define SERVER_H

#include <qserversocket.h>
#include <qstring.h>
#include <qsocket.h>
#include <qsocketdevice.h>
#include <map>
#include <string>
#include <vector>
#include "connectionobject.h"
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"

//#define EXPERIMENTAL_SERVER

class Server : public QServerSocket
{
  Q_OBJECT

    public:
  Server(const char* db, int port, const char* dataTable="data", QObject* parent=0, const char* name=0 ); 
  ~Server();
  // hope that works !! .. 

 signals:
  void statusMessage(QString);

  private slots:
    void discardClient(QSocketDevice* s);
  void newConnection(int socket);

 private:
  map<QSocketDevice*, ConnectionObject*> connections;
  ProbeSetSet2 pSet;
  //vector<probe_data> data;
  //const char* conninfo;      // for connection to database.. 
  string conninfo;            // this is a kludge perhaps. 
  int commandCounter;        // total number of commands served.. update on connection object destruction.. 
  int clientCounter;         // number of clients that have connected.. increment on destruction
  protected :
    virtual void customEvent(QCustomEvent* e);
};

#endif
