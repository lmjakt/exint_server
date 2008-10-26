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

#include "server.h"
#include "connectionobject.h"
#include "../raw/probe_set.h"
#include "../raw/probeSetSet2.h"
#include <qserversocket.h>
#include <qsocketdevice.h>
#include <qstring.h>
#include <qapplication.h>            // for qApp ?? 
#include <iostream>
#include <qtextstream.h>
#include <stdlib.h>
#include <map>
#include <unistd.h>
#include <crypt.h>
#include <string>
#include <vector>
// for shared memory. 
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

using namespace std;                   // USE port 8090 for experimental.. 

//const ProbeSetSet2 pSet("expression");

Server::Server(const char* db, int port, const char* dataTable, QObject* parent, const char* name) :
  QServerSocket(port, 0, parent, name)
{
  // hope that works, I've seens something like this somewhere else, maybe in the tutorial 
  if( !ok() ){
    cerr << "Failed to bind to port 8108" << endl;
    exit(1);       
  }
  cout << "Seems that the socket got created after all, how strange!! " << endl;

  clientCounter = 0;
  commandCounter = 0;

  conninfo = "dbname=";
  conninfo.append(db);
 
  //#ifdef EXPERIMENTAL_SERVER 
  //conninfo = "dbname=expression_2";
  //cout << "Experimental server is defined : " << endl;
  //#endif
  cout << "conninfo is " << conninfo << endl;
  
  //  const char* cinfo = conninfo.c_str();
  pSet = ProbeSetSet2(conninfo.c_str(), dataTable);
  cout << "size of probe set is " << sizeof(pSet) << endl;
  cout << "Sizeo of a pointer to probe set is : " << sizeof(&pSet) << endl;

  // let's try to create a shared memory segment.. hmm
  //void* sharedSet;      // but make this a pointer to the pSet.. and see if we can access it..
 //  char c;
//   int shmid;
//   //int shsize = 27;
//   int shsize = sizeof(ProbeSetSet2*);
//   //char *shm, *s;
//   ProbeSetSet2* shm;
//   //  int shsize = sizeof(&pSet);
//   //  int shsize = sizeof(ProbeSetSet2);
//   int key = 56789;
//   if ((shmid = shmget(key, shsize, IPC_CREAT | 0666)) < 0) {
//     cerr << "couldn't create shmid .. bugger " << endl ;
//     exit(1);
//   }
//   cout << "smid is : " << shmid << endl;
//   // and try to attach it.. hmm 
//   //sharedSet = shmat(shmid, NULL, 0);
//   shm = (ProbeSetSet2*)shmat(shmid, NULL, 0);
  
//   if((int)shm == -1){
//     //  if((int)sharedSet == -1){
//     cerr << "couldn't attach memory, what's going on " << endl;
//     //cerr << "shared set is : " << (int)sharedSet << endl;
//     exit(1);
//   }
//   //sharedSet = (void*)&pSet;
//   //s = shm;       // just a pointer...
//   //for(c = 'a'; c <= 'z'; c++){
//   //  *s++ = c;
//   //}
//   //(*s) = NULL; 
//   shm = &pSet;
//   cout << "Address of pSet is   : " << &pSet << endl
//        << "Value of sharedSet : " << shm << endl
//        << "Address of the pSet is : " << &pSet << endl
//        << "Location of ShardSet : " << &shm << endl;


  cout << "pSet loaded" << endl;
  cout << "Socket created and data loaded" << endl;;
}

Server::~Server(){
  // close the current connection,, -but how do I know if there is one??
  // map<QSocket*, ConnectionObject*>::iterator it;
  //for(it = connections.begin(); it != connections.end(); it++){
  //  delete connections[it->first];      // this should cause the socket to be closed by the destructor function.
  //  delete it->first;                   // might take some time as I'm writing to each socket.. 
  //}
}

void Server::newConnection(int socket){
  cout << "New connection received" << endl;
  QSocketDevice* s = new QSocketDevice();
  s->setSocket(socket, QSocketDevice::Stream);
  connections[s] = new ConnectionObject(s, &pSet, (QWidget*)this);
  //connect(connections[s], SIGNAL(deleteMe(QSocket*)), this, SLOT(discardClient(QSocket*)) );
  //connect(connections[s], SIGNAL(statusMessage(QString)), this, SIGNAL(statusMessage(QString)) );
  //connections[s]->start();
  //cout << "new connection object created and connections set up" << endl;;
}

void Server::discardClient(QSocketDevice* socket){
  cout << "deleting socket from: " << socket->peerAddress().toString() << endl;
  cout << "         socket   no: " << socket->socket() << endl;
  cout << "memory   address    : " << socket << endl;
  clientCounter++;

  //delete connections[socket->socket()];
  if(connections.find(socket) != connections.end()){
    // wait one second for the thread to terminate..
    commandCounter += connections[socket]->commandCounter;
    if(connections[socket]->wait(1000)){
      cout << "deleting data from connectionData in the discardClient bit" << endl;
      delete connections[socket];
    }else{
      cout << "connection object not returning from thread ,, try to terminate it" << endl;
      connections[socket]->terminate();
      if(connections[socket]->wait(1000)){
	delete connections[socket];
      }else{
	cerr << "can't seem to get the connectionbject to stop, try something else" << endl;
	delete connections[socket];
      }
    }
  }
  // erase it..
  connections.erase(socket);   
  delete socket;    // do NOT delete in the connectionObject
  cout << "\tCONNECTION OBJECT DELETED" << endl
       << "\tTotal number of clients served : " << clientCounter << endl
       << "\tTotal number of commands served : " << commandCounter << endl;
}

void Server::customEvent(QCustomEvent* e){
  if(e->type() == 1001){
    cout << "Update User Event information received.. should do something" << endl;;
    map<QSocketDevice*, ConnectionObject*>::iterator it;
    for(it = connections.begin(); it != connections.end(); it++){
      qApp->lock();
      (*it).second->updateUsers = true;
      qApp->unlock();           // the aApp lock and unlock should actually be unnecessary as there is only one other place where this is called.. but better safe than sorry.. 
    }
    return;
  }
  if(e->type() == 9999){
    // we want to delete the socket.. 
    QSocketDevice* s = (QSocketDevice*)e->data();
    discardClient(s);
    return;
  }
  cout << "Unknown User Event received. Am deeply confused" << endl;
}
