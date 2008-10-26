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

#include "blastClient.h"
#include <iostream>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <string.h>   // for memset .. ?
#include <unistd.h>   // for close ?? i.e. file descriptors calls I suppose
#include <errno.h>

using namespace std;

BlastClient::BlastClient(QString seq){
  sequence = seq;
  // create the socket and everything in the run function..
  cout << "Blast Client Constructed : " << endl;
}

BlastClient::~BlastClient(){
  cout << "Destroying Blast Client.. " << endl;
  for(int i=0; i < blastOut.size(); i++){
    delete []blastOut[i];
  }
  cout << "\t\tBlast Client data deleted " << endl;
}

void BlastClient::run(){                     // create the socket, send the sequence and wait for the return..
  int sock;
  sockaddr_in blastServer;
  unsigned short serverPort = 12345;
  const char* ip = "127.0.0.1";         // i.e. local host..
  int bufferSize = 100;           // a fairly small buffer size .. 

  if((sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0){
    cerr << "Blast Client failed to create a socket .. " << endl;
    return;   // nothing needs to be done really..
  }
  //  cout << "Blast Client created the socket .. " << endl;
  // the server address structure..
  memset(&blastServer, 0, sizeof(blastServer));
  blastServer.sin_family = AF_INET;
  blastServer.sin_addr.s_addr = inet_addr(ip);
  blastServer.sin_port = htons(serverPort);
  //cout << "Blast client set the socket address : " << endl;
  // see if we can establish a connection..
  if(connect(sock, (sockaddr*)&blastServer, sizeof(blastServer)) < 0){
    cerr << "Blast Client failed to connect to server.. " << endl;
    return;    // set some error I suppose, but..
  }
  //cout << "Blast client connected to server " << endl;
  // then see if we can send something to the socket.. -- just the sequence in latin1..
  sequence.append("|..|\n");           // this is the end character thingy.. maybe not brilliant, but should work.. 
  size_t nleft = sequence.length();
  ssize_t nwritten;
  const char* ptr = sequence.latin1();
  //cout << "begninning to write to the blast server .. " << endl;
  while(nleft > 0){
    if( (nwritten = write(sock, ptr, nleft)) <= 0) {
      if( errno == EINTR ){
	nwritten = 0;
	cout << "\t\twriten received an interrupt, socket number is " << sock << endl;
      }else{
	cout << "\t\tblast Client received an error number : " << errno << endl;
	return;
      }
    }
    nleft -= nwritten;
    ptr += nwritten;
  }
  //cout << "Finished writing to the blast server am about to start listening.. " << endl;
  // at this point we should have sent the whole thingy, now, let's just listen to the socket, until it closes
  // at which point we will do something else..
  char* data = new char[bufferSize + 1];
  int charsRead;
  bool keepOnGoing = true;
  while(keepOnGoing){
    while((charsRead = read(sock, data, bufferSize)) > 0){
      //cout << "Read " << charsRead << "  of output from the server .. " << endl;
      data[charsRead] = 0;            // allows me to pass it as a string.. 
      //  cout << data;           // 
      blastOut.push_back(data);
      lineLengths.push_back(charsRead);
      data = new char[bufferSize+1];
      
    }
    // here charsRead is less than 0, if equel to EINTR, then keep on going, otherwise, just end the reading.. and close
    if(charsRead != EINTR){
      keepOnGoing = false;
      cout << "Socket either closed, or we received an error of some sort " << endl;
    }
  }
  cout << endl;
  // ok, we are here, let's close the socket,, -- and get out of here.
  close(sock);
}
    
  
