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

#include <iostream>
//#include <sys/socket.h>
//#include <netinet/in.h>
//#include <sys/types.h>
//#include <arpa/inet.h>
//#include <netdb.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/select.h>
#include <netinet/in.h>
#include "netArray.h"


using namespace std;

void NetArray::empty(){
  curSize = 0;
  readPos = -1;                    // so I can increment before returning,, (then I can return in one statement)... 
}

NetArray::NetArray(){
  // assume a 1000 bytes is long enough.. 
  memLength = 1000;
  data = new char[memLength];
  curSize = 0;    // no values allocated.
  readPos = -1;
}

NetArray::NetArray(int inSize, int i, const char* nm){
  memLength = inSize;
  data = new char[memLength];
  curSize = 0;
  readPos = -1;
  id = i;
  name = nm;
}

NetArray::~NetArray(){
  delete []data;
}

void NetArray::app(char c){
  if(curSize >= memLength){
    grow();
  }
  data[curSize] = c;
  curSize++;
}

void NetArray::sapp(string s){
  iapp(s.size());
  for(int i=0; i < s.size(); i++){
    app(s[i]);
  }
}

void NetArray::iapp(int v){
  int q = htonl(v);
  char* c = (char*)&q;
  for(int i=0; i < 4; i++){
    app(*(c+i));
  }
}

void NetArray::fapp(float f){
  char* c = (char*)&f;
  for(int i=0; i < 4; i++){
    app(*(c+i));
  }
}

void NetArray::dapp(double d){
  char* c = (char*)&d;
  for(int i=0; i < 8; i++){            // is a double really 8 bytes ?? should be using sizeof.. ahh, well we can change later.. 
    app(*(c+i));
  }
}

int NetArray::appendFromSocket(int socketNumber, int length){
  //  cout << "beginning of appendFromSocket" << endl;
  while(curSize + length > memLength){
    grow();
  }
  //cout << "errNo is now : " << errno << endl;
  //  alarm(2);

  ///////////////// DON'T use alarm, this seems to cause trouble with multiple threads requiring a static function
  ///////////////// that is not bound to an object. Try to use select instead.. 
  
  // loop on a select.. 
  //while(1){
  
//   // set up the select parameters..
  fd_set readSet;
  FD_ZERO(&readSet);
  FD_SET(socketNumber, &readSet);
  timeval timeout;
  timeout.tv_sec = 2;
  timeout.tv_usec = 0;
  int selreturn;
  if((selreturn = select(socketNumber+1, &readSet, NULL, NULL, &timeout)) > 0){
    int charsRead = read(socketNumber, data+curSize, length);
    //alarm(0);
    if(charsRead < 0){   // some error..
      return(-errno);    // bit dodgy, but hmm, maybe ok..  
    }
    //cout << "NetArray : charsRead is " << charsRead << endl;
    //cout << "and errno is now : " << errno << endl;
    curSize += charsRead;
    return(charsRead);
  }
  // if we get here.. 
  // we either got interrupted, some error occured or we timed out..
  // --- This is bad but,, if we timed out,, let's return -EINTR as this
  // will cause the connectionobject to check flags and eventually return here... 
  if(selreturn == 0){
    return(-EINTR);
  }
  return(-errno);
}

char* NetArray::readBinaryFromSocket(int socketNumber, int length){
  // read the number of bytes specified by the argument into a char array,,,
  // first take the immediate bytes in the data array and copy up to length 
  // bytes.. -if need more than read directly into the new char*,, -then empty the data array (for completeness)..
  char* rdata = new char[length];
  int i;
  for(i=0; i < length; i++){
    readPos++;
    if(readPos < curSize){
      rdata[i] = *(data + readPos);
    }else{
      break;
    }
  }
  if(i == length){     // 
    return(rdata);     // and we don't have to do anything more.. -the user has to take care of deleting the data. 
  }
  // if I am here, the readPos is now larger than the currentPosition.. which means I have run out of data.. 
  // lets empty the data.. and then use a socket function to read directly into the rdata..
  empty();
  /// I now need to read length - i char's into the array.. 
  int remaining = length - i;
  while(remaining){
    int charsRead = read(socketNumber, rdata + i, remaining);
    // if error do again..
    if(charsRead < 0){
      continue;     // just go to the top again..
    }
    if(charsRead == 0){    // socket is broken.. 
      delete rdata;      
      return(0);
    }
    remaining -= charsRead;
    i += charsRead;
  }
  // remaining should now be 0.. and i should be equal to length.. let's print somethinng out and see if that's correct.
  cout << "NetArray readBinaryFromSocket remaining is " << remaining << "  length : " << length << "  and i is : " << i << endl;
  return(rdata);
}


char* NetArray::readChar(){
  readPos++;
  if(readPos < curSize){
    return(data + readPos);
  }
  return(0);
}

bool NetArray::hasMoreData(){
  return(curSize > (readPos + 1));
}

void NetArray::resetReadPointer(){
  readPos = -1;
}

void NetArray::grow(){
  memLength = memLength*2;
  cout << "name : " << name << " id: " << id << "  NetArray trying to grow to size :  " << memLength << endl;
  char* temp = new char[memLength];
  cout << "name : " << name << " id: " << id << "  NetArray grew to : " << memLength << endl;
  //// now the values from data have to be copied over..
  for(int i=0; i < curSize; i++){
    temp[i] = data[i];
  }
  // delete data..
  delete []data;
  data = temp;
}

int NetArray::size(){
  return(curSize);
}

int NetArray::capacity(){
  return(memLength - curSize);
}
