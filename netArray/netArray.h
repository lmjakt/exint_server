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

#ifndef NETARRAY_H
#define NETARRAY_H

#include <string>
using namespace std;

class NetArray {
  public :
    NetArray();      // no indication of size
  NetArray(int inSize, int i=0, const char* nm=0);  // initial size.. 
  ~NetArray();
  void empty();          // just 0's the counter, but it keeps the memory,,
  void app(char c);     // append a char..
  void sapp(string s);
  void iapp(int v);
  void fapp(float f);
  void dapp(double d);

  char* readBinaryFromSocket(int socketNumber, int length);   // returns a null value on failure (and deletes the data).. 
                                                              // on sucess the user has to take care of deleting the data.. 
  int appendFromSocket(int socketNumber, int length);
  char* readChar();
  void resetReadPointer();
  int size();            // returns the size..
  int capacity();        // the capacity.. 
  char* data;            // dangerous but make it open.. 
  bool hasMoreData();        // if I can read one more.. char
  int id;
  const char* name;         // some identifiers.. 


  private :
    int memLength;            // amount of memory allocated..
  int curSize;                // the currentSize;
  int readPos;
  void grow();                // double our current size.. 
};

#endif

