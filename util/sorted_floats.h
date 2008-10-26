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

#ifndef SORTED_FLOATS_H
#define SORTED_FLOATS_H

// a small struct that maintains a set of floats in order.
// supports inserting a new number, traversing the numbers and deleting the set.
// -- use to get medians..

/// THIS IS A DANGEROUS FUNCTION WHICH CAN CAUSE MEMORY SEGMENTATION AS IT MAKES NO CHECKS
/// WHATSOEVER IN ORDER TO GET FASTER SPEED

#include <string.h>

//#include <iostream>
//using namespace std;

struct fpoint {
  float value;
  fpoint* next;
  fpoint* previous;
  fpoint(float v){
    value = v;
    next = 0;
    previous = 0;
  }
  fpoint(){
    value = 0;
    next = previous = 0;
  }
};

struct sorted_floats {
  sorted_floats(int cap){      // cap is the capacity
    capacity = cap;
    number = 0;
    values = new fpoint[cap];
    startPoint = values;
    endPoint = values;
  }
  sorted_floats(int cap, float iv){         // the initial value
    capacity = cap;
    number = 1;
    values = new fpoint[cap];
    values[0].value = iv;   // ofcourse there is no next... nor any previous
    startPoint = values;
    endPoint = values;
  }
  ~sorted_floats(){
    delete []values;
  }

  void erase(){
    number = 0;
    memset((void*)values, 0, capacity * sizeof(fpoint));
  }

  void insert(float v){
    values[number].value = v;    // hmm to check or not to check...
    if(!number){
      values[0].value = v;
      number++;
      return;
    }
    
    // and then work out if we need to change any pointers..
    fpoint* it = startPoint;
    while(it){
      if(it->value > v){
	values[number].next = it;
	values[number].previous = it->previous;
	if(it->previous){
	  it->previous->next = values + number;
	}else{
	  startPoint = values + number;
	}
	it->previous = values + number;
	break;
      }
      it = it->next;
    }
    // if values[number].next is 0, then our value must be larger than everything..
    if(!values[number].next){
      values[number].previous = endPoint;
      endPoint->next = values + number;
      endPoint = values + number;
    }
    number++;
  }

	
  fpoint* begin(){
    return(startPoint);
  }

  float median(){
    if(!number){
      return(0);
    }
    if(number == 1){
      return(values[0].value);
    }
    int midpoint = number / 2;
    fpoint* it = startPoint;
    int c = 0;
    float med = 0;
    while(it){
      if(c >= midpoint){    // if only one value, well we won't be here anyway..
	if(number % 2){   // odd number , return this value..
	  return(it->value);
	}else{
	  return((it->value + it->previous->value)/2);
	}
      }
      it = it->next;
      c++;
    }
    return(med);   // shouldn't happen, but..
  }
	
  private :
  int capacity;
  int number;   // the current number..
  fpoint* values;
  fpoint* startPoint;
  fpoint* endPoint;
};
  
#endif
