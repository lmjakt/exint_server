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

#ifndef PATHTRACER_H
#define PATHTRACER_H

#include <qthread.h>
#include <qmutex.h>
#include <set>

// some structs and an object inheriting qthread to implement a method
// for tracing a path through a set of points in N-dimensions
// The method may be referred to as something like a nearest neighbourhood 
// method, or perhaps a weighted nearest neighbour.
// Basically follows the following steps ..
// 1. find the two closest points, define one as the start, and one as the next one..
// 2. find a mean position for the points included in the path, weighted on the last included point..
// 3. find the nearest point not already included in the path to above weighted mean..
// 4. repeat from 2 until all points are included in the path..

// optionally, include some function for assigning 2d positions for the path based on the distances of the point being
// added to the weighted mean and the last point added. --> unfortunately, this will either have 2 good solutions, or
// no good solutions, --> so we will have some decisions to make for that.

// ok, first a struct for a point. 

struct npoint {
  float* p;     // p for position.. 
  uint n;       // number of dimensions and the size of p;
  uint id;      // some identifier.. 
  float x;
  float y;     // in case we calculate these positions.. 
  float d;     // distance to the next point.. 
  npoint(uint dims, float* pos, uint i){
    p = pos;
    n = dims;
    id = i;
    x = y = d = 0;     // default, but can be set by external user.. 
  }
  ~npoint(){
    delete []p;
  }
};

struct pointLink {    // a link in a doubly linked list..
  npoint* point;
  pointLink* next;
  pointLink* prior;
  pointLink(npoint* p, pointLink* pr, pointLink* nxt){
    point = p;
    prior = pr;
    next = nxt;
  }
  ~pointLink(){
    //    delete point;  // can't do this as there may be multiple things pointing to this.. 
  }
};

// so I know I should make things private, and have accessor functions. but am afraid of any performance penalties for 
// these operations, so am trying to keep it simple..

class PathTracer : public QThread
{
  public :
    PathTracer(float** values, uint* id, uint pn, uint dn, std::set<PathTracer*>* trcrs, QMutex* mut, float sm);     // just get the user to use the bool finished() function to determine stuff.. 
                                                                    // id == point ids.. usually same as experiment ids.. bummer.. 
  ~PathTracer();

  // we will also need some sort of function to return the data to the caller, but work that out later.
  pointLink* chain(){
    return(startPoint);
  }

  protected :
    void run();
  
  private :
    float distance(npoint* p1, npoint* p2);
  void findFirstPoint();
  void findNextPoint();
  void weightedCenter(pointLink* pl);   // rather than returning anything these functions set state variables in the object (as we need to define these)..
  void removePoint(pointLink* pl);      // remove from startChain,, make sure no problem with it.. 
  void set_x_y();                       // set the x and y position of next point, depending on the distances to central point and trace.. 
  void update_center(pointLink* pl);    // relies on knowing the size of the trace, which must be the size of the chain after incorporating the link.

  // and some of these variables that I need to know..
  int pointNumber;          // the starting point Number.. (perhaps not strictly necessary, but may come in handy
  int dimNumber;                // the number of dimensions in the space.. 
  float sigmaMultiplier;     // a value for choosing the sigma mutliplier.. 
  pointLink* startChain;    // The starting chain defined by the constructor
  pointLink* trace;         // the final trace. start with none, and then input stuff into this.. as we go along..
  int traceSize;            // keep track of the trace size so that we can 
  
  pointLink* startPoint;       // set initially as one of the members of the smallest distance..
  pointLink* nextPoint;        // set after finding each closest point each time..
  npoint* centerPoint;         // the center set by weightedCenter.. and used to create the x-y coordinates.. 
  npoint* realCenter;          // the mean center, not the weighted mean... -use this for determining stuff.. 
                               // needs to have an n-dimensional point, and also a 2d point.. -- update every time something is added to the trace.. 
  float* addPoints;            // used to calculate the centerPoint. Just so we don't have to do anything to it.. 

  // stuff for communicationg with my owner, or whoever might want to care..
  std::set<PathTracer*>* tracers;
  QMutex* mutex;

};


#endif
