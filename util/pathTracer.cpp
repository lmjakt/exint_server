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

#include "pathTracer.h"
#include <iostream>
#include <math.h>   // for calculating weighted centers using a gaussian (?) distribution.. 

using namespace std;    // just for iostream so I can do some reasonable checking.. or something like that..

PathTracer::PathTracer(float** values, uint* id, uint pn, uint dn, set<PathTracer*>* trcrs, QMutex* mut, float sm){       // we don't need to pass any arguments to QThread..
  // set all pointers to 0, and then check through..
  sigmaMultiplier = sm;
  pointNumber = pn;
  dimNumber = dn;
  startChain = 0;
  trace = 0;
  startPoint = 0;
  nextPoint = 0;
  centerPoint = 0;
  tracers = trcrs;
  traceSize = 0;
  mutex = mut;
  if(pointNumber < 3){
    cerr << "PathTracer constrcutor point number is too small : " << pointNumber << endl
	 << "Will return without doing anything" << endl;
    return;
  }

  npoint* point = new npoint(dimNumber, values[0], id[0]);
  startChain = new pointLink(point, 0, 0);       // set the next one later.. 
  pointLink* lastPoint = startChain;
  for(uint i=1; i < pointNumber; i++){
    point = new npoint(dimNumber, values[i], id[i]);
    pointLink* thisPoint = new pointLink(point, lastPoint, 0);
    lastPoint->next = thisPoint;
    lastPoint = thisPoint;
  }
  // and we should set up the weighted center.. so that we don't have to do any defines of it at a later stage.. 
  centerPoint = new npoint(dimNumber, new(float[dimNumber]), 0);
  realCenter = new npoint(dimNumber, new(float[dimNumber]), 0);
  addPoints = new float[dimNumber];   // which will be worked out each time.. 
  // which has now set up the whole chain... don't delete anything at the moment as that would tend to really screw things up... 
  // we have to be very careful with our pointers here, or we'll end up getting into trouble...
}

PathTracer::~PathTracer(){
  // haven't decided yet,, as we may want to tell someone about the thingies..
  // perhaps we can say if we have passed the pointer, then leave it, otherwise 
  // delete it, as noone can know it..
  // but have to wait until I've written accessor functions.. 
  cout << "deleting path tracer .. " << endl;
  while(startChain){
    pointLink* link = startChain;
    startChain = startChain->next;
    delete link->point;
    delete link;
  }
  while(trace){          // I should add something to make sure that it is not possible for something to be in both trace and start..
    cout << "deleting point with id : " << trace->point->id << "  x: " << trace->point->x << "\ty: " << trace->point->y << "\tdistance : " << trace->point->d << endl;
    pointLink* link = trace;
    trace = trace->prior;         // trace has the last thingy in it.. 
    delete link->point;
    delete link;
  }
  cout << "trace deleted" << endl;
  delete centerPoint;
  cout << "centerPoint deleted " << endl;
  delete realCenter;
  cout << "realCenter deleted " << endl;
  delete []addPoints;
  cout << "addPoints deleted " << endl;
  //delete nextPoint;
  cout << "nextPoint deleted" << endl;
  //delete startPoint;
  cout << "done deleting path tracer.. " << endl;
}

float PathTracer::distance(npoint* p1, npoint* p2){
  // simply define the euclidean distance..
  float distance = 0;
  for(uint i=0; i < dimNumber; i++){
    distance += ((p1->p[i] - p2->p[i]) * (p1->p[i] - p2->p[i]));
  }
  return(sqrt(distance));
}

void PathTracer::findFirstPoint(){
  // go through all against all and find the smallest distance..
  pointLink* pt1;
  pointLink* pt2;
  
  // traverse the chain in two waves..
  pt1 = startChain;
  pt2 = startChain->next;
  float smallestDistance = distance(startChain->point, startChain->next->point);
  cout << "At start of finding smallest distance distance : " << smallestDistance << "  for " << pt1->point->id << "  and   " << pt2->point->id << endl;
  startPoint = pt1;
  nextPoint = pt2;
  while(pt1->next){
    // do something.. 
    pt2 = pt1->next;    // first time round this doesn't actually do anything, but what the hell. 
    while(pt2->next){
      float dist = distance(pt1->point, pt2->point);
      //cout << " | " << dist << " | " << endl;
      if(dist < smallestDistance && dist > 0){       // repeats itself the first time,, but no really good way of doing this.
	//cout << "*";
	startPoint = pt1;
	nextPoint = pt2;
	nextPoint->point->x = dist;      // start with first point at origin, and the second on the y axis at an appropriate distance.. 
	startPoint->point->d = dist;
	smallestDistance = dist;
      }
      if(dist == 0){
	pointLink* deadLink = pt2;
	pt2 = pt2->next;
	cerr << "Removing a link with id " << deadLink->point->id << "  as it seems to be identical to : " << pt1->point->id << endl;
	removePoint(deadLink);
	delete deadLink->point;
	delete deadLink;
      }else{
	pt2 = pt2->next;
      }
    }
    pt1 = pt1->next;
  }
  cout << endl << "Smallest distance is now : " << smallestDistance << "  for " << startPoint->point->id << "   and   " << nextPoint->point->id << endl;
  // and at this point actually everything should be set, unless I'm very much mistaken..
  startPoint->point->x = 0;
  startPoint->point->y = 0;     // I shouldn't really have to set this, but seems strange.. 
  removePoint(startPoint);
  removePoint(nextPoint);       // they exist only as the links to the thingies.. 

  // ok, bit of a kludge, but define the real center point, in this function.. 
  for(uint i=0; i < dimNumber; i++){
    realCenter->p[i] = (startPoint->point->p[i] + nextPoint->point->p[i])/(float)2;
  }
  realCenter->x = nextPoint->point->x/(float)2;   // 
  realCenter->y = 0;                              // dependent on the above procedure, so be careful.. 

}

void PathTracer::weightedCenter(pointLink* pl){
  // the biggest issue at this point is what we can use for a sigma.. 
  // how about arbitrarily setting sigma to be equal to 2x the distance of the previous connection,
  // which has at least some reasonable grounds to it. This will mean the weighing, in terms of multiples
  // the preceding connection will be..
  // 1x => 0.6, 2x => 0.38, 3x => 0.2, 4x => 0.14    (approximate numbers eyeballed from a plot of 2.7^(-x/2)
  // this is a pretty localised thing, but if there's a good cluster we may pick it up.. we'll see...
  // ofcourse, one woould rather not have to set any parameters at all, but that is always going to be difficult..
  
  // first 0 the addPoints...
  for(uint i=0; i < dimNumber; i++){
    addPoints[i] = 0;
  }
  // and then work out what we'll use for sigma..
  // well first make sure that there is a prior for the thingy..
  if(!pl->prior){
    for(uint i=0; i < dimNumber; i++){
      centerPoint->p[i] = pl->point->p[i];
    }
    return;
  }
  float sigma = (sigmaMultiplier * distance(pl->point, pl->prior->point));   // which we should make changeable at some point..
  if(sigma <= 0){
    cerr << "sigma has invalid value... " << sigma << "setting to 1. " << endl;
    sigma = 1;
  }
  float count = 1;
  pointLink* pt = pl;       // as we are determining the distance to the current point !!
  while(pt->prior){
    float d = distance(pl->point, pt->prior->point);
    float f = exp(-d/sigma);      // so if distance is equal to the closest distance, the the value is exp(-0.5) which is about 0.6
    for(uint i=0; i < dimNumber; i++){
      addPoints[i] += f * (pt->prior->point->p[i] - pl->point->p[i]);
    }
    count += 1;
    pt = pt->prior;
  }
  
  // and divide by the count...
  for(uint i=0; i < dimNumber; i++){
    // and define center point..
    centerPoint->p[i] = pl->point->p[i] + addPoints[i]/count;
  }
  /// which is a weighted average of a a set of n-dimensional points. kind of neat, but I wonder if it will actually work.. ok la.. 
}

void PathTracer::findNextPoint(){
  // this is simply a case of going through the remaining points, and finding the one closest to the 
  // weighted Center.. (note that it is the responsibility of the run function to remove points which 
  // have been incorporated into the new chain... 
  float minDistance = distance(centerPoint, startChain->point);
  pointLink* pl = startChain;
  nextPoint = pl;
  while(pl->next){
    float d = distance(centerPoint, pl->point);
    if(d < minDistance){
      minDistance = d;
      nextPoint = pl;
    }
    pl = pl->next;
  }
  // and set the distance to the trace..
  trace->point->d = distance(nextPoint->point, trace->point);
  //cout << "set trace distance to : " << trace->point->d << endl;
  // remove the nextPoint from the startChain, carefully making sure to reallocate the 
  // startChain if the nextPoint is equal to the startChain.
  removePoint(nextPoint);
}

void PathTracer::removePoint(pointLink* pl){
  if(!pl->next && !pl->prior){      // we're done.
    startChain = 0;
    pl->next = 0;
    pl->prior = 0;
    return;
  }
  if(pl == startChain){    // then also pl->prior = 0 and so pl->next must be true.. at least if I didn't screw up somewhere along the lines.. 
    startChain = pl->next;
    pl->next->prior = 0;
    pl->next = 0;
    pl->prior = 0;
    return;
  }
  if(!pl->next){   // end of the line..
    pl->prior->next = 0;
    pl->next = 0;
    pl->prior = 0;
    return;
  }
  pl->next->prior = pl->prior;
  pl->prior->next = pl->next;
  pl->next = 0;
  pl->prior = 0;
}



void PathTracer::set_x_y(){
  // ok, though I'm pretty sure that there's a way of doing this directly invovlving just simple
  // algebra, and a little geometry, (i.e. no need for the sins and cosines which I'm afraid may be
  // a little slow, I can't seem to solve the appropriate equation for this. (I get there halfway, but
  // get scared off by the number of terms, can't see myself doing it without screwing up somewhere).
  //
  // so, do it the hard way. 
  // Pretend that the trace->point (i.e. current point) is at the origin (0, 0).
  // create a transformed point of the previous point (or the centerpoint) which is 
  // 1. translated by -Tx, -Ty          (where T is the trace positions). call this point P for previous..
  // 2. rotate around the origin (i.e. the trace point transformed) by -a 
  //  where a is defined as :  atan(-Px/-Py), where Px are the coordinates of the transformed point Px. 
  // 3. find the two possible points for the Next point (Nx, Ny). Nx will be the same for the two transformed points..
  // 4. rotate these two points arount the origin by a degrees. 
  // 5. translate these points by Tx, Ty. Determine by some criteria which is the better point and assign the position
  //    of the nextPoint...
  // its a long way around, but I think it should always work.. 

  // use a triangulation between the centerpoint and the prior point as these have been defined in 2d space. Centerpoint, though
  // it is useful doesn't get defined in 2d space (its no quite clear how to define it, and it might take a little more time to work out..
  cout << "trace->prior address : " << trace->prior << endl;
  cout << "trace prior point x " << trace->prior->point->x << "\ttrace point x " << trace->point->x << endl;

  double px_o = (double)(trace->prior->point->x - trace->point->x);
  double py_o = (double)(trace->prior->point->y - trace->point->y);   // original.. value.. 
  double p = (double)distance(trace->prior->point, nextPoint->point);
  double t = (double)distance(trace->point, nextPoint->point);

  double a = atan(py_o/px_o);

  cout << "px_o: " << px_o << "\tpy_o: " << py_o << "\tp: " << p << "\tt: " << t << "\ta: " << a << endl;
  
  double px = px_o * cos(-a) - py_o * sin(-a);
  double py = px_o * sin(-a) + py_o * cos(-a);       // this is the equation for rotation around a certain angle.. 
  cout << "rot px: " << px << "\tpy: " << py ;

  double nx = ((px*px) + (t*t) - (p*p)) / (2 * px);   // next positive x.. 
  double d = (t*t) - (nx*nx);
  cout << "\td: " << d << endl;
  if(d < 0){ d = 0;}        // this shouldn't be possible, but rounding errors mean it could happen... 
  double npy = sqrt(d);     // next positive y, next negative y.. 
  double nny = -sqrt(d);     // there are two options eh.. -- actually we need two nx's as well as these will change after rotation..
  double npx, nnx;           // the rotated values. need extra identifiers to make sure one doesn't change during the transformation.. 

  cout << "nx: " << nx << "\tnpy: " << npy << "\tnny: " << nny << endl;
  
  // rotate the two pairs.. of figures.. by a
  npx = nx * cos(a) - npy * sin(a);
  npy = nx * sin(a) + npy * cos(a);
  // and .. 
  nnx = nx * cos(a) - nny * sin(a);
  nny = nx * sin(a) + nny * cos(a);

  // and translate these in an appropriate manner..
  npx = npx + trace->point->x;
  npy = npy + trace->point->y;
  
  nnx = nnx + trace->point->x;
  nny = nny + trace->point->y;
  
  cout << "npx: " << npx << "\tnpy: " << npy << "\tnnx: " << nnx << "\tnny: " << nny << endl;

  // calculate the global error for the two potential points..
  pointLink* link = trace;
  float pError = 0;
  float nError = 0;   // the error for the positive and negative points respectively.. 
  while(link){
    float dst = distance(nextPoint->point, link->point);
    float pdst = sqrt((npx - link->point->x) * (npx - link->point->x) + (npy - link->point->y)*(npy - link->point->y));
    float ndst = sqrt((nnx - link->point->x) * (nnx - link->point->x) + (nny - link->point->y)*(nny - link->point->y));
    // that's ugly, but should at least be reasonably fast.. 
    pError += fabs(log(dst/pdst));
    nError += fabs(log(dst/ndst));
    link = link->prior;
  }
  //  cout << "pError: " << pError << "\tnError: " << nError << endl;
  if(pError < nError){
    nextPoint->point->x = npx;
    nextPoint->point->y = npy;
    return;
  }
  nextPoint->point->x = nnx;
  nextPoint->point->y = nny;
  // make some reporting as to the different values.. at some point.. 
}
    
void PathTracer::update_center(pointLink* pl){
  /// relies on the correct value of traceSize.. so be careful..
  float r1 = ((float)traceSize-1)/((float)traceSize);
  float r2 = 1.0 / float(traceSize);
  for(uint i=0; i < dimNumber; i++){
    realCenter->p[i] = realCenter->p[i] * r1 + pl->point->p[i] * r2;
  }
  realCenter->x = realCenter->x * r1 + pl->point->x * r2;
  realCenter->y = realCenter->y * r1 + pl->point->y * r2;
}
  
void PathTracer::run(){
  // pretty simple really..
  cout << "PathTracer run funtion " << endl;
  findFirstPoint();
  // which actually defines both the first point and the nextpoint.. 
  // then actually we need to incorporate the next point in,, as we won't need to do anything directly with it..
  cout << "Got first point " << endl;

  // startpoint is already set, but we should set..
  startPoint->next = nextPoint;
  trace = nextPoint;       // so we know where we are..
  trace->prior = startPoint;
  traceSize = 2;            // 
  //cout << "got a few more points " << endl;
  // so what is the trace for.. we just need startPoint,, as we can go through everything through it..
  while(startChain){
    cout << ".";
    weightedCenter(trace);
    findNextPoint();      // this sets nextPoint to something, so I need to remember the last point.. call it the trace.. 
    set_x_y();            // whoaaah.. dodgy perhaps.. 
    nextPoint->prior = trace;
    trace->next = nextPoint;
    trace = nextPoint;
    traceSize++;                // hmmm, 
    update_center(nextPoint);   // which is incorporated into the chain, and for which traceSize has already been updated.. 
  }
  cout << endl;
  // done... ??? but do we really need the trace ?? 
  mutex->lock();
  tracers->insert(this);
  mutex->unlock();
  ///cout << "inserted into thingy, and will do some more stuff later on.. " << endl;
}

// that's awfully simple, only 190 lines.. plus 86 of header info stuff.. and looks to me like it should be pretty fast.. 
