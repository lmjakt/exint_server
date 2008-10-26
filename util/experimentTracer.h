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

#ifndef EXPERIMENTTRACER_H
#define EXPERIMENTTRACER_H

#include "pathTracer.h"
#include "../raw/probe_set.h"
#include <vector>       // for the constructor
#include <set>
#include <qmutex.h>

// a class that takes a vector of probe_set points, and a vector of experiments, -- then using a global normalisation
// (well that can be changed at later stages, creates a set of mean expression patterns for those probe sets and those
// experiments. That matrix is then fed into a pathTracer which attempts to plot a path through the n-dimensional points
// in some sort of reasonable manner. The resulting path, along with its distances can then be used for various types of
// clustering approaches. 
//
// PathTracer runs in a separate thread once set up (though it in itself uses only one thread, so we will have to decide on
// some way of informing the caller when it is done, so that the resulting data can be sent back to the user. Path-tracer, might
// be possible to use on genes as well, though, I'm not sure that it will scale well to thousands of genes, as it uses a rather 
// tiresome approach to finding the stuff (sort of around O^3 or O^2, not sure at the moment).

// basically get the data, hand it off to a pathTracer, start the pathTracer, and then disappear.. 
// tracers is a pointer to a parent structure into which the pathTracer inserts itself at the end of execution.
// the parent can thus check. The mutex is to make sure the parent isn't doing anything with the mutex at that time.. 

class ExperimentTracer {
  public :
    ExperimentTracer(std::vector<probe_set*> ps, std::vector<uint> expt, std::set<PathTracer*>* tracers, QMutex* mutex, float sm);
  ~ExperimentTracer();    // but actually we are not going to do much with this..
};

// and actually there is not much need for the thingy.. 

#endif
