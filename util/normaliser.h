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

#ifndef NORMALISER_H
#define NORMALISER_H

#include <math.h>     // well you never know..

class Normaliser {          // a class containing some convenient normalising functions.. 
  public :
    Normaliser(){}
  ~Normaliser(){}

  void zScore(float* v, unsigned int s);
  void zScore(float** v, unsigned int ps, unsigned int es);    // as above but for a whole load in one go.. 

  void mScore(float** v, unsigned int ps, unsigned int s);     // don't know how to describe this, but base the normalisation 
                                               // on the mean of the std deviations and means of the individual probe series
};

class Euclid {
  public :
    Euclid(){}

  float euclidean(float* v1, float* v2, unsigned int l);
  float sqEuclidean(float* v1, float* v2, unsigned int l);   // returns the square euclidean distance.. 
};

#endif
