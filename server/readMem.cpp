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

/*
 * shm-client - client program to demonstrate shared memory.
 */
#include "../raw/probeSetSet2.h"
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
    int shmid;
    key_t key;
    //    void* shm;   // the shared memory.. 
    //char *shm, *s;
    ProbeSetSet2** shm;
    /*
     * We need to get the segment named
     * "5678", created by the server.
     */
    key = 56789;
    int shsize = sizeof(ProbeSetSet2**);
    //int shsize = 27;
    /*
     * Locate the segment.
     */
    if ((shmid = shmget(key, shsize, 0666)) < 0) {
      cerr << "couldn't get access to the memoery.. ";
      return(1);
    }

    /*
     * Now we attach the segment to our data space.
     */
    shm = (ProbeSetSet2**)shmat(shmid, NULL, 0);
    if ((int)shm == -1) {
      cerr << "couldn't attach the memory " << endl;
      return(1);
    }
    // but just checking so..
    //for(s = shm; *s != NULL; s++){
    //  cout << *s;
    //}
    cout << endl;
    // if I get here, screw it let's just print something out..
    cout << "and shmid is : " << shmid << endl;
    cout << "and the value of the shm is : " << shm << endl;
    cout << "and the location of shm is : " << &shm << endl;
    cout << "And what if I deref the shm : " << *shm << endl;
    // can I cast the value of shm to a raw probe set .. 
    //ProbeSetSet2* pSet = (ProbeSetSet2*)shm;    // I doubt this is going to work though...
    cout << "and the size of data is " << (*shm)->data.size() << endl;
}
