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

#include "server/server.h"
#include "version.h"
#include <qstring.h>
#include <qapplication.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char** argv){

  #ifndef EXPERIMENTAL_SERVER
  int dem = daemon(1, 1);
  cout << "deamon returned " << dem << endl;
  #endif

  QApplication app(argc, argv, false);

  // parse the commandline for options for the database name and the port to use..
  int port = 8109;    // the default port for this version.. 

  const char* db = "expression";  // the default..
  const char* dataTable = "data";


  char* db_env = getenv("EXPRESSION_DATABASE");  // if environment variable is defined override the option
  if(db_env){
    db = db_env;
  }

  int c;
  bool report_version = false;
  while((c = getopt(argc, argv, "d:p:t:v")) != -1){
    switch(c){
    case 'd':
      db = optarg;
      break;
    case 'p':
      port = atoi(optarg);
      break;
    case 't':
      dataTable = optarg;
      break;
    case 'v':
      report_version = true;
      break;
    case '?':
      cerr << "Unknown option : " << c << endl;
    default :
      abort();
    }
  }

  if(report_version){
    cout << "mtserver version " << server_version << endl << version_info << endl;
    exit(0);
  }

  cout << "port is : " << port << "  database is : " << db << endl;

  Server server(db, port, dataTable);

  return app.exec();
}

