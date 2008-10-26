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

#ifndef SERVERWINDOW_H
#define SERVERWINDOW_H

//#include <qtextview.h>
#include <qstring.h>
#include <qobject.h>
#include <string>
#include <vector>
#include "server/server.h"

class ServerWindow : public QObject
{
  Q_OBJECT

    public:
  ServerWindow(QObject* parent=0, const char* name=0);
  ~ServerWindow();

  private slots:

 private:
  Server* server;

};

#endif
