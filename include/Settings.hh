#ifndef __SETTINGS_HH
#define __SETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TSystem.h"
#include "TEnv.h"

#include "TRexSettings.hh"

class Settings : public TRexSettings {
public:
  Settings();//default ctor
  Settings(const char*);
  Settings(const TRexSettings&);
  ~Settings();

  void ReadSettings(const char* filename = nullptr);
  void PrintSettings();

  double FBMaxStripPos(int quadr)    { if(0 <= quadr && quadr < 4) return fFBMaxStripPos[quadr];    return 0.; }
  double FBMinStripPos(int quadr)    { if(0 <= quadr && quadr < 4) return fFBMinStripPos[quadr];    return 0.; }
  double BBMaxStripPos(int quadr)    { if(0 <= quadr && quadr < 4) return fBBMaxStripPos[quadr];    return 0.; }
  double BBMinStripPos(int quadr)    { if(0 <= quadr && quadr < 4) return fBBMinStripPos[quadr];    return 0.; }

  bool SmearStrip() { return fSmearStrip; }
  bool IncludeDeadLayers() { return fDeadLayers; }

protected:

  int fVerboseLevel;
  
  double fFBMaxStripPos[4];
  double fFBMinStripPos[4];
  double fBBMaxStripPos[4];
  double fBBMinStripPos[4];

  bool fSmearStrip;
  bool fDeadLayers;
};
#endif
