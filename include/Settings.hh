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

  bool SmearStrip() { return fSmearStrip; }
  bool IncludeDeadLayers() { return fDeadLayers; }

protected:
  int fVerboseLevel;
  
  bool fSmearStrip;
  bool fDeadLayers;
};
#endif
