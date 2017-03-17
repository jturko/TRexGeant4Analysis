#ifndef __HIT_HH
#define __HIT_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TFile.h"
#include "TVector3.h"

#include "ParticleMC.hh"
#include "Settings.hh"

#ifndef PI
#define PI                       (TMath::Pi())
#endif

class HitSim {
public:
  HitSim(Settings* settings);
 
  void Clear();
  void InitBarrel(ParticleMC* barrel, std::string direction);
  //secondbarrel
  void InitSecondBarrel(ParticleMC* secondbarrel, std::string direction);                    //#B.Wach

  TVector3 BPosition(bool smear);
  //SecondBPosition
  TVector3 SecondBPosition(bool smear);                                                //#B.Wach

private:
  Settings* fSett;
  TRandom* fRand;
  ParticleMC* fBarrel;
  ParticleMC* fSecondBarrel;                                                    //#B.Wach
  std::string fDirection;
  
  double fMaxStrip[4];
  double fMinStrip[4];
};
#endif
