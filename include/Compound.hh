#ifndef __COMPOUND_HH
#define __COMPOUND_HH

#include <iostream>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

#include "Nucleus.hh"

class Compound{
 public:
  Compound();
  Compound(const char*);
  Compound(Nucleus*);
  ~Compound();
  //Compound(int nofelements, Nucleus* nuclei, double* fracs);
  void SetNofElements(int nofelements){
    fNofElements = nofelements;
  };
  int GetNofElements(){
    return fNofElements;
  };
  double GetMass(){
    return fMass;
  };
  const char* GetSymbol(){
    return fSymbol;
  };
  Nucleus* GetNucleus(int);
  double GetFrac(int);
 private:
  Nucleus** fNuclei;
  double* fFrac;
  int fNofElements;
  double fMass;
  const char* fSymbol;
};
#endif
