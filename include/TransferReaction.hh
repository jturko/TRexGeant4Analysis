/****************************************************************************
 *
 * Transfer reaction class:
 *  - contains the projectile, the target, the recoil and the ejectile
 *    of the 1n transfer experiment with a dPE target
 *  - kinematic splines
 *
 ***************************************************************************/

#ifndef TRANSFERREACTION_H_
#define TRANSFERREACTION_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Settings.hh"
#include "Nucleus.hh"
#include "Compound.hh"
#include "Kinematics.hh"
#include "Reconstruction.hh"

class TFile;

class TransferReaction {
private:
  // define the nuclei of the reaction
  void DefineReactionNuclei();

  // define target compound 
  // (needed for the calculation of the energy loss in the target)
  void DefineTargetCompound();

  // kinematics splines for elastic p, d and t
  void DoKinematicCalculationsForLightElastics();
  void WriteKinematicSplinesOfLightElastics(TFile &file);

  // transfer splines
  void DoKinematicCalculationsForTransferParticles();
  void WriteKinematicSplinesOfTransferParticles(TFile &file);


  //////////////////////////////
  // private member variables
  //////////////////////////////

  // file containing the experimental settings
  std::string fSettingFile; 

  // class to read the settings file
  Settings fSett;

  // nuclei of the reaction: projectile, target, ejectile and recoil
  Nucleus fProj, fTarg, fEjec, fReco;

  // target compound
  Compound fTarget;

  // target thickness
  double fTargetThickness;

  // kinematics
  Kinematics fElaP_beforeTarget, fElaP_middleTarget, fElaP_afterTarget;
  Kinematics fElaD_beforeTarget, fElaD_middleTarget, fElaD_afterTarget;
  
  Kinematics fTransferP_beforeTarget, fTransferP_middleTarget, fTransferP_afterTarget ;
 
  std::vector<Kinematics> fTransferP_Target;

public:
  // constructor
  TransferReaction(std::string SettingFile);
  TransferReaction(Settings&);
  
  // destructor
  virtual ~TransferReaction();


  ////////////
  // getter
  ////////////
  Settings* GetSettings() {return &fSett; };

  Nucleus* GetElasticEjectile() {return &fProj; };
  Nucleus* Get1nTransferEjectile() {return &fEjec; };

  double GetProjMass(){return fProj.GetMass(); }; // in MeV
  
  // spaeter loeschen
  double GetProjectileMass(){return fProj.GetMass(); }; // in MeV
  double GetEjectileMass(){return fEjec.GetMass(); }; // in MeV
  double GetRecoilMass(){return fReco.GetMass(); }; // in MeV
  double GetTargetMass(){return fTarg.GetMass(); }; // in MeV
  
  double GetElasticEjecMass() {return fProj.GetMass(); };
  double Get1nTransferEjecMass() {return fEjec.GetMass(); };

  Compound GetTargetMaterial(){return fTarget; };
 
  // getter for kinematics in the middle of the target
  Kinematics* GetElaP() {return &fElaP_middleTarget; };
  Kinematics* GetElaD() {return &fElaD_middleTarget; };
  Kinematics* GetGroundStateTransferP() {return &fTransferP_middleTarget; };
 

  //////////Beam energy along target for Q-value calculation/////////////
  std::vector<Kinematics>* GetGroundStateTransferP_V() {return &fTransferP_Target; };


  // do kinematic calculations 
  void DoKinematicCalculations(); 

  // write kinematic calculations into the Root file
  void WriteKinematicCalculationsToFile(TFile &file); 
};

#endif /* TRANSFERREACTIONDEUTERIUM_H_ */

