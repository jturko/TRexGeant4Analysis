#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <stdlib.h>     // for (rand, srand)
#include <time.h>       // to add (time)
#include <string>

#include "CommandLineInterface.hh"
#include "TRexSettings.hh"
#include "HitSim.hh"
#include "ParticleMC.hh"
#include "Compound.hh"
#include "Reconstruction.hh"
#include "Kinematics.hh"
#include "Particle.hh"
#include "Germanium.hh"

using namespace TMath;
using namespace std;

ClassImp(ParticleMC);
ClassImp(Particle);

//Ekin_cm = m_d*Ekin_lab/(m_d+m_132Sn)
 
/**********************************************************************************************************************************
 * main function
 *********************************************************************************************************************************/
int main(int argc, char* argv[]) {
	char* InputFile;
	char* OutputFile = nullptr;
	string particleType;
	bool verbose = false;
	Long64_t maxEntries = -1;
	bool writeTree = false;
	bool dontSmear = false;
	bool doubleSidedFirstLayer = false;
	// command line interface
	CommandLineInterface* interface = new CommandLineInterface();

	interface->Add("-i", "inputfile (result Root file from simulation)", &InputFile);
	interface->Add("-o", "outputfile", &OutputFile);
	interface->Add("-particleType", "p or d", &particleType);
	interface->Add("-v", "verbose mode", &verbose);
	interface->Add("-t", "write output tree", &writeTree);
	interface->Add("-nosmear", "don't smear strips", &dontSmear);
	interface->Add("-n", "max. # of entries processed", &maxEntries);
	interface->Add("-d", "double sided first layer", &doubleSidedFirstLayer);
	interface->CheckFlags(argc, argv);

	if(InputFile == nullptr || OutputFile == nullptr) {
		std::cerr<<"You have to provide at least one input file and the output file!"<<std::endl;
		return 1;
	}
	std::cout<<"input file: "<<InputFile<<std::endl;

	// open input file
	TFile infile(InputFile);

	// load settings file
	TRexSettings* sett = static_cast<TRexSettings*>(infile.Get("settings"));//new Settings(SettingFile);
	if(sett == nullptr) {
		std::cerr<<"Failed to find \"settings\" in \""<<InputFile<<"\": "<<infile.Get("settings")<<std::endl;
		return 1;
	}
	if(verbose) sett->Print();

	std::cout<<"beam N = "<<sett->GetProjectileA()-sett->GetProjectileZ()<<", Z = "<<sett->GetProjectileZ() 
		<< " on target N = "<<sett->GetTargetA()-sett->GetTargetZ()<<", Z = "<<sett->GetTargetZ()<<std::endl;

	double beamEnergy = sett->GetBeamEnergy(); //initial beam energy (total) in MeV

	bool isSolid = (sett->GetGasTargetLength() == 0.);
	if (verbose) {
		if(isSolid) cout<<"Using a solid target!"<<endl;
		else        cout<<"Using a gas taget!"<<endl;
	}
	// nStrips are only used for the second layer (which is double-sided), and we assume that forward and backward detectors are the same
	int nStripsX = static_cast<int>(sett->GetSecondFBarrelDeltaESingleLengthX()/sett->GetSecondFBarrelDeltaESingleStripWidthPer());
	int nStripsY = static_cast<int>(sett->GetSecondFBarrelDeltaESingleLengthY()/sett->GetSecondFBarrelDeltaESingleStripWidthPar());

	std::cout<<"number of strips in x "<<nStripsX<<" ("<<sett->GetSecondFBarrelDeltaESingleLengthX()<<"/"<<sett->GetSecondFBarrelDeltaESingleStripWidthPer()<<")"<<std::endl;
	std::cout<<"number of strips in y "<<nStripsY<<" ("<<sett->GetSecondFBarrelDeltaESingleLengthY()<<"/"<<sett->GetSecondFBarrelDeltaESingleStripWidthPar()<<")"<<std::endl;

	////////////////////////////////////
	// prepare input and output trees
	////////////////////////////////////
	// prepare input tree
	TTree* trGen = (TTree*) infile.Get("treeGen");
	TTree* tr = (TTree*) infile.Get("treeDet");
	if(tr == nullptr) {
		cout<<"could not find tree tr in file "<<infile.GetName()<<endl;
		return 3;
	}

	vector<ParticleMC>* firstDeltaE[2] =  {nullptr, nullptr}; // new vector<ParticleMC>, new vector<ParticleMC>};
	vector<ParticleMC>* secondDeltaE[2] = {nullptr, nullptr}; // new vector<ParticleMC>, new vector<ParticleMC>};
	vector<ParticleMC>* pad[2] =          {nullptr, nullptr}; // new vector<ParticleMC>, new vector<ParticleMC>};
	vector<Germanium>* mBall = {nullptr};

	tr->SetBranchAddress("FBarrelDeltaEMC",&firstDeltaE[0]);
	tr->SetBranchAddress("SecondFBarrelDeltaEMC",&secondDeltaE[0]);
	tr->SetBranchAddress("FBarrelErestMC",&pad[0]);
	tr->SetBranchAddress("BBarrelDeltaEMC",&firstDeltaE[1]);
	tr->SetBranchAddress("SecondBBarrelDeltaEMC",&secondDeltaE[1]);
	tr->SetBranchAddress("BBarrelErestMC",&pad[1]);
	tr->SetBranchAddress("Miniball", &mBall);
	
	double reactionEnergyBeam;
	trGen->SetBranchAddress("reactionEnergy", &reactionEnergyBeam);
	
	double reactionEnergyBeamCM;
	trGen->SetBranchAddress("reactionEnergyCM", &reactionEnergyBeamCM);

	double reactionXSim;
	trGen->SetBranchAddress("reactionX", &reactionXSim);

	double reactionYSim;
	trGen->SetBranchAddress("reactionY", &reactionYSim);

	double reactionZSim;
	trGen->SetBranchAddress("reactionZ", &reactionZSim);

	double recoilThetaSim;                 //in radiants
	trGen->SetBranchAddress("recoilTheta", &recoilThetaSim);

	double recoilPhiSim;
	trGen->SetBranchAddress("recoilPhi", &recoilPhiSim);

	double recoilEnergySim;
	trGen->SetBranchAddress("recoilEnergy", &recoilEnergySim);
	
	std::vector<double>* gammaThetaSim = new std::vector<double>;
    trGen->SetBranchAddress("gammaTheta", &gammaThetaSim);
    
    std::vector<double>* gammaEnergySim = new std::vector<double>;
    trGen->SetBranchAddress("gammaEnergy", &gammaEnergySim);
    
    double ejectileEnergySim;
    trGen->SetBranchAddress("ejectileEnergy", &ejectileEnergySim);

	UInt_t reactionSim;
	trGen->SetBranchAddress("reaction", &reactionSim);

	// create output file and output tree
	TFile outfile(OutputFile,"recreate");  

	TTree* rectr = new TTree("rectr","reconstructed events");
	vector<Particle>* ParticleBranch = new vector<Particle>;
	rectr->Branch("Particle",&ParticleBranch);

	// initialize hit class
	HitSim* hit = new HitSim(sett);

	// set massfile
	const char* massfile = sett->GetMassFile().c_str();
	std::cout<<"Mass file "<<massfile<<std::endl;

	//set Compounds for now we always assume density of H2 = 0.180e-3 g/cm3 at 1 bar
	if(verbose) std::cout<<"creating compound \""<<sett->GetTargetMaterialName().c_str()<<"\""<<std::endl;
	Compound* targetMat = new Compound(sett->GetTargetMaterialName().c_str());
	//if(!isSolid) targetMat->SetDensity(0.180e-3*sett->GetTargetPressure()/6.24151e+08); // internal geant4 units are such that 1 bar = 6.24151e+08 whatever it is MeV/mm3:-)))
	if(!isSolid) targetMat->SetDensity(targetMat->GetDensity()*sett->GetTargetPressure()/6.24151e+08); // 4He gas density in Geant Material Lib.
	std::cout<<"\n set target material density to "<<0.180e-3*sett->GetTargetPressure()/6.24151e+08<<" = "<<targetMat->GetDensity()<< "g/cm3"<<std::endl;
	
    if(verbose) std::cout<<"creating compound \"MY\""<<std::endl;
	Compound* foilMat = new Compound("MY");
	if(verbose) std::cout<<"creating compound \""<<sett->GetVacuumChamberGas().c_str()<<"\""<<std::endl;
	Compound* chamberGasMat = new Compound(sett->GetVacuumChamberGas().c_str());
	chamberGasMat->SetDensity(chamberGasMat->GetDensity()*sett->GetVacuumChamberGasPressure()/6.24151e+08); // internal geant4 units are such that 1 bar = 6.24151e+08 whatevers *** changed by Leila
	std::cout<<"set chamber gas density to "<<0.18e-3*sett->GetVacuumChamberGasPressure()/6.24151e+08<<" = "<<chamberGasMat->GetDensity()<<std::endl;// is 0.180e-3 correct? yes, He gas density is 0.179e-3 g/cm3 at stp (1 bar & 0 deg). Leila!
	
    Compound* layerMat = new Compound("silicon");

	// 1n transfer
	Nucleus* projectile = new Nucleus(sett->GetProjectileZ(), sett->GetProjectileA() - sett->GetProjectileZ(),   massfile);
	Nucleus* target     = new Nucleus(sett->GetTargetZ(),     sett->GetTargetA()-sett->GetTargetZ(),             massfile);
	Nucleus* ejectile   = new Nucleus(sett->GetProjectileZ(), sett->GetProjectileA()-sett->GetProjectileZ() + 1, massfile);
	Nucleus* recoil     = new Nucleus(sett->GetTargetZ(),     sett->GetTargetA()-sett->GetTargetZ() - 1,         massfile);

	std::cout<<"projectile "<<projectile->GetSymbol()<<" ("<<projectile->GetA()<<", "<<projectile->GetZ()<<"; "<<projectile->GetMass()<<")"<<std::endl;
	std::cout<<"target "<<target->GetSymbol()<<" ("<<target->GetA()<<", "<<target->GetZ()<<"; "<<target->GetMass()<<")"<<std::endl;
	std::cout<<"ejectile "<<ejectile->GetSymbol()<<" ("<<ejectile->GetA()<<", "<<ejectile->GetZ()<<"; "<<ejectile->GetMass()<<")"<<std::endl;
	std::cout<<"recoil "<<recoil->GetSymbol()<<" ("<<recoil->GetA()<<", "<<recoil->GetZ()<<"; "<<recoil->GetMass()<<")"<<std::endl;

	// transfer reaction and kinematics class for Q-value calculation
	Kinematics* transferP = new Kinematics(projectile, target, recoil, ejectile, beamEnergy, 0.); //reaction.GetGroundStateTransferP();

	//Kinematics* transferP = new Kinematics(projectile, target, beamEnergy); //reaction.GetGroundStateTransferP();

	// variables for reconstruction
	Reconstruction* beamTarget = new Reconstruction(projectile, targetMat);
	double beamEnergyRec;      //beam energy at reaction, reconstructed 
	double recoilEnergyRec;    //recoil energy, reconstructed
	double recoilEnergyRecdE;
	double recoilEnergyRecErest;
	double recoilThetaRec;
	double recoilPhiRec;

	double targetThickness = sett->GetTargetThicknessMgPerCm2();
	double targetLength    = sett->GetTargetPhysicalLength();
	double targetForwardZ  =   targetLength/2.;
	double targetBackwardZ = - targetLength/2.;
	if(verbose) { 
		cout <<"Target Thickness from Input File: "<< targetThickness<<endl;
		cout <<"Target ForwardZ from Input File: "<< targetForwardZ<<endl;
		cout <<"Target BackwardZ from Input File: "<< targetBackwardZ<<endl;
		cout <<"Target Length from Input File: "<< targetLength<<endl;
	}
	TSpline3* energyInTarget = beamTarget->Thickness2EnergyAfter(beamEnergy, targetThickness, targetThickness/1000., true);
	energyInTarget->Write("energyInTarget");

	// variables for recoil energy loss reconstruction (assuming max. recoil energy of 100 MeV!!!)
	Reconstruction* recoilTarget = new Reconstruction(recoil, targetMat);
	Reconstruction* recoilFoil = new Reconstruction(recoil, foilMat);
	Reconstruction* recoilLayer = new Reconstruction(recoil, layerMat);
	Reconstruction* recoilChamberGas = new Reconstruction(recoil, chamberGasMat);
	TSpline3* recoilTargetRange  = recoilTarget->Energy2Range(100., 0.1, !isSolid);
	TSpline3* recoilTargetEnergy = recoilTarget->Range2Energy(100., 0.1, !isSolid);
	TSpline3* recoilFoilRange  = recoilFoil->Energy2Range(100., 0.1, false);
	TSpline3* recoilFoilEnergy = recoilFoil->Range2Energy(100., 0.1, false);
	TSpline3* recoilLayerRange  = recoilLayer->Energy2Range(100., 0.1, false);
	TSpline3* recoilLayerEnergy = recoilLayer->Range2Energy(100., 0.1, false);
	TSpline3* recoilChamberGasRange  = recoilChamberGas->Energy2Range(100., 0.1, true);
	TSpline3* recoilChamberGasEnergy = recoilChamberGas->Range2Energy(100., 0.1, true);
	
	double foilDistance = sett->GetTargetDiameter()/2.; 
	
	double firstLayerDistance = sett->GetFBarrelDeltaESingleDistanceToBeam()[0];
	double secondLayerDistance = sett->GetSecondFBarrelDeltaESingleDistanceToBeam()[0];
	double padDistance = sett->GetFBarrelErestSingleDistanceToBeam()[0];
	
	double targetWidthMgCm2 = foilDistance*100.*targetMat->GetDensity();// convert from mm to cm, and from there to mg/cm2 using density
    double foilThicknessMgCm2 = sett->GetFBarrelDeltaESingleFoilThickness()*0.1*foilMat->GetDensity();// convert from um to cm, and from there to mg/cm2 using the solid mylar density     
    double firstGasLayerThicknessMgCm2 = (firstLayerDistance - foilDistance)*100.*chamberGasMat->GetDensity(); 
    double secondGasLayerThicknessMgCm2 = (secondLayerDistance - firstLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 1. and 2. layer in cm and converted to mg/cm2 using density
    double thirdGasLayerThicknessMgCm2 = (padDistance - secondLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 2. layer and pad in cm and converted to mg/cm2 using density
    double secondLayerThicknessMgCm2 = sett->GetSecondFBarrelDeltaESingleThickness()[0]*0.1*layerMat->GetDensity(); // density Silicon = 2.328 g/cm3
    double firstLayerThicknessMgCm2 = sett->GetFBarrelDeltaESingleThickness()[0]*0.1*layerMat->GetDensity(); // density Silicon = 2.328 g/cm3	

	std::cout<<"\n foil distance "<<foilDistance<<" mm => target width = "<<targetWidthMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"first layer distance "<<firstLayerDistance<<" mm and second layer: "<<secondLayerDistance<<" mm"<<std::endl;
	std::cout<<"foil thickness "<<sett->GetFBarrelDeltaESingleFoilThickness()<<" mm => foilThicknessMgCm2 "<<foilThicknessMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"1. layer thickness "<<sett->GetFBarrelDeltaESingleThickness()[0]<<" mm => 1. layerThicknessMgCm2 "<<firstLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"2. layer thickness "<<sett->GetSecondFBarrelDeltaESingleThickness()[0]<<" mm => 2. layerThicknessMgCm2 "<<secondLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"distance foil - 1. layer "<<firstLayerDistance - foilDistance<<" mm => 1. gas layer thickness = "<<firstGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl; // Leila
	std::cout<<"distance 1. layer - 2. layer "<<secondLayerDistance - firstLayerDistance<<" mm and chambergasdensity "<<chamberGasMat->GetDensity()<<" => 2. gas layer thickness = "<<secondGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"distance 2. layer - pad "<<padDistance - secondLayerDistance<<" mm => 3. gas layer thickness = "<<thirdGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;

	// Define Histograms
	TList list;
	TH2F* originXY = new TH2F("originXY", "y vs. x of reconstructed origin", 200, -10., 10., 200, -10., 10.); list.Add(originXY);
	TH2F* originXYErr = new TH2F("originXYErr", "Error y vs. error x of reconstructed origin - simulated origin", 200, -10., 10., 200, -10., 10.); list.Add(originXYErr);
	TH2F* errorOrigin = new TH2F("errorOrigin", "Error between reconstructed and true origin vs. true origin", 200, -100., 100., 1000, -5, 5); list.Add(errorOrigin);
	TH2F* errorThetaPhi = new TH2F("errorThetaPhi", "Error between reconstructed and true phi vs. error in theta", 600, -30, 30, 720, -360., 360.); list.Add(errorThetaPhi);
	TH1F* excEnProton = new TH1F("excEnProton", "Excitaiton Energy Spectrum from reconstructed Protons", 5000, -20000, 20000); list.Add(excEnProton);
	TH1F* reaction = new TH1F("reaction", "Simulated reaction/level", 10, -0.5, 9.5); list.Add(reaction);
	TH2F* phiErrorVsPhi = new TH2F("phiErrorVsPhi","Error in reconstructed #varphi vs. simulated #varphi", 360, -180., 180., 720, -360., 360.); list.Add(phiErrorVsPhi);
	TH2F* phiErrorVsPhiF0 = new TH2F("phiErrorVsPhiF0","Error in reconstructed #varphi vs. simulated #varphi, forward  detector #0 only", 360, -180., 180., 720, -360., 360.); list.Add(phiErrorVsPhiF0);
	TH2F* phiErrorVsPhiB0 = new TH2F("phiErrorVsPhiB0","Error in reconstructed #varphi vs. simulated #varphi, backward detector #0 only", 360, -180., 180., 720, -360., 360.); list.Add(phiErrorVsPhiB0);

	TH2F* dE12VsPad = new TH2F("dE12VsPad", "energy loss first+second layer vs. pad energy", 2000, 0, 50000, 1000, 0, 10000); list.Add(dE12VsPad);
	TH2F* dE12VsE = new TH2F("dE12VsE", "energy loss first+second layer vs. total energy", 2000, 0, 50000, 1000, 0, 10000); list.Add(dE12VsE);
	TH2F* dE1VsE = new TH2F("dE1VsE", "energy loss first layer vs. total energy", 5000, 0, 50000, 1000, 0, 10000); list.Add(dE1VsE);
	TH2F* dE2VsE = new TH2F("dE2VsE", "energy loss second layer vs. total energy", 5000, 0, 50000, 1500, 0, 15000); list.Add(dE2VsE);
	TH2F* dE1VsdE2 = new TH2F("dE1VsdE2", "energy loss second layer vs. energy loss first layer", 1000, 0, 5000, 2000, 0, 10000); list.Add(dE1VsdE2);
	TH2F* eVsTheta = new TH2F("eVsTheta", "recoil energy vs. theta (lab)", 360, 0, 180, 1000, 0, 35000); list.Add(eVsTheta);
	TH2F* eVsZ = new TH2F("eVsZ", "recoil energy vs. z", 1000, -100., 100., 1000, 0, 35000); list.Add(eVsZ);
	TH2F* eVsZThetaCut = new TH2F("eVsZThetaCut", "recoil energy vs. z with 30<#Theta<50", 1000, -100., 100., 1000, 0, 35000); list.Add(eVsZThetaCut);
	TH2F* eVsZSame = new TH2F("eVsZSame", "recoil energy vs. z, first and second layer both forward or both backward", 1000, -100., 100., 1000, 0, 35000); list.Add(eVsZSame);
	TH2F* eVsZCross = new TH2F("eVsZCross", "recoil energy vs. z, first and second layer over cross", 1000, -100., 100., 1000, 0, 35000); list.Add(eVsZCross);
	
	// ******************************** my histograms start ************************************
    
	TH2F* dE2VsZ = new TH2F("dE2VsZ","energy loss second layer vs. z",1000,-100,100,2000,0,10000);list.Add(dE2VsZ);	
	TH2F* dE2VsY = new TH2F("dE2VsY","energy loss second layer vs. y",1000,-100,100,2000,0,10000);list.Add(dE2VsY);	
	TH2F* dE2VsX = new TH2F("dE2VsX","energy loss second layer vs. x",1000,-100,100,2000,0,10000);list.Add(dE2VsX);	
	TH2F* dE2VsR = new TH2F("dE2VsR","energy loss second layer vs. R",1000,-100,100,2000,0,10000);list.Add(dE2VsR);	
	
	TH1F* ErecoilSimMindE12  = new TH1F("ErecoilSimMindE12","recoil energy sim - (dE1 + dE2)",2000,-500,39500);list.Add(ErecoilSimMindE12);	
	TH1F* ErecoilSimMindE12Epad0  = new TH1F("ErecoilSimMindE12Epad0","recoil energy sim - (dE1 + dE2) for Epad0",2000,-500,39500);list.Add(ErecoilSimMindE12Epad0);
	TH2F* ErecoilSimMindE12VsdE2  = new TH2F("ErecoilSimMindE12VsdE2","recoil energy sim - (dE1 + dE2) vs dE2 (y)",2000,-500,39500,2000,0,10000);list.Add(ErecoilSimMindE12VsdE2);
	TH2F* ErecoilSimMindE1VsdE2  = new TH2F("ErecoilSimMindE1VsdE2","recoil energy sim - dE1 vs dE2 (y)",2000,-500,39500,2000,0,10000);list.Add(ErecoilSimMindE1VsdE2);
	TH2F* ErecoilSimMindE1VsdE2Epad0  = new TH2F("ErecoilSimMindE1VsdE2Epad0","recoil energy sim - dE1 vs dE2 (y) for Epad=0",2000,-500,39500,2000,0,10000);list.Add(ErecoilSimMindE1VsdE2Epad0);
	TH2F* ErecoilSimMindE12VsdE2Epad0  = new TH2F("ErecoilSimMindE12VsdE2Epad0","recoil energy sim - (dE1 + dE2) vs dE2 (y) for Epad=0",2000,-500,39500,2000,0,10000);list.Add(ErecoilSimMindE12VsdE2Epad0);
	TH2F* ErecoilSimMindE12VsdE1  = new TH2F("ErecoilSimMindE12VsdE1","recoil energy sim - (dE1 + dE2) vs dE1 (y)",2000,-500,39500,1000,0,5000);list.Add(ErecoilSimMindE12VsdE1);	
	TH2F* ErecoilSimMindE12VsdE1Epad0  = new TH2F("ErecoilSimMindE12VsdE1Epad0","recoil energy sim - (dE1 + dE2) vs dE1 (y) for Epad=0",2000,-500,39500,1000,0,5000);list.Add(ErecoilSimMindE12VsdE1Epad0);
	TH2F* ErecoilSimMindE12VsThetaLab  = new TH2F("ErecoilSimMindE12VsThetaLab","recoil energy sim - (dE1 + dE2) vs #vartheta_{lab}",1000,-20,180,2000,-500,39500);list.Add(ErecoilSimMindE12VsThetaLab);
	TH2F* ErecoilSimMindE1VsThetaLab  = new TH2F("ErecoilSimMindE1VsThetaLab","recoil energy sim - dE1 vs #vartheta_{lab}",1000,-20,180,2000,-500,39500);list.Add(ErecoilSimMindE1VsThetaLab);
	TH2F* ErecoilSimMindE2VsThetaLab  = new TH2F("ErecoilSimMindE2VsThetaLab","recoil energy sim - dE2 vs #vartheta_{lab}",1000,-20,180,2000,-500,39500);list.Add(ErecoilSimMindE2VsThetaLab);
	TH2F* ErecoilSimMindE12VsThetaLabEpad0  = new TH2F("ErecoilSimMindE12VsThetaLabEpad0","recoil energy sim - (dE1 + dE2) vs #vartheta_{lab} for Epad=0",1000,-20,180,2000,-500,39500);list.Add(ErecoilSimMindE12VsThetaLabEpad0);
	TH2F* ErecoilSimMindE12VsdThetaLab  = new TH2F("ErecoilSimMindE12VsdThetaLab","recoil energy sim - (dE1 + dE2) vs #Delta#vartheta_{lab}",900,-30,30,2000,-500,39500);list.Add(ErecoilSimMindE12VsdThetaLab);
	TH2F* ErecoilSimMindE12VsdThetaLabEpad0  = new TH2F("ErecoilSimMindE12VsdThetaLabEpad0","recoil energy sim - (dE1 + dE2) vs #Delta#vartheta_{lab} for Epad=0",900,-30,30,2000,-500,39500);list.Add(ErecoilSimMindE12VsdThetaLabEpad0);	
	TH2F* ErecoilSimVsThetaLab  = new TH2F("ErecoilSimVsThetaLab","recoil energy sim vs #vartheta_{lab}",1000,-20,180,2000,-500,39500);list.Add(ErecoilSimVsThetaLab);
	TH2F* ErecoilSimVsZ  = new TH2F("ErecoilSimVsZ","recoil energy sim vs z",1000,-100,100,2000,-500,39500);list.Add(ErecoilSimVsZ);
	TH2F* ErecoilSimMindE12VsZ  = new TH2F("ErecoilSimMindE12VsZ","recoil energy sim - (dE1 + dE2) vs z",1000,-100,100,2000,-500,39500);list.Add(ErecoilSimMindE12VsZ);
	TH2F* ErecoilSimVsRec = new TH2F("ErecoilSimVsRec","recoil energy simulated (y) vs reconstructed (x)",2000,-500,39500,2000,-500,39500);list.Add(ErecoilSimVsRec);
	TH2F* ErecoilSimVsRecEpad0 = new TH2F("ErecoilSimVsRecEpad0","recoil energy simulated (y) vs reconstructed (x) for Epad=0",2000,-500,39500,2000,-500,39500);list.Add(ErecoilSimVsRecEpad0);	
	TH1F* ErecoilSim  = new TH1F("ErecoilSim","recoil energy simulated",2000,-500,39500);list.Add(ErecoilSim);
	TH1F* ErecoilRec  = new TH1F("ErecoilRec","recoil energy reconstructed",2000,-500,39500);list.Add(ErecoilRec);	
	TH2F* recoil2posZVsX = new TH2F("recoil2posZVsX", "z (x-axis) vs x positions on the second layer", 1000, -100., 100., 1000, -100., 100.); list.Add(recoil2posZVsX);
	TH2F* recoil1posZVsX = new TH2F("recoil1posZVsX", "z (x-axis) vs x positions on the first layer", 1000, -100., 100., 1000, -100., 100.); list.Add(recoil1posZVsX);
	TH2F* recoilVertexZVsX = new TH2F("recoilVertexZVsX", "z (x-axis) vs x  vertex positions", 1000, -100., 100., 1000, -100., 100.); list.Add(recoilVertexZVsX);
	
	TH1F* hERec = new TH1F("hERec", "reconstructed energy of recoil", 1000, -5000., 35000.); list.Add(hERec);
	TH1F* hERecdE = new TH1F("hERecdE", "reconstructed dE energy of recoil", 1000, -5000., 35000.); list.Add(hERecdE);
	TH1F* hERecErest = new TH1F("hERecErest", "reconstructed pad energy of recoil", 10000, -5000., 35000.); list.Add(hERecErest);
	TH2F* hERecErestVsTheta = new TH2F("hERecErestVsTheta", "Erest vs Theta", 1000, -200., 200., 1000, -20000., 40000.);list.Add(hERecErestVsTheta);
	TH2F* hERecdEVsTheta = new TH2F("hERecdEVsTheta", "dE vs Theta", 1000, -200., 200., 1000, -20000., 20000.); list.Add(hERecdEVsTheta);
	TH2F* hrecoilERecErestElossVsThetaLeila = new TH2F("hrecoilERecErestElossVsThetaLeila", "Erest loss vs Theta leila", 1000, -200., 200., 1000, -20000., 40000.); list.Add(hrecoilERecErestElossVsThetaLeila);
	TH2F* hrecoilERecdE2ElossVsThetaLeila = new TH2F("hrecoilERecdE2ElossVsThetaLeila", "dE2 loss vs Theta leila", 1000, -200., 200., 1000, -10000., 40000.); list.Add(hrecoilERecdE2ElossVsThetaLeila);
	TH2F* hrecoilERecdE1ElossVsThetaLeila = new TH2F("hrecoilERecdE1ElossVsThetaLeila", "dE1 loss vs Theta leila", 1000, -200., 200., 1000, -10000., 40000.); list.Add(hrecoilERecdE1ElossVsThetaLeila);
	TH2F* hrecoilERecFoilElossVsThetaLeila = new TH2F("hrecoilERecFoilElossVsThetaLeila", "Foil energy loss vs Theta leila", 1000, -200., 200., 1000, -10000., 40000.); list.Add(hrecoilERecFoilElossVsThetaLeila);
	TH2F* hrecoilERecTargetElossVsThetaLeila = new TH2F("hrecoilERecTargetElossVsThetaLeila", "Target energy loss vs Theta leila", 1000, -200., 200., 1000, -10000., 40000.); list.Add(hrecoilERecTargetElossVsThetaLeila);
	TH1F* excEnProtonCorr = new TH1F("excEnProtonCorr", "Excitaiton Energy after E loss correction for reconstructed protons", 1000, -1000, 11000); list.Add(excEnProtonCorr);
	TH1F* excEnProtonCorrEpadCut = new TH1F("excEnProtonCorrEpadCut", "Excitaiton Energy after E loss correction for reconstructed protons with Epad>0", 5000, -20000, 20000); list.Add(excEnProtonCorrEpadCut);
	TH2F* excEnProtonCorrVsX = new TH2F("excEnProtonCorrVsX", "Excitation Energy vs Vertex X after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); list.Add(excEnProtonCorrVsX);
	TH2F* excEnProtonCorrVsY = new TH2F("excEnProtonCorrVsY", "Excitation Energy vs Vertex Y after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); list.Add(excEnProtonCorrVsY);
	TH2F* excEnProtonCorrVsZ = new TH2F("excEnProtonCorrVsZ", "Excitation Energy vs Vertex Z after E loss correction for reconstructed protons", 5000,-100,100,5000, -20000, 20000); list.Add(excEnProtonCorrVsZ);
	TH2F* excEnProtonCorrVsT = new TH2F("excEnProtonCorrVsT", "Excitation Energy vs Vertex T after E loss correction for reconstructed protons", 5000,-10,10,5000, -20000, 20000); list.Add(excEnProtonCorrVsT);
	TH2F* excEnProtonCorrVsR = new TH2F("excEnProtonCorrVsR", "Excitation Energy vs Vertex R after E loss correction for reconstructed protons", 4400,-10,100,5000, -20000, 20000); list.Add(excEnProtonCorrVsR);
	TH1F* excEnProtonCorrdE1Sigma1 = new TH1F("excEnProtonCorrdE1Sigma1", "Excitation Energy for dE1 1#sigma range after Eloss correction for protons", 5000, -5000, 5000); list.Add(excEnProtonCorrdE1Sigma1);
	TH1F* excEnProtonCorrdE1Sigma2 = new TH1F("excEnProtonCorrdE1Sigma2", "Excitation Energy for dE1 2#sigma range after Eloss correction for protons", 5000, -5000, 5000); list.Add(excEnProtonCorrdE1Sigma2);
	
	TH2F* hrecoilERecErestElossRange1VsTheta = new TH2F("hrecoilERecErestElossRange1VsTheta", "Erest loss range1 vs Theta", 1000, -200., 200., 1000, -20000., 10000.); list.Add(hrecoilERecErestElossRange1VsTheta);
	TH2F* hrecoilERecErestElossRange2VsTheta = new TH2F("hrecoilERecErestElossRange2VsTheta", "Erest loss range2 vs Theta", 1000, -200., 200., 1000, -20000., 10000.); list.Add(hrecoilERecErestElossRange2VsTheta);
	TH2F* hrecoilERecdE2ElossRangeVsTheta = new TH2F("hrecoilERecdE2ElossRangeVsTheta", "dE2 loss range vs Theta", 1000, -200., 200., 1000, -10000., 10000.); list.Add(hrecoilERecdE2ElossRangeVsTheta);
	TH2F* hrecoilERecdE1ElossRangeVsTheta = new TH2F("hrecoilERecdE1ElossRangeVsTheta", "dE1 loss range vs Theta", 1000, -200., 200., 1000, -10000., 10000.); list.Add(hrecoilERecdE1ElossRangeVsTheta);
	TH2F* hrecoilERecFoilElossRangeVsTheta = new TH2F("hrecoilERecFoilElossRangeVsTheta", "Foil energy loss range vs Theta", 1000, -200., 200., 1000, -10000., 10000.); list.Add(hrecoilERecFoilElossRangeVsTheta);
	TH2F* hrecoilERecTargetElossRangeVsTheta = new TH2F("hrecoilERecTargetElossRangeVsTheta", "Target energy loss range vs Theta", 1000, -200., 200., 1000, -10000., 10000.); list.Add(hrecoilERecTargetElossRangeVsTheta);
	
	TH1F* hdE1ElossRange = new TH1F("hdE1ElossRange", "corrected energy range in the first layer", 1000, -100., 3900.); list.Add(hdE1ElossRange);
	TH1F* hdE2ElossRange = new TH1F("hdE2ElossRange", "corrected energy range in the second layer", 2000, -1000., 9000.); list.Add(hdE2ElossRange);
	TH1F* hdE2ElossRangeWoEpad0 = new TH1F("hdE2ElossRangeWoEpad0", "corrected energy range in the second layer with Epad>0", 2000, -1000., 9000.); list.Add(hdE2ElossRangeWoEpad0);
	TH1F* hdE1Eloss = new TH1F("hdE1Eloss", "corrected energy loss in the first layer", 1000, -100., 9900.); list.Add(hdE1Eloss);
	TH1F* hdE2Eloss = new TH1F("hdE2Eloss", "corrected energy loss in the second layer", 2000, -100., 18900.); list.Add(hdE2Eloss);
    TH1F* hdE1Measured = new TH1F("hdE1Measured", "measured energy loss in the first layer", 1000, -100., 9900.); list.Add(hdE1Measured);
	TH1F* hdE2Measured = new TH1F("hdE2Measured", "measured energy loss in the second layer", 2000, -100., 18900.); list.Add(hdE2Measured);	
	TH1F* hErestMeasured = new TH1F("hErestMeasured", "measured energy loss in the pad", 2000, -1000., 39000.); list.Add(hErestMeasured);	
	TH2F* hdE1ElossVsMeasured = new TH2F("hdE1ElossVsMeasured", "corrected energy loss in the first layer vs measured", 1000, -100., 2400., 1000,-100.,2400.); list.Add(hdE1ElossVsMeasured);
	TH2F* hdE2ElossVsMeasuredWoEpad0 = new TH2F("hdE2ElossVsMeasuredWoEpad0", "corrected energy loss in the second layer vs measured for Epad>0", 2000, -100., 5400., 2000,-100.,5400.); list.Add(hdE2ElossVsMeasuredWoEpad0);	
	TH2F* hdE2ElossVsMeasured = new TH2F("hdE2ElossVsMeasured", "corrected energy loss in the second layer vs measured", 2000, -100., 9400., 2000,-100.,9400.); list.Add(hdE2ElossVsMeasured);	
    TH2F* dE2VsdE2Pad = new TH2F("dE2VsdE2Pad", "energy loss in second layer vs. dE2+pad energy", 2000, 0, 50000, 1000, 0, 10000); list.Add(dE2VsdE2Pad);	
    TH2F* EPadVsThetaLab = new TH2F("EPadVsThetaLab", "pad energy vs theta lab", 720,0,180.,2000, -1000, 49000); list.Add(EPadVsThetaLab);
    TH2F* EPadVsZ = new TH2F("EPadVsZ", "pad energy vs z", 1000,-100.,100.,2000, -1000, 49000); list.Add(EPadVsZ);			
    TH2F* dE2VsThetaLabEpadCut = new TH2F("dE2VsThetaLabEpadCut", "second layer energy vs theta lab for Epad>0", 720,0,180.,1000, -100, 4900); list.Add(dE2VsThetaLabEpadCut);
    TH2F* dE2VsEPadThetaCut = new TH2F("dE2VsEPadThetaCut", "energy loss in second layer vs. pad energy for Theta=40", 2000, -1000, 49000, 1000, -100, 9900); list.Add(dE2VsEPadThetaCut);
    TH2F* dE1VsThetaLab = new TH2F("dE1VsThetaLab", "first layer energy vs theta lab", 720,0,180.,2000, -100, 4900); list.Add(dE1VsThetaLab);
    TH2F* dE2VsThetaLab = new TH2F("dE2VsThetaLab", "second layer energy vs theta lab", 1800,0,180.,3000, -100, 14900); list.Add(dE2VsThetaLab);
	TH2F* dE12VsThetaLab = new TH2F("dE12VsThetaLab", "first+second layer energy vs theta lab", 720,0,180.,2000, -100, 7900); list.Add(dE12VsThetaLab);	
	TH2F* dE12VsPhiLab = new TH2F("dE12VsPhiLab", "first+second layer energy vs phi lab", 720,-180,180.,2000, -100, 7900); list.Add(dE12VsPhiLab);	
	TH2F* dE1VsPhiLab = new TH2F("dE1VsPhiLab", "first layer energy vs phi lab", 720,-180,180.,1000, -100, 5000); list.Add(dE1VsPhiLab);
	TH2F* dE2VsPhiLab = new TH2F("dE2VsPhiLab", "second layer energy vs phi lab", 720,-180,180.,2000, -100, 9900); list.Add(dE2VsPhiLab);
	TH2F* dE1EpadVsThetaLab = new TH2F("dE1EpadVsThetaLab", "first layer + pad energy vs theta lab", 720,0,180.,2000, -1000, 49900); list.Add(dE1EpadVsThetaLab);
	TH2F* dE2EpadVsThetaLab = new TH2F("dE2EpadVsThetaLab", "second layer + pad energy vs theta lab", 720,0,180.,2000, -1000, 49900); list.Add(dE2EpadVsThetaLab);	
	TH2F* hdE1ElossVsMeasuredEpad0 = new TH2F("hdE1ElossVsMeasuredEpad0", "corrected energy loss in the first layer vs measured for Epad=0", 1000, -100., 2400., 1000,-100.,2400.); list.Add(hdE1ElossVsMeasuredEpad0);
	TH2F* hdE1ElossVsMeasuredEpadWo0 = new TH2F("hdE1ElossVsMeasuredEpadWo0", "corrected energy loss in the first layer vs measured for Epad>0", 1000, -100., 2400., 1000,-100.,2400.); list.Add(hdE1ElossVsMeasuredEpadWo0);
	TH1F* hdE1MeasMinRec = new TH1F("hdE1MeasMinRec", "(measured - reconstructed) energy loss in the first layer", 1000, -500., 500.); list.Add(hdE1MeasMinRec);
	TH1F* hdE2MeasMinRec = new TH1F("hdE2MeasMinRec", "(measured - reconstructed) energy loss in the second layer", 1000, -500., 9500.); list.Add(hdE2MeasMinRec);
	
	TH1F* hrecoilEnergyRecEloss = new TH1F("hrecoilEnergyRecEloss","reconstructed energy loss of recoils", 1000, -5000., 35000.); list.Add(hrecoilEnergyRecEloss);
	TH1F* hrecoilThetaRec = new TH1F("hrecoilThetaRec", "reconstructed theta of recoil", 1000, -60., 300.); list.Add(hrecoilThetaRec);
	TH1F* hrecoilThetaRecDiffSim = new TH1F("hrecoilThetaRecDiffSim", "difference reconstructed and simulated theta of recoil", 1000, -60., 300.); list.Add(hrecoilThetaRecDiffSim);
	TH1F* hrecoilEnergyRecDiffSim = new TH1F("hrecoilEnergyRecDiffSim", "difference reconstructed and simulated energy of recoil", 1000, -25000., 25000.); list.Add(hrecoilEnergyRecDiffSim);
	
	TH2F* thetaVsZ2 = new TH2F("thetaVsZ2","#vartheta_{lab} vs. z in 2. layer", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZ2);
	TH2F* thetaVsZ1 = new TH2F("thetaVsZ1","#vartheta_{lab} vs. z in 1. layer", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZ1);
	TH2F* thetaVsX2 = new TH2F("thetaVsX2","#vartheta_{lab} vs. x in 2. layer", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsX2);
	TH2F* thetaVsX1 = new TH2F("thetaVsX1","#vartheta_{lab} vs. x in 1. layer", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsX1);
	TH2F* phiVsZ2 = new TH2F("phiVsZ2","#varphi_{lab} vs. z in 2. layer", 200, -100., 100., 180, 0., 180.); list.Add(phiVsZ2);
	TH2F* phiVsZ1 = new TH2F("phiVsZ1","#varphi_{lab} vs. z in 1. layer", 200, -100., 100., 180, 0., 180.); list.Add(phiVsZ1);
	TH2F* phiVsX2 = new TH2F("phiVsX2","#varphi_{lab} vs. x in 2. layer", 200, -100., 100., 180, 0., 180.); list.Add(phiVsX2);
	TH2F* phiVsX1 = new TH2F("phiVsX1","#varphi_{lab} vs. x in 1. layer", 200, -100., 100., 180, 0., 180.); list.Add(phiVsX1);
	TH2F* thetaVsX = new TH2F("thetaVsX","#vartheta_{lab} vs. vertex x ", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsX);
	TH2F* phiVsX = new TH2F("phiVsX","#varphi_{lab} vs. vertex x ", 200, -100., 100., 180, 0., 180.); list.Add(phiVsX);
	
	TH1F* firstLayerX = new TH1F("firstLayerX", "simulated x of the 1. layer", 1000, -100., 100.); list.Add(firstLayerX);
	TH1F* firstLayerY = new TH1F("firstLayerY", "simulated y of the 1. layer", 1000, -100., 100.); list.Add(firstLayerY);
	TH1F* firstLayerZ = new TH1F("firstLayerZ", "simulated z of the 1. layer", 1000, -100., 100.); list.Add(firstLayerZ);
	TH1F* secondLayerX = new TH1F("secondLayerX", "simulated x of the 2. layer", 1000, -100., 100.); list.Add(secondLayerX);
	TH1F* secondLayerY = new TH1F("secondLayerY", "simulated y of the 2. layer", 1000, -100., 100.); list.Add(secondLayerY);
	TH1F* secondLayerZ = new TH1F("secondLayerZ", "simulated z of the 2. layer", 1000, -100., 100.); list.Add(secondLayerZ);
	TH2F* posSimVsRecZ1 = new TH2F("posSimVsRecZ1", " simulated (y) vs reconstructed (x) z for 1. layer", 1000, -100., 100., 1000, -100., 100); list.Add(posSimVsRecZ1);
	TH2F* posSimVsRecX1 = new TH2F("posSimVsRecX1", " simulated (y) vs reconstructed (x) x for 1. layer", 1000, -100., 100., 1000, -100., 100); list.Add(posSimVsRecX1);
	TH2F* posSimVsRecZ2 = new TH2F("posSimVsRecZ2", " simulated (y) vs reconstructed (x) z for 2. layer", 1000, -100., 100., 1000, -100., 100); list.Add(posSimVsRecZ2);
	TH2F* posSimVsRecX2 = new TH2F("posSimVsRecX2", " simulated (y) vs reconstructed (x) x for 2. layer", 1000, -100., 100., 1000, -100., 100); list.Add(posSimVsRecX2);
	TH1F* posSimMinRecX1 = new TH1F("posSimMinRecX1", "simulated - reconstructed x for 1. layer", 1000, -50., 50.); list.Add(posSimMinRecX1);
	TH1F* posSimMinRecZ1 = new TH1F("posSimMinRecZ1", "simulated - reconstructed z for 1. layer", 1000, -50., 50.); list.Add(posSimMinRecZ1);
	TH1F* posSimMinRecX2 = new TH1F("posSimMinRecX2", "simulated - reconstructed x for 2. layer", 1000, -50., 50.); list.Add(posSimMinRecX2);
	TH1F* posSimMinRecZ2 = new TH1F("posSimMinRecZ2", "simulated - reconstructed z for 2. layer", 1000, -50., 50.); list.Add(posSimMinRecZ2);
	TH2F* posSimVsSimMinRecX1 = new TH2F("posSimVsSimMinRecX1", "rec Vs sim - rec x for 1. layer", 1000, -50., 50.,1000, -50., 50.);list.Add(posSimVsSimMinRecX1);
	TH2F* posSimVsSimMinRecZ1 = new TH2F("posSimVsSimMinRecZ1", "rec Vs sim - rec z for 1. layer", 1000, -50., 50.,1000, -50., 50.);list.Add(posSimVsSimMinRecZ1);
	TH2F* posSimVsSimMinRecX2 = new TH2F("posSimVsSimMinRecX2", "rec Vs sim - rec x for 2. layer", 1000, -50., 50.,1000, -50., 50.);list.Add(posSimVsSimMinRecX2);
	TH2F* posSimVsSimMinRecZ2 = new TH2F("posSimVsSimMinRecZ2", "rec Vs sim - rec z for 2. layer", 1000, -50., 50.,1000, -50., 50.);list.Add(posSimVsSimMinRecZ2);
	
	TH2F* dE1VsdE2EpadHit = new TH2F("dE1VsdE2EpadHit", "energy loss second layer vs. energy loss first layer for Epad>0", 1000, 0, 5000, 2000, 0, 10000); list.Add(dE1VsdE2EpadHit);
	
	// ******************************** my histograms end *************************************
	
	TH2F* eRecErrVsESim = new TH2F("eRecErrVsESim", "error of reconstructed energy vs. simulated energy of recoil", 1000, 0., 25000., 1000, -5000., 5000.); list.Add(eRecErrVsESim);
	TH2F* thetaErrorVsZ = new TH2F("thetaErrorVsZ", "Error in #vartheta_{lab} reconstruction vs. simulated z-position;z [mm];#Delta#vartheta_{lab} [^{o}]", 200, -100, 100, 100, -15, 15); list.Add(thetaErrorVsZ);
	TH2F* thetaErrorVsTheta = new TH2F("thetaErrorVsTheta", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}]", 180, 0, 180, 100, -15, 15); list.Add(thetaErrorVsTheta);
	TH2F* thetaErrorVsThetaEpadCut = new TH2F("thetaErrorVsThetaEpadCut", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}] with Epad>0", 180, 0, 180, 100, -15, 15); list.Add(thetaErrorVsThetaEpadCut);
	TH2F* zReactionEnergy = new TH2F("zReactionEnergy", "z position of reaction vs. Beam energy (rec.)", 200, -100, 100, 1000, 0, 1.1*beamEnergy); list.Add(zReactionEnergy);
	TH2F* excEnProtonVsTheta = new TH2F("excEnProtonVsTheta", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsTheta);
	TH2F* excEnProtonVsPhi = new TH2F("excEnProtonVsPhi", "Excitation Energy Spectrum from reconstructed Protons;#varphi_{lab}[^{o}];E_{exc} [keV]", 360, -180., 180., 5000, -20000, 20000); list.Add(excEnProtonVsPhi);
	TH2F* excEnProtonVsZ = new TH2F("excEnProtonVsZ", "Excitation Energy Spectrum from reconstructed Protons;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); list.Add(excEnProtonVsZ);
	TH2F* excEnProtonVsThetaGS = new TH2F("excEnProtonVsThetaGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsThetaGS);
	TH2F* excEnProtonVsZGS = new TH2F("excEnProtonVsZGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); list.Add(excEnProtonVsZGS);
	TH2F* excEnProtonVsThetaCm = new TH2F("excEnProtonVsThetaCm", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{cm}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsThetaCm);
	TH2F* thetaVsZ = new TH2F("thetaVsZ","#vartheta_{lab} vs. vertex z", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZ);
	TH2F* thetaVsZSame = new TH2F("thetaVsZSame","#vartheta_{lab} vs. z, first and second layer both forward or both backward", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZSame);
	TH2F* thetaVsZCross = new TH2F("thetaVsZCross","#vartheta_{lab} vs. z, first and second layer over cross", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZCross);
	TH2F* phiVsZ = new TH2F("phiVsZ","#varphi_{lab} vs. vertex z", 200, -100., 100., 360, -180., 180.); list.Add(phiVsZ);
	TH2F* hitpattern = new TH2F("hitpattern","detector # of second layer vs. detector # of first layer", 2, -0.5, 1.5, 2, -0.5, 1.5); list.Add(hitpattern);
	TH2F* eCmVsZ = new TH2F("eCmVsZ","energy of cm-system vs. z;z [mm];e-cm [GeV]", 200, -100., 100., 2000, transferP->GetCmEnergy(0.)/1000., transferP->GetCmEnergy(beamEnergy)/1000.); list.Add(eCmVsZ);
	TH2F* betaCmVsZ = new TH2F("betaCmVsZ","#beta of cm-system vs. z", 200, -100., 100., 2000, 0., 0.2); list.Add(betaCmVsZ);
	TH2F* stripPattern = new TH2F("stripPattern","Parallel strip # (#varphi) vs. perpendicular strip # (#vartheta)", 2*nStripsY, 0., 2.*nStripsY, 4*nStripsX, 0., 4.*nStripsX); list.Add(stripPattern);
	TH2F* recBeamEnergyErrVsZ = new TH2F("recBeamEnergyErrVsZ","Error in reconstructed beam energy vs. z", 200, -100., 100., 1000, -50., 50.); list.Add(recBeamEnergyErrVsZ);
	TH2F* thetaCmVsThetaLab = new TH2F("thetaCmVsThetaLab", "#vartheta_{cm} vs. #vartheta_{lab};#vartheta_{cm} [^{o}];#vartheta_{lab} [^{o}]", 180,0.,180., 180,0.,180.); list.Add(thetaCmVsThetaLab);
	TH2F* zErrorVsThetaSim = new TH2F("zErrorVsThetaSim", "Error between reconstructed and true z-position vs. simulated #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); list.Add(zErrorVsThetaSim);
	TH2F* zErrorVsThetaRec = new TH2F("zErrorVsThetaRec", "Error between reconstructed and true z-position vs. recon structed #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); list.Add(zErrorVsThetaRec);
	TH2F* zErrorVsthetaError = new TH2F("zErrorVstheaError", "z Erorr vs error in theta", 200, -10., 10., 1000, -5., 5.); list.Add(zErrorVsthetaError);
	TH2F* elossVsTheta = new TH2F("elossVsTheta", "reconstructed energy loss vs. #vartheta;#vartheta_lab [^{o}];energy loss [keV]", 360, -180., 180., 10000, -100000., 100000.); list.Add(elossVsTheta);
	TH2F* elossVsPhi = new TH2F("elossVsPhi", "reconstructed energy loss vs. #varphi;#varphi_lab [^{o}];energy loss [keV]", 360, -180., 180., 1000, -10000., 10000.); list.Add(elossVsPhi);
	TH2F* excEnElossVsTheta = new TH2F("excEnElossVsTheta", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV]", 190, -10., 180., 1000,-1500, 9500); list.Add(excEnElossVsTheta);
	TH2F* excEnElossVsThetaPhiGated = new TH2F("excEnElossVsThetaPhiGated", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV] Phi gated", 190, -10., 180., 1000,-1500, 9500); list.Add(excEnElossVsThetaPhiGated);
	TH2F* excEnElossVsThetaEpadCut = new TH2F("excEnElossVsThetaEpadCut", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV] with Epad>0", 180, 0., 180., 5000, -2000, 20000); list.Add(excEnElossVsThetaEpadCut);
	UInt_t nofLevels = trGen->GetMaximum("reaction")+1;
	if(nofLevels < 1) nofLevels = 1;
	if(nofLevels > 10) nofLevels = 10;
	std::vector<TH2F*> excEnElossVsThetaLevel(nofLevels);
	for(size_t r = 0; r < nofLevels; ++r) {
		excEnElossVsThetaLevel[r] = new TH2F(Form("excEnElossVsThetaLevel_%d",static_cast<int>(r)), Form("Excitation Energy Spectrum from reconstructed energy loss for level %d;#vartheta_{lab}[^{o}];E_{exc} [keV]",static_cast<int>(r)), 180, 0., 180., 5000, -20000, 20000); list.Add(excEnElossVsThetaLevel[r]);
	}
	
	// ******************************** Miniball Histograms *************************************
	
	TH1F* betaHisto = new TH1F ("betaHisto", "Beta", 1000, 0, .5); list.Add(betaHisto);	
	TH1F* gammaEnergyDet = new TH1F ("gammaEnergyDet", "gamma core energy (keV) detected in MB",1000,-1000,4000); list.Add(gammaEnergyDet);	
	TH1F* gammaTimeDet = new TH1F ("gammaTimeDet", "gamma time (ns) detected in MB",1000,-100,100); list.Add(gammaTimeDet);		
	//TH1F* noComptSpec = new TH1F("noComptSpec", "Miniball gamma sum spectrum of all clusters and crystals no compton ", 5000, 0, 15000); list.Add(noComptSpec);
    TH1F* gatedSpec = new TH1F("gatedSpec", "Miniball gated gamma spectrum", 5000, 0, 2500); list.Add(gatedSpec);
    TH1F* gatedSpecShift = new TH1F("gatedSpecShift", "Miniball gated gamma spectrum post shift (doppler corrected)", 5000, 0, 2500); list.Add(gatedSpecShift);    
    TH2F* gammaVsExc = new TH2F("gammaVsExc", "Miniball gamma sum energy vs excitaiton energy", 1000, 0., 10000., 4000, 0, 3000); list.Add(gammaVsExc);
    TH2F* gammaCluID = new TH2F("gammaCluID", "Miniball gamma sum energy vs Cluster ID", 1000, 0., 3000., 1000, 0, 40); list.Add(gammaCluID);
    TH2F* gammaCryID = new TH2F("gammaCryID", "Miniball gamma sum energy vs Crystal ID", 1000, 0., 3000.,1000, 0, 40); list.Add(gammaCryID);
	TH1F* preShiftSum = new TH1F("preShiftSum", "Miniball gamma sum spectrum of all clusters and crystals ", 5000, 0, 15000); list.Add(preShiftSum);	
	TH2F* realThetavsEnergy = new TH2F("realThetavsEnergy", "Miniball gamma sum energy vs real Theta", 180, 0.,180,4000, 0, 2600); list.Add(realThetavsEnergy);
	TH2F* realCorrTheta = new TH2F("realCorrTheta", "Miniball corrected gamma energy vs Real Theta", 180, 0.,180,4000, 0, 11000); list.Add(realCorrTheta);
    TH2F* gammaSimThetavsEnergy = new TH2F("gammaSimThetavsEnergy", "Miniball gamma energy vs Simulated Gamma Theta", 20, 0., 3.5,4000, 0, 2600); list.Add(gammaSimThetavsEnergy);
    TH1F* cluCorr = new TH1F("cluCorr", "Cluster Theta Doppler correction ", 5000, 0, 2500); list.Add(cluCorr);
    TH1F* realCorr = new TH1F("realCorr", "Real Theta Doppler correction ", 5000, 0, 15000); list.Add(realCorr);
    TH2F* corrVsDist = new TH2F("corrVsDist", "Correction amount as a function of distance", 1000, -100, 100, 1000, 0, 200);list.Add(corrVsDist);
    TH2F* realEnergyVsZ = new TH2F("realEnergyVsZ", "Real Doppler Shifted Energy VS Z", 1000, -100, 100, 1000, 0, 2500);list.Add(realEnergyVsZ);
    TH1F* idealCorr = new TH1F("idealCorr", "Ideal Doppler correction ", 5000, 0, 15000); list.Add(idealCorr);
    TH1F* recBeamEnergy = new TH1F ("recBeamEnergy", "Reconstructed Beam Energy", 1000, 0, 10000); list.Add(recBeamEnergy);
    TH2F* recBeamVsRealBeam = new TH2F ("recBeamVsRealBeam", "Reconstructed Beam Energy vs Real Beam Energy", 1000, 0, 1000, 1000, 0, 1000); list.Add(recBeamVsRealBeam);
    TH2F* recErrVsZ= new TH2F("recErrVsZ", "recoil error vs z", 400, -400, 400, 400, -1000, 1000);list.Add(recErrVsZ);
    TH2F* ejcErrVsZ= new TH2F("ejcErrVsZ", "ejectile error vs z", 400, -400, 400, 400, -1000, 100);list.Add(ejcErrVsZ);
    TH2F* ejcVsThetaRec= new TH2F("ejcVsThetaRec", "energy vs theta of ejectile", 180, 0, 180, 500, 0, 2500);list.Add(ejcVsThetaRec);
    TH2F* idealEnergyVsZ = new TH2F("idealEnergyVsZ", "Ideal Doppler Shifted Energy VS Z", 1000, -100, 100, 500, 0, 2500);list.Add(idealEnergyVsZ);
    TH2F* ejectileRecVsSim = new TH2F("ejectileRecVsSim", "ejectile energy in MeV sim (y) vs rec (x)",2000,0,10000,2000,0,10000);list.Add(ejectileRecVsSim);
    TH2F* gammaESimVsRec = new TH2F("gammaESimVsRec", "gamma sum energy sim (y) vs rec (x)",500,0,2500,500,0,2500);list.Add(gammaESimVsRec);
    TH2F* gammaErealSimVsRec = new TH2F("gammaErealSimVsRec", "gamma sum energy doppler corrected sim (y) vs rec (x) [keV]",500,0,2500,500,0,2500);list.Add(gammaErealSimVsRec);
    TH1F* coreEnergy = new TH1F("coreEnergy", "detected core energy", 5000, 0, 5000); list.Add(coreEnergy);
    
    TH1F* clu0 = new TH1F("clu0", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu0);
    TH1F* clu1 = new TH1F("clu1", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu1);
    TH1F* clu2 = new TH1F("clu2", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu2);
    TH1F* clu3 = new TH1F("clu3", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu3);
    TH1F* clu4 = new TH1F("clu4", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu4);
    TH1F* clu5 = new TH1F("clu5", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu5);
    TH1F* clu6 = new TH1F("clu6", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu6);
    TH1F* clu7 = new TH1F("clu7", "Reconstructed gamma sum energy no doppler", 1000, 0, 2500); list.Add(clu7);

    TH2F* etclu0 = new TH2F("etclu0", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu0);
    TH2F* etclu1 = new TH2F("etclu1", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu1);
    TH2F* etclu2 = new TH2F("etclu2", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu2);
    TH2F* etclu3 = new TH2F("etclu3", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu3);
    TH2F* etclu4 = new TH2F("etclu4", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu4);
    TH2F* etclu5 = new TH2F("etclu5", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,25000); list.Add(etclu5);
    TH2F* etclu6 = new TH2F("etclu6", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu6);
    TH2F* etclu7 = new TH2F("etclu7", "reconstructed gamma theta (x) vs doppler corrected energy (y)",180,0,180,500,0,2500); list.Add(etclu7);
    
    TH2F* etforw = new TH2F("etforw", "reconstructed gamma theta (x) vs doppler corrected energy (y) for z > 0", 180, 0, 180, 500, 0, 2500); list.Add(etforw);
    TH2F* etback = new TH2F("etback", "reconstructed gamma theta (x) vs doppler corrected energy (y) for z < 0", 180, 0, 180, 500, 0, 2500); list.Add(etback);
    TH2F* eiforw = new TH2F("eiforw", "ideal gamma theta (x) vs ideal doppler corrected energy (y) for z > 0", 180, 0, 180, 500, 0, 2500); list.Add(eiforw);
    TH2F* eiback = new TH2F("eiback", "ideal gamma theta (x) vs ideal doppler corrected energy (y) for z < 0", 180, 0, 180, 500, 0, 2500); list.Add(eiback);
    TH2F* thetaIdVsSim = new TH2F("thetaIdVsSim", "ideal #vartheta_{#gamma} (y) vs. simulated #vartheta_{#gamma} (x)", 180, 0, 180, 180, 0, 180); list.Add(thetaIdVsSim);
    TH2F* thetaReVsSim = new TH2F("thetaReVsSim", "real #vartheta_{#gamma} (y) vs. simulated #vartheta_{#gamma} (x)", 180, 0, 180, 180, 0, 180); list.Add(thetaReVsSim);
    TH2F* thetaReVsId  = new TH2F("thetaReVsId" , "real #vartheta_{#gamma} (y) vs. ideal #vartheta_{#gamma} (x)", 180, 0, 180, 180, 0, 180); list.Add(thetaReVsId);
    
	// ******************************** End of Miniball Histograms *************************************

	Particle part; 

	// create and save energy vs. theta-lab splines for reaction at front/middle/back of target
	//std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy/sett->GetProjectileA()<<"/";
	std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy<<"/";
	//transferP->SetEBeam(beamEnergy/sett->GetProjectileA());
	transferP->SetEBeam(beamEnergy);
	TSpline3* front = transferP->Evslab(0., 180., 1.);
	front->Write("RecoilEVsThetaLabFront");
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<"/";
	std::cout<<energyInTarget->Eval(targetThickness/2.)/1000.<<"/";
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(energyInTarget->Eval(targetThickness/2.)/1000.);
	TSpline3* middle = transferP->Evslab(0., 180., 1.);
	middle->Write("RecoilEVsThetaLabMiddle");
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<std::endl;
	std::cout<<energyInTarget->Eval(targetThickness)/1000.<<std::endl;
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(energyInTarget->Eval(targetThickness)/1000.);
	
	TSpline3* back = transferP->Evslab(0., 180., 1.);
	back->Write("RecoilEVsThetaLabBack");

	/************************************************************************************************************
	 * loop over all entries
	 ************************************************************************************************************/

     double dE1Eloss=0;
     double dE2Eloss=0;
     double dE1ElossRange=0;
     double dE2ElossRange=0;
     double dE1MeasMinRec;
     double dE2MeasMinRec;

	Long64_t nEntries = tr->GetEntries();
	if(maxEntries > 0 && maxEntries < nEntries) nEntries = maxEntries;

	cout<<"Nr of Events "<<nEntries<<endl;

	if(verbose) {
		std::cout<<"***********************************************************************************************"<<std::endl
			<<"*    Row   * Instance * FBarrelDe * SecondFBa * FBarrelEr * FBarrelDe * SecondFBa * FBarrelEr *"<<std::endl
			<<"************************************************************************************************"<<std::endl;
	}

	for(Long64_t i = 0; i < nEntries; ++i) {
		if(verbose) cout <<"Loop over entry Nr "<<i<<endl; 
		hit->Clear();
		ParticleBranch->clear();
		tr->GetEntry(i);
		trGen->GetEntry(i);
		Int_t silicon_mult_first = firstDeltaE[0]->size()+ firstDeltaE[1]->size();				
		Int_t silicon_mult_second = secondDeltaE[0]->size()+ secondDeltaE[1]->size();
		TVector3 firstposition;
		TVector3 secondposition;
		TVector3 firstLayerPosSim;
		TVector3 secondLayerPosSim;

		if(verbose) {
			//std::cout<<"*  what is setw(8)?: * "<<setw(8)<<i<<" *";
			size_t d;
			for(d = 0; d < firstDeltaE[0]->size() || d < secondDeltaE[0]->size(); ++d) {
				std::cout<<" "<<setw(8)<<d<<" *";
				if(d < firstDeltaE[0]->size()) std::cout<<" "<<setw(8)<<firstDeltaE[0]->at(d).GetID()<<" *";
				else                           std::cout<<"          *";
				if(d < secondDeltaE[0]->size()) std::cout<<" "<<setw(8)<<secondDeltaE[0]->at(d).GetID()<<" *";
				else                            std::cout<<"          *";
				if(d < pad[0]->size()) std::cout<<" "<<setw(8)<<pad[0]->at(d).GetID()<<" *";
				else                   std::cout<<"          *";
				if(d < firstDeltaE[0]->size() && firstDeltaE[0]->at(d).GetStripEnergy().size() > 0) std::cout<<" "<<setw(8)<<firstDeltaE[0]->at(d).GetStripEnergy()[0]<<" *";
				else                                                                                std::cout<<"          *";
				if(d < secondDeltaE[0]->size() && secondDeltaE[0]->at(d).GetStripEnergy().size() > 0) std::cout<<" "<<setw(8)<<secondDeltaE[0]->at(d).GetStripEnergy()[0]<<" *";
				else                                                                                  std::cout<<"          *";
				if(d < pad[0]->size()) std::cout<<" "<<setw(8)<<pad[0]->at(d).GetEdet()<<" *"<<std::endl;
				else                   std::cout<<"          *"<<std::endl;
				//std::cout<<" simulated forward x: "<<firstDeltaE[0]->at(d).GetPosGlobalX()[0]<<"	y: "<<firstDeltaE[0]->at(d).GetPosGlobalY()[0]<<"	z: "<<firstDeltaE[0]->at(d).GetPosGlobalZ()[0]<<std::endl;
			}
			if(d == 0) std::cout<<std::endl;
		}

		Int_t index_first = 0;
		Int_t index_second = 0;

		// exactly one hit in any first layer and one hit in any second layer ?
		//TODO: take into account multiple hits

		if(silicon_mult_first > 1 || silicon_mult_second > 1) {
			std::cout<<"Warning: Multiple hits in Silicon Tracker! "<<std::endl
			         <<"First layer:  "<<silicon_mult_first<<" ( ";
			for(auto dir : firstDeltaE) {
				std::cout<<dir->size()<<" ";
			}
			std::cout<<")  Second layer:  "<<silicon_mult_second<<" ( ";
			for(auto dir : secondDeltaE) {
				std::cout<<dir->size()<<" ";
			}
			std::cout<<")"<<std::endl;
		}  
		
		
		    if(silicon_mult_first == 1 && (silicon_mult_second == 1 || isSolid)) 
		    { // because solid target is not a tracker we do not care the mult
				
			//if(isSolid && silicon_mult_second < 2) { // begin of second layer mult<2 loop for solid target, comment in this for the solid target 

			if(firstDeltaE[0]->size() == 1 ) {
				hit->SetFirstDeltaE(firstDeltaE[0]->at(0), kForward);
				index_first = 0;
				firstLayerX->Fill(firstDeltaE[0]->at(0).GetPosGlobalX()[0]);
				firstLayerY->Fill(firstDeltaE[0]->at(0).GetPosGlobalY()[0]);
				firstLayerZ->Fill(firstDeltaE[0]->at(0).GetPosGlobalZ()[0]);
				
				firstLayerPosSim.SetXYZ(firstDeltaE[0]->at(0).GetPosGlobalX()[0],firstDeltaE[0]->at(0).GetPosGlobalY()[0],firstDeltaE[0]->at(0).GetPosGlobalZ()[0]);
				
			} else {
				hit->SetFirstDeltaE(firstDeltaE[1]->at(0), kBackward);
				index_first = 1;
				firstLayerX->Fill(firstDeltaE[1]->at(0).GetPosGlobalX()[0]);
				firstLayerY->Fill(firstDeltaE[1]->at(0).GetPosGlobalY()[0]);
				firstLayerZ->Fill(firstDeltaE[1]->at(0).GetPosGlobalZ()[0]);
				
				firstLayerPosSim.SetXYZ(firstDeltaE[1]->at(0).GetPosGlobalX()[0],firstDeltaE[1]->at(0).GetPosGlobalY()[0],firstDeltaE[1]->at(0).GetPosGlobalZ()[0]);
			}

			if(secondDeltaE[0]->size() == 1 ) {
				hit->SetSecondDeltaE(secondDeltaE[0]->at(0), kForward);
				index_second = 0;
				secondLayerX->Fill(secondDeltaE[0]->at(0).GetPosGlobalX()[0]);
				secondLayerY->Fill(secondDeltaE[0]->at(0).GetPosGlobalY()[0]);
				secondLayerZ->Fill(secondDeltaE[0]->at(0).GetPosGlobalZ()[0]);
				
				secondLayerPosSim.SetXYZ(secondDeltaE[0]->at(0).GetPosGlobalX()[0],secondDeltaE[0]->at(0).GetPosGlobalY()[0],secondDeltaE[0]->at(0).GetPosGlobalZ()[0]);
				
			} else if(secondDeltaE[1]->size() == 1) {
				hit->SetSecondDeltaE(secondDeltaE[1]->at(0), kBackward);
				index_second = 1;
				secondLayerX->Fill(secondDeltaE[1]->at(0).GetPosGlobalX()[0]);
				secondLayerY->Fill(secondDeltaE[1]->at(0).GetPosGlobalY()[0]);
				secondLayerZ->Fill(secondDeltaE[1]->at(0).GetPosGlobalZ()[0]);
				
				secondLayerPosSim.SetXYZ(secondDeltaE[1]->at(0).GetPosGlobalX()[0],secondDeltaE[1]->at(0).GetPosGlobalY()[0],secondDeltaE[1]->at(0).GetPosGlobalZ()[0]);
				
			}

			if(silicon_mult_second == 1) {
				if(pad[index_second]->size() == 1 ) {
					hit->SetPad(pad[index_second]->at(0));
				}

				if(verbose) {
					std::cout<<"Using pad "<<index_second<<" with "<<pad[index_second]->size()<<" detectors"<<std::endl;
					for(int p = 0; p < 2; ++p) {
						for(size_t d = 0; d < pad[p]->size(); ++d) {
							std::cout<<p<<": pad "<<pad[p]->at(d).GetID()<<" = "<<pad[p]->at(d).GetEdet()<<" keV / "<<pad[p]->at(d).GetRear()<<" keV"<<std::endl;
						}
					}
				}
			}
						
			//get position of hit in first layer
			firstposition = hit->FirstPosition(doubleSidedFirstLayer, !dontSmear); 
			firstposition.SetX(firstposition.X() + 13.0/3);//- 13.0/3 wrong direction

			// get position of hit in second layer
			if(!isSolid) {secondposition = hit->SecondPosition(!dontSmear);
				          secondposition.SetX(secondposition.X() + 13.0); // - 13.0 is wrong direction
					  }
			
			else        secondposition.SetXYZ(13.0, 0., 0.); // original --> secondposition.SetXYZ(0., 0., 0.);

			part.Clear();

			// vector between two hits in Siliocn Tracker
			if(isSolid) part.SetPosition(firstposition);
			else        part.SetPosition(secondposition - firstposition); 
			if(verbose) {
				cout<<"Position to first hit: "<< firstposition.X()<<"  "<<firstposition.Y()<<"   "<<firstposition.Z()<<endl;
				cout<<"Position to second hit: "<< secondposition.X()<<"  "<<secondposition.Y()<<"   "<<secondposition.Z()<<endl;
				cout<<"Position of relative vector: "<< part.GetPosition().X()<<"  "<<part.GetPosition().Y()<<"  "<<part.GetPosition().Z()<<endl;
			}

			// reaction angles
			recoilThetaSim = recoilThetaSim*180./TMath::Pi();
			recoilThetaRec = part.GetPosition().Theta()*180./TMath::Pi(); 
			recoilPhiSim = recoilPhiSim*180./TMath::Pi();
			recoilPhiRec = part.GetPosition().Phi()*180./TMath::Pi(); 
			if(verbose) std::cout<<"reaction phi from position: "<<recoilPhiRec<<" - "<<recoilPhiSim<<" = "<<(recoilPhiRec - recoilPhiSim)<<std::endl;
			if(verbose) std::cout<<"reaction theta from position: "<<recoilThetaRec<<" - "<<recoilThetaSim<<" = "<<(recoilThetaRec - recoilThetaSim)<<std::endl;


			TVector3 vertex;                   //reconstructed vertex
			if(!isSolid) {
				//find the closest point between beam axis and vector of the two hits in the silicon tracker
				TVector3 r = part.GetPosition();  //relative vector from first hit to second hit
				TVector3 r2 = secondposition;     // vector to second hit
				double t = 0;                          //line parameter to calculate vertex; temp use only
				if((r*r - r.Z()*r.Z()) != 0 ) t = (r2*r - (r2.Z()*r.Z()))/(r*r - r.Z()*r.Z());
				vertex = r2 -( t*r); 
			} else {
				vertex.SetXYZ(0., 0., (targetForwardZ + targetBackwardZ)/2.); // middle of target
			}
			
			//vertex.SetXYZ(reactionXSim, reactionYSim, reactionZSim); ///****************** for testing!
			
			if(verbose) {
				cout <<"Vertex: "<< vertex.X() <<"  "<<vertex.Y() <<"   "<< vertex.Z()<<endl;
				cout <<"Z from simu: "<< (reactionZSim)<<endl;
			}
			//update particle information
			if(vertex.Z() > targetForwardZ) {
				if(verbose) std::cout<<"Correcting vertex z from "<<vertex.Z();
				vertex.SetZ(targetForwardZ);
				if(verbose) std::cout<<" to "<<vertex.Z()<<std::endl;
			}
			
			if(vertex.Z() < targetBackwardZ) {
				if(verbose) std::cout<<"Correcting vertex z from "<<vertex.Z();
				vertex.SetZ(targetBackwardZ);
				if(verbose) std::cout<<" to "<<vertex.Z()<<std::endl;
			}

			// target length at reaction
			double targetThickEvent;
			if(isSolid) targetThickEvent = targetThickness/2.;
			else        targetThickEvent = targetThickness * ( vertex.Z() - targetBackwardZ ) / targetLength; 
			if(verbose) cout<<"Target Thickness at reaction: "<<targetThickEvent<<" = "<<targetThickness<<" * ( "<<vertex.Z()<<" - "<<targetBackwardZ<<" ) / "<<targetLength<<endl;


			//calculate target thickness for reconstruction of          
			if(targetThickEvent > 0) beamEnergyRec = energyInTarget->Eval(targetThickEvent)/1000.;
			else                     beamEnergyRec = beamEnergy;

			// reconstruct energy of recoil
			recoilEnergyRecdE    =  hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose);
			recoilEnergyRecErest =  hit->GetPadEnergy();
			recoilEnergyRec = recoilEnergyRecdE + recoilEnergyRecErest;
			//recoilEnergyRec =  recoilEnergyRecErest;
			// if(verbose && index_first == 0 && index_second == 0)  std::cout<<" 1. layer energy: "<<hit->GetFirstDeltaEEnergy()<<" .. "<<endl;
			if(verbose) std::cout<<" 1. layer energy: "<<hit->GetFirstDeltaEEnergy()<<" 2. layer energy: "<<hit->GetSecondDeltaEEnergy()<<" pad energy: "<<hit->GetPadEnergy()<<" => dE: "<<recoilEnergyRecdE<<", Erest: "<<recoilEnergyRecErest<<" => Erec: "<<recoilEnergyRec<<std::endl;
			// reconstruct energy loss in gas and foil
			double sinTheta = TMath::Sin(recoilThetaRec/180.*TMath::Pi());
			double cosTheta = TMath::Cos(recoilThetaRec/180.*TMath::Pi());
			double tmpPhi = recoilPhiRec; while(tmpPhi > 45.) tmpPhi -= 90.; while(tmpPhi < -45.) tmpPhi += 90.;
			double cosPhi   = TMath::Cos(tmpPhi/180.*TMath::Pi());
			double recoilEnergyRecEloss;

			if(isSolid) {
				//for the solid target we only need to reconstruct the energy loss in the foil and the target (no chamber gas)
				double range = recoilFoilRange->Eval(recoilEnergyRec);
				//recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta*cosPhi));// original
				recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta)); // changed bei Dennis
				range = recoilTargetRange->Eval(recoilEnergyRecEloss);
				recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetThickness/2./TMath::Abs(cosTheta));
								
				// due to changing of the target foil from box to cylinder the ernergy loss is corrected by ommitting cosphi. Leila & Dennis
				
			} else {
				// need to reconstruct energy loss in chamber gas, foil, and target
				double range;
				if(verbose) std::cout<<"theta "<<recoilThetaRec<<", phi "<<recoilPhiRec<<" (sinTheta "<<sinTheta<<", cosTheta "<<cosTheta<<", cosPhi "<<cosPhi<<"): ";
				
				//*** energy loss through the pad ***
				if(recoilEnergyRecErest > 0. ) {					
					if(verbose) std::cout<<"\n\n *** energy loss through the pad *** "<<std::endl;
					
					if(verbose) std::cout<<"from pad "<<recoilEnergyRecErest<<" through "<<(padDistance - secondLayerDistance)/(sinTheta*cosPhi)<<" mm gas \n";
					range = recoilChamberGasRange->Eval(recoilEnergyRecErest);
					hrecoilERecErestElossRange1VsTheta->Fill(recoilThetaRec,range);
					if(verbose) std::cout<<"1. range from the pad Erest "<<range<<std::endl;					
					
					dE2ElossRange = recoilLayerRange->Eval(recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi)));
				    if(verbose) std::cout<<"5.a second layer thickness: "<<sett->GetSecondFBarrelDeltaESingleThickness()[0]<<" range of the second layer: "<<dE2ElossRange<<std::endl;
				    dE2Eloss = recoilLayerEnergy->Eval(dE2ElossRange + secondLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi));
				    dE2MeasMinRec = TMath::Abs(hit->GetSecondDeltaEEnergy(verbose) - dE2Eloss);
				    hdE2MeasMinRec->Fill(hit->GetSecondDeltaEEnergy(verbose) - dE2Eloss);
				    hdE2ElossRangeWoEpad0->Fill(dE2ElossRange);								
					hdE2Eloss->Fill(dE2Eloss);
					hdE2Measured->Fill(hit->GetSecondDeltaEEnergy(verbose));
					hdE2ElossVsMeasuredWoEpad0->Fill(dE2Eloss,hit->GetSecondDeltaEEnergy(verbose));		
					
					recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta*cosPhi)) + hit->GetSecondDeltaEEnergy(verbose);
					hrecoilERecErestElossVsThetaLeila->Fill(recoilThetaRec,recoilEnergyRecEloss);
					if(verbose) std::cout<<"2. energy loss from the pad + de2 "<<recoilEnergyRecEloss<<std::endl;
					range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);		
					if(verbose) std::cout<<"3. range from the pad Eloss "<<range<<std::endl;	
					hrecoilERecErestElossRange2VsTheta->Fill(recoilThetaRec,range);	
					hErestMeasured->Fill(hit->GetPadEnergy());				
					
				} else {
					range = recoilChamberGasRange->Eval(hit->GetSecondDeltaEEnergy(verbose));
					if(verbose) std::cout<<"\n range from dE2 for Erest=0 "<<range<<std::endl;
					dE2ElossRange = recoilLayerRange->Eval(range);
					dE2Eloss = recoilLayerEnergy->Eval(dE2ElossRange + secondLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range); // unshifted anticorr 
					hdE2ElossVsMeasured->Fill(dE2Eloss,hit->GetSecondDeltaEEnergy(verbose));
				}
				
				//*** energy loss through the second layer ***	
								
				//if(recoilEnergyRecErest == 0. )
				{dE1ElossRange = recoilLayerRange->Eval(recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi)));
				dE1Eloss = recoilLayerEnergy->Eval(dE1ElossRange + firstLayerThicknessMgCm2/(sinTheta*cosPhi)) - recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi));
				dE1MeasMinRec = TMath::Abs(hit->GetFirstDeltaEEnergy(verbose) - dE1Eloss);
				hdE1ElossRange->Fill(dE1ElossRange);
				hdE1Eloss->Fill(dE1Eloss);
				hdE1Measured->Fill(hit->GetFirstDeltaEEnergy(verbose));
				hdE1MeasMinRec->Fill(hit->GetFirstDeltaEEnergy(verbose) - dE1Eloss);
				hdE1ElossVsMeasured->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(verbose));
				if(recoilEnergyRecErest == 0.) hdE1ElossVsMeasuredEpad0->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(verbose));
				if(recoilEnergyRecErest > 0.) hdE1ElossVsMeasuredEpadWo0->Fill(dE1Eloss,hit->GetFirstDeltaEEnergy(verbose));
				}		
				
				if(verbose) std::cout<<"\n\n *** energy loss through the second layer *** "<<std::endl;
				if(verbose) std::cout<<" with 2. layer "<<hit->GetSecondDeltaEEnergy(verbose)<<" ("<<recoilEnergyRecEloss<<") through "<<(secondLayerDistance - firstLayerDistance)/(sinTheta*cosPhi)<<" mm gas \n";
				hrecoilERecdE2ElossRangeVsTheta->Fill(recoilThetaRec,range);
				recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi)) + hit->GetFirstDeltaEEnergy(verbose);
				hrecoilERecdE2ElossVsThetaLeila->Fill(recoilThetaRec,recoilEnergyRecEloss);			
				
				// for now the foil is box-shaped as well, so we can just continue the same way
				
				//*** energy loss through the first layer ***
				
				if(verbose) std::cout<<"\n\n *** energy loss through the first layer *** "<<std::endl;
				
				if(verbose) std::cout<<" with 1. layer "<<hit->GetFirstDeltaEEnergy(verbose)<<" ("<<recoilEnergyRecEloss<<") through "<<(firstLayerDistance - foilDistance)/(sinTheta*cosPhi)<<" mm gas \n";
				range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);
				hrecoilERecdE1ElossRangeVsTheta->Fill(recoilThetaRec,range);
				recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + (targetWidthMgCm2*(1 - cosPhi) + 100.*chamberGasMat->GetDensity()*(firstLayerDistance - foilDistance)/(sinTheta*cosPhi)));
				//effectiveLength(firstLayer - vertex)/sinTheta*cosPhi - targetRadii/sinTheta = a0/sinTheta*cosPhi - Rt/sinTheta 
				hrecoilERecdE1ElossVsThetaLeila->Fill(recoilThetaRec,recoilEnergyRecEloss);
				
				//*** energy loss through the foil ***
				if(verbose) std::cout<<"\n\n *** energy loss through the foil *** "<<std::endl;
				if(verbose) std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<foilThicknessMgCm2/(sinTheta*cosPhi)<<" mm foil \n";
				range = recoilFoilRange->Eval(recoilEnergyRecEloss);
				hrecoilERecFoilElossRangeVsTheta->Fill(recoilThetaRec,range);
				if(verbose) std::cout<<"8. range from the foil "<<range<<std::endl;
				recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta*cosPhi));// original 
				//recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta));// changed by Leila
				hrecoilERecFoilElossVsThetaLeila->Fill(recoilThetaRec,recoilEnergyRecEloss);
				if(verbose) std::cout<<"9. energy loss from the foil "<<recoilEnergyRecEloss<<std::endl;
				if(verbose) std::cout<<"\n recoilEnergyRec: "<<recoilEnergyRec<<", range: "<<range<<" + foilThicknessMgCm2/(sinTheta*cosPhi) "<<foilThicknessMgCm2/(sinTheta)<<" => recoilEnergyRecEloss: "<<recoilEnergyRecEloss<<std::endl;
				
				//*** energy loss through the target ***
				// for now assume that the "box" inside the foil is filled with target gas. Not box any more. It is a cylinder --> phi is ommitted!
				if(verbose) std::cout<<"\n\n *** energy loss through the target *** "<<std::endl;
				if(verbose) std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<targetWidthMgCm2/(sinTheta*cosPhi)<<" mm gas \n";
				range = recoilTargetRange->Eval(recoilEnergyRecEloss);
				hrecoilERecTargetElossRangeVsTheta->Fill(recoilThetaRec,range);
				recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetWidthMgCm2/(sinTheta));// changed by Leila
				hrecoilERecTargetElossVsThetaLeila->Fill(recoilThetaRec,recoilEnergyRecEloss);
				//std::cout<<" => "<<recoilEnergyRecEloss<<std::endl;
				if(verbose) std::cout<<" => "<<recoilEnergyRecEloss<<std::endl; 
			}

			// position has already been set above
			part.SetRecEnergy(recoilEnergyRec);
			part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
			part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 


			//////////////////////////
			// Q-value reconstruction
			//////////////////////////

			transferP->SetEBeam(beamEnergyRec);
			transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
			transferP->SetAngles(recoilThetaRec/180.*TMath::Pi(), 2, true);
			double excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose); 
                        
			double recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.; 

			///////////////////////
			// Fill some histograms
			///////////////////////

			reaction->Fill(reactionSim);
			hitpattern->Fill(index_first, index_second);
			originXY->Fill(vertex.X(), vertex.Y());
			originXYErr->Fill(vertex.X() - reactionXSim, vertex.Y() - reactionYSim);
			errorOrigin->Fill(vertex.Z(),  vertex.Z()-reactionZSim );
			errorThetaPhi->Fill(recoilThetaRec - recoilThetaSim, recoilPhiRec - recoilPhiSim);
			dE12VsPad->Fill(recoilEnergyRecErest, recoilEnergyRecdE );
			dE12VsE->Fill(recoilEnergyRec, recoilEnergyRecdE );
			dE1VsE->Fill(recoilEnergyRec, hit->GetFirstDeltaEEnergy(verbose));//(firstDeltaE[index_first]->at(0)).GetRear() );
			dE2VsE->Fill(recoilEnergyRec, hit->GetSecondDeltaEEnergy(verbose));//(secondDeltaE[index_second]->at(0)).GetRear() );
			dE1VsdE2->Fill(hit->GetFirstDeltaEEnergy(verbose), hit->GetSecondDeltaEEnergy(verbose));
			eVsTheta->Fill(recoilThetaRec, recoilEnergyRec);
			eVsZ->Fill(vertex.Z(), recoilEnergyRec);
			if(recoilThetaRec<50. && recoilThetaRec>30.) eVsZThetaCut->Fill(vertex.Z(), recoilEnergyRec);
			
			//if(recoilEnergyRecErest == 0.0)
			{
			dE2VsZ->Fill(secondposition.Z(), hit->GetSecondDeltaEEnergy(verbose));
			dE2VsY->Fill(secondposition.Y(), hit->GetSecondDeltaEEnergy(verbose));
			dE2VsX->Fill(secondposition.X(), hit->GetSecondDeltaEEnergy(verbose));
			dE2VsR->Fill(sqrt(secondposition.X()*secondposition.X()+secondposition.Y()*secondposition.Y()+secondposition.Z()*secondposition.Z()), hit->GetSecondDeltaEEnergy(verbose));
		    }
		    		    
		    ErecoilSimMindE12->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    ErecoilSimMindE12VsdE2->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)),hit->GetSecondDeltaEEnergy(verbose));
		    ErecoilSimMindE1VsdE2->Fill(recoilEnergySim-hit->GetFirstDeltaEEnergy(verbose),hit->GetSecondDeltaEEnergy(verbose));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE1VsdE2Epad0->Fill(recoilEnergySim-hit->GetFirstDeltaEEnergy(verbose),hit->GetSecondDeltaEEnergy(verbose));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE12VsdE2Epad0->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)),hit->GetSecondDeltaEEnergy(verbose));
		    ErecoilSimMindE12VsdE1->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)),hit->GetFirstDeltaEEnergy(verbose));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE12VsdE1Epad0->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)),hit->GetFirstDeltaEEnergy(verbose));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE12Epad0->Fill(recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));		    		    
		    ErecoilSimMindE12VsThetaLab->Fill(recoilThetaRec,recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    ErecoilSimMindE1VsThetaLab->Fill(recoilThetaRec,recoilEnergySim-hit->GetFirstDeltaEEnergy(verbose));
		    ErecoilSimMindE2VsThetaLab->Fill(recoilThetaRec,recoilEnergySim-hit->GetSecondDeltaEEnergy(verbose));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE12VsThetaLabEpad0->Fill(recoilThetaRec,recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    ErecoilSimMindE12VsdThetaLab->Fill(recoilThetaRec-recoilThetaSim,recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    if(recoilEnergyRecErest == 0.0) ErecoilSimMindE12VsdThetaLabEpad0->Fill(recoilThetaRec-recoilThetaSim,recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    		    
		    ErecoilSimVsThetaLab->Fill(recoilThetaRec,recoilEnergySim);
		    ErecoilSimVsZ->Fill(vertex.Z(),recoilEnergySim);
		    ErecoilSimMindE12VsZ->Fill(vertex.Z(),recoilEnergySim-(hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose)));
		    ErecoilSimVsRec->Fill(recoilEnergyRec,recoilEnergySim);
		    if(recoilEnergyRecErest == 0.0)  ErecoilSimVsRecEpad0->Fill(recoilEnergyRec,recoilEnergySim);
		    ErecoilSim->Fill(recoilEnergySim);
		    ErecoilRec->Fill(recoilEnergyRec);
		    recoil2posZVsX->Fill(secondposition.Z(),secondposition.X());
		    recoil1posZVsX->Fill(firstposition.Z(),firstposition.X());
		    recoilVertexZVsX->Fill(vertex.Z(),vertex.X());
		   
			eRecErrVsESim->Fill(recoilEnergySim, recoilEnergyRec - recoilEnergySim);
			thetaErrorVsZ->Fill(vertex.Z(), recoilThetaRec - recoilThetaSim);
			//if(hit->GetPadEnergy()>1.00) thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
			thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim); 
			thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
			if(recoilEnergyRecErest > 0.) {thetaErrorVsThetaEpadCut->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);}
			if(reactionEnergyBeamCM > 0.0) zReactionEnergy->Fill(vertex.Z(), beamEnergyRec);
			excEnProton->Fill(excEnergy);
			excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);
			excEnProtonVsPhi->Fill(recoilPhiRec, excEnergy);
			excEnProtonVsThetaCm->Fill(recoilThetaCmRec, excEnergy);
			excEnProtonVsZ->Fill(vertex.Z(), excEnergy);
			if(reactionSim == 0) {
				excEnProtonVsThetaGS->Fill(recoilThetaRec, excEnergy);
				excEnProtonVsZGS->Fill(vertex.Z(), excEnergy);
			}
			
			thetaVsZ->Fill(vertex.Z(), recoilThetaRec);
			thetaVsX->Fill(vertex.X(), recoilThetaRec);
			phiVsZ->Fill(vertex.Z(), recoilPhiRec);
			phiVsX->Fill(vertex.X(), recoilPhiRec);
			
			if(index_first == index_second) {
				eVsZSame->Fill(vertex.Z(), recoilEnergyRec);
				thetaVsZSame->Fill(vertex.Z(), recoilThetaRec);
			} else {
				eVsZCross->Fill(vertex.Z(), recoilEnergyRec);
				thetaVsZCross->Fill(vertex.Z(), recoilThetaRec);
			}
			phiVsZ->Fill(vertex.Z(), recoilPhiRec);
			phiErrorVsPhi->Fill(recoilPhiRec, recoilPhiRec - recoilPhiSim);
			if(firstDeltaE[index_first]->at(0).GetID() == 0) {
				if(index_first == 0) phiErrorVsPhiF0->Fill(recoilPhiRec, recoilPhiRec - recoilPhiSim);
				else                 phiErrorVsPhiB0->Fill(recoilPhiRec, recoilPhiRec - recoilPhiSim);
			}
			betaCmVsZ->Fill(vertex.Z(), transferP->GetBetacm());
			eCmVsZ->Fill(vertex.Z(), transferP->GetCmEnergy()/1000.);
			if(silicon_mult_second > 0) stripPattern->Fill(index_second*nStripsY + secondDeltaE[index_second]->at(0).GetStripNr()[0], secondDeltaE[index_second]->at(0).GetID()*nStripsX + secondDeltaE[index_second]->at(0).GetRingNr()[0]);
			recBeamEnergyErrVsZ->Fill(vertex.Z(), beamEnergyRec - reactionEnergyBeam);
			thetaCmVsThetaLab->Fill(recoilThetaRec, recoilThetaCmRec);
			zErrorVsthetaError->Fill(recoilThetaRec - recoilThetaSim, vertex.Z() - reactionZSim);
			elossVsTheta->Fill(recoilThetaRec, recoilEnergyRecEloss - recoilEnergyRec);
			elossVsPhi->Fill(recoilPhiRec, recoilEnergyRecEloss - recoilEnergyRec);
			
			hERec->Fill(recoilEnergyRec);
			hERecdE->Fill(recoilEnergyRecdE);
			hERecErest->Fill(recoilEnergyRecErest);
			hrecoilEnergyRecEloss->Fill(recoilEnergyRecEloss);
			hERecErestVsTheta->Fill(recoilThetaRec,recoilEnergyRecErest);
			hERecdEVsTheta->Fill(recoilThetaRec,recoilEnergyRecdE);
			
			hrecoilThetaRec->Fill(recoilThetaRec);
            hrecoilThetaRecDiffSim->Fill(recoilThetaRec - recoilThetaSim);
            hrecoilEnergyRecDiffSim->Fill(recoilEnergyRec - recoilEnergySim);        
            
			dE2VsdE2Pad->Fill(hit->GetSecondDeltaEEnergy(verbose)+hit->GetPadEnergy(),hit->GetSecondDeltaEEnergy(verbose));    
			EPadVsThetaLab->Fill(recoilThetaRec,hit->GetPadEnergy());
            EPadVsZ->Fill(vertex.Z(),hit->GetPadEnergy());
            if(recoilEnergyRecErest > 0.0) {dE2VsThetaLabEpadCut->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(verbose));}
            if(recoilEnergyRecErest > 0.0) {dE1VsdE2EpadHit->Fill(hit->GetFirstDeltaEEnergy(verbose),hit->GetSecondDeltaEEnergy(verbose));}
            if(recoilThetaRec>0.0 && recoilThetaRec<180.0) {dE2VsEPadThetaCut->Fill(hit->GetPadEnergy(),hit->GetSecondDeltaEEnergy(verbose));}
            dE1VsThetaLab->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(verbose));
            dE2VsThetaLab->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(verbose));
            dE12VsThetaLab->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(verbose)/(cosPhi)+hit->GetSecondDeltaEEnergy(verbose)/(cosPhi));
            dE12VsPhiLab->Fill(recoilPhiRec,hit->GetFirstDeltaEEnergy(verbose)/(firstposition.Mag()*cosPhi)+hit->GetSecondDeltaEEnergy(verbose)/(secondposition.Mag()*cosPhi));
            dE1VsPhiLab->Fill(recoilPhiRec,hit->GetFirstDeltaEEnergy(verbose)/(firstposition.Mag()*cosPhi));
            dE2VsPhiLab->Fill(recoilPhiRec,hit->GetSecondDeltaEEnergy(verbose)/(secondposition.Mag()*(sinTheta*cosPhi)));
            dE1EpadVsThetaLab->Fill(recoilThetaRec,hit->GetFirstDeltaEEnergy(verbose)+hit->GetPadEnergy());
            dE2EpadVsThetaLab->Fill(recoilThetaRec,hit->GetSecondDeltaEEnergy(verbose)+hit->GetPadEnergy());
            
            thetaVsZ2->Fill(secondposition.Z(),recoilThetaRec);
            thetaVsX2->Fill(secondposition.X(),recoilThetaRec);
            thetaVsZ1->Fill(firstposition.Z(),recoilThetaRec);
            thetaVsX1->Fill(firstposition.X(),recoilThetaRec);
            phiVsZ2->Fill(secondposition.Z(),recoilPhiRec);
            phiVsX2->Fill(secondposition.X(),recoilPhiRec);
            phiVsZ1->Fill(firstposition.Z(),recoilPhiRec);
            phiVsX1->Fill(firstposition.X(),recoilPhiRec);
                     
                     
            posSimVsRecX1->Fill(firstposition.X(),firstLayerPosSim.X());
            posSimVsRecZ1->Fill(firstposition.Z(),firstLayerPosSim.Z());
            posSimVsRecX2->Fill(secondposition.X(),secondLayerPosSim.X());
            posSimVsRecZ2->Fill(secondposition.Z(),secondLayerPosSim.Z());
            posSimMinRecX1->Fill(firstLayerPosSim.X()-firstposition.X());
            posSimMinRecZ1->Fill(firstLayerPosSim.Z()-firstposition.Z());
            posSimMinRecX2->Fill(secondLayerPosSim.X()-secondposition.X());
            posSimMinRecZ2->Fill(secondLayerPosSim.Z()-secondposition.Z());
            posSimVsSimMinRecX1->Fill(firstposition.X(),firstLayerPosSim.X()-firstposition.X());
            posSimVsSimMinRecZ1->Fill(firstposition.Z(),firstLayerPosSim.Z()-firstposition.Z());
            posSimVsSimMinRecX2->Fill(secondposition.X(),secondLayerPosSim.X()-secondposition.X());
            posSimVsSimMinRecZ2->Fill(secondposition.Z(),secondLayerPosSim.Z()-secondposition.Z());
                 
            // ************************* Q-value using the reconstructed energy loss **************************8

			// now we reconstruct the q-value using the reconstructed energy loss
			// position has already been set above
			part.SetRecEnergy(recoilEnergyRecEloss); // 1 --> projectile, 2 --> recoil, 3 --> ejectile
			part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
			part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 
			transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
			transferP->SetAngles(recoilThetaRec/180.*TMath::Pi(), 2, true);
			excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose); 
                       
			recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;
			
			double ejectileEnergyRec = transferP->GetTlab(3);//Changed from TLab to ELab
			
			excEnElossVsTheta->Fill(recoilThetaRec, excEnergy);
		    if(firstposition.X()>-5.0 && firstposition.X()<5.0) excEnElossVsThetaPhiGated->Fill(recoilThetaRec,excEnergy); 
			if(recoilEnergyRecErest > 0. ) {excEnElossVsThetaEpadCut->Fill(recoilThetaRec, excEnergy);}
			excEnProtonCorr->Fill(excEnergy);
			if(recoilEnergyRecErest > 0. ) {excEnProtonCorrEpadCut->Fill(excEnergy);}
			excEnProtonCorrVsX->Fill(vertex.X(),excEnergy);
			excEnProtonCorrVsY->Fill(vertex.Y(),excEnergy);
			excEnProtonCorrVsZ->Fill(vertex.Z(),excEnergy);
			excEnProtonCorrVsT->Fill(sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y()),excEnergy);
			excEnProtonCorrVsR->Fill(sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y()+vertex.Z()*vertex.Z()),excEnergy);
			if(dE1MeasMinRec<20.0 && dE2MeasMinRec<40.0) excEnProtonCorrdE1Sigma1->Fill(excEnergy);
			if(dE1MeasMinRec<40.0 && dE2MeasMinRec<80.0) excEnProtonCorrdE1Sigma2->Fill(excEnergy);
						
			if(0 <= reactionSim && reactionSim < nofLevels) excEnElossVsThetaLevel[reactionSim]->Fill(recoilThetaRec, excEnergy);
			//if(0 <= reactionSim && reactionSim < nofLevels-1) excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);//leila	
			
			
			///************************** Loop over Miniball **************************///

            for(auto cluster : *mBall) // loop over MB
            {
              int cluId = cluster.GetCluID();
              auto crystals = cluster.GetCrystal();              
              
              for(auto crystal : crystals)
              {
                double idealTheta, realTheta = 0.;
                double gamTime = crystal.GetTime();
                auto gamPos = cluster.GetPosition();
                //double gamNoCompt = cluster.GetnoComptE();

                double eLab = crystal.GetCore();//Core Energy keV
                coreEnergy->Fill(eLab);
                double massMeV = ejectile->GetMass();// mass in MeV
                double enMeV = massMeV + ejectileEnergyRec;//Using reconstructed energy beamEnergyRec in MeV
                double idealMeV = massMeV + ejectileEnergySim/1000.;//Using  true reaction energy
                double beta = sqrt(1 - pow((massMeV/enMeV), 2));//Calculation of relativistic beam velocity (real) reconstructed
                double idealBeta = sqrt(1 - pow((massMeV/idealMeV), 2));//Calculation of relativistic velocity (ideal)
                betaHisto->Fill(beta);
                
                gammaEnergyDet->Fill(eLab);
                gammaTimeDet->Fill(gamTime);
                
                if(verbose) std::cout<<" core energy keV: "<<eLab<<" rec beta: "<<beta<<" sim beta: "<<idealBeta<<" ejectile rest mass MeV: "<<massMeV<<" ejectile rec energy MeV: "<<ejectileEnergyRec<<" ejectile total energy MeV: "<<enMeV<<std::endl;

                double tigDist = 11.;//Disance from center of TI-STAR to TIGRESS (radius) (cm)
                double tigRes = .1;//Spatial Resolution of TIGRESS (cm)
                double thetaErr = TMath::ATan((tigRes/tigDist));//The max possible error in theta = arctan (tigRes/tigDist)
                
                idealTheta = TMath::ACos((gamPos[2]-reactionZSim)/TMath::Sqrt((gamPos[0]-reactionXSim)*(gamPos[0]-reactionXSim) + (gamPos[1]-reactionYSim)*(gamPos[1]-reactionYSim) + (gamPos[2]-reactionZSim)*(gamPos[2]-reactionZSim))); // in radian
                gamPos[2] = gamPos[2] - vertex.Z();
                realTheta = TMath::ACos(gamPos[2]/TMath::Sqrt(gamPos[0]*gamPos[0] + gamPos[1]*gamPos[1] + gamPos[2]*gamPos[2])); // in radian
                
                if(verbose) std::cout<<"idealTheta [deg]: "<<idealTheta/TMath::Pi()*180.<<" gamPosZ: "<<gamPos[2]<<" realTheta [deg]: "<<realTheta/TMath::Pi()*180.<<std::endl;
                
                std::random_device rd; // obtain a random number from hardware
                std::mt19937 eng(rd()); // seed the generator
                std::uniform_real_distribution<> distr(-1, 1); // define the range (-1, 1)

                double thetaSmear = thetaErr * distr(eng);//Multiply max theta error by double between -1 and 1 (+- random error)
                realTheta += thetaSmear;//Add theta error to our realTheta
                
                if(verbose) std::cout<<"idealTheta deg: "<<idealTheta/TMath::Pi()*180.<<" gamPosZ: "<<gamPos[2]<<" realTheta deg: "<<realTheta/TMath::Pi()*180.<<" real theta smeared out: "<<realTheta/TMath::Pi()*180.<<std::endl;
                                
                double realShift;//Real theta, rec energy
                double idealShift;//True shift

                realShift = (eLab/sqrt(1 - beta*beta)) * (1 - beta*TMath::Cos(realTheta));//Calculating the doppler shifted energies using real theta and ideal theta
                idealShift = (eLab/sqrt(1 - idealBeta*idealBeta)) * (1 - idealBeta*TMath::Cos(idealTheta));
                //gamNoCompt = (gamNoCompt/sqrt(1 - beta*beta)) * (1 - beta*TMath::Cos(realTheta));
                //noComptSpec->Fill(gamNoCompt);
                
                if (gamTime > 20) // in ns
                {
                        gatedSpec->Fill(crystal.GetCore());
                        gatedSpecShift->Fill(realShift);
                }
                
                gammaVsExc->Fill(excEnergy, crystal.GetCore());
                preShiftSum->Fill(crystal.GetCore());
                gammaCryID->Fill(crystal.GetCore(), (cluId*3+crystal.GetCryID()));
                gammaCluID->Fill(crystal.GetCore(), cluId);
                
                realThetavsEnergy->Fill(realTheta/TMath::Pi()*180., crystal.GetCore());                
                ejectileRecVsSim->Fill(ejectileEnergyRec, ejectileEnergySim);
                recErrVsZ->Fill(vertex.Z(), recoilEnergySim-recoilEnergyRec);
                ejcErrVsZ->Fill(vertex.Z(), ejectileEnergySim-ejectileEnergyRec);
                ejcVsThetaRec->Fill(vertex.Z(), ejectileEnergyRec);
                idealEnergyVsZ->Fill(vertex.Z(), idealShift);
                realEnergyVsZ->Fill(vertex.Z(), realShift);
                
                realCorr->Fill(realShift);
                idealCorr->Fill(idealShift);
                realCorrTheta->Fill(realTheta/TMath::Pi()*180.,realShift);

                corrVsDist->Fill(vertex.Z(), idealShift-eLab);
                
                 if (gamTime > 20){switch (cluId)
                {
                        case 0:
                                clu0->Fill(crystal.GetCore());
                                etclu0->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 1:
                                clu1->Fill(crystal.GetCore());
                                etclu1->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 2:
                                clu2->Fill(crystal.GetCore());
                                etclu2->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 3:
                                clu3->Fill(crystal.GetCore());
                                etclu3->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 4:
                                clu4->Fill(crystal.GetCore());
                                etclu4->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 5:
                                clu5->Fill(crystal.GetCore());
                                etclu5->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 6:
                                clu6->Fill(crystal.GetCore());
                                etclu6->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                        case 7:
                                clu7->Fill(crystal.GetCore());
                                etclu7->Fill(realTheta/TMath::Pi()*180., realShift);
                                break;
                } // loop over clusterID
				} // loop over gamma gate
								
	            if(reactionZSim > 0)
				{
					etforw->Fill(realTheta/TMath::Pi()*180., realShift);
					eiforw->Fill(idealTheta/TMath::Pi()*180., idealShift);
				}
				else
				{
					etback->Fill(realTheta/TMath::Pi()*180., realShift);
					eiback->Fill(idealTheta/TMath::Pi()*180., idealShift);
				}
				
				for(auto simTh : *gammaThetaSim) 
				{
					thetaIdVsSim->Fill(simTh/TMath::Pi()*180., idealTheta/TMath::Pi()*180.);
					thetaReVsSim->Fill(simTh/TMath::Pi()*180., realTheta/TMath::Pi()*180.);
					thetaReVsId->Fill(idealTheta/TMath::Pi()*180., realTheta/TMath::Pi()*180.);
				}
				
				for(auto simEgamma : *gammaEnergySim)
				{
					gammaESimVsRec->Fill(crystal.GetCore(),simEgamma);
					gammaErealSimVsRec->Fill(realShift,simEgamma);
				}
                             
			  } // loop over crystals             
			} // loop over clusters end of the MB loop 
			
	//}   //end of second layer mult<2 loop for solid target
			
	}   //end of mult = 1 events
		
		if(i%1000 == 0){
			cout<<setw(5)<<std::fixed<<setprecision(1)<<(100.*i)/nEntries<<setprecision(3)<<" % done\r"<<flush;
		}

} // end of loop over all events
	
	

	if(verbose) std::cout<<"***********************************************************************************************"<<std::endl;


	//////////////////
	// write results
	//////////////////
	outfile.cd();
	if(writeTree) rectr->Write("",TObject::kOverwrite);

	list.Sort();
	list.Write();

	cout<<"closing file ..."<<endl;
	infile.Close();
	outfile.Close();

	cout<<endl<<"Done!"<<endl;

	return 0;

}

