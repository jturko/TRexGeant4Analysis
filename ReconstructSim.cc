#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "CommandLineInterface.hh"
#include "TRexSettings.hh"
#include "HitSim.hh"
#include "ParticleMC.hh"
#include "Compound.hh"
#include "Reconstruction.hh"
#include "Kinematics.hh"
#include "Particle.hh"

using namespace TMath;
using namespace std;

ClassImp(ParticleMC);
ClassImp(Particle);

// * changes done by Leila for elastic scattering lines: 528-530, 596-597
// * for the beam energy 5-6 Mev/u (or 660 and 792 Mev total) the beam energy in the middle and after the target negative. But for the total beam energy of 1340 MeB the beam energy is only after the target negative. For 2640 MeV all are positive. for gas and foil target both. 
//
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

	tr->SetBranchAddress("FBarrelDeltaEMC",&firstDeltaE[0]);
	tr->SetBranchAddress("SecondFBarrelDeltaEMC",&secondDeltaE[0]);
	tr->SetBranchAddress("FBarrelErestMC",&pad[0]);
	tr->SetBranchAddress("BBarrelDeltaEMC",&firstDeltaE[1]);
	tr->SetBranchAddress("SecondBBarrelDeltaEMC",&secondDeltaE[1]);
	tr->SetBranchAddress("BBarrelErestMC",&pad[1]);

	double reactionEnergyBeam;
	trGen->SetBranchAddress("reactionEnergy", &reactionEnergyBeam);

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

	//set Compounds for now we always assume density of H2 = 0.180e-3 at 1 bar
	if(verbose) std::cout<<"creating compound \""<<sett->GetTargetMaterialName().c_str()<<"\""<<std::endl;
	Compound* targetMat = new Compound(sett->GetTargetMaterialName().c_str());
	if(!isSolid) targetMat->SetDensity(0.180e-3*sett->GetTargetPressure()/6.24151e+08); // internal geant4 units are such that 1 bar = 6.24151e+08 whatever it is MeV/mm3:-)))
	std::cout<<"set target material density to "<<0.180e-3*sett->GetTargetPressure()/6.24151e+08<<" = "<<targetMat->GetDensity()<<std::endl;
	if(verbose) std::cout<<"creating compound \"MY\""<<std::endl;
	Compound* foilMat = new Compound("MY");
	if(verbose) std::cout<<"creating compound \""<<sett->GetVacuumChamberGas().c_str()<<"\""<<std::endl;
	Compound* chamberGasMat = new Compound(sett->GetVacuumChamberGas().c_str());
	chamberGasMat->SetDensity(0.180e-3*sett->GetVacuumChamberGasPressure()/6.24151e+08); // internal geant4 units are such that 1 bar = 6.24151e+08 whatevers
	std::cout<<"set chamber gas density to "<<0.180e-3*sett->GetVacuumChamberGasPressure()/6.24151e+08<<" = "<<chamberGasMat->GetDensity()<<std::endl;// is 0.180e-3 correct? He das density is 0.164e-3 g/cm3. Leila!

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

        // Elastic Scattering reaction and kinematics class for Q-value = 0 and no outgoing articles *** Leila *** 
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
	double targetLength    = sett->GetTargetPhysicalLength();// original
	//double targetLength    = sett->GetTargetPhysicalLength()/2;// changed by Leila
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
	Reconstruction* recoilChamberGas = new Reconstruction(recoil, chamberGasMat);
	TSpline3* recoilTargetRange  = recoilTarget->Energy2Range(100., 0.1, !isSolid);
	TSpline3* recoilTargetEnergy = recoilTarget->Range2Energy(100., 0.1, !isSolid);
	TSpline3* recoilFoilRange  = recoilFoil->Energy2Range(100., 0.1, false);
	TSpline3* recoilFoilEnergy = recoilFoil->Range2Energy(100., 0.1, false);
	TSpline3* recoilChamberGasRange  = recoilChamberGas->Energy2Range(100., 0.1, true);
	TSpline3* recoilChamberGasEnergy = recoilChamberGas->Range2Energy(100., 0.1, true);

	//double foilDistance = sett->GetFBarrelDeltaESingleDistanceToBeam()[0] - 3.; // there seems to be a hard-coded distance of three mm between foil and detector ??????????????? check this out!!!!!!!!!!!!1 original
	
	double foilDistance = 5.0; // there seems to be a hard-coded distance of three mm between foil and detector ??????????????? check this out!!!!!!!!!!!!1 changed by Leila
	
	double firstLayerDistance = sett->GetFBarrelDeltaESingleDistanceToBeam()[0];
	double secondLayerDistance = sett->GetSecondFBarrelDeltaESingleDistanceToBeam()[0];
	double padDistance = sett->GetFBarrelErestSingleDistanceToBeam()[0];

	double targetWidthMgCm2 = foilDistance*100.*targetMat->GetDensity();// convert from mm to 10 um, and from there to mg/cm2 using density
	double foilThicknessMgCm2 = sett->GetFBarrelDeltaESingleFoilThickness()*100.*1.4;// convert from mm to 10 um, and from there to mg/cm2 using density 1.4 g/cm3. it is 5 mu from setting file
	double firstGasLayerThicknessMgCm2 = 300.*chamberGasMat->GetDensity(); // distance between foil and 1. layer hard-coded to 3 mm
	double secondGasLayerThicknessMgCm2 = (secondLayerDistance - firstLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 1. and 2. layer in 10 um and converted to mg/cm2 using density
	double thirdGasLayerThicknessMgCm2 = (padDistance - secondLayerDistance)*100.*chamberGasMat->GetDensity(); // distance 2. layer and pad in 10 um and converted to mg/cm2 using density

	std::cout<<"foil distance "<<foilDistance<<" mm => target width = "<<targetWidthMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"first layer distance "<<firstLayerDistance<<" mm and second layer: "<<secondLayerDistance<<" mm"<<std::endl;
	std::cout<<"foil thickness "<<sett->GetFBarrelDeltaESingleFoilThickness()<<" mm = "<<targetWidthMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"distance foil - 1. layer 3. mm => 1. gas layer thickness = "<<firstGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
	std::cout<<"distance 1. layer - 2. layer "<<secondLayerDistance - firstLayerDistance<<" mm => 2. gas layer thickness = "<<secondGasLayerThicknessMgCm2<<" mg/cm^2"<<std::endl;
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

	TH2F* dE12VsPad = new TH2F("dE12VsPad", "energy loss first+second layer vs. pad energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE12VsPad);
	TH2F* dE12VsE = new TH2F("dE12VsE", "energy loss first+second layer vs. total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE12VsE);
	TH2F* dE1VsE = new TH2F("dE1VsE", "energy loss first layer vs. total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE1VsE);
	TH2F* dE2VsE = new TH2F("dE2VsE", "energy loss second layer vs. total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE2VsE);
	TH2F* dE1VsdE2 = new TH2F("dE1VsdE2", "energy loss second layer vs. energy loss first layer", 200, 0, 10000, 200, 0, 10000); list.Add(dE1VsdE2);
	TH2F* eVsTheta = new TH2F("eVsTheta", "recoil energy vs. theta (lab)", 180, 0, 180, 200, 0, 25000); list.Add(eVsTheta);
	TH2F* eVsZ = new TH2F("eVsZ", "recoil energy vs. z", 200, -100., 100., 200, 0, 25000); list.Add(eVsZ);
	TH2F* eVsZSame = new TH2F("eVsZSame", "recoil energy vs. z, first and second layer both forward or both backward", 200, -100., 100., 200, 0, 25000); list.Add(eVsZSame);
	TH2F* eVsZCross = new TH2F("eVsZCross", "recoil energy vs. z, first and second layer over cross", 200, -100., 100., 200, 0, 25000); list.Add(eVsZCross);
	
	TH1F* hERecLeila = new TH1F("hERecLeila", "reconstructed energy of recoil leila", 1000, -5000., 35000.); list.Add(hERecLeila);
	TH1F* hrecoilEnergyRecElossLeila = new TH1F("hrecoilEnergyRecElossLeila", "reconstructed energy loss of recoil leila", 1000, -5000., 35000.); list.Add(hrecoilEnergyRecElossLeila);
	TH1F* hrecoilThetaRecLeila = new TH1F("hrecoilThetaRecLeila", "reconstructed theta of recoil leila", 1000, -60., 300.); list.Add(hrecoilThetaRecLeila);
	TH1F* hrecoilThetaRecDiffSimLeila = new TH1F("hrecoilThetaRecDiffSimLeila", "difference reconstructed and simulated theta of recoil leila", 1000, -60., 300.); list.Add(hrecoilThetaRecDiffSimLeila);
	TH1F* hrecoilEnergyRecDiffSimLeila = new TH1F("hrecoilEnergyRecDiffSimLeila", "difference reconstructed and simulated energy of recoil leila", 1000, -25000., 25000.); list.Add(hrecoilEnergyRecDiffSimLeila);
	TH1F* hpos1xLeila = new TH1F("hpos1xLeila", "x positoin on the first layer leila", 1000, -50., 50.); list.Add(hpos1xLeila);
	TH1F* hpos1yLeila = new TH1F("hpos1yLeila", "y positoin on the first layer leila", 1000, -50., 50.); list.Add(hpos1yLeila);
	TH1F* hpos1zLeila = new TH1F("hpos1zLeila", "z positoin on the first layer leila", 1000, -250., 250.); list.Add(hpos1zLeila);
	TH1F* hpos2xLeila = new TH1F("hpos2xLeila", "x positoin on the second layer leila", 1000, -100., 100.); list.Add(hpos2xLeila);
	TH1F* hpos2yLeila = new TH1F("hpos2yLeila", "y positoin on the second layer leila", 1000, -50., 50.); list.Add(hpos2yLeila);
	TH1F* hpos2zLeila = new TH1F("hpos2zLeila", "z positoin on the second layer leila", 1000, -250., 250.); list.Add(hpos2zLeila);
	TH2F* hpos1xyLeila = new TH2F("hpos1xyLeila", "x vs y positoin on the first layer leila", 1000, -50., 50., 1000, -50., 50.); list.Add(hpos1xyLeila);
	TH2F* hpos1ztLeila = new TH2F("hpos1ztLeila", "z vs t positoin on the first layer leila", 1000, -500., 500., 1000, -50., 50.); list.Add(hpos1ztLeila);
	TH2F* hpos2xyLeila = new TH2F("hpos2xyLeila", "x vs y positoin on the second layer leila", 1000, -100., 100., 1000, -100., 100.); list.Add(hpos2xyLeila);
	TH2F* hpos2ztLeila = new TH2F("hpos2ztLeila", "z vs t positoin on the second layer leila", 1000, -500., 500., 1000, -50., 150.); list.Add(hpos2ztLeila);
	TH2F* hexcEnProtonVsZLeila = new TH2F("hexcEnProtonVsZLeila", "recoil exc. energy vs z leila", 1000, -150., 150., 1000, -20000., 20000.); list.Add(hexcEnProtonVsZLeila);
	TH2F* hexcEnProtonVsXLeila = new TH2F("hexcEnProtonVsXLeila", "recoil exc. energy vs x leila", 1000, -100., 100., 1000, -20000., 20000.); list.Add(hexcEnProtonVsXLeila);
	TH2F* hexcEnProtonVsYLeila = new TH2F("hexcEnProtonVsYLeila", "recoil exc. energy vs y leila", 1000, -150., 150., 1000, -20000., 20000.); list.Add(hexcEnProtonVsYLeila);
	TH1F* hposdiff12xLeila = new TH1F("hposdiff12xLeila", "x positoin difference between layers leila", 1000, -50., 50.); list.Add(hposdiff12xLeila);
	TH1F* hposdiff12yLeila = new TH1F("hposdiff12yLeila", "y positoin difference between layers leila", 1000, -50., 50.); list.Add(hposdiff12yLeila);
	TH1F* hposdiff12zLeila = new TH1F("hposdiff12zLeila", "z positoin difference between layers leila", 1000, -50., 50.); list.Add(hposdiff12zLeila);
	
	TH2F* eRecErrVsESim = new TH2F("eRecErrVsESim", "error of reconstructed energy vs. simulated energy of recoil", 1000, 0., 25000., 1000, -5000., 5000.); list.Add(eRecErrVsESim);
	TH2F* thetaErrorVsZ = new TH2F("thetaErrorVsZ", "Error in #vartheta_{lab} reconstruction vs. simulated z-position;z [mm];#Delta#vartheta_{lab} [^{o}]", 200, -100, 100, 100, -15, 15); list.Add(thetaErrorVsZ);
	TH2F* thetaErrorVsTheta = new TH2F("thetaErrorVsTheta", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}]", 180, 0, 180, 100, -15, 15); list.Add(thetaErrorVsTheta);
	TH2F* zReactionEnergy = new TH2F("zReactionEnergy", "z position of reaction vs. Beam energy (rec.)", 200, -100, 100, 1000, 0, 1.1*beamEnergy); list.Add(zReactionEnergy);
	TH2F* excEnProtonVsTheta = new TH2F("excEnProtonVsTheta", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsTheta);
	TH2F* excEnProtonVsPhi = new TH2F("excEnProtonVsPhi", "Excitation Energy Spectrum from reconstructed Protons;#varphi_{lab}[^{o}];E_{exc} [keV]", 360, -180., 180., 5000, -20000, 20000); list.Add(excEnProtonVsPhi);
	TH2F* excEnProtonVsZ = new TH2F("excEnProtonVsZ", "Excitation Energy Spectrum from reconstructed Protons;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); list.Add(excEnProtonVsZ);
	TH2F* excEnProtonVsThetaGS = new TH2F("excEnProtonVsThetaGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsThetaGS);
	TH2F* excEnProtonVsZGS = new TH2F("excEnProtonVsZGS", "Excitation Energy Spectrum from reconstructed Protons, ground state only;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); list.Add(excEnProtonVsZGS);
	TH2F* excEnProtonVsThetaCm = new TH2F("excEnProtonVsThetaCm", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{cm}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsThetaCm);
	TH2F* thetaVsZ = new TH2F("thetaVsZ","#vartheta_{lab} vs. z", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZ);
	TH2F* thetaVsZSame = new TH2F("thetaVsZSame","#vartheta_{lab} vs. z, first and second layer both forward or both backward", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZSame);
	TH2F* thetaVsZCross = new TH2F("thetaVsZCross","#vartheta_{lab} vs. z, first and second layer over cross", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZCross);
	TH2F* phiVsZ = new TH2F("phiVsZ","#varphi_{lab} vs. z", 200, -100., 100., 360, -180., 180.); list.Add(phiVsZ);
	TH2F* hitpattern = new TH2F("hitpattern","detector # of second layer vs. detector # of first layer", 2, -0.5, 1.5, 2, -0.5, 1.5); list.Add(hitpattern);
	TH2F* eCmVsZ = new TH2F("eCmVsZ","energy of cm-system vs. z;z [mm];e-cm [GeV]", 200, -100., 100., 2000, transferP->GetCmEnergy(0.)/1000., transferP->GetCmEnergy(beamEnergy)/1000.); list.Add(eCmVsZ);
	TH2F* betaCmVsZ = new TH2F("betaCmVsZ","#beta of cm-system vs. z", 200, -100., 100., 2000, 0., 0.2); list.Add(betaCmVsZ);
	TH2F* stripPattern = new TH2F("stripPattern","Parallel strip # (#varphi) vs. perpendicular strip # (#vartheta)", 2*nStripsY, 0., 2.*nStripsY, 4*nStripsX, 0., 4.*nStripsX); list.Add(stripPattern);
	TH2F* recBeamEnergyErrVsZ = new TH2F("recBeamEnergyErrVsZ","Error in reconstructed beam energy vs. z", 200, -100., 100., 1000, -50., 50.); list.Add(recBeamEnergyErrVsZ);
	TH2F* thetaCmVsThetaLab = new TH2F("thetaCmVsThetaLab", "#vartheta_{cm} vs. #vartheta_{lab};#vartheta_{cm} [^{o}];#vartheta_{lab} [^{o}]", 180,0.,180., 180,0.,180.); list.Add(thetaCmVsThetaLab);
	TH2F* zErrorVsThetaSim = new TH2F("zErrorVsThetaSim", "Error between reconstructed and true z-position vs. simulated #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); list.Add(zErrorVsThetaSim);
	TH2F* zErrorVsThetaRec = new TH2F("zErrorVsThetaRec", "Error between reconstructed and true z-position vs. reconstructed #vartheta_{lab}", 180, 0., 180., 1000, -5, 5); list.Add(zErrorVsThetaRec);
	TH2F* zErrorVsthetaError = new TH2F("zErrorVstheaError", "z Erorr vs error in theta", 200, -10., 10., 1000, -5., 5.); list.Add(zErrorVsthetaError);
	TH2F* elossVsTheta = new TH2F("elossVsTheta", "reconstructed energy loss vs. #vartheta;#vartheta_lab [^{o}];energy loss [keV]", 180, 0., 180., 1000, -10000., 10000.); list.Add(elossVsTheta);
	TH2F* elossVsPhi = new TH2F("elossVsPhi", "reconstructed energy loss vs. #varphi;#varphi_lab [^{o}];energy loss [keV]", 360, -180., 180., 1000, -10000., 10000.); list.Add(elossVsPhi);
	TH2F* excEnElossVsTheta = new TH2F("excEnElossVsTheta", "Excitation Energy Spectrum from reconstructed energy loss;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnElossVsTheta);
	UInt_t nofLevels = trGen->GetMaximum("reaction")+1;
	if(nofLevels < 1) nofLevels = 1;
	if(nofLevels > 10) nofLevels = 10;
	std::vector<TH2F*> excEnElossVsThetaLevel(nofLevels);
	for(size_t r = 0; r < nofLevels; ++r) {
		excEnElossVsThetaLevel[r] = new TH2F(Form("excEnElossVsThetaLevel_%d",static_cast<int>(r)), Form("Excitation Energy Spectrum from reconstructed energy loss for level %d;#vartheta_{lab}[^{o}];E_{exc} [keV]",static_cast<int>(r)), 180, 0., 180., 5000, -20000, 20000); list.Add(excEnElossVsThetaLevel[r]);
	}

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

		if(verbose) {
			std::cout<<"*  what is setw(8)?: * "<<setw(8)<<i<<" *";
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

		if(silicon_mult_first == 1 && (silicon_mult_second == 1 || isSolid)) { 
			if(firstDeltaE[0]->size() == 1 ) {
				hit->SetFirstDeltaE(firstDeltaE[0]->at(0), kForward);
				index_first = 0;
			} else {
				hit->SetFirstDeltaE(firstDeltaE[1]->at(0), kBackward);
				index_first = 1;
			}

			if(secondDeltaE[0]->size() == 1 ) {
				hit->SetSecondDeltaE(secondDeltaE[0]->at(0), kForward);
				index_second = 0;
			} else if(secondDeltaE[1]->size() == 1) {
				hit->SetSecondDeltaE(secondDeltaE[1]->at(0), kBackward);
				index_second = 1;
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

			// get position of hit in second layer
			if(!isSolid) secondposition = hit->SecondPosition(!dontSmear);
			else         secondposition.SetXYZ(0., 0., 0.);

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
			if(verbose) std::cout<<"Beam Energy at Reaction: "<<beamEnergyRec<<" MeV"<<std::endl;


			// reconstruct energy of recoil
			recoilEnergyRecdE    =  hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose);
			recoilEnergyRecErest =  hit->GetPadEnergy();
			recoilEnergyRec = recoilEnergyRecdE + recoilEnergyRecErest;
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
			} else {
				// need to reconstruct energy loss in chamber gas, foil, and target
				double range;
				if(verbose) std::cout<<"theta "<<recoilThetaRec<<", phi "<<recoilPhiRec<<" (sinTheta "<<sinTheta<<", cosTheta "<<cosTheta<<", cosPhi "<<cosPhi<<"): ";
				if(recoilEnergyRecErest > 0. ) {
					//if(verbose) std::cout<<"from pad "<<recoilEnergyRecErest<<" through "<<(padDistance - secondLayerDistance)/(sinTheta*cosPhi)<<" mm gas ";
					range = recoilChamberGasRange->Eval(recoilEnergyRecErest);// original
					if(verbose) std::cout<<"from pad "<<recoilEnergyRecErest<<" through "<<(padDistance - secondLayerDistance)/(sinTheta)<<" mm gas ";
					range = recoilChamberGasRange->Eval(recoilEnergyRecErest);// changed bei Leila
					recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + thirdGasLayerThicknessMgCm2/(sinTheta)) + hit->GetSecondDeltaEEnergy(verbose);
					range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);
					
					// due to changing of the target foil from box to sylinder the ernergy loss is corrected by ommitting cosphi. Leila & Dennis
					
				} else {
					range = recoilChamberGasRange->Eval(hit->GetSecondDeltaEEnergy(verbose));
				}
				if(verbose) std::cout<<" with 2. layer "<<hit->GetSecondDeltaEEnergy(verbose)<<" ("<<recoilEnergyRecEloss<<") through "<<(secondLayerDistance - firstLayerDistance)/(sinTheta*cosPhi)<<" mm gas ";
				recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + secondGasLayerThicknessMgCm2/(sinTheta*cosPhi)) + hit->GetFirstDeltaEEnergy(verbose);
				// for now the foil is box-shaped as well, so we can just continue the same way
				if(verbose) std::cout<<" with 1. layer "<<hit->GetFirstDeltaEEnergy(verbose)<<" ("<<recoilEnergyRecEloss<<") through "<<(firstLayerDistance - foilDistance)/(sinTheta*cosPhi)<<" mm gas ";
				range = recoilChamberGasRange->Eval(recoilEnergyRecEloss);
				recoilEnergyRecEloss = recoilChamberGasEnergy->Eval(range + firstGasLayerThicknessMgCm2/(sinTheta*cosPhi));
				if(verbose) std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<foilThicknessMgCm2/(sinTheta*cosPhi)<<" mm foil ";
				range = recoilFoilRange->Eval(recoilEnergyRec);
				recoilEnergyRecEloss = recoilFoilEnergy->Eval(range + foilThicknessMgCm2/(sinTheta*cosPhi));
				if(verbose) std::cout<<"/n recoilEnergyRec: "<<recoilEnergyRec<<", range: "<<range<<" + foilThicknessMgCm2/(sinTheta*cosPhi) "<<foilThicknessMgCm2/(sinTheta*cosPhi)<<" => recoilEnergyRecEloss: "<<recoilEnergyRecEloss<<std::endl;
				// for now assume that the "box" inside the foil is filled with target gas
				if(verbose) std::cout<<" at "<<recoilEnergyRecEloss<<" through "<<targetWidthMgCm2/(sinTheta*cosPhi)<<" mm gas ";
				range = recoilTargetRange->Eval(recoilEnergyRecEloss);
				//std::cout<<recoilEnergyRecEloss<<", "<<range<<" + "<<targetWidthMgCm2/(sinTheta*cosPhi);
				recoilEnergyRecEloss = recoilTargetEnergy->Eval(range + targetWidthMgCm2/(sinTheta*cosPhi));
				//std::cout<<" => "<<recoilEnergyRecEloss<<std::endl;
				if(verbose) std::cout<<" => "<<recoilEnergyRecEloss<<std::endl;
			}

			//std::cout<<"theta "<<recoilThetaRec<<", phi "<<recoilPhiRec<<" (sinTheta "<<sinTheta<<", tmpPhi "<<tmpPhi<<", cosPhi "<<cosPhi<<"): "<<foilThicknessMgCm2/(sinTheta*cosPhi)<<" mg/cm^2, "<<recoilEnergyRec<<" => "<<recoilEnergyRecEloss<<", diff "<<recoilEnergyRecEloss-recoilEnergyRec<<std::endl;

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
			double excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose); // commented out by **** Leila *****
                        //double excEnergy = 0.0; //transferP->GetExcEnergy(part.GetReconstructed(), verbose); // added by **** Leila ****
			double recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;  // commented out by **** Leila *****
			if(verbose) {
				std::cout<<"beamEnergyRec "<<beamEnergyRec<<" => eex = "<<excEnergy<<" (middle spline at "<<recoilThetaRec<<" = "<<middle->Eval(recoilThetaRec)<<", recoilEnergyRec = "<<recoilEnergyRec<<")"<<std::endl;
				std::cout<<"recoilThetaCmRec = "<<recoilThetaCmRec<<", "<<transferP->GetThetacm(3)<<", "<<transferP->GetThetacm(2)<<", "<<transferP->GetThetacm(1)<<", "<<transferP->GetThetacm(0)<<std::endl;
			
				if(excEnergy>5000) std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1********************************1 excEnergy: "<<excEnergy<<" 1. layer E "<< hit->GetFirstDeltaEEnergy(verbose)<<" 2. layer E "<< hit->GetSecondDeltaEEnergy(verbose)<<" Epad: "<<recoilEnergyRecErest<<" Epad from hit: "<<hit->GetPadEnergy()<<std::endl;
			
			
			}

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
			dE1VsdE2->Fill(hit->GetFirstDeltaEEnergy(verbose), hit->GetSecondDeltaEEnergy(verbose));//(firstDeltaE[index_first]->at(0)).GetRear() );
			eVsTheta->Fill(recoilThetaRec, recoilEnergyRec);
			eVsZ->Fill(vertex.Z(), recoilEnergyRec);
			eRecErrVsESim->Fill(recoilEnergySim, recoilEnergyRec - recoilEnergySim);
			thetaErrorVsZ->Fill(vertex.Z(), recoilThetaRec - recoilThetaSim);
			//if(hit->GetPadEnergy()>1.00) thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim); 
			thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
			zReactionEnergy->Fill(vertex.Z(), beamEnergyRec);
			excEnProton->Fill(excEnergy);
			//if(recoilEnergyRecErest>1.00) excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy); ???????????????????????????
			excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);
			excEnProtonVsPhi->Fill(recoilPhiRec, excEnergy);
			excEnProtonVsThetaCm->Fill(recoilThetaCmRec, excEnergy);
			excEnProtonVsZ->Fill(vertex.Z(), excEnergy);
			if(reactionSim == 0) {
				excEnProtonVsThetaGS->Fill(recoilThetaRec, excEnergy);
				excEnProtonVsZGS->Fill(vertex.Z(), excEnergy);
			}
			thetaVsZ->Fill(vertex.Z(), recoilThetaRec);
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
			
			hERecLeila->Fill(recoilEnergyRec);
			hrecoilEnergyRecElossLeila->Fill(recoilEnergyRecEloss);
			hrecoilThetaRecLeila->Fill(recoilThetaRec);
			hrecoilThetaRecLeila->SetLineColor(kRed);
			//hrecoilThetaRecLeila->Fill(recoilThetaSim/180.*TMath::Pi());
            //if(part.GetPosition().Z()>0.) hrecoilThetaRecDiffSimLeila->Fill(recoilThetaRec - recoilThetaSim);
            //if(secondposition.Z()<0.0 || secondposition.Z()>-7.0) continue;
            hrecoilThetaRecDiffSimLeila->Fill(recoilThetaRec - recoilThetaSim);
            hrecoilEnergyRecDiffSimLeila->Fill(recoilEnergyRec - recoilEnergySim);
            hpos1xLeila->Fill(firstposition.X());
            hpos1yLeila->Fill(firstposition.Y());
            hpos1zLeila->Fill(firstposition.Z());
            hpos2xLeila->Fill(secondposition.X());
            hpos2yLeila->Fill(secondposition.Y());
            hpos2zLeila->Fill(secondposition.Z());
            hpos1xyLeila->Fill(firstposition.X(), firstposition.Y());
            hpos1ztLeila->Fill(firstposition.Z(), sqrt(firstposition.X()*firstposition.X() + firstposition.Y()*firstposition.Y()));
            hpos2xyLeila->Fill(secondposition.X(), secondposition.Y());
            hpos2ztLeila->Fill(secondposition.Z(), sqrt(secondposition.X()*secondposition.X() + secondposition.Y()*secondposition.Y()));
            hexcEnProtonVsXLeila->Fill(secondposition.X(),excEnergy);
            hexcEnProtonVsYLeila->Fill(secondposition.Y(),excEnergy);
            hexcEnProtonVsZLeila->Fill(secondposition.Z(),excEnergy);
            hposdiff12xLeila->Fill(firstposition.X() - secondposition.X());
            hposdiff12yLeila->Fill(firstposition.Y() - secondposition.Y());
            hposdiff12zLeila->Fill(firstposition.Z() - secondposition.Z());

			// now we reconstruct the q-value using the reconstructed energy loss
			// position has already been set above
			part.SetRecEnergy(recoilEnergyRecEloss);
			part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
			part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 
			transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
			transferP->SetAngles(recoilThetaRec/180.*TMath::Pi(), 2, true);
			excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose); // commented out by **** Leila ****
                        //excEnergy = 0.0; // added by **** Leila **** 
			recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;
			excEnElossVsTheta->Fill(recoilThetaRec, excEnergy);
			if(0 <= reactionSim && reactionSim < nofLevels) excEnElossVsThetaLevel[reactionSim]->Fill(recoilThetaRec, excEnergy);
			//if(0 <= reactionSim && reactionSim < nofLevels-1) excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);//leila
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


