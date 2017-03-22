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
	// command line interface
	CommandLineInterface* interface = new CommandLineInterface();

	interface->Add("-i", "inputfile (result Root file from simulation)", &InputFile);
	interface->Add("-o", "outputfile", &OutputFile);
	interface->Add("-particleType", "p or d", &particleType);
	interface->Add("-v", "verbose mode", &verbose);
	interface->Add("-t", "write output tree", &writeTree);
	interface->Add("-nosmear", "don't smear strips", &dontSmear);
	interface->Add("-n", "max. # of entries processed", &maxEntries);
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
		if (isSolid) cout <<"Using a solid target!"<<endl;
		else cout <<"Using a gas taget!"<<endl;
	}
	// nStrips are only used for the second layer (which is double-sided), and we assume that forward and backward detectors are the same
	int nStripsX = static_cast<int>(sett->GetSecondFBarrelDeltaESingleLengthX()/sett->GetSecondFBarrelDeltaESingleStripWidth());
	int nStripsY = static_cast<int>(sett->GetSecondFBarrelDeltaESingleLengthY()/sett->GetSecondFBarrelDeltaESingleStripWidth());

	std::cout<<"number of strips in x "<<nStripsX<<" ("<<sett->GetSecondFBarrelDeltaESingleLengthX()<<"/"<<sett->GetSecondFBarrelDeltaESingleStripWidth()<<")"<<std::endl;
	std::cout<<"number of strips in y "<<nStripsY<<" ("<<sett->GetSecondFBarrelDeltaESingleLengthY()<<"/"<<sett->GetSecondFBarrelDeltaESingleStripWidth()<<")"<<std::endl;

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
	//set Compounds

	if(verbose) std::cout<<"creating compound \""<<sett->GetTargetMaterialName().c_str()<<"\""<<std::endl;
	Compound* targetMat = new Compound(sett->GetTargetMaterialName().c_str());

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

	// variables for reconstruction
	Reconstruction* beamTarget = new Reconstruction(projectile, targetMat);
	double beamEnergyRec;      //beam energy at reaction, reconstructed 
	double recoilEnergyRec;    //recoil energy, reconstructed
	double recoilEnergyRecdE;
	double recoilEnergyRecErest;
	double recoilThetaRec;
	double recoilPhiRec;

	double targetThick     = sett->GetTargetThicknessMgPerCm2();
	double targetLength    = sett->GetTargetPhysicalLength();
	double targetForwardZ  =   targetLength/2.;
	double targetBackwardZ = - targetLength/2.;
	if(verbose) { 
		cout <<"Target Thickness from Input File: "<< targetThick<<endl;
		cout <<"Target ForwardZ from Input File: "<< targetForwardZ<<endl;
		cout <<"Target BackwardZ from Input File: "<< targetBackwardZ<<endl;
		cout <<"Target Length from Input File: "<< targetLength<<endl;
	}
	TSpline3* energyInTarget = beamTarget->Thickness2EnergyAfter(beamEnergy, targetThick, targetThick/1000., true);
	energyInTarget->Write("energyInTarget");

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
	TH2F* eRecErrVsESim = new TH2F("eRecErrVsESim", "error of reconstructed energy vs. simulated energy of recoil", 1000, 0., 25000., 1000, -5000., 5000.); list.Add(eRecErrVsESim);
	TH2F* thetaErrorVsZ = new TH2F("thetaErrorVsZ", "Error in #vartheta_{lab} reconstruction vs. simulated z-position;z [mm];#Delta#vartheta_{lab} [^{o}]", 200, -100, 100, 100, -15, 15); list.Add(thetaErrorVsZ);
	TH2F* thetaErrorVsTheta = new TH2F("thetaErrorVsTheta", "Error in #vartheta_{lab} reconstruction vs. simulated #vartheta_{lab};#vartheta_{lab} [^{o}];#Delta#vartheta_{lab} [^{o}]", 180, 0, 180, 100, -15, 15); list.Add(thetaErrorVsTheta);
	TH2F* zReactionEnergy = new TH2F("zReactionEnergy", "z position of reaction vs. Beam energy (rec.)", 200, -100, 100, 1000, 0, 1.1*beamEnergy); list.Add(zReactionEnergy);
	TH2F* excEnProtonVsTheta = new TH2F("excEnProtonVsTheta", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(excEnProtonVsTheta);
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

	Particle part; 

	// create and save energy vs. theta-lab splines for reaction at front/middle/back of target
	//std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy/sett->GetProjectileA()<<"/";
	std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy<<"/";
	//transferP->SetEBeam(beamEnergy/sett->GetProjectileA());
	transferP->SetEBeam(beamEnergy);
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabFront");
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<"/";
	std::cout<<energyInTarget->Eval(targetThick/2.)/1000.<<"/";
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(energyInTarget->Eval(targetThick/2.)/1000.);
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabMiddle");
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<std::endl;
	std::cout<<energyInTarget->Eval(targetThick)/1000.<<std::endl;
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(energyInTarget->Eval(targetThick)/1000.);
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabBack");

	/************************************************************************************************************
	 * loop over all entries
	 ************************************************************************************************************/

	Long64_t nEntries = tr->GetEntries();
	if(maxEntries > 0 && maxEntries < nEntries) nEntries = maxEntries;

	cout<<"Nr of Events "<<nEntries<<endl;

	if(verbose) {
		std::cout<<"***********************************************************************************************"<<std::endl
			<<"*    Row   * Instance * FBarrelDe * SecondFBa * FBarrelEr * FBarrelDe * SecondFBa * FBarrelEr *"<<std::endl
			<<"*    ***********************************************************************************************"<<std::endl;
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
			std::cout<<"*    * "<<setw(8)<<i<<" *";
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
			firstposition = hit->FirstPosition(!dontSmear); 

			// get position of hit in second layer
			if(isSolid) secondposition = hit->SecondPosition(!dontSmear);
			else        secondposition.SetXYZ(0., 0., 0.);

			part.Clear();

			// vector between two hits in Siliocn Tracker
			if(isSolid) part.SetPosition(firstposition);
			else        part.SetPosition(secondposition - firstposition); 
			if(verbose) {
				cout<<"Position to first hit: "<< firstposition.X()<<"  "<<firstposition.Y()<<"   "<<firstposition.Z()<<endl;
				cout<<"Position to second hit: "<< secondposition.X()<<"  "<<secondposition.Y()<<"   "<<secondposition.Z()<<endl;
				cout<<"Position of relative vextor: "<< part.GetPosition().X()<<"  "<<part.GetPosition().Y()<<"  "<<part.GetPosition().Z()<<endl;
			}


			// reaction angles
			recoilThetaSim = recoilThetaSim*180./TMath::Pi();
			recoilThetaRec = part.GetPosition().Theta()*180./TMath::Pi(); 
			recoilPhiSim = recoilPhiSim*180./TMath::Pi();
			recoilPhiRec = part.GetPosition().Phi()*180./TMath::Pi(); 
			if(verbose) std::cout<<recoilPhiRec<<" - "<<recoilPhiSim<<" = "<<(recoilPhiRec - recoilPhiSim)<<std::endl;


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
			if(isSolid) targetThickEvent = targetThick/2.;
			else        targetThickEvent = targetThick * ( vertex.Z() - targetBackwardZ ) / targetLength; 
			if(verbose) cout<<"Target Thickness at reaction: "<<targetThickEvent<<" = "<<targetThick<<" * ( "<<vertex.Z()<<" - "<<targetBackwardZ<<" ) / "<<targetLength<<endl;


			//calculate target thickness for reconstruction of beam energy
			if(targetThickEvent > 0) beamEnergyRec = energyInTarget->Eval(targetThickEvent)/1000.;
			else                     beamEnergyRec = beamEnergy;
			if(verbose) std::cout<<"Beam Energy at Reaction: "<<beamEnergyRec<<" keV"<<std::endl;


			// reconstruct energy of recoil
			if(silicon_mult_second == 1) {
				recoilEnergyRecdE    =  hit->GetFirstDeltaEEnergy(verbose) + hit->GetSecondDeltaEEnergy(verbose);
				recoilEnergyRecErest =  hit->GetPadEnergy();
			} else {
				recoilEnergyRecdE    =  hit->GetFirstDeltaEEnergy(verbose);
				recoilEnergyRecErest =  0.;
			}
			recoilEnergyRec = recoilEnergyRecdE + recoilEnergyRecErest;
			if(verbose && index_first == 0 && index_second == 0) std::cout<<" "<<hit->GetFirstDeltaEEnergy()<<" , "<<hit->GetSecondDeltaEEnergy()<<", "<<hit->GetPadEnergy()<<" => "<<recoilEnergyRecdE<<", "<<recoilEnergyRecErest<<" => "<<recoilEnergyRec<<std::endl;
			//update particle information
			// TODO: reconstruct energy loss in gas and foil !!

			// position has already been set above
			part.SetRecEnergy(recoilEnergyRec);
			part.SetType(2); //proton; this is for one-neutron transfer, only; this sets the mass of the particle
			part.SetReconstructed(); // set TLorentzVector using mass, rec. energy, and position 


			//////////////////////////
			// Q-value reconstruction
			//////////////////////////

			//transferP->SetEBeam(beamEnergyRec/sett->GetProjectileA());
			transferP->SetEBeam(beamEnergyRec);
			transferP->Final(recoilThetaRec/180.*TMath::Pi(), 2, true);
			double excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose);
			double recoilThetaCmRec = transferP->GetThetacm(3)/TMath::Pi()*180.;
			if(verbose) {
				std::cout<<"beamEnergyRec "<<beamEnergyRec<<" => eex = "<<excEnergy<<" (spline at "<<recoilThetaRec<<" = "<<transferP->Evslab(0., 180., 1.)->Eval(recoilThetaRec)<<", recoilEnergyRec = "<<recoilEnergyRec<<")"<<std::endl;
				std::cout<<"recoilThetaCmRec = "<<recoilThetaCmRec<<", "<<transferP->GetThetacm(3)<<", "<<transferP->GetThetacm(2)<<", "<<transferP->GetThetacm(1)<<", "<<transferP->GetThetacm(0)<<std::endl;
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
			thetaErrorVsTheta->Fill(recoilThetaSim , recoilThetaRec - recoilThetaSim);
			zReactionEnergy->Fill(vertex.Z(), beamEnergyRec);
			excEnProton->Fill(excEnergy);
			excEnProtonVsTheta->Fill(recoilThetaRec, excEnergy);
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
			if(silicon_mult_second) stripPattern->Fill(index_second*nStripsY + secondDeltaE[index_second]->at(0).GetStripNr()[0], secondDeltaE[index_second]->at(0).GetID()*nStripsX + secondDeltaE[index_second]->at(0).GetRingNr()[0]);
			recBeamEnergyErrVsZ->Fill(vertex.Z(), beamEnergyRec - reactionEnergyBeam);
			thetaCmVsThetaLab->Fill(recoilThetaRec, recoilThetaCmRec);
			zErrorVsthetaError->Fill(recoilThetaRec - recoilThetaSim, vertex.Z() - reactionZSim);
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


