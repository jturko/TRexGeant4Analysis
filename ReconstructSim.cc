#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "CommandLineInterface.hh"
#include "Settings.hh"
#include "HitSim.hh"
#include "ParticleMC.hh"
#include "Compound.hh"
#include "Reconstruction.hh"
#include "Kinematics.hh"
#include "Particle.hh"
#include "TransferReaction.hh"

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
	char* SettingFile = nullptr;
	string particleType;
	bool verbose = false;
	bool test = false;
	bool writeTree = false;
	// command line interface
	CommandLineInterface* interface = new CommandLineInterface();

	interface->Add("-s", "settings file", &SettingFile);
	interface->Add("-i", "inputfile (result Root file from simulation)", &InputFile);
	interface->Add("-o", "outputfile", &OutputFile);
	interface->Add("-particleType", "p or d", &particleType);
	interface->Add("-v", "verbose mode", &verbose);
	interface->Add("-t", "write output tree", &writeTree);
	interface->Add("-test", "testing: just 10 entries", &test);
	interface->CheckFlags(argc, argv);

	if(InputFile == nullptr || OutputFile == nullptr) {
		std::cerr<<"You have to provide at least one input file and the output file!"<<std::endl;
		return 1;
	}
	if(SettingFile == nullptr) {
		std::cerr<<"You have to provide a settings file!"<<std::endl;
		return 1;
	}
	std::cout<<"input file: "<<InputFile<<std::endl;

	// open input file
	TFile infile(InputFile);

	// load settings file
	TRexSettings* sim_sett = static_cast<TRexSettings*>(infile.Get("settings"));//new Settings(SettingFile);
	if(sim_sett == nullptr) {
		std::cerr<<"Failed to find \"settings\" in \""<<InputFile<<"\": "<<infile.Get("settings")<<std::endl;
		return 1;
	}
	if(verbose) sim_sett->Print();

	Settings* sett = new Settings(*sim_sett);
	if(verbose) sett->PrintSettings();
	sett->ReadSettings(SettingFile);
	if(verbose) sett->PrintSettings();

	std::cout << "beam N = " << sett->GetProjectileA()-sett->GetProjectileZ() << ", Z = " << sett->GetProjectileZ() 
		<< " on target N = " << sett->GetTargetA()-sett->GetTargetZ() << ", Z = " << sett->GetTargetZ() << std::endl;

	double beamEnergy = sett->GetBeamEnergy(); //initial beam energy (total) in MeV

	////////////////////////////////////
	// prepare input and output trees
	////////////////////////////////////
	// prepare input tree
	TTree* trGen = (TTree*) infile.Get("treeGen");
	TTree* tr = (TTree*) infile.Get("treeDet");
	if(tr == nullptr) {
		cout << "could not find tree tr in file " << infile.GetName() << endl;
		return 3;
	}

	vector<ParticleMC> *FirstBarrel[2] = {new vector<ParticleMC>, new vector<ParticleMC>} ;
	vector<ParticleMC> *SecondBarrel[2] = { new vector<ParticleMC>, new vector<ParticleMC> };

	tr->SetBranchAddress("FBarrelMC",&FirstBarrel[0]);
	tr->SetBranchAddress("SecondFBarrelDeltaEMC",&SecondBarrel[0]);
	tr->SetBranchAddress("BBarrelMC",&FirstBarrel[1]);
	tr->SetBranchAddress("SecondBBarrelDeltaEMC",&SecondBarrel[1]);                  

	double reactionEnergyBeam;
	trGen->SetBranchAddress("reactionEnergy", &reactionEnergyBeam);

	double reactionZSim;
	trGen->SetBranchAddress("reactionZ", &reactionZSim);

	double recoilThetaSim;                 //in radiants
	trGen->SetBranchAddress("recoilTheta", &recoilThetaSim);

	double recoilPhiSim;
	trGen->SetBranchAddress("recoilPhi", &recoilPhiSim);

	double recoilEnergySim;
	trGen->SetBranchAddress("recoilEnergy", &recoilEnergySim);

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
	//TransferReaction reaction(*sett);
	//Kinematics* transferP = new Kinematics(projectile, target, recoil, ejectile, beamEnergy/sett->GetProjectileA(), 0.); //reaction.GetGroundStateTransferP();
	Kinematics* transferP = new Kinematics(projectile, target, recoil, ejectile, beamEnergy, 0.); //reaction.GetGroundStateTransferP();

	// variables for reconstruction
	Reconstruction* beamTarget = new Reconstruction(projectile, targetMat);
	double beamEnergyRec;      //beam energy at reaction, reconstructed 
	double recoilEnergyRec;    //recoil energy, reconstructed
	double recoilEnergyRecdE;
	double recoilEnergyRecErest;
	double recoilThetaRec;

	double targetThick     = sett->GetTargetThicknessMgPerCm2();
	double targetLength    = sett->GetTargetPhysicalLength();
	double targetForwardZ  =   targetLength/2.;
	double targetBackwardZ = - targetLength/2.;
	if(verbose) { 
		cout <<"Target Thickness from Input File: "<< targetThick << endl;
		cout <<"Target ForwardZ from Input File: "<< targetForwardZ << endl;
		cout <<"Target BackwardZ from Input File: "<< targetBackwardZ << endl;
		cout <<"Target Length from Input File: "<< targetLength << endl;
	}

	// Define Histograms
	TList list;
	TH2F* errorOrigin = new TH2F("errorOrigin", "Error between reconstructed and true origin vs. true origin", 200, -100., 100., 1000, -5, 5); list.Add(errorOrigin);
	TH1F* errorTheta = new TH1F("errorTheta", "Error between reconstructed and true theta", 600, -30, 30); list.Add(errorTheta);
	TH1F* qValueProton = new TH1F("qValueProton", "Excitaiton Energy Spectrum from reconstructed Protons", 5000, -20000, 20000); list.Add(qValueProton);

	TH2F* dE12E = new TH2F("dE12E", "energy loss first +second layer vs total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE12E);
	TH2F* dE1E = new TH2F("dE1E", "energy loss first layer vs total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE1E);
	TH2F* dE2E = new TH2F("dE2E", "energy loss second layer vs total energy", 200, 0, 25000, 200, 0, 10000); list.Add(dE2E);
	TH2F* eVsTheta = new TH2F("eVsTheta", "recoil energy vs theta (lab)", 180, 0, 180, 200, 0, 25000); list.Add(eVsTheta);
	TH2F* eVsZ = new TH2F("eVsZ", "recoil energy vs z", 200, -100., 100., 200, 0, 25000); list.Add(eVsZ);
	TH2F* eVsZSame = new TH2F("eVsZSame", "recoil energy vs z, first and second layer both forward or both backward", 200, -100., 100., 200, 0, 25000); list.Add(eVsZSame);
	TH2F* eVsZCross = new TH2F("eVsZCross", "recoil energy vs z, first and second layer over cross", 200, -100., 100., 200, 0, 25000); list.Add(eVsZCross);
	TH2F* eRecESim = new TH2F("eRecESim", "reconstructed energy vs simulated energy of recoil", 1000, 0, 25000, 1000,0, 25000); list.Add(eRecESim);
	TH2F* zErrorTheta = new TH2F("zErrorTheta", "z position of reaction vs Error (deg.) in theta reconstruction", 200, -100, 100, 100, -15, 15); list.Add(zErrorTheta);
	TH2F* thetaErrorTheta = new TH2F("thetaErrorTheta", "Theta (recoil) vs Error (deg.) in theta reconstruction", 180, 0, 180, 100, -15, 15); list.Add(thetaErrorTheta);
	TH2F* zReactionEnergy = new TH2F("zReactionEnergy", "z position of reaction vs Beam energy (rec.)", 200, -100, 100, 1000, 0, 1.1*beamEnergy); list.Add(zReactionEnergy);
	TH2F* qValueProtonVsTheta = new TH2F("qValueProtonVsTheta", "Excitation Energy Spectrum from reconstructed Protons;#vartheta_{lab}[^{o}];E_{exc} [keV]", 180, 0., 180., 5000, -20000, 20000); list.Add(qValueProtonVsTheta);
	TH2F* qValueProtonVsZ = new TH2F("qValueProtonVsZ", "Excitation Energy Spectrum from reconstructed Protons;z [mm];E_{exc} [keV]", 200, -100., 100., 5000, -20000, 20000); list.Add(qValueProtonVsZ);
	TH2F* thetaVsZ = new TH2F("thetaVsZ","#vartheta_{lab} vs. z", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZ);
	TH2F* thetaVsZSame = new TH2F("thetaVsZSame","#vartheta_{lab} vs. z, first and second layer both forward or both backward", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZSame);
	TH2F* thetaVsZCross = new TH2F("thetaVsZCross","#vartheta_{lab} vs. z, first and second layer over cross", 200, -100., 100., 180, 0., 180.); list.Add(thetaVsZCross);
	TH2F* hitpattern = new TH2F("hitpattern","detector # of second layer vs. detector # of first layer", 2, -0.5, 1.5, 2, -0.5, 1.5); list.Add(hitpattern);
	TH2F* eCmVsZ = new TH2F("eCmVsZ","energy of cm-system vs. z;z [mm];e-cm [GeV]", 200, -100., 100., 2000, transferP->GetCmEnergy(0.)/1000., transferP->GetCmEnergy(beamEnergy)/1000.); list.Add(eCmVsZ);
	TH2F* betaCmVsZ = new TH2F("betaCmVsZ","#beta of cm-system vs. z", 200, -100., 100., 2000, 0., 0.2); list.Add(betaCmVsZ);

	Particle part; 

	// create and save energy vs theta-lab splines for reaction at front/middle/back of target
	//std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy/sett->GetProjectileA()<<"/";
	std::cout<<"beam energy at front/middle/back of target: "<<beamEnergy<<"/";
	//transferP->SetEBeam(beamEnergy/sett->GetProjectileA());
	transferP->SetEBeam(beamEnergy);
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabFront");
	beamTarget->SetTargetThickness(targetThick/2.);
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<"/";
	std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)<<"/";
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true));
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabMiddle");
	beamTarget->SetTargetThickness(targetThick);
	//std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA()<<std::endl;
	std::cout<<beamTarget->EnergyAfter(beamEnergy, -3, true)<<std::endl;
	//transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true)/sett->GetProjectileA());
	transferP->SetEBeam(beamTarget->EnergyAfter(beamEnergy, -3, true));
	transferP->Evslab(0., 180., 1.)->Write("RecoilEVsThetaLabBack");

	/************************************************************************************************************
	 * loop over all entries
	 ************************************************************************************************************/

	Double_t nentries = tr->GetEntries();
	if(test)
		nentries = 10;

	cout << "Nr of Events " << nentries << endl;

	for(int i=0; i<nentries;i++){
		if(verbose) cout <<"Loop over entry Nr "<<i<<endl;
		ParticleBranch->clear();
		tr->GetEntry(i);
		trGen->GetEntry(i);
		Int_t silicon_mult_first = FirstBarrel[0]->size()+ FirstBarrel[1]->size();
		Int_t silicon_mult_second = SecondBarrel[0]->size()+ SecondBarrel[1]->size();
		TVector3 firstposition;
		TVector3 secondposition;

		Int_t index_first = 0;
		Int_t index_second = 0;

		// exactly one hit in any first layer and one hit in any second layer ?
		//TODO: take into account multiple hits

		if(silicon_mult_first > 1 || silicon_mult_second > 1) {
			cout <<"Warning: Multiple hits in Silicon Tracker! "<< endl;
			cout <<"First layer:  "<< silicon_mult_first<< "  Second layer:  "<<silicon_mult_second<< endl;
		}

		if(silicon_mult_first == 1 && silicon_mult_second == 1) { 
			if(FirstBarrel[0]->size() == 1 ) {
				hit->InitBarrel(&(FirstBarrel[0]->at(0)), "forward");
				index_first = 0;
			} else {
				hit->InitBarrel(&(FirstBarrel[1]->at(0)), "backward");
				index_first = 1;
			}

			//get position of hit in first layer
			firstposition = hit->BPosition(sett->SmearStrip()); 

			if(SecondBarrel[0]->size() == 1 ) {
				hit->InitSecondBarrel(&(SecondBarrel[0]->at(0)), "forward");
				index_second = 0;
			} else {
				hit->InitSecondBarrel(&(SecondBarrel[1]->at(0)), "backward");
				index_second = 1;
			}
			// get position of hit in second layer
			secondposition = hit->SecondBPosition(sett->SmearStrip());

			part.Clear();

			// vector between two hits in Siliocn Tracker
			part.SetPosition(secondposition - firstposition); 
			if(verbose) {
				cout << "Position to first hit: "<< firstposition.X() << "  "<<firstposition.Y()<<"   "<<firstposition.Z()<<endl;
				cout << "Position to second hit: "<< secondposition.X() << "  "<<secondposition.Y()<<"   "<<secondposition.Z()<<endl;
				cout << "Position of relative vextor: "<< part.GetPosition().X() << "  " << part.GetPosition().Y() << "  " << part.GetPosition().Z() << endl;
			}


			// reaction angles
			recoilThetaSim = recoilThetaSim*180./TMath::Pi();
			recoilThetaRec = part.GetPosition().Theta()*180./TMath::Pi(); 


			//find the closest point between beam axis and vector of the two hits in the silicon tracker
			TVector3 r = part.GetPosition();  //relative vector from first hit to second hit
			TVector3 r2 = secondposition;     // vector to second hit
			TVector3 vertex;                   //reconstructed vertex
			double t = 0;                          //line parameter to calculate vertex; temp use only
			if((r*r - r.Z()*r.Z()) != 0 ) t = (r2*r - (r2.Z()*r.Z()))/(r*r - r.Z()*r.Z());
			vertex = r2 -( t*r); 
			if(verbose) {
				cout <<"Vertex: "<< vertex.X() <<"  "<<vertex.Y() <<"   "<< vertex.Z()<<endl;
				cout <<"Z from simu: "<< (reactionZSim)<<endl;
			}
			//update particle information

			// target length at reaction
			double targetThickEvent = targetThick * ( vertex.Z() - targetBackwardZ ) / targetLength; 
			if(verbose) cout<<"Target Thickness at reaction: "<<targetThickEvent<<" = "<<targetThick<<" * ( "<<vertex.Z()<<" - "<<targetBackwardZ<<" ) / "<<targetLength<<endl;


			//calculate target thickness for reconstruction of beam energy
			beamTarget->SetTargetThickness(targetThickEvent);
			beamEnergyRec = beamTarget->EnergyAfter(beamEnergy, -3, true);
			if(verbose) cout <<"Beam Energy at Reaction: "<< beamEnergyRec << endl;


			// reconstruct energy of recoil
			recoilEnergyRecErest =  (FirstBarrel[index_first]->at(0)).GetEdet() ;
			recoilEnergyRecdE    =  (FirstBarrel[index_first]->at(0)).GetRear()+ (SecondBarrel[index_second]->at(0)).GetRear() ;
			recoilEnergyRec = recoilEnergyRecdE + recoilEnergyRecErest;

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
			double excEnergy = transferP->GetExcEnergy(part.GetReconstructed(), verbose);
			if(verbose) {
				std::cout<<"beamEnergyRec "<<beamEnergyRec<<" => eex = "<<excEnergy<<" (spline at "<<recoilThetaRec<<" = "<<transferP->Evslab(0., 180., 1.)->Eval(recoilThetaRec)<<", recoilEnergyRec = "<<recoilEnergyRec<<")"<<std::endl;
			}

			///////////////////////
			// Fill some histograms
			///////////////////////

			hitpattern->Fill(index_first, index_second);
			errorOrigin->Fill(reactionZSim,  vertex.Z()-reactionZSim );
			errorTheta->Fill( recoilThetaSim - recoilThetaRec );
			dE12E->Fill( recoilEnergyRec, recoilEnergyRecdE );
			dE1E->Fill( recoilEnergyRec, (FirstBarrel[index_first]->at(0)).GetRear() );
			dE2E->Fill( recoilEnergyRec, (SecondBarrel[index_second]->at(0)).GetRear() );
			eVsTheta->Fill( recoilThetaRec, recoilEnergyRec );
			eVsZ->Fill(vertex.Z(), recoilEnergyRec);
			if(index_first == index_second) eVsZSame->Fill(vertex.Z(), recoilEnergyRec);
			else                            eVsZCross->Fill(vertex.Z(), recoilEnergyRec);
			eRecESim->Fill( recoilEnergySim, recoilEnergyRec );
			zErrorTheta->Fill( reactionZSim, recoilThetaSim - recoilThetaRec );
			thetaErrorTheta->Fill( recoilThetaSim , recoilThetaSim - recoilThetaRec );
			zReactionEnergy->Fill(  reactionZSim, beamEnergyRec );
			qValueProton->Fill( excEnergy );
			qValueProtonVsTheta->Fill(recoilThetaRec, excEnergy);
			qValueProtonVsZ->Fill(vertex.Z(), excEnergy);
			thetaVsZ->Fill(vertex.Z(), recoilThetaRec);
			if(index_first == index_second) thetaVsZSame->Fill(vertex.Z(), recoilThetaRec);
			else                            thetaVsZCross->Fill(vertex.Z(), recoilThetaRec);
			betaCmVsZ->Fill(vertex.Z(), transferP->GetBetacm());
			eCmVsZ->Fill(vertex.Z(), transferP->GetCmEnergy()/1000.);
		}   //end of mult = 1 events
		if(i%1000 == 0){
			cout<<setw(5)<<std::fixed<<setprecision(1)<<(100.*i)/nentries<<setprecision(3)<<" % done\r"<<flush;
		}

	} // end of loop over all events



	//////////////////
	// write results
	//////////////////
	outfile.cd();
	if(writeTree) rectr->Write("",TObject::kOverwrite);

	list.Sort();
	list.Write();

	cout << "closing file ..." << endl;
	infile.Close();
	outfile.Close();

	cout<<endl<<"Done!"<<endl;

	return 0;

}


