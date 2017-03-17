#include "TransferReaction.hh"
#include "TFile.h"
#include "TSpline.h"
#include "TH1F.h"

/******************************************************************************
 *
 * Constructor
 *
 *****************************************************************************/
TransferReaction::TransferReaction(std::string SettingFile) {
	fTransferP_Target.resize(1600);

	// load settings file
	fSett = Settings(SettingFile.c_str());
	std::cout << "setting file: " << SettingFile  << std::endl;

	// define reaction nuclei (do this only once as it is quite slow)
	DefineReactionNuclei();

	// define target compound
	DefineTargetCompound();

	// do kinematic calculations
	DoKinematicCalculations(); 
}

TransferReaction::TransferReaction(Settings& sett) {
	fTransferP_Target.resize(1600);

	// load settings file
	fSett = sett;

	// define reaction nuclei (do this only once as it is quite slow)
	DefineReactionNuclei();

	// define target compound
	DefineTargetCompound();

	// do kinematic calculations
	DoKinematicCalculations(); 
}


/******************************************************************************
 *
 * Destructor 
 *
 *****************************************************************************/
TransferReaction::~TransferReaction(){
	// TODO Auto-generated destructor stub
}


/******************************************************************************
 *
 * Define the nuclei of the reaction
 *
 *****************************************************************************/
void TransferReaction::DefineReactionNuclei(){
	// load mass file
	const char* massfile = fSett.GetMassFile().c_str();

	std::cout << "massfile = " << massfile << std::endl;

	// define nuclei 
	std::cout << "Z = " << fSett.GetProjectileZ() << " , N = " << fSett.GetProjectileA()-fSett.GetProjectileZ() << std::endl;
	std::cout << "Z = " << fSett.GetTargetZ() << " , N = " << fSett.GetTargetA()-fSett.GetTargetZ() << std::endl;

	fProj = Nucleus(fSett.GetProjectileZ(), fSett.GetProjectileA()-fSett.GetProjectileZ(), massfile);
	fTarg = Nucleus(fSett.GetTargetZ(), fSett.GetTargetA()-fSett.GetTargetZ(), massfile);

	// 1n transfer
	fEjec = Nucleus(fSett.GetProjectileZ(), fSett.GetProjectileA()-fSett.GetProjectileZ() + 1, massfile);
	fReco = Nucleus(fSett.GetTargetZ(), fSett.GetTargetA()-fSett.GetTargetZ() - 1, massfile);

	std::cout << "Ebeam = " << fSett.GetBeamEnergy()<< std::endl;
}


/******************************************************************************
 *
 * Define target compound
 *
 *****************************************************************************/
void TransferReaction::DefineTargetCompound(){
	fTarget = Compound(fSett.GetTargetMaterialName().c_str());

	fTargetThickness = fSett.GetTargetThicknessMgPerCm2();

	std::cout<<"target: "<<fTarget.GetSymbol()<<", Thickness: "<<fTargetThickness<<" mg/cm^2 "<<std::endl;
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for elastic p, d and t
 *
 *****************************************************************************/
void TransferReaction::DoKinematicCalculations(){
	// elastic p and d 
	DoKinematicCalculationsForLightElastics();

	// transfer p
	DoKinematicCalculationsForTransferParticles();
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for elastic p, d and t
 *
 *****************************************************************************/
void TransferReaction::DoKinematicCalculationsForLightElastics(){
	std::cout << "Doing kinematic calculation for light elastic particles ..." << std::endl;

	// for energy loss calulations of the projectile in the target
	Reconstruction targetbeam(&fProj, &fTarget);


	//////////////////////////////////
	// reaction is before the target
	//////////////////////////////////
	double Ebeam_beforeTarget = fSett.GetBeamEnergy();
	std::cout << "beam energy before target: " << Ebeam_beforeTarget << std::endl;

	// kinematic before target: elastic 
	fElaP_beforeTarget = Kinematics(&fProj, &fReco, &fReco, &fProj, Ebeam_beforeTarget, 0.); 
	fElaD_beforeTarget = Kinematics(&fProj, &fTarg, &fTarg, &fProj, Ebeam_beforeTarget, 0.);


	////////////////////////////////////////////
	// reaction is in the middle of the target
	////////////////////////////////////////////
	targetbeam.SetTargetThickness(fTargetThickness/2.);
	double Ebeam_middleTarget = targetbeam.EnergyAfter(fSett.GetBeamEnergy(),-5);
	std::cout << "beam energy middle target: " << Ebeam_middleTarget << std::endl;

	// kinematic middle target: elastic 
	fElaP_middleTarget = Kinematics(&fProj, &fReco, &fReco, &fProj, Ebeam_middleTarget, 0.); 
	fElaD_middleTarget = Kinematics(&fProj, &fTarg, &fTarg, &fProj, Ebeam_middleTarget, 0.);


	/////////////////////////////////
	// reaction is after the target
	/////////////////////////////////
	targetbeam.SetTargetThickness(fTargetThickness);
	double Ebeam_afterTarget = targetbeam.EnergyAfter(fSett.GetBeamEnergy(),-5);
	std::cout << "beam energy after target: " << Ebeam_afterTarget << std::endl;

	// kinematic after target: elastic 
	fElaP_afterTarget = Kinematics(&fProj, &fReco, &fReco, &fProj, Ebeam_afterTarget, 0.); 
	fElaD_afterTarget = Kinematics(&fProj, &fTarg, &fTarg, &fProj, Ebeam_afterTarget, 0.);
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for transfer p, d, alphas
 *
 *****************************************************************************/
void TransferReaction::DoKinematicCalculationsForTransferParticles(){
	std::cout << "Doing kinematic calculation for transfer particles ..." << std::endl;

	// for energy loss calulations of the projectile in the target
	Reconstruction targetbeam(&fProj, &fTarget);


	//////////////////////////////////
	// reaction is before the target
	//////////////////////////////////
	double Ebeam_beforeTarget = fSett.GetBeamEnergy();
	std::cout << "beam energy before target: " << Ebeam_beforeTarget << std::endl;

	// 1n transfer
	fTransferP_beforeTarget = Kinematics(&fProj, &fTarg, &fReco, &fEjec, Ebeam_beforeTarget, 0.);


	////////////////////////////////////////////
	// reaction is in the middle of the target
	////////////////////////////////////////////
	targetbeam.SetTargetThickness(fTargetThickness/2.);
	double Ebeam_middleTarget = targetbeam.EnergyAfter(fSett.GetBeamEnergy(),-5);
	std::cout << "beam energy middle target: " << Ebeam_middleTarget << std::endl;

	// 1n transfer
	fTransferP_middleTarget = Kinematics(&fProj, &fTarg, &fReco, &fEjec, Ebeam_middleTarget, 0.);


	///////////////Beam energy along extended target for Q-value calculation////////////////////////

	double Ebeam_Target;
	double thickness;
	for(unsigned int i = 0; i < fTransferP_Target.size(); i++){
		thickness = fTargetThickness/fTransferP_Target.size()*i + fTargetThickness/fTransferP_Target.size() * 0.5; 
		targetbeam.SetTargetThickness(thickness);
		Ebeam_Target = targetbeam.EnergyAfter(fSett.GetBeamEnergy(),-5);
		fTransferP_Target[i] = Kinematics(&fProj, &fTarg, &fReco, &fEjec,Ebeam_Target, 0.);
		// std::cout << "index  " << i <<"  " << Ebeam_Target << std::endl;

	}



	/////////////////////////////////
	// reaction is after the target
	/////////////////////////////////
	targetbeam.SetTargetThickness(fTargetThickness);
	double Ebeam_afterTarget = targetbeam.EnergyAfter(fSett.GetBeamEnergy(),-5);
	std::cout << "beam energy after target: " << Ebeam_afterTarget << std::endl;

	// 1n transfer
	fTransferP_afterTarget = Kinematics(&fProj, &fTarg, &fReco, &fEjec, Ebeam_afterTarget, 0.);
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for elastic p, d and t
 *
 *****************************************************************************/
void TransferReaction::WriteKinematicCalculationsToFile(TFile &file){
	// elastic p, d and t
	WriteKinematicSplinesOfLightElastics(file);

	// transfer p
	WriteKinematicSplinesOfTransferParticles(file);
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for elastic p, d and t
 *
 *****************************************************************************/
void TransferReaction::WriteKinematicSplinesOfLightElastics(TFile &file){
	std::cout << "Write kinematic splines of light elastic particles ..." << std::endl;

	TSpline3* spline = NULL;

	// open directory
	file.cd();

	// create a subdirectory to store the splines
	TDirectory* subDir = file.mkdir("KinSplinesLightElastics", "KinSplinesLightElastics");

	// write splines to file
	file.cd();
	subDir->cd();

	// reaction before the target
	spline = fElaP_beforeTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaP_beforeTarget");
	spline->Write("",TObject::kOverwrite);

	spline = fElaD_beforeTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaD_beforeTarget");
	spline->Write("",TObject::kOverwrite);


	// reaction in the middle of the target
	spline = fElaP_middleTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaP_middleTarget");
	spline->Write("",TObject::kOverwrite);

	spline = fElaD_middleTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaD_middleTarget");
	spline->Write("",TObject::kOverwrite);


	// reaction after the target
	spline = fElaP_afterTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaP_afterTarget");
	spline->Write("",TObject::kOverwrite);

	spline = fElaD_afterTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("ElaD_afterTarget");
	spline->Write("",TObject::kOverwrite);


	// wenn man das macht, dann werden die Splines nicht in das Root file geschrieben
	//delete spline;
	//delete subDir;

	file.cd();
}


/******************************************************************************
 *
 * Do kinematic calculations: kinematics splines for transfer reactions
 *
 *****************************************************************************/
void TransferReaction:: WriteKinematicSplinesOfTransferParticles(TFile &file){
	std::cout << "Write kinematic splines of transfer particles ..." << std::endl;

	TSpline3* spline = NULL;

	// open directory
	file.cd();

	// create a subdirectory to store the splines
	TDirectory* subDir = file.mkdir("KinSplinesTransfer", "KinSplinesTransfer");

	// write splines to file
	file.cd();
	subDir->cd();


	// reaction before the target
	// 1n transfer
	spline = fTransferP_beforeTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("TransferP_beforeTarget_0");
	spline->Write("",TObject::kOverwrite);

	// reaction in the middle of the target
	// 1n transfer
	spline = fTransferP_middleTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("TransferP_middleTarget_0");
	spline->Write("",TObject::kOverwrite);



	///////////////Beam energy along extended target for Q-value calculation////////////////////////

	for(unsigned int i = 0; i < fTransferP_Target.size(); i++){
		spline = fTransferP_Target[i].Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps 
		spline->SetName(Form("TransferP_Target_%d_0", i));
		spline->Write("",TObject::kOverwrite);
	}



	// reaction after the target
	// 1n transfer
	spline = fTransferP_afterTarget.Evslab(1,179,1);  // from 1deg to 179deg in 1deg steps
	spline->SetName("TransferP_afterTarget_0");
	spline->Write("",TObject::kOverwrite);

	file.cd();
}
