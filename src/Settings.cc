#include "Settings.hh"

Settings::Settings() {
}

Settings::Settings(const char* filename) {
	SetSettingsFile(filename);
	ReadSettings();
	if(fVerboseLevel>0) PrintSettings();
}

Settings::Settings(const TRexSettings& sett) : TRexSettings(sett) {
}

Settings::~Settings() {
}

void Settings::ReadSettings(const char* filename) {
	if(filename != nullptr) SetSettingsFile(filename);
	std::cout<<"Reading settings file ..."<<std::endl; 

	TEnv* sett = new TEnv(GetSettingsFile().c_str());
	fVerboseLevel = sett->GetValue("VerboseLevel",1);

	for(int i = 0; i < 4; i++) {
		fFBMaxStripPos[i] = sett->GetValue(Form("Forward.Max.Pos.%d",i),1.0);
		fFBMinStripPos[i] = sett->GetValue(Form("Forward.Min.Pos.%d",i),0.0);
		fBBMaxStripPos[i] = sett->GetValue(Form("Backward.Max.Pos.%d",i),1.0);
		fBBMinStripPos[i] = sett->GetValue(Form("Backward.Min.Pos.%d",i),0.0);
	}
	fSmearStrip = sett->GetValue("SmearStrip", true);
	fDeadLayers = sett->GetValue("IncludeDeadLayers",1);
}

void Settings::PrintSettings() {
	std::cout<<"MassFile\t"<<GetMassFile()<<std::endl;

	std::cout<<"BeamEnergy.Nucl\t"<<GetBeamEnergy()<<std::endl;

	std::cout<<"Target.Thickness in mg/cm2 \t"<<GetTargetThicknessMgPerCm2()<<std::endl;
	std::cout<<"Target.Material\t"<<GetTargetMaterialName()<<std::endl;
	std::cout<<"Target.Ratio\t"<<GetTargetAtomicRatio()<<std::endl;
	std::cout<<"Target.N\t"<<GetTargetA()-GetTargetZ()<<std::endl;
	std::cout<<"Target.Z\t"<<GetTargetZ()<<std::endl;
	std::cout<<"Target.Diameter\t"<<GetTargetDiameter()<<std::endl;
	std::cout<<"Projectile.N\t"<<GetProjectileA()-GetProjectileZ()<<std::endl;
	std::cout<<"Projectile.Z\t"<<GetProjectileZ()<<std::endl;

	std::cout<<"Foil.FBarrel.Thickness\t"<<GetFBarrelDeltaESingleFoilThickness()<<std::endl;
	std::cout<<"Foil.BBarrel.Thickness\t"<<GetBBarrelDeltaESingleFoilThickness()<<std::endl;

	for(int i = 0; i < 4; i++) {
		std::cout<<Form("FBarrel.Delta.Thick.%d\t",i)<<GetFBarrelDeltaESingleThickness()[i]<<std::endl;
		std::cout<<Form("FBarrel.E.Thick.%d\t",i)<<GetFBarrelErestSingleThickness()[i]<<std::endl;
		std::cout<<Form("BBarrel.Delta.Thick.%d\t",i)<<GetBBarrelDeltaESingleThickness()[i]<<std::endl;
		std::cout<<Form("BBarrel.E.Thick.%d\t",i)<<GetBBarrelErestSingleThickness()[i]<<std::endl;
		std::cout<<Form("Forward.Max.Pos.%d\t",i)<<fFBMaxStripPos[i]<<std::endl;
		std::cout<<Form("Forward.Min.Pos.%d\t",i)<<fFBMinStripPos[i]<<std::endl;
		std::cout<<Form("Backward.Max.Pos.%d\t",i)<<fBBMaxStripPos[i]<<std::endl;
		std::cout<<Form("Backward.Min.Pos.%d\t",i)<<fBBMinStripPos[i]<<std::endl;
	}

	std::cout<<"SmearStripNr\t"<<fSmearStrip<<std::endl;
	std::cout<<"Barrel.DistBeam\t"<<GetFBarrelDeltaESingleDistanceToBeam()[0]<<std::endl;
	std::cout<<"Barrel.Length\t"<<GetFBarrelDeltaESingleLengthX()<<std::endl;
	std::cout<<"Forward.Barrel.Z\t"<<GetFBarrelDeltaESinglePosZ()[0]<<std::endl;
	std::cout<<"Backward.Barrel.Z\t"<<GetBBarrelDeltaESinglePosZ()[0]<<std::endl;
	std::cout<<"Second.Barrel.DistBeam\t"<<GetSecondFBarrelDeltaESingleDistanceToBeam()[0]<<std::endl;                 //#B.Wach
	std::cout<<"Second.Forward.Barrel.Z\t"<<GetSecondFBarrelDeltaESinglePosZ()[0]<<std::endl;
	std::cout<<"Second.Backward.Barrel.Z\t"<<GetSecondBBarrelDeltaESinglePosZ()[0]<<std::endl;
	std::cout<<"Barrel.Pitch\t"<<GetFBarrelDeltaESingleStripWidth()<<std::endl;
	std::cout<<"Second.Barrel.Pitch\t"<<GetSecondFBarrelDeltaESingleStripWidth()<<std::endl;
}
