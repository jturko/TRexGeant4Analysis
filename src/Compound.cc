#include "Compound.hh"

#define debug

Compound::Compound() {
}

Compound::Compound(const char* symbol) {
	std::string massfile = SIMULATION_PATH;
	massfile.append("/../MassFile.dat");
	int length = strlen(symbol);
	if(length == 0) {
		std::cerr<<"Error, length of Material string is zero!"<<std::endl;
		exit(1);
	}
	fSymbol = symbol;
	if(strstr(symbol,"PE")) {
		std::cout << "Polyethylene!" << std::endl;
		//H4C2
		SetNofElements(2);
		fNuclei = new Nucleus*[2];
		fFrac = new double[2];

		fNuclei[0] = new Nucleus(1,0,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());

		fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
		fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

		fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
	} else if(strstr(symbol,"DPE")) {
		std::cout << "deuterated Polyethylene!" << std::endl;
		//D4C2
		SetNofElements(2);
		fNuclei = new Nucleus*[2];
		fFrac = new double[2];

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());

		fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
		fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

		fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
	} else if(strstr(symbol,"MY")) {
		std::cout << "Mylar!" << std::endl;
		//H8C10O4
		SetNofElements(3);
		fNuclei = new Nucleus*[3];
		fFrac = new double[3];

		fNuclei[0] = new Nucleus(1,0,massfile.c_str());
		fNuclei[1] = new Nucleus(6,6,massfile.c_str());
		fNuclei[2] = new Nucleus(8,8,massfile.c_str());

		fFrac[0] = 8.*fNuclei[0]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
		fFrac[1] = 10.*fNuclei[1]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
		fFrac[2] = 4.*fNuclei[2]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());

		fMass = fNuclei[0]->GetMass()*8. + fNuclei[1]->GetMass()*10. + fNuclei[2]->GetMass()*4.;
	} else if(strstr(symbol,"TTI")) {
		std::cout << "Tritiated Titanium Target!" << std::endl;
		// ratioTTI = atomic ratio Tritium/Titanium
		if(isalpha(symbol[0])) {
			std::cerr<<"give atomic ratio of Tritium to Titanium!"<< std::endl;
			exit(1);
		} else{
			double ratio = atof(symbol);
			//std::cout << ratio << std::endl;
			SetNofElements(2);
			fNuclei = new Nucleus*[2];
			fFrac = new double[2];
			//     std::cout << "fFrac " << fFrac << std::endl;
			fNuclei[0] = new Nucleus(1,2,massfile.c_str());
			fNuclei[1] = new Nucleus(18,26,massfile.c_str());
			fFrac[0] = ratio*fNuclei[0]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
			fFrac[1] = fNuclei[1]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
			fMass = fNuclei[0]->GetMass()*ratio + fNuclei[1]->GetMass();
		}

	} else if(strstr(symbol,"DTI")) {
		//std::cout << "Deuterated Titanium Target!" << std::endl;
		if(isalpha(symbol[0])) {
			std::cerr<<"give atomic ratio of Deuterium to Titanium!"<< std::endl;
			exit(1);
		} else{
			double ratio = atof(symbol);
			std::cout << ratio << std::endl;
			SetNofElements(2);
			fNuclei = new Nucleus*[2];
			fFrac = new double[2];

			fNuclei[0] = new Nucleus(1,1,massfile.c_str());
			fFrac[0] = ratio/(1+ratio);
			fNuclei[1] = new Nucleus(22,26,massfile.c_str());
			fFrac[1] = 1/(1+ratio);

			fMass = fNuclei[0]->GetMass()*fFrac[0] + fNuclei[1]->GetMass()*fFrac[1];
		}
	} else if(strstr(symbol,"2H")) {
		std::cout << "DeuteriumGas!" << std::endl;
		SetNofElements(1);
		fNuclei = new Nucleus*[1];
		fFrac = new double[1];

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 2.0;              //*2.0 da Molekuel
	} else if(strstr(symbol,"SolidDeuterium")) {
		std::cout << "Solid Deuterium!" << std::endl;
		SetNofElements(1);
		fNuclei = new Nucleus*[1];
		fFrac = new double[1];

		fNuclei[0] = new Nucleus(1,1,massfile.c_str());

		fFrac[0] = 1;

		fMass = fNuclei[0]->GetMass() * 2.0;              //*2.0 da Molekuel
	} 

	else {
		std::cerr<<"Compound \""<<symbol<<"\" not implemented yet!"<< std::endl;
		exit(1);
	}
}

Compound::Compound(Nucleus* target) {
	SetNofElements(1);
	fNuclei = new Nucleus*[1];
	fFrac = new double[1];
	fNuclei[0] = target;
	fFrac[0] = 1;

}

Compound::~Compound() {
	// if(fFrac!=NULL)
	//   delete[] fFrac;
	// if(fNuclei!=NULL) {
	//   for(int i=0;i<GetNofElements();i++) {
	//     if(fNuclei[i]!=NULL)
	// 	delete fNuclei[i];
	//   }
	//   delete[] fNuclei;
	// }
	// if(fSymbol!=NULL)
	//   delete[] fSymbol;
}

Nucleus* Compound::GetNucleus(int i) {
	if(i<GetNofElements())
		return fNuclei[i];
	else
		return NULL;
}

double Compound::GetFrac(int i) {
	// std::cout << "i" << i << "GetNofElements()" << GetNofElements() << std::endl;
	if(i<GetNofElements())
		return fFrac[i];
	else
		return 0;
}
