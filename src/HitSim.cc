#include "HitSim.hh"

HitSim::HitSim(Settings* setting) {
	fSett = setting;
	fRand = new TRandom();
}; // end constructor


void HitSim::Clear() {
	fBarrel = NULL;
	fSecondBarrel = NULL;             //#B.Wach
}


void HitSim::InitBarrel(ParticleMC* barrel, std::string direction) {
	Clear();
	fBarrel = barrel;
	fDirection = direction;

	//limits for positions

	for(int i=0;i<4;i++) {
		if(fDirection == "forward") {
			fMaxStrip[i] = fSett->FBMaxStripPos(i);
			fMinStrip[i] = fSett->FBMinStripPos(i);
		} else {
			fMaxStrip[i] = fSett->BBMaxStripPos(i);
			fMinStrip[i] = fSett->BBMinStripPos(i);
		}
	}
}



void HitSim::InitSecondBarrel(ParticleMC* secondbarrel, std::string direction) {                     
	Clear();
	fSecondBarrel = secondbarrel;
	fDirection = direction;

	//limits for positions

	for(int i=0;i<4;i++) {
		if(fDirection == "forward") {
			fMaxStrip[i] = fSett->FBMaxStripPos(i);
			fMinStrip[i] = fSett->FBMinStripPos(i);
		} else {
			fMaxStrip[i] = fSett->BBMaxStripPos(i);
			fMinStrip[i] = fSett->BBMinStripPos(i);
		}
	}
}


/******************************************************************
 * Barrel hit position 
 * 
 * returns TVector3(0,0,0) if:
 *                        -> two not neighboring strips are hit
 *                        -> broken strip is hit (not possible in current simulation)
 *                        -> if strip position is not in the allowed range (typically: 0-1)
 *
 *****************************************************************/
TVector3 HitSim::BPosition(bool smear) {
	// variables for final hit position
	double x,y,z;

	// quadrant
	int quadr = fBarrel->GetID();

	// position along the strip (should be between 0 and 1)
	double along = 0;

	// strip number
	double strip = 0;

	// two neighboring strips hit: calculate mean strip number and mean strip position
	if(fBarrel->GetNeighbor()) { 
		for(unsigned int i = 0; i < fBarrel->GetStripNr().size(); i++) {
			if(fBarrel->GetStripPos()[i] > fMinStrip[quadr]) {
				along += fBarrel->GetStripPos()[i];
				strip += fBarrel->GetStripNr()[i];
			}
		}
		along /= fBarrel->GetStripNr().size();
		strip /= fBarrel->GetStripNr().size();
	} else if(fBarrel->GetStripNr().size() > 1) { 
		// two not neighbooring strips: ignore them
		std::cerr << "found strips " << fBarrel->GetStripNr()[0] << " and " << fBarrel->GetStripNr()[1] << "but not neighboring!" << std::endl; 
		return TVector3(0,0,0);
	} else if(fBarrel->GetStripNr().size() == 1) { 
		// one hit only
		along = fBarrel->GetStripPos()[0];
		strip = fBarrel->GetStripNr()[0];
	} else { 
		// no hit
		std::cerr << "can not find any hit " << std::endl;
		return TVector3(0,0,0);
	}

	along-=fMinStrip[quadr];
	along/=(1-fMinStrip[quadr]);

	if( (along < 0) || (along > fMaxStrip[quadr]) ) { 
		return TVector3(0,0,0);
	}

	// smear strip position
	if(smear) {
		strip+=(1-fRand->Uniform());  //uniform (0,1] - > 1-uniform [0,1)
	} else { // use mean strip position
		strip+=0.5;
	}

	// strip 0 closest to target
	if(fDirection == "forward") {
		z = fSett->GetFBarrelDeltaESinglePosZ()[quadr] - fSett->GetFBarrelDeltaESingleLengthY()/2. + strip*fSett->GetFBarrelDeltaESingleStripWidth();
		//quadr   0 top   1 lef   2 bot   3 rig
		//x       +pos    +dtb    -pos    -dtb
		//y       +dtb    -pos    -dtb    +pos

		switch(quadr) {
			case 0:
				x = 0;
				y = fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 1:
				x = fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = 0;
				break;
			case 2:
				x = 0;
				y = -fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 3:
				x = -fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = 0;
				break;
			default:
				break;
		}
	} else { // backward
		z = fSett->GetBBarrelDeltaESinglePosZ()[quadr] + fSett->GetFBarrelDeltaESingleLengthY()/2. - strip*fSett->GetFBarrelDeltaESingleStripWidth();

		switch(quadr) {
			case 0:
				x = 0;
				y = fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 1:
				x = fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = 0;
				break;
			case 2:
				x = 0;
				y = -fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 3:
				x = -fSett->GetFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = 0;
				break;
			default:
				break;
		}
	}

	// final position in lab frame
	TVector3 pos;
	pos.SetXYZ(x, y, z);

	return pos;
}

//#B.Wach

/******************************************************************
 * Second Barrel hit position 
 * 
 * returns TVector3(0,0,0) if:
 *                        -> two not neighboring strips are hit
 *                        -> broken strip is hit (not possible in current simulation)
 *                        -> if strip position is not in the allowed range (typically: 0-1)
 *
 *****************************************************************/ 
TVector3 HitSim::SecondBPosition(bool smear) {                             
	// variables for final hit position
	double x,y,z;

	// quadrant
	int quadr = fSecondBarrel->GetID();                                                   

	// position along the strip (should be between 0 and 1)
	double along = 0;

	// strip number
	double strip = 0;

	// two neighboring strips hit: calculate mean strip number and mean strip position
	if(fSecondBarrel->GetNeighbor()) {
		for(unsigned int i = 0; i < fSecondBarrel->GetStripNr().size(); i++) { 
			if(fSecondBarrel->GetStripPos()[i] > fMinStrip[quadr]) {                                    
				along += fSecondBarrel->GetStripPos()[i];                                    
				strip += fSecondBarrel->GetStripNr()[i];
			}
		}
		along /= fSecondBarrel->GetStripNr().size();
		strip /= fSecondBarrel->GetStripNr().size();
	} else if(fSecondBarrel->GetStripNr().size() > 1) { 
		// two not neighbooring strips: ignore them
		std::cerr << "second : found strips " << fSecondBarrel->GetStripNr()[0] << " and " << fSecondBarrel->GetStripNr()[1] << "but not neighboring!" << std::endl; 
		return TVector3(0,0,0);
	} else if(fSecondBarrel->GetStripNr().size() == 1) { 
		// one hit only
		along = fSecondBarrel->GetStripPos()[0];
		strip = fSecondBarrel->GetStripNr()[0];
	} else { 
		// no hit
		std::cerr << "can not find any hit " << std::endl;
		return TVector3(0,0,0);
	}

	along-=fMinStrip[quadr];
	along/=(1-fMinStrip[quadr]);

	if( (along < 0) || (along > fMaxStrip[quadr]) ) { 
		return TVector3(0,0,0);
	}

	// smear strip position
	if(smear) {
		strip+=(1-fRand->Uniform());  //uniform (0,1] - > 1-uniform [0,1)
	} else { // use mean strip position
		strip+=0.5;
	}

	// strip 0 closest to target
	if(fDirection == "forward") {
		z = fSett->GetSecondFBarrelDeltaESinglePosZ()[quadr] - fSett->GetSecondFBarrelDeltaESingleLengthY()/2. + strip*fSett->GetSecondFBarrelDeltaESingleStripWidth();
		//quadr   0 top   1 lef   2 bot   3 rig
		//x       +pos    +dtb    -pos    -dtb
		//y       +dtb    -pos    -dtb    +pos

		switch(quadr) {
			case 0:
				x = (along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				y = fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 1:
				x = fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = -(along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				break;
			case 2:
				x = -(along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				y = -fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 3:
				x = -fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = (along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				break;
			default:
				break;
		}
	} else { // backward
		z = fSett->GetBBarrelDeltaESinglePosZ()[quadr] + fSett->GetSecondFBarrelDeltaESingleLengthY()/2. - strip*fSett->GetSecondFBarrelDeltaESingleStripWidth();

		switch(quadr) {
			case 0:
				x = -(along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				y = fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 1:
				x = fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = +(along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				break;
			case 2:
				x = (along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				y = -fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				break;
			case 3:
				x = -fSett->GetSecondFBarrelDeltaESingleDistanceToBeam()[quadr];
				y = -(along-0.5)*fSett->GetSecondFBarrelDeltaESingleLengthX();
				break;
			default:
				break;
		}
	}

	// final position in lab frame
	TVector3 pos;
	pos.SetXYZ(x, y, z);

	return pos;
}
