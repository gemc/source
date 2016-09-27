// G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// gemc headers
#include "MPrimaryGeneratorAction.h"
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// C++ headers
#include <iostream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

MPrimaryGeneratorAction::MPrimaryGeneratorAction(goptions *opts)
{
	gemcOpt = opts;
	hd_msg         = gemcOpt->optMap["LOG_MSG"].args + " Beam Settings >> " ;
	input_gen      = gemcOpt->optMap["INPUT_GEN_FILE"].args;
	background_gen = gemcOpt->optMap["MERGE_LUND_BG"].args;
	cosmics        = gemcOpt->optMap["COSMICRAYS"].args;
	GEN_VERBOSITY  = gemcOpt->optMap["GEN_VERBOSITY"].arg;

	particleTable = G4ParticleTable::GetParticleTable();

	beamPol  = 0;

	setBeam();

	particleGun = new G4ParticleGun(1);

	if(input_gen == "gemc_internal")
	{
		vector<string> cvalues = get_info(gemcOpt->optMap["COSMICAREA"].args, string(",\""));

		if(cvalues.size() < 4)
		cout << "  !!!  Warning:  COSMICAREA flag not set correctly. It should be 4 numbers: x,y,z and R." << endl;

		cosmicTarget = G4ThreeVector(get_number(cvalues[0]), get_number(cvalues[1]), get_number(cvalues[2]));
		cosmicRadius = get_number(cvalues[3]);
		if(cvalues.size() == 5){
			cosmicGeo = cvalues[4];
		}else{
			cosmicGeo = "sph";
		}


		if(cosmics == "no")
		{
			cout << endl << hd_msg << " Beam Type: "      << Particle->GetParticleName() << endl;
			cout << hd_msg << " Beam Momentum: "    << G4BestUnit(mom, "Energy") ;
			if(dmom > 0) cout << " +- " << G4BestUnit(dmom, "Energy") ;
			cout << endl;
			cout << hd_msg << " Beam Direction: (theta, phi) = (" << theta/deg << ", " << phi/deg << ") deg" ;
			if(dtheta > 0 || dphi > 0) cout << " +- (" << dtheta/deg << ", " << dphi/deg << ") deg ";
			cout << endl << hd_msg << " Beam Vertex: (" << vx/cm << ", " << vy/cm << ", " << vz/cm << ") cm" ;
			if(drdzOrdxdydz == 0) cout << " (radius, z-spread) = (" << dvr/cm << ", " << dvz/cm << ") cm, " ;
			else                  cout << " (dvx, dvy, dvz) = (" << dvx/cm << ", " << dvy/cm << ", " << dvz/cm << ") cm, " ;
			cout << (gaussOrFlatV ?  "gaussian spread" : "flat spread");
			cout << endl << hd_msg << " Beam polarization: "    << polDeg << "%" << endl ;
			cout << hd_msg << " Polarization Direction: (theta, phi) = (" << polTheta/deg << ", " << polPhi/deg << ")" ;
			cout << endl;
		}
		else
		{
			cout << endl << hd_msg << " Beam Type: Cosmic rays."   << endl;
			cout << hd_msg << " Beam Parameters :" << cosmics << endl;
			cout << hd_msg << " a =  :" << cosmicA << endl;
			cout << hd_msg << " b =  :" << cosmicB << endl;
			cout << hd_msg << " c =  :" << cosmicC << endl;
			cout << hd_msg << " Momentum Range: [" << cminp/GeV << " - " << cmaxp/GeV << "] GeV" << endl ;
			cout << hd_msg << " Cosmic Area :" << cosmicTarget << endl;
			cout << hd_msg << " Cosmic Radius :" << cosmicRadius/cm << " cm " << endl;
			cout << hd_msg << " Cosmic Surface Type: " << cosmicGeo << endl;
			cout << hd_msg << " Cosmic Particle Type: " << cosmicParticle << endl;
		}
	}


	if(NP>0)
	{
		cout << endl << hd_msg << " Luminosity Particle Type: "      << L_Particle->GetParticleName() << endl;
		cout << hd_msg << " Luminosity Particle Momentum: "    << G4BestUnit(L_mom, "Energy") ;
		if(L_dmom > 0) cout << " +- " << G4BestUnit(L_dmom, "Energy") ;
		cout << endl;
		cout << hd_msg << " Luminosity Particle Direction: (theta, phi) = (" << L_theta/deg << ", " << L_phi/deg << ") deg" ;
		if(L_dtheta > 0 || L_dphi > 0) cout << " +- (" << L_dtheta/deg << ", " << L_dphi/deg << ") deg" ;
		cout << endl << hd_msg << " Luminosity Particle Vertex: (" << L_vx/cm << ", " << L_vy/cm << ", " << L_vz/cm << ") cm" ;
		if(L_dvr + L_dvz > 0) cout << " (radius, z-spread) = (" << L_dvr/cm << ", " << L_dvz/cm << ")" ;
		cout << endl << hd_msg << " Number of Luminosity Particles: " << NP << endl;
		cout << hd_msg << " Luminosity Time Window: " << TWINDOW/ns << " nanoseconds." << endl ;
		cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH/ns << " nanoseconds." << endl;
	}

	if(NP2>0)
	{
		cout << endl << hd_msg << " Luminosity Particle 2 Type: "      << L2_Particle->GetParticleName() << endl;
		cout << hd_msg << " Luminosity Particle 2 Momentum: "    << G4BestUnit(L2_mom, "Energy") ;
		if(L2_dmom > 0) cout << " +- " << G4BestUnit(L2_dmom, "Energy") ;
		cout << endl;
		cout << hd_msg << " Luminosity Particle 2 Direction: (theta, phi) = (" << L2_theta/deg << ", " << L2_phi/deg << ") deg" ;
		if(L2_dtheta > 0 || L2_dphi > 0) cout << " +- (" << L2_dtheta/deg << ", " << L2_dphi/deg << ") deg" ;
		cout << endl << hd_msg << " Luminosity Particle Vertex: (" << L2_vx/cm << ", " << L2_vy/cm << ", " << L2_vz/cm << ") cm" ;
		if(L2_dvr + L2_dvz > 0) cout << " (radius, z-spread) = (" << L2_dvr/cm << ", " << L2_dvz/cm << ") cm" ;
		cout << endl << hd_msg << " Number of Luminosity Particles 2: " << NP2 << endl;
		cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH2/ns << " nanoseconds." << endl;
	}

}


void MPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// internal generator. Particle defined by command line
	if(input_gen == "gemc_internal")
	{
		lundUserDefined.clear();
		for(unsigned i=0; i<8; i++) lundUserDefined.push_back(0);

		if(cosmics == "no")
		{
			// redefining particle if in graphic mode
			if( gemcOpt->optMap["USE_GUI"].arg > 0)
			setBeam();

			// Primary Particle
			particleGun->SetParticleDefinition(Particle);

			G4ThreeVector beam_dir;

			// 4-momenta
			double Mom   = mom/MeV   + (2.0*G4UniformRand()-1.0)*dmom/MeV;
			double Theta = acos(G4UniformRand()*(cos(theta/rad-dtheta/rad)-cos(theta/rad+dtheta/rad)) + cos(theta/rad+dtheta/rad))/rad;
			if(primaryFlat)
			Theta = theta/rad + (2.0*G4UniformRand()-1.0)*dtheta/rad;

			double Phi   = phi/rad   + (2.0*G4UniformRand()-1.0)*dphi/rad;
			double mass = Particle->GetPDGMass();
			double akine = sqrt(Mom*Mom + mass*mass) - mass ;
			if(gemcOpt->optMap["ALIGN_ZAXIS"].args == "no")
			beam_dir = G4ThreeVector(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad));
			else if(gemcOpt->optMap["ALIGN_ZAXIS"].args == "beamp")
			{
				beam_dir = G4ThreeVector(cos(phi/rad)*sin(theta/rad), sin(phi/rad)*sin(theta/rad), cos(theta/rad));
				const G4ThreeVector beam_axis(cos(phi/rad)*sin(theta/rad), sin(phi/rad)*sin(theta/rad), cos(theta/rad));
				const G4ThreeVector rotx1(cos(phi/rad), -sin(phi/rad), 0);
				beam_dir.rotate((2.0*G4UniformRand()-1.0)*dtheta/rad, rotx1);
				beam_dir.rotate((2.0*G4UniformRand()-1.0)*dphi/rad, beam_axis);
			}
			else
			{
				beam_dir = G4ThreeVector(cos(cphi/rad)*sin(ctheta/rad), sin(cphi/rad)*sin(ctheta/rad), cos(ctheta/rad));
				//		const G4ThreeVector beam_axis(cos(cphi/rad)*sin(ctheta/rad), sin(cphi/rad)*sin(ctheta/rad), cos(ctheta/rad));
				const G4ThreeVector beam_axis(beam_dir);
				const G4ThreeVector rotx1(cos(cphi/rad), -sin(cphi/rad), 0);
				beam_dir.rotate(Theta, rotx1);
				beam_dir.rotate(Phi, beam_axis);
			}

			particleGun->SetParticleEnergy(akine);
			particleGun->SetParticleMomentumDirection(beam_dir);
			G4ThreeVector beam_vrt;

			// vertex
			if(drdzOrdxdydz == 0)
			{

				double VR;

				if(gaussOrFlatV == 1) {
					// r is gaussian
					VR  = G4RandGauss::shoot(sqrt(vx*vx + vy*vy), dvr);
				} else {
					// r2 is flat
					VR  = sqrt(G4UniformRand())*dvr;
				}


				double PHI = 2.0*pi*G4UniformRand();

				double Vx = vx/mm + VR*cos(PHI);
				double Vy = vy/mm + VR*sin(PHI);
				double Vz = vz/mm + (2.0*G4UniformRand()-1.0)*dvz/mm;


				beam_vrt = G4ThreeVector(Vx, Vy, Vz);
			}
			else {

				double Vx, Vy, Vz;

				if(gaussOrFlatV == 1) {
					Vx  = G4RandGauss::shoot(vx/mm, dvx/mm);
					Vy  = G4RandGauss::shoot(vy/mm, dvy/mm);
					Vz  = G4RandGauss::shoot(vz/mm, dvz/mm);
				} else {
					Vx = vx/mm + (2.0*G4UniformRand()-1.0)*dvx/mm;
					Vy = vy/mm + (2.0*G4UniformRand()-1.0)*dvy/mm;
					Vz = vz/mm + (2.0*G4UniformRand()-1.0)*dvz/mm;
				}
				beam_vrt = G4ThreeVector(Vx, Vy, Vz);

			}
			particleGun->SetParticlePosition(beam_vrt);


			// polarization
			double partPol = 0.0;
			double polCast = 100.0 * G4UniformRand();
			if( polCast <= polDeg ) partPol = 1;
			double polX = partPol * sin( polTheta/rad ) * cos( polPhi/rad );
			double polY = partPol * sin( polTheta/rad ) * sin( polPhi/rad );
			double polZ = partPol * cos( polTheta/rad );
			particleGun->SetParticlePolarization(G4ThreeVector( polX, polY, polZ ));

			// Primary particle generated int the middle of Time window
			particleGun->SetParticleTime(TWINDOW/2);
			particleGun->SetNumberOfParticles(1);
			particleGun->GeneratePrimaryVertex(anEvent);
			if(GEN_VERBOSITY > 3)
			{
				cout << hd_msg << " Particle id=" <<  Particle->GetParticleName()
				<< "  Vertex=" << beam_vrt/cm << "cm,  momentum=" << Mom/GeV << " GeV, theta="
				<< Theta/deg <<  " degrees,   phi=" << Phi/deg << " degrees" << endl;
				if( partPol > 0 )
				cout << hd_msg << "   with polarization  angles: polar - " << polTheta/deg << " degrees, "
				<< "azimuthal - " << polPhi/deg << " degrees " ;
				cout << endl;
			}
		}
		// cosmic model
		// paper: A. Dar, Phys.Rev.Lett, 51,3,p.227 (1983)
		else
		{
			bool cosmicNeutrons=false;
			if(cosmicParticle != "muon") cosmicNeutrons=true;
			double thisMom;
			double thisthe;
			double thisPhi;
			double akine;

			// first randomly pick a number inside the sphere
			double cosmicVX = 100000;
			double cosmicVY = 100000;
			double cosmicVZ = 100000;

			while( (cosmicVX - cosmicTarget.x() )*(cosmicVX - cosmicTarget.x() ) +
					 (cosmicVY - cosmicTarget.y() )*(cosmicVY - cosmicTarget.y() ) +
					 (cosmicVZ - cosmicTarget.z() )*(cosmicVZ - cosmicTarget.z() ) >= cosmicRadius*cosmicRadius )
			{

				// point generated inside spherical or cylindrical volume
				if(cosmicGeo == "sph" || cosmicGeo == "sphere"){
			  // point inside spherical volume
			  cosmicVX = cosmicTarget.x() - cosmicRadius + 2*cosmicRadius*G4UniformRand();
			  cosmicVY = cosmicTarget.y() - cosmicRadius + 2*cosmicRadius*G4UniformRand();
			  cosmicVZ = cosmicTarget.z() - cosmicRadius + 2*cosmicRadius*G4UniformRand();
				}else{
			  // point inside a cylinder, height of the cylinder = cosmicRadius/2.
			  double h = cosmicRadius/2.;
			  cosmicVX = -cosmicRadius + 2*cosmicRadius*G4UniformRand();
			  double sig=1.;
			  if((2.*G4UniformRand()-1)<0) sig =-1.;
			  cosmicVY = h*(2.*G4UniformRand()-1)*sig;
			  cosmicVZ = -cosmicRadius + 2*cosmicRadius*G4UniformRand();
				}
			}

			if(cosmicNeutrons) {
				if (cminp<0.1 || cmaxp>10000) cout <<"WARNING !!!! COSMIC NEUTRONS E (MeV) is OUT OF THE VALID RANGE !!!"<<endl;
				// Model by Ashton (1973)
				// now generating random momentum, cos(theta)
				// the maximum of the distribution is the lowest momentum and 0 theta
				// normalizing by that number
				double cosmicProb = G4UniformRand()*cosmicNeutBeam(0., cminp/GeV);
				thisMom = cminp + (cmaxp-cminp)*G4UniformRand(); // momentum in MeV/c
				thisthe = pi*G4UniformRand()/2.0; // [0,pi/2] zenith angle
				while(cosmicNeutBeam(thisthe, thisMom/GeV) < cosmicProb){
					thisMom = cminp + (cmaxp-cminp)*G4UniformRand();
					thisthe = pi*G4UniformRand()/2.0;
				}
			}else{
				// muons
				// now generating random momentum, cos(theta)
				// the maximum of the distribution is the lowest momentum and 0 theta
				// normalizing by that number
				double cosmicProb = G4UniformRand()*cosmicMuBeam(0, cminp/GeV);
				thisMom = (cminp + (cmaxp-cminp)*G4UniformRand());
				thisthe = pi*G4UniformRand()/2.0;
				while(cosmicMuBeam(thisthe, thisMom/GeV) < cosmicProb)
				{
			  thisMom = (cminp + (cmaxp-cminp)*G4UniformRand());
			  thisthe = pi*G4UniformRand()/2.0;
				}
			}

			// isotropic in phi
			thisPhi = -pi + 2*pi*G4UniformRand();

			// now finding the vertex. Assuming twice the radius as starting point
			// attention:
			// axis transformation, z <> y,  x <> -x
			// cause the cosmics come from the sky
			double pvx = cosmicVX - 2*cosmicRadius*sin(thisthe)*cos(thisPhi);
			double pvy = cosmicVY + 2*cosmicRadius*cos(thisthe);
			double pvz = cosmicVZ + 2*cosmicRadius*sin(thisthe)*sin(thisPhi);
			particleGun->SetParticlePosition(G4ThreeVector(pvx, pvy, pvz));

			if(cosmicNeutrons) {
				Particle= particleTable->FindParticle("neutron");
			}else{
				// choosing charge of the muons
				string muonType = "mu+";
				if(G4UniformRand() <= 0.5)
				muonType = "mu-";
				Particle = particleTable->FindParticle(muonType);
			}
			double mass = Particle->GetPDGMass();
			akine = sqrt(thisMom*thisMom + mass*mass) - mass ;

			particleGun->SetParticleDefinition(Particle);

			// when assigning momentum the direction is reversed
			G4ThreeVector beam_dir(cos(thisPhi)*sin(thisthe), -cos(thisthe), -sin(thisPhi)*sin(thisthe));

			if(GEN_VERBOSITY > 3)
			{
				cout << hd_msg << " Particle id=" <<  Particle->GetParticleName()
				<< "  Vertex=" << G4ThreeVector(pvx, pvy, pvz)/cm << "cm,  momentum=" << thisMom/GeV << " GeV, theta="
				<< thisthe/deg <<  " degrees,   phi=" << thisPhi/deg << " degrees" << endl;
				cout << endl;
			}


			particleGun->SetParticleEnergy(akine);
			particleGun->SetParticleMomentumDirection(beam_dir);
			particleGun->SetNumberOfParticles(1);
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}
	else
	// external generator: input file
	{
		// LUND format:
		// Header (Event Info):
		// These are the original LUND variables, however after # particles, and except beam polarization, these can be user defined.
		// 1               2                     3                    4               5                 6  7  8   9   10
		// # of Particles, # of Target Nucleons, # of Target Protons, Pol. of Target, Pol. of Electron, x, y, W2, Q2, nu
		//
		// Body (Particle Info):
		// 1       2      3     4            5       6         7    8    9    10   11    12        13        14
		// index, charge, type, particle id, parent, daughter, p_x, p_y, p_z, p_t, mass, x vertex, y vertex, z vertex
		// type is 1 for particles in the detector
		if((gformat == "LUND" || gformat == "lund") && !gif.eof())
		{
			lundUserDefined.clear();
			int nparticles;

			gif >> nparticles ;
			for(unsigned i=0; i<9; i++)
			{
				double tmp;

				if(i==3)
				{
					gif >> beamPol;
					if(beamPol>1)
					beamPol = 1;
				}
				else
				{
					gif >> tmp;
					lundUserDefined.push_back(tmp);
				}
			}

			for(int p=0; p<nparticles; p++)
			{
				double tmp, px, py, pz;
				int pdef, type, parent, daughter, pindex;
				double Vx, Vy, Vz;
				gif >> pindex >> tmp >> type >> pdef >> parent >> daughter >> px >> py >> pz >> tmp >> tmp >> Vx >> Vy >> Vz;
				if(type == 1 && pindex == p+1)
				{
					// Primary Particle
					Particle = particleTable->FindParticle(pdef);
					if(!Particle)
					{
						cout << hd_msg << " Particle id " << pdef << " not found in G4 table." << endl << endl;

						return;
					}
					particleGun->SetParticleDefinition(Particle);

					// 4-momenta
					G4ThreeVector pmom(px*GeV, py*GeV, pz*GeV);
					double Mom = pmom.mag();
					double Phi   = pmom.getPhi();
					double Theta = pmom.getTheta();
					double mass = Particle->GetPDGMass();
					double akine = sqrt(Mom*Mom + mass*mass) - mass ;

					particleGun->SetParticleEnergy(akine);
					particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad)));

					// vertex
					G4ThreeVector beam_vrt(Vx*cm, Vy*cm, Vz*cm);
					particleGun->SetParticlePosition(beam_vrt);


					// beam polarization only along the beam
					// only for the first particle
					if(p==0)
					{
						particleGun->SetParticlePolarization(G4ThreeVector( 0, 0, beamPol ));
					}

					// Primary particle generated int the middle of Time window
					particleGun->SetParticleTime(TWINDOW/2);
					particleGun->SetNumberOfParticles(1);
					particleGun->GeneratePrimaryVertex(anEvent);
					if(GEN_VERBOSITY > 3)
					cout << hd_msg << " Particle Number:  " << p+1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
					<< "  Vertex=" << beam_vrt/cm << "cm,  momentum=" << pmom/GeV << " GeV" << endl;
				}
				else if(pindex != p+1)
				if(GEN_VERBOSITY > 3)
				cout << hd_msg << " Warning: file particle index " << tmp << " does not match read particle index " << p+1 << endl;

			}
		}
		else if((gformat == "stdhep" || gformat == "STDHEP" || gformat == "StdHep" || gformat == "StdHEP"))
		{
			//
			// StdHep is an (old like LUND) MC generator format in binary form.
			//

			long lerr=stdhep_reader->readEvent();  // Read the next event from the file.

			if( lerr ==  LSH_ENDOFFILE)
			{
				return;
			}
			else if(lerr != LSH_SUCCESS)
			{
				cout << hd_msg << " -- Error reading stdhep file: " << lerr << "\n";
				return;
			}

			int NPART = stdhep_reader->nTracks();

			for(int p=0;p<NPART;p++)
			{
				if( stdhep_reader->daughter1(p)>0)
				{
					// The particle has daughters, so we do not want to generate this one.
					continue;
				}
				else
				{

					Particle = particleTable->FindParticle(stdhep_reader->pid(p));
					if(!Particle)
					{
						cout << hd_msg << " Particle id " << stdhep_reader->pid(p) << " not found in G4 table." << endl << endl;

						return;
					}

					particleGun->SetParticleDefinition(Particle);

					// 4-momenta
					G4ThreeVector pmom(stdhep_reader->Px(p)*GeV,stdhep_reader->Py(p)*GeV, stdhep_reader->Pz(p)*GeV);
					double Mom = pmom.mag();
					double Phi   = pmom.getPhi();
					double Theta = pmom.getTheta();
					double mass = Particle->GetPDGMass();
					double akine = sqrt(Mom*Mom + mass*mass) - mass ;

					G4ThreeVector beam_dir(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad));

					if(gemcOpt->optMap["STEER_BEAM"].arg != 0)
					{
						beam_dir.rotateY(theta);
						beam_dir.rotateZ(phi);
					}

					particleGun->SetParticleEnergy(akine);
					particleGun->SetParticleMomentumDirection(beam_dir);

					G4ThreeVector beam_vrt;

					// vertex
					if(gemcOpt->optMap["STEER_BEAM"].arg == 0)
					{
						beam_vrt = G4ThreeVector(stdhep_reader->X(p)*cm,
														 stdhep_reader->Y(p)*cm,
														 stdhep_reader->Z(p)*cm);
					}
					else
					{
						// vertex smear and offset
						double VR  = sqrt(G4UniformRand())*dvr/mm;
						double PHI = 2.0*pi*G4UniformRand();

						beam_vrt = G4ThreeVector(stdhep_reader->X(p)*cm + vx/mm + VR*cos(PHI),
														 stdhep_reader->Y(p)*cm + vy/mm + VR*sin(PHI),
														 stdhep_reader->Z(p)*cm + vz/mm +  (2.0*G4UniformRand()-1.0)*dvz/mm);

					}
					particleGun->SetParticlePosition(beam_vrt);


					// beam polarization only along the beam
					// only for the first particle
					if(p==0)
					{
						particleGun->SetParticlePolarization(G4ThreeVector( 0, 0, beamPol ));
					}

					// Primary particle generated int the middle of Time window
					particleGun->SetParticleTime(TWINDOW/2);
					particleGun->GeneratePrimaryVertex(anEvent);
					if(GEN_VERBOSITY > 3)
					cout << hd_msg << " Particle Number:  " << p << ", id=" << stdhep_reader->pid(p) << "(" << Particle->GetParticleName() << ")"
					<< "  Vertex=" << beam_vrt << "cm,  momentum=" << pmom/GeV << " GeV" << endl;
				}
			}
		}
	}

	// merging (background) events from LUND format
	if(background_gen != "no")
	{
		int nparticles;
		double tmp;

		bgif >> nparticles ;
		for(unsigned i=0; i<9; i++)
		bgif >> tmp;

		for(int p=0; p<nparticles; p++)
		{
			double px, py, pz;
			int pdef, pindex;
			double time;
			double Vx, Vy, Vz;
			bgif >> pindex >> tmp >> tmp >> pdef >> tmp >> tmp >> px >> py >> pz >> time >> tmp >> Vx >> Vy >> Vz;
			if(pindex == p+1)
			{
				// Primary Particle
				Particle = particleTable->FindParticle(pdef);
				if(!Particle)
				{
					cout << hd_msg << " Particle id " << pdef << " not found in G4 table." << endl << endl;

					return;
				}
				particleGun->SetParticleDefinition(Particle);

				// 4-momenta
				G4ThreeVector pmom(px*GeV, py*GeV, pz*GeV);
				double Mom = pmom.mag();
				double Phi   = pmom.getPhi();
				double Theta = pmom.getTheta();
				double mass = Particle->GetPDGMass();
				double akine = sqrt(Mom*Mom + mass*mass) - mass ;

				particleGun->SetParticleEnergy(akine);
				particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad)));

				// vertex
				G4ThreeVector beam_vrt(Vx*cm, Vy*cm, Vz*cm);
				particleGun->SetParticlePosition(beam_vrt);


				// Primary particle generated int the middle of Time window
				particleGun->SetParticleTime(time);
				particleGun->GeneratePrimaryVertex(anEvent);
				if(GEN_VERBOSITY > 3)
				cout << hd_msg << " Merged Particle Number:  " << p+1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
				<< "  Vertex=" << beam_vrt << "cm,  momentum=" << pmom/GeV << " GeV" << endl;
			}
			else if(pindex != p+1)
			if(GEN_VERBOSITY > 3)
			cout << hd_msg << " Warning: file particle index " << tmp << " does not match read particle index " << p+1 << endl;

		}

	}

	// Luminosity Particles
	int NBUNCHES   = (int) floor(TWINDOW/TBUNCH);
	int PBUNCH     = (int) floor((double)NP/NBUNCHES) - 1;


	// there will be some remaining particles, these will be added at the last bunch
	int NREMAINING = NP - NBUNCHES*PBUNCH;

	// cout << PBUNCH << " " << NBUNCHES <<  " " << NBUNCHES*PBUNCH << endl;

	if(PBUNCH > 0)
	{
		particleGun->SetParticleDefinition(L_Particle);

		// getting kinematics
		double L_mass  = L_Particle->GetPDGMass();
		double L_Mom   = L_mom;
		double L_Theta = L_theta;
		double L_Phi   = L_phi;

		// all particles in a bunch are identical
		for(int b=0; b<NBUNCHES; b++)
		{
			// spread momentum if requested
			if(L_dmom > 0)
			{
				L_Mom   = L_mom + (2.0*G4UniformRand()-1.0)*L_dmom;
				L_Theta = acos(G4UniformRand()*(cos(L_theta/rad-L_dtheta/rad)-cos(L_theta/rad+L_dtheta/rad)) + cos(L_theta/rad+L_dtheta/rad))/rad;
				if(lumiFlat)
				L_Theta = L_theta + (2.0*G4UniformRand()-1.0)*L_dtheta;

				L_Phi = L_phi + (2.0*G4UniformRand()-1.0)*L_dphi;
			}
			double L_akine = sqrt(L_Mom*L_Mom + L_mass*L_mass) - L_mass ;
			particleGun->SetParticleEnergy(L_akine);
			particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(L_Phi/rad)*sin(L_Theta/rad), sin(L_Phi/rad)*sin(L_Theta/rad), cos(L_Theta/rad)));

			// luminosity vertex
			// vertex has uniform density across the cilinder
			double lvx = L_vx;
			double lvy = L_vy;

			do {
				lvx = L_vx - L_dvr + 2*G4UniformRand()*L_dvr;
				lvy = L_vy - L_dvr + 2*G4UniformRand()*L_dvr;
			} while (sqrt(lvx*lvx + lvy*lvy) > L_dvr);

			//			double L_VR  = G4UniformRand()*L_dvr;
			//			double L_PHI = 2.0*pi*G4UniformRand();
			//			double lvx = L_vx + L_VR*cos(L_PHI);
			//			double lvy = L_vy + L_VR*sin(L_PHI);
			double lvz = L_vz;

			// spread vertex if requested
			if(L_dvz > 0)
			{
				lvz = L_vz + (2.0*G4UniformRand()-1.0)*L_dvz;
			}

			particleGun->SetNumberOfParticles(PBUNCH);

			if(b == NBUNCHES-1)
			particleGun->SetNumberOfParticles(PBUNCH + NREMAINING);


			// cout << " bunch " << b << " " << PBUNCH << endl;

			particleGun->SetParticleTime(TBUNCH*b);
			particleGun->SetParticlePosition(G4ThreeVector(lvx, lvy, lvz));
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}

	// Luminosity Particles2
	int NBUNCHES2   = (int) floor(TWINDOW/TBUNCH2);
	int PBUNCH2     = (int) floor((double)NP2/NBUNCHES2);

	if(PBUNCH2 > 0)
	{
		particleGun->SetParticleDefinition(L2_Particle);
		particleGun->SetNumberOfParticles(PBUNCH);

		// getting kinematics
		double L2_mass  = L2_Particle->GetPDGMass();
		double L2_Mom   = L2_mom;
		double L2_Theta = L2_theta;
		double L2_Phi   = L2_phi;


		// all particles in a bunch are identical
		for(int b=0; b<NBUNCHES2; b++)
		{
			particleGun->SetParticleTime(TBUNCH2*b);
			// spread momentum if requested
			if(L2_dmom > 0)
			{
				L2_Mom   += (2.0*G4UniformRand()-1.0)*L2_dmom;
				L2_Theta = acos(G4UniformRand()*(cos(L2_theta/rad-L2_dtheta/rad)-cos(L2_theta/rad+L2_dtheta/rad)) + cos(L2_theta/rad+L2_dtheta/rad))/rad;
				if(lumi2Flat)
				L2_Theta += (2.0*G4UniformRand()-1.0)*L2_dtheta;
				L2_Phi   += (2.0*G4UniformRand()-1.0)*L2_dphi;
			}
			double L2_akine = sqrt(L2_Mom*L2_Mom + L2_mass*L2_mass) - L2_mass ;
			particleGun->SetParticleEnergy(L2_akine);
			particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(L2_Phi/rad)*sin(L2_Theta/rad), sin(L2_Phi/rad)*sin(L2_Theta/rad), cos(L2_Theta/rad)));

			// luminosity vertex 2
			double L2_VR  = G4UniformRand()*L2_dvr/mm;
			double L2_PHI = 2.0*pi*G4UniformRand();
			L2_vx += L2_VR*cos(L2_PHI);
			L2_vy += L2_VR*sin(L2_PHI);

			// spread vertex if requested
			if(L2_dvz > 0)
			{
				L2_vz += (2.0*G4UniformRand()-1.0)*L2_dvz;
			}
			particleGun->SetParticlePosition(G4ThreeVector(L2_vx, L2_vy, L2_vz));

			particleGun->GeneratePrimaryVertex(anEvent);
		}

	}




	if(GEN_VERBOSITY > 5)
	cout << " Generation done " << endl;
}



void MPrimaryGeneratorAction::setBeam()
{
	string hd_msg    = gemcOpt->optMap["LOG_MSG"].args + " Beam Settings >> " ;

	// vector of string - filled from the various option
	vector<string> values;
	string units;

	if(input_gen == "gemc_internal")
	{
		if(cosmics == "no")
		{
			// Getting particle name,  momentum from option value
			values       = get_info(gemcOpt->optMap["BEAM_P"].args, string(",\""));
			string pname = trimSpacesFromString(values[0]);

			if(values.size() == 4)
			{
				mom          = get_number(values[1]);
				theta        = get_number(values[2]);
				phi          = get_number(values[3]);
			}

			// making sure the particle exists
			Particle = particleTable->FindParticle(pname);
			if(!Particle)
			{
				// it may be the "show_all" option. In this case print all available particle names
				if(pname == "show_all")
				{
					for(int i=0; i<particleTable->entries(); i++)
					cout << hd_msg << " g4 particle: "  << particleTable->GetParticleName(i)
					<< " pdg encoding: " << particleTable->GetParticle(i)->GetPDGEncoding() << endl;
 			}
				// otherwise it's not found. Need to exit here.
				else
				cout << hd_msg << " Particle " << pname << " not found in G4 table. Exiting" << endl << endl;

				exit(0);
			}

			// Getting custom beam direction if it's set
			values = get_info(gemcOpt->optMap["ALIGN_ZAXIS"].args);
			string align = trimSpacesFromString(values[0]);
			if(align == "custom")
			{
				ctheta = get_number(values[1]);
				cphi   = get_number(values[2]);
			}

			// Getting momentum spread from option value
			primaryFlat = 0;
			values = get_info(gemcOpt->optMap["SPREAD_P"].args);
			dmom   = get_number(values[0]);
			dtheta = get_number(values[1]);
			dphi   = get_number(values[2]);
			if(values.size() == 4)
			if(trimSpacesFromString(values[3]) == "flat")
			primaryFlat = 1;

			// Getting vertex from option value
			values = get_info(gemcOpt->optMap["BEAM_V"].args);
			units = trimSpacesFromString(values[3]);
			vx = get_number(values[0] + "*" + units);
			vy = get_number(values[1] + "*" + units);
			vz = get_number(values[2] + "*" + units);

			// Getting vertex spread from option value
			values = get_info(gemcOpt->optMap["SPREAD_V"].args);
			// distribution type is not given. defaults to "flat"
			// number of argument can be 3 (drdz) or 4 (dxdydz)
			if(values.back().find("gauss") == string::npos && values.back().find("flat") == string::npos)
			{
				gaussOrFlatV = 0;

				// check if it's (dr, dz) or (dx, dy, dz)
				if(values.size() == 3) {
					drdzOrdxdydz = 0;
					units = trimSpacesFromString(values[2]);
					dvr = get_number(values[0] + "*" + units);
					dvz = get_number(values[1] + "*" + units);
				}
				else if(values.size() == 4) {
					drdzOrdxdydz = 1;
					units = trimSpacesFromString(values[3]);
					dvx = get_number(values[0] + "*" + units);
					dvy = get_number(values[1] + "*" + units);
					dvz = get_number(values[2] + "*" + units);
				}
			} else {
				if(values.back().find("gauss") == string::npos) gaussOrFlatV = 0;
				else gaussOrFlatV = 1;

				// check if it's (dr, dz) or (dx, dy, dz)
				if(values.size() == 4) {
					drdzOrdxdydz = 0;
					units = trimSpacesFromString(values[2]);
					dvr = get_number(values[0] + "*" + units);
					dvz = get_number(values[1] + "*" + units);
				}
				else if(values.size() == 5) {
					drdzOrdxdydz = 1;
					units = trimSpacesFromString(values[3]);
					dvx = get_number(values[0] + "*" + units);
					dvy = get_number(values[1] + "*" + units);
					dvz = get_number(values[2] + "*" + units);
				}
			}


			// Getting polarization from option value
			values = get_info(gemcOpt->optMap["POLAR"].args);
			polDeg   = get_number(values[0]);
			polTheta = get_number(values[1]);
			polPhi   = get_number(values[2]);
		}
		else
		{
			vector<string> csettings = get_info(cosmics, string(",\""));
			int len = csettings.size();
			// parsing information for COSMIC RAYS option
			if(csettings[0] == "default")
			{
				cosmicA = 55.6;
				cosmicB = 1.04;
				cosmicC = 64;

				cminp = get_number(csettings[1], 0)*GeV;
				cmaxp = get_number(csettings[2], 0)*GeV;

				// model is valid only starting at 1 GeV for now
				if(cminp < 1) cminp = 1;

				// select cosmic ray particle from data card
				if(len>3){
					cosmicParticle = csettings[3];
				}else{
					cosmicParticle = "muon";
				}
			}
			else
			{
				cosmicA = get_number(csettings[0], 0);
				cosmicB = get_number(csettings[1], 0);
				cosmicC = get_number(csettings[2], 0);

				cminp = get_number(csettings[3], 0)*GeV;
				cmaxp = get_number(csettings[4], 0)*GeV;

				// model is valid only starting at 1 GeV for now
				if(cminp < 1) cminp = 1;

				// select cosmic ray particle from data card
				if(len>5){
					cosmicParticle = csettings[5];
				}else{
					cosmicParticle = "muon";
				}
			}
		}
	}

	else if( input_gen.compare(0,4,"LUND")==0 || input_gen.compare(0,4,"lund")==0 )
	{
		gformat.assign(  input_gen, 0, input_gen.find(",")) ;
		gfilename.assign(input_gen,    input_gen.find(",") + 1, input_gen.size()) ;
		cout << hd_msg << "LUND: Opening  " << gformat << " file: " << trimSpacesFromString(gfilename).c_str() << endl;
		gif.open(trimSpacesFromString(gfilename).c_str());
		if(!gif)
		{
			cerr << hd_msg << " Can't open input file " << trimSpacesFromString(gfilename).c_str() << ". Exiting. " << endl;
			exit(1);
		}
	}

	else if( input_gen.compare(0,6,"stdhep")==0 || input_gen.compare(0,6,"STDHEP")==0 ||
			  input_gen.compare(0,6,"StdHep")==0 || input_gen.compare(0,6,"StdHEP")==0 )
	{
		// StdHep is an (old like LUND) MC generator format in binary form.
		gformat.assign(  input_gen, 0, input_gen.find(",")) ;
		gfilename.assign(input_gen,    input_gen.find(",") + 1, input_gen.size()) ;
		cout << hd_msg << "StdHEP: Opening  " << gformat << " file: " << trimSpacesFromString(gfilename).c_str() << endl;
		stdhep_reader = new lStdHep(trimSpacesFromString(gfilename).c_str());

		if(!stdhep_reader)
		{
			cerr << hd_msg << " Can't open input file " << trimSpacesFromString(gfilename).c_str() << ". Exiting. " << endl;
			exit(1);
		}

		// For the STEER_BEAM option, we need to have the angles and vertex of the GCARD in BEAM_P and BEAM_V, SPREAD_V
		// Getting particle name,  momentum from option value
		values       = get_info(gemcOpt->optMap["BEAM_P"].args);
		string pname = trimSpacesFromString(values[0]);
		mom          = get_number(values[1]);
		theta        = get_number(values[2]);
		phi          = get_number(values[3]);

		// Getting vertex from option value
		values = get_info(gemcOpt->optMap["BEAM_V"].args);
		units = trimSpacesFromString(values[3]);
		vx = get_number(values[0] + "*" + units);
		vy = get_number(values[1] + "*" + units);
		vz = get_number(values[2] + "*" + units);

		// Getting vertex spread from option value
		values = get_info(gemcOpt->optMap["SPREAD_V"].args);
		units = trimSpacesFromString(values[2]);
		dvr = get_number(values[0] + "*" + units);
		dvz = get_number(values[1] + "*" + units);

	}


	// merging (background) events from LUND format
	if(background_gen != "no")
	{
		// file may be already opened cause setBeam is called again in graphic mode
		if(!bgif.is_open() )
		{
			bgif.open(trimSpacesFromString(background_gen).c_str());
			if(!bgif)
			{
				cerr << hd_msg << " Can't open background input file >" << trimSpacesFromString(background_gen).c_str() << "<. Exiting. " << endl;
				exit(1);
			}
		}
	}


	// %%%%%%%%%%%%%%%
	// Luminosity Beam
	// %%%%%%%%%%%%%%%

	// Getting particle name,  momentum from option value
	values         = get_info(gemcOpt->optMap["LUMI_P"].args);
	string L_pname = trimSpacesFromString(values[0]);
	L_mom          = get_number(values[1]);
	L_theta        = get_number(values[2]);
	L_phi          = get_number(values[3]);


	// Getting momentum spread from option value
	lumiFlat = 0;
	values = get_info(gemcOpt->optMap["LUMI_SPREAD_P"].args);
	L_dmom   = get_number(values[0]);
	L_dtheta = get_number(values[1]);
	L_dphi   = get_number(values[2]);
	if(values.size() == 4)
	if(trimSpacesFromString(values[3]) == "flat")
	lumiFlat = 1;


	// making sure the particle exists
	L_Particle = particleTable->FindParticle(L_pname);
	if(!L_Particle)
	{
		// it may be the "show_all" option. In this case print all available particle names
		if(L_pname == "show_all")
		{
			for(int i=0; i<particleTable->entries(); i++)
			cout << hd_msg << " g4 particle: " << particleTable->GetParticleName(i) << endl;
		}
		// otherwise it's not found. Need to exit here.
		else
		cout << hd_msg << " Particle " << L_pname << " not found in G4 table. Exiting" << endl << endl;

		exit(0);
	}

	// Getting vertex from option value
	values = get_info(gemcOpt->optMap["LUMI_V"].args);
	units = trimSpacesFromString(values[3]);
	L_vx = get_number(values[0] + "*" + units);
	L_vy = get_number(values[1] + "*" + units);
	L_vz = get_number(values[2] + "*" + units);

	// Getting vertex spread from option value
	values = get_info(gemcOpt->optMap["LUMI_SPREAD_V"].args);
	units = trimSpacesFromString(values[2]);
	L_dvr = get_number(values[0] + "*" + units);
	L_dvz = get_number(values[1] + "*" + units);

	// Getting parameters from option value
	values   = get_info(gemcOpt->optMap["LUMI_EVENT"].args);
	NP       = (int) get_number(values[0]);
	TWINDOW  = get_number(values[1]);
	TBUNCH   = get_number(values[2]);



	// %%%%%%%%%%%%%%%%%
	// Luminosity Beam 2
	// %%%%%%%%%%%%%%%%%

	// Getting particle name,  momentum from option value
	values          = get_info(gemcOpt->optMap["LUMI2_P"].args);
	string L2_pname = trimSpacesFromString(values[0]);
	L2_mom          = get_number(values[1]);
	L2_theta        = get_number(values[2]);
	L2_phi          = get_number(values[3]);

	// Getting momentum spread from option value
	lumi2Flat = 0;
	values = get_info(gemcOpt->optMap["LUMI2_SPREAD_P"].args);
	L2_dmom   = get_number(values[0]);
	L2_dtheta = get_number(values[1]);
	L2_dphi   = get_number(values[2]);
	if(values.size() == 4)
	if(trimSpacesFromString(values[3]) == "flat")
	lumi2Flat = 1;


	// making sure the particle exists
	L2_Particle = particleTable->FindParticle(L2_pname);
	if(!L2_Particle)
	{
		// it may be the "show_all" option. In this case print all available particle names
		if(L_pname == "show_all")
		{
			for(int i=0; i<particleTable->entries(); i++)
			cout << hd_msg << " g4 particle: " << particleTable->GetParticleName(i) << endl;
		}
		// otherwise it's not found. Need to exit here.
		else
		cout << hd_msg << " Particle " << L2_pname << " not found in G4 table. Exiting" << endl << endl;

		exit(0);
	}

	// Getting vertex from option value
	values = get_info(gemcOpt->optMap["LUMI2_V"].args);
	units = trimSpacesFromString(values[3]);
	L2_vx = get_number(values[0] + "*" + units);
	L2_vy = get_number(values[1] + "*" + units);
	L2_vz = get_number(values[2] + "*" + units);


	// Getting vertex spread from option value
	values = get_info(gemcOpt->optMap["LUMI2_SPREAD_V"].args);
	units = trimSpacesFromString(values[2]);
	L2_dvr = get_number(values[0] + "*" + units);
	L2_dvz = get_number(values[1] + "*" + units);

	// Getting parameters from option value
	values    = get_info(gemcOpt->optMap["LUMI2_EVENT"].args);
	NP2       = (int) get_number(values[0]);
	TBUNCH2   = get_number(values[1]);

}


MPrimaryGeneratorAction::~MPrimaryGeneratorAction()
{
	delete particleGun;
	gif.close();
	bgif.close();
}


double MPrimaryGeneratorAction::cosmicMuBeam(double t, double p)
{
	return pow(cosmicA, cosmicB*cos(t))/(cosmicC*p*p);
}

double MPrimaryGeneratorAction::cosmicNeutBeam(double t, double p)
{
	// cosmic neutrons spectrum as a function of kinetic energy (GeV) and
	// zenith angle
	double massNeut = particleTable->FindParticle("neutron")->GetPDGMass()/GeV;
	double En = sqrt(p*p+massNeut*massNeut)-massNeut;
	double I0 = pow(En, -2.95);
	return I0*pow(cos(t), 3.5);
}






