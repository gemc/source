// G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// gemc headers
#include "MPrimaryGeneratorAction.h"
#include "detector.h"

// C++ headers
#include <iostream>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

MPrimaryGeneratorAction::MPrimaryGeneratorAction(goptions *opts)
{
	gemcOpt = opts;
	hd_msg        = gemcOpt->optMap["LOG_MSG"].args + " Beam Settings >> " ;
	input_gen     = gemcOpt->optMap["INPUT_GEN_FILE"].args;
	GEN_VERBOSITY = gemcOpt->optMap["GEN_VERBOSITY"].arg;
	
	particleTable = G4ParticleTable::GetParticleTable();
	
	beamPol  = 0;
	beam_pol = G4ThreeVector(0.0, 0.0, 0.0);
	
	
	setBeam();
	
	particleGun = new G4ParticleGun(1);
	
	if(input_gen == "gemc_internal")
	{
		cout << endl << hd_msg << " Beam Type: "      << Particle->GetParticleName() << endl;
		cout << hd_msg << " Beam Momentum: "    << G4BestUnit(mom, "Energy") ;
		if(dmom > 0) cout << " +- " << G4BestUnit(dmom, "Energy") ;
		cout << endl;
		cout << hd_msg << " Beam Direction: (theta, phi) = (" << theta/deg << ", " << phi/deg << ")" ;
		if(dtheta > 0 || dphi > 0) cout << " +- (" << dtheta/deg << ", " << dphi/deg << ")" ;
		cout << " deg " << endl;
		cout << hd_msg << " Beam Vertex: (" << vx/cm << ", " << vy/cm << ", " << vz/cm << ")" ;
		if(dvr + dvz > 0) cout << " (radius, z-spread) = (" << dvr/cm << ", " << dvz/cm << ")" ;
		cout << " cm " << endl;
		cout << hd_msg << " Beam polarization: "    << polDeg << "%" ;
		cout << endl ;
		cout << hd_msg << " Polarization Direction: (theta, phi) = (" << polTheta/deg << ", " << polPhi/deg << ")" ;
		cout << endl;
	}
	
	
	if(NP>0)
	{
		cout << endl << hd_msg << " Luminosity Particle Type: "      << L_Particle->GetParticleName() << endl;
		cout << hd_msg << " Luminosity Particle Momentum: "    << G4BestUnit(L_Mom, "Energy") ;
		cout << endl;
		cout << hd_msg << " Luminosity Particle Direction: (theta, phi) = (" << L_Theta/deg << ", " << L_Phi/deg << ")" ;
		cout << " deg " << endl;
		cout << hd_msg << " Luminosity Particle Vertex: (" << L_vx/cm << ", " << L_vy/cm << ", " << L_vz/cm << ")" ;
		if(L_dvr + L_dvz > 0) cout << " (radius, z-spread) = (" << L_dvr/cm << ", " << L_dvz/cm << ")" ;
		cout << " cm " << endl;
		cout << hd_msg << " Number of Luminosity Particles: " << NP << endl;
		cout << hd_msg << " Luminosity Time Window: " << TWINDOW/ns << " nanoseconds." << endl ;
		cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH/ns << " nanoseconds." << endl;
	}
	
	if(NP2>0)
	{
		cout << endl << hd_msg << " Luminosity Particle 2 Type: "      << L2_Particle->GetParticleName() << endl;
		cout << hd_msg << " Luminosity Particle 2 Momentum: "    << G4BestUnit(L2_Mom, "Energy") ;
		cout << endl;
		cout << hd_msg << " Luminosity Particle 2 Direction: (theta, phi) = (" << L2_Theta/deg << ", " << L2_Phi/deg << ")" ;
		cout << " deg " << endl;
		cout << hd_msg << " Luminosity Particle Vertex: (" << L2_vx/cm << ", " << L2_vy/cm << ", " << L2_vz/cm << ")" ;
		if(L2_dvr + L2_dvz > 0) cout << " (radius, z-spread) = (" << L2_dvr/cm << ", " << L2_dvz/cm << ")" ;
		cout << " cm " << endl;
		cout << hd_msg << " Number of Luminosity Particles 2: " << NP2 << endl;
		cout << hd_msg << " Luminosity Time Between Bunches: " << TBUNCH2/ns << " nanoseconds." << endl;
	}
	
}


void MPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// internal generator. Particle defined by command line
	if(input_gen == "gemc_internal")
	{
		// redefining particle if in graphic mode
		if( gemcOpt->optMap["USE_GUI"].arg > 0)
		setBeam();
		
		lundUserDefined.clear();
		for(unsigned i=0; i<8; i++) lundUserDefined.push_back(0);
		
		// Primary Particle
		particleGun->SetParticleDefinition(Particle);
		
		// 4-momenta
		Mom   = mom/MeV   + (2.0*G4UniformRand()-1.0)*dmom/MeV;
		Theta = theta/rad + (2.0*G4UniformRand()-1.0)*dtheta/rad;
		Phi   = phi/rad   + (2.0*G4UniformRand()-1.0)*dphi/rad;
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
		
		// vertex
		double VR  = sqrt(G4UniformRand())*dvr/mm;
		double PHI = 2.0*pi*G4UniformRand();
		Vx = vx/mm + VR*cos(PHI);
		Vy = vy/mm + VR*sin(PHI);
		Vz = vz/mm + (2.0*G4UniformRand()-1.0)*dvz/mm;
		beam_vrt = G4ThreeVector(Vx, Vy, Vz);
		particleGun->SetParticlePosition(beam_vrt);
		
		// polarization
		double partPol = 0.0;
		double polCast = 100.0 * G4UniformRand();
		if( polCast <= polDeg ) partPol = 1;
		double polX = partPol * sin( polTheta/rad ) * cos( polPhi/rad );
		double polY = partPol * sin( polTheta/rad ) * sin( polPhi/rad );
		double polZ = partPol * cos( polTheta/rad );
		beam_pol = G4ThreeVector( polX, polY, polZ );
		particleGun->SetParticlePolarization(beam_pol);
		
		// Primary particle generated int the middle of Time window
		particleGun->SetParticleTime(TWINDOW/2);
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
					Mom = pmom.mag();
					Phi   = pmom.getPhi();
					Theta = pmom.getTheta();
					double mass = Particle->GetPDGMass();
					double akine = sqrt(Mom*Mom + mass*mass) - mass ;
					
					beam_dir = G4ThreeVector(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad));
					particleGun->SetParticleEnergy(akine);
					particleGun->SetParticleMomentumDirection(beam_dir);
					
					// vertex
					beam_vrt = G4ThreeVector(Vx*cm, Vy*cm, Vz*cm);
					particleGun->SetParticlePosition(beam_vrt);
					
					
					// beam polarization only along the beam
					// only for the first particle
					if(p==0)
					{
						beam_pol = G4ThreeVector( 0, 0, beamPol );
						particleGun->SetParticlePolarization(beam_pol);
					}
					
					// Primary particle generated int the middle of Time window
					particleGun->SetParticleTime(TWINDOW/2);
					particleGun->GeneratePrimaryVertex(anEvent);
					if(GEN_VERBOSITY > 3)
						cout << hd_msg << " Particle Number:  " << p+1 << ", id=" << pdef << " (" << Particle->GetParticleName() << ")"
					         << "  Vertex=" << beam_vrt << "cm,  momentum=" << pmom/GeV << " GeV" << endl;
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
					Mom = pmom.mag();
					Phi   = pmom.getPhi();
					Theta = pmom.getTheta();
					double mass = Particle->GetPDGMass();
					double akine = sqrt(Mom*Mom + mass*mass) - mass ;
					//double akine = stdhep_reader->E(p)*GeV;
					
					beam_dir = G4ThreeVector(cos(Phi/rad)*sin(Theta/rad), sin(Phi/rad)*sin(Theta/rad), cos(Theta/rad));
					
					if(gemcOpt->optMap["STEER_BEAM"].arg != 0){
						beam_dir.rotateY(theta);
						beam_dir.rotateZ(phi);
					}
					
					particleGun->SetParticleEnergy(akine);
					particleGun->SetParticleMomentumDirection(beam_dir);
					
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
						beam_pol = G4ThreeVector( 0, 0, beamPol );
						particleGun->SetParticlePolarization(beam_pol);
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
	
	
	// Luminosity Particles
	int NBUNCHES   = (int) floor(TWINDOW/TBUNCH);
	int PBUNCH     = (int) floor((double)NP/NBUNCHES);
	
	particleGun->SetParticleDefinition(L_Particle);
	double L_mass = L_Particle->GetPDGMass();

	for(int b=0; b<NBUNCHES; b++)
	{
		for(int p=0; p<PBUNCH; p++)
		{
			
			// luminosity momentum
			L_Mom   = L_mom/MeV   + (2.0*G4UniformRand()-1.0)*L_dmom/MeV;
			L_Theta = L_theta/rad + (2.0*G4UniformRand()-1.0)*L_dtheta/rad;
			L_Phi   = L_phi/rad   + (2.0*G4UniformRand()-1.0)*L_dphi/rad;
			
			double L_akine = sqrt(L_Mom*L_Mom + L_mass*L_mass) - L_mass ;
			particleGun->SetParticleEnergy(L_akine);

			L_beam_dir     = G4ThreeVector(cos(L_Phi/rad)*sin(L_Theta/rad), sin(L_Phi/rad)*sin(L_Theta/rad), cos(L_Theta/rad));
			particleGun->SetParticleMomentumDirection(L_beam_dir);

			
			// luminosity vertex
			double L_VR  = G4UniformRand()*L_dvr/mm;
			double L_PHI = 2.0*pi*G4UniformRand();
			L_vx = L_beam_vrt.x()/mm + L_VR*cos(L_PHI);
			L_vy = L_beam_vrt.y()/mm + L_VR*sin(L_PHI);
			L_vz = L_beam_vrt.z()/mm + (2.0*G4UniformRand()-1.0)*L_dvz/mm;
			G4ThreeVector LUMI_V(L_vx, L_vy, L_vz);
			particleGun->SetParticlePosition(LUMI_V);
			particleGun->  SetParticleTime(TBUNCH*b);
			particleGun->GeneratePrimaryVertex(anEvent);
			
			if(GEN_VERBOSITY > 5)
				cout << " Bunch " << b << "  particle no. " << p << endl;

		}
	}
	
	
	// Luminosity Particles2
	int NBUNCHES2   = (int) floor(TWINDOW/TBUNCH2);
	int PBUNCH2     = (int) floor((double)NP2/NBUNCHES2);
	
	particleGun->SetParticleDefinition(L2_Particle);
	double L2_mass  = L2_Particle->GetPDGMass();
	
	for(int b=0; b<NBUNCHES2; b++)
	{
		for(int p=0; p<PBUNCH2; p++)
		{
			// luminosity momentum
			L2_Mom   = L2_mom/MeV   + (2.0*G4UniformRand()-1.0)*L2_dmom/MeV;
			L2_Theta = L2_theta/rad + (2.0*G4UniformRand()-1.0)*L2_dtheta/rad;
			L2_Phi   = L2_phi/rad   + (2.0*G4UniformRand()-1.0)*L2_dphi/rad;
			
			
			double L2_akine = sqrt(L2_Mom*L2_Mom + L2_mass*L2_mass) - L2_mass ;
			particleGun->SetParticleEnergy(L2_akine);
			
			L2_beam_dir     = G4ThreeVector(cos(L2_Phi/rad)*sin(L2_Theta/rad), sin(L2_Phi/rad)*sin(L2_Theta/rad), cos(L2_Theta/rad));
			particleGun->SetParticleMomentumDirection(L2_beam_dir);

			
			// luminosity vertex 2
			double L2_VR  = G4UniformRand()*L2_dvr/mm;
			double L2_PHI = 2.0*pi*G4UniformRand();
			L2_vx = L2_beam_vrt.x()/mm + L2_VR*cos(L2_PHI);
			L2_vy = L2_beam_vrt.y()/mm + L2_VR*sin(L2_PHI);
			L2_vz = L2_beam_vrt.z()/mm + (2.0*G4UniformRand()-1.0)*L2_dvz/mm;
			G4ThreeVector LUMI2_V(L2_vx, L2_vy, L2_vz);
			particleGun->SetParticlePosition(LUMI2_V);
			
			particleGun->  SetParticleTime(TBUNCH2*b);
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
		
		// Primary Beam
		
		// Getting particle name,  momentum from option value
		values       = get_info(gemcOpt->optMap["BEAM_P"].args, string(",\""));
		string pname = TrimSpaces(values[0]);
		
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
		string align = TrimSpaces(values[0]);
		if(align == "custom")
		{
			ctheta = get_number(values[1]);
			cphi   = get_number(values[2]);
		}
		
		// Getting momentum spread from option value
		values = get_info(gemcOpt->optMap["SPREAD_P"].args);
		dmom   = get_number(values[0]);
		dtheta = get_number(values[1]);
		dphi   = get_number(values[2]);
		
		// Getting vertex from option value
		values = get_info(gemcOpt->optMap["BEAM_V"].args);
		units = TrimSpaces(values[3]);
		vx = get_number(values[0] + "*" + units);
		vy = get_number(values[1] + "*" + units);
		vz = get_number(values[2] + "*" + units);
		
		// Getting vertex spread from option value
		values = get_info(gemcOpt->optMap["SPREAD_V"].args);
		units = TrimSpaces(values[2]);
		dvr = get_number(values[0] + "*" + units);
		dvz = get_number(values[1] + "*" + units);
		
		// Getting polarization from option value
		values = get_info(gemcOpt->optMap["POLAR"].args);
		polDeg   = get_number(values[0]);
		polTheta = get_number(values[1]);
		polPhi   = get_number(values[2]);
		
		
	}
	else  if( input_gen.compare(0,4,"LUND")==0 || input_gen.compare(0,4,"lund")==0 )
	{
		gformat.assign(  input_gen, 0, input_gen.find(",")) ;
		gfilename.assign(input_gen,    input_gen.find(",") + 1, input_gen.size()) ;
		cout << hd_msg << "LUND: Opening  " << gformat << " file: " << TrimSpaces(gfilename).c_str() << endl;
		gif.open(TrimSpaces(gfilename).c_str());
		if(!gif)
		{
			cerr << hd_msg << " Can't open input file " << TrimSpaces(gfilename).c_str() << ". Exiting. " << endl;
			exit(1);
		}
	}
	else if( input_gen.compare(0,6,"stdhep")==0 || input_gen.compare(0,6,"STDHEP")==0 ||
			 input_gen.compare(0,6,"StdHep")==0 || input_gen.compare(0,6,"StdHEP")==0 )
	{
		// StdHep is an (old like LUND) MC generator format in binary form.
		//
		gformat.assign(  input_gen, 0, input_gen.find(",")) ;
		gfilename.assign(input_gen,    input_gen.find(",") + 1, input_gen.size()) ;
		cout << hd_msg << "StdHEP: Opening  " << gformat << " file: " << TrimSpaces(gfilename).c_str() << endl;
		stdhep_reader = new lStdHep(TrimSpaces(gfilename).c_str());
		
		if(!stdhep_reader)
		{
			cerr << hd_msg << " Can't open input file " << TrimSpaces(gfilename).c_str() << ". Exiting. " << endl;
			exit(1);
		}
    
    	// For the STEER_BEAM option, we need to have the angles and vertex of the GCARD in BEAM_P and BEAM_V, SPREAD_V
    	// Getting particle name,  momentum from option value
		values       = get_info(gemcOpt->optMap["BEAM_P"].args);
		string pname = TrimSpaces(values[0]);
		mom          = get_number(values[1]);
		theta        = get_number(values[2]);
		phi          = get_number(values[3]);
    
    	// Getting vertex from option value
		values = get_info(gemcOpt->optMap["BEAM_V"].args);
		units = TrimSpaces(values[3]);
		vx = get_number(values[0] + "*" + units);
		vy = get_number(values[1] + "*" + units);
		vz = get_number(values[2] + "*" + units);
		
		// Getting vertex spread from option value
		values = get_info(gemcOpt->optMap["SPREAD_V"].args);
		units = TrimSpaces(values[2]);
		dvr = get_number(values[0] + "*" + units);
		dvz = get_number(values[1] + "*" + units);


	}
	
	
	// %%%%%%%%%%%%%%%
	// Luminosity Beam
	// %%%%%%%%%%%%%%%
	
	// Getting particle name,  momentum from option value
	values         = get_info(gemcOpt->optMap["LUMI_P"].args);
	string L_pname = TrimSpaces(values[0]);
	L_mom          = get_number(values[1]);
	L_theta        = get_number(values[2]);
	L_phi          = get_number(values[3]);
	
	
	// Getting momentum spread from option value
	values = get_info(gemcOpt->optMap["LUMI_SPREAD_P"].args);
	L_dmom   = get_number(values[0]);
	L_dtheta = get_number(values[1]);
	L_dphi   = get_number(values[2]);

	
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
	units = TrimSpaces(values[3]);
	L_vx = get_number(values[0] + "*" + units);
	L_vy = get_number(values[1] + "*" + units);
	L_vz = get_number(values[2] + "*" + units);
	L_beam_vrt = G4ThreeVector(L_vx, L_vy, L_vz);
	
	// Getting vertex spread from option value
	values = get_info(gemcOpt->optMap["LUMI_SPREAD_V"].args);
	units = TrimSpaces(values[2]);
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
	string L2_pname = TrimSpaces(values[0]);
	L2_mom          = get_number(values[1]);
	L2_theta        = get_number(values[2]);
	L2_phi          = get_number(values[3]);
	
	// Getting momentum spread from option value
	values = get_info(gemcOpt->optMap["LUMI2_SPREAD_P"].args);
	L2_dmom   = get_number(values[0]);
	L2_dtheta = get_number(values[1]);
	L2_dphi   = get_number(values[2]);

	
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
	units = TrimSpaces(values[3]);
	L2_vx = get_number(values[0] + "*" + units);
	L2_vy = get_number(values[1] + "*" + units);
	L2_vz = get_number(values[2] + "*" + units);
	L2_beam_vrt = G4ThreeVector(L2_vx, L2_vy, L2_vz);
	
	
	// Getting vertex spread from option value
	values = get_info(gemcOpt->optMap["LUMI2_SPREAD_V"].args);
	units = TrimSpaces(values[2]);
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
}










