// gemc headers
#include "HitProcess.h"
#include "HitProcess_MapRegister.h"
#include "flux_hitprocess.h"         ///< flux hit process common to all

// CLAS12
#include "clas12/bonus_hitprocess.h"          ///< Barrel Silicon Tracker (bst)
#include "clas12/svt/bst_hitprocess.h"          ///< Barrel Silicon Tracker (bst)
#include "clas12/cnd_hitprocess.h"              ///< Central Neutron Detector
#include "clas12/ctof_hitprocess.h"             ///< Central TOF
#include "clas12/dc_hitprocess.h"               ///< Drift Chambers
#include "clas12/ec_hitprocess.h"               ///< Forward Electromagnetic Calorimeter EC
#include "clas12/ftof_hitprocess.h"             ///< Forward TOF
#include "clas12/ft_cal_hitprocess.h"           ///< Forward Tagger Calorimeter
#include "clas12/ft_hodo_hitprocess.h"          ///< Forward Tagger Hodoscope
#include "clas12/micromegas/ftm_hitprocess.h"   ///< Forward Tagger Micromegas
#include "clas12/htcc_hitprocess.h"             ///< High Threshold Cherenkov Counter
#include "clas12/ltcc_hitprocess.h"             ///< Low Threshold Cherenkov Counter
#include "clas12/micromegas/FMT_hitprocess.h"   ///< forward micromegas
#include "clas12/micromegas/BMT_hitprocess.h"   ///< barrel micromegas
#include "clas12/pcal_hitprocess.h"             ///< Pre-shower calorimeter
#include "clas12/rich_hitprocess.h"             ///< Pre-shower calorimeter

// Beam Dump eXperiment
#include "bdx/cormo_hitprocess.h"               ///< Cormorino detector
#include "bdx/veto_hitprocess.h"                ///< Cormorino vetos
#include "bdx/crs_hitprocess.h"                ///< Cormorino vetos

// APrime
#include "HPS/ECAL_hitprocess.h"       ///< Calorimeter Crystals
#include "HPS/SVT_hitprocess.h"        ///< Silicon Vertex Trackers.
#include "HPS/muon_hodo_hitprocess.h"  ///< HPS Muon Hodoscopes

// Cebaf
#include "injector/bubble_hitprocess.h"       ///< Calorimeter Crystals


map<string, HitProcess_Factory> HitProcess_Map(string experiments)
{
	
	map<string, HitProcess_Factory> hitMap;
	
	stringstream exps(experiments);
	string EXP;
	
	cout << endl;
	while(!exps.eof())
	{
		exps >> EXP;
		cout << "  >> Registering experiment \"" << EXP << "\" hit processes " << endl;
		
		// flux is independent of experiment
		hitMap["flux"]           = &flux_HitProcess::createHitClass;
		// mirror is also a flux detector
		hitMap["mirror"]         = &flux_HitProcess::createHitClass;

		// CLAS12
		if(EXP == "clas12")
		{
			hitMap["bonus"]    = &bonus_HitProcess::createHitClass;
			hitMap["bst"]      = &bst_HitProcess::createHitClass;
			hitMap["cnd"]      = &cnd_HitProcess::createHitClass;
			hitMap["ctof"]     = &ctof_HitProcess::createHitClass;
			hitMap["dc"]       = &dc_HitProcess::createHitClass;
			hitMap["ec"]       = &ec_HitProcess::createHitClass;
			hitMap["ecs"]      = &ec_HitProcess::createHitClass;
			hitMap["ftof_p1a"] = &ftof_HitProcess::createHitClass;
			hitMap["ftof_p1b"] = &ftof_HitProcess::createHitClass;
			hitMap["ftof_p2"]  = &ftof_HitProcess::createHitClass;
			hitMap["ft_cal"]   = &ft_cal_HitProcess::createHitClass;
			hitMap["ft_hodo"]  = &ft_hodo_HitProcess::createHitClass;
			hitMap["htcc"]     = &htcc_HitProcess::createHitClass;
			hitMap["ltcc"]     = &ltcc_HitProcess::createHitClass;
			hitMap["fmt"]      = &FMT_HitProcess::createHitClass;
			hitMap["bmt"]      = &BMT_HitProcess::createHitClass;
			hitMap["ftm"]      = &FTM_HitProcess::createHitClass;
			hitMap["pcal"]     = &pcal_HitProcess::createHitClass;
			hitMap["rich"]     = &rich_HitProcess::createHitClass;
		}
		// Aprime
		else if(EXP == "HPS")
		{
			hitMap["SVT"]        = &SVT_HitProcess::createHitClass;
			hitMap["ECAL"]       = &ECAL_HitProcess::createHitClass;
			hitMap["muon_hodo"]  = &muon_hodo_HitProcess::createHitClass;
		}
		// GlueX
		else if( EXP == "gluex" )
		{
		}
		// SoLID
		else if( EXP == "solid" )
		{
		}
		else if( EXP == "BDX" )
		{
			hitMap["cormo"]      = &cormo_HitProcess::createHitClass;
			hitMap["veto"]       = &veto_HitProcess::createHitClass;
			hitMap["crs"]        = &crs_HitProcess::createHitClass;
		}
		else if( EXP == "injector" )
		{
			hitMap["bubble"]     = &bubble_HitProcess::createHitClass;
		}
		
	}
	
	cout << endl;
	return hitMap;
	
}
