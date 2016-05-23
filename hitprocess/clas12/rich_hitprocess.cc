// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// gemc headers
#include "rich_hitprocess.h"

map<string, double> rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
  map<string, double> dgtz;
  vector<identifier> identity = aHit->GetId();
	
  int sector = identity[0].id;
  int pmt    = identity[1].id;
  int pixel  = identity[2].id;
	
	
  // Ahmed El Alaoui, May, 2010
  dgtz["sector"] = sector;
  dgtz["pmt"]    = pmt;
  dgtz["pixel"]  = pixel;
  dgtz["hitn"]   = hitn;

  return dgtz;
}


// this routine needs to be modified
// no point drawing should be made here, but in MHit
// finding the PMT should be in a different function,
// with parameters coming from DB

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"

vector<identifier> rich_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{

  vector<identifier> id2 = id;

  G4StepPoint *prestep = aStep->GetPreStepPoint();
//  G4StepPoint *poststep = aStep->GetPostStepPoint();
  G4ThreeVector  xyz   = aStep->GetPostStepPoint()->GetPosition();
  G4ThreeVector Lxyz = prestep->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz);
  const double xpos = Lxyz.x();
  const double ypos = Lxyz.y();
  const double Energy = aStep->GetTrack()->GetTotalEnergy()/eV;
  G4ThreeVector vtx = aStep->GetTrack()->GetVertexPosition();
  int pad = id2[1].id;
  int i;
  int iPMT=-1;
  int jPMT=-1;
  double lx=0.,ly=0.;

  //cout << " Pixalization begin " << id2[0].id << " " << id2[1].id << endl;
  // H8500 specs:
  /*static const int nbins=15;
    static const double p[nbins]={1.86,   1.96,  2.03,  2.13, 2.28, 2.37, 2.48, 2.59, 2.71, 3.11, 3.37, 3.73, 4.13, 4.39, 4.64};
    static const double q[nbins]={0.0002, 0.002, 0.007, 0.02, 0.04, 0.08, 0.13, 0.17, 0.20, 0.26, 0.27, 0.28, 0.22, 0.18, 0.10};*/
  
  static const int nbins=85;
  static const double p[nbins]={
    4.624,      4.558,      4.510,      4.432,      4.372,      4.285,      4.201,      4.081,      3.980,      3.815, 
    3.726,      3.591,      3.503,      3.410,      3.314,      3.214,      3.128,      3.054,      2.983,      2.902, 
    2.845,      2.783,      2.724,      2.668,      2.624,      2.572,      2.531,      2.487,      2.454,      2.431, 
    2.408,      2.381,      2.364,      2.342,      2.325,      2.305,      2.284,      2.260,      2.240,      2.221, 
    2.202,      2.180,      2.165,      2.147,      2.136,      2.119,      2.108,      2.091,      2.084,      2.071, 
    2.058,      2.048,      2.041,      2.032,      2.022,      2.013,      2.003,      1.997,      1.991,      1.982, 
    1.975,      1.969,      1.963,      1.957,      1.951,      1.946,      1.940,      1.934,      1.931,      1.922, 
    1.919,      1.914,      1.908,      1.902,      1.897,      1.894,      1.888,      1.883,      1.877,      1.872, 
    1.867,      1.861,      1.858,      1.853,      1.848}; 
  static const double q[nbins]={
    0.11241,    0.12738,    0.14272,    0.16358,    0.18120,    0.20301,    0.22234,    0.24351,    0.25775,    0.27282, 
    0.27594,    0.27594,    0.27594,    0.27282,    0.26974,    0.26669,    0.26070,    0.25484,    0.24629,    0.23534, 
    0.22745,    0.21489,    0.20301,    0.18963,    0.17713,    0.16358,    0.15279,    0.13951,    0.12738,    0.11500, 
    0.10381,    0.09161,    0.08177,    0.07298,    0.06514,    0.05881,    0.05249,    0.04685,    0.04181,    0.03775, 
    0.03369,    0.03007,    0.02684,    0.02396,    0.02163,    0.01887,    0.01703,    0.01503,    0.01342,    0.01184, 
    0.01045,    0.00933,    0.00832,    0.00726,    0.00656,    0.00572,    0.00511,    0.00456,    0.00402,    0.00351, 
    0.00317,    0.00276,    0.00247,    0.00215,    0.00192,    0.00169,    0.00150,    0.00130,    0.00118,    0.00103, 
    0.00093,    0.00082,    0.00072,    0.00063,    0.00056,    0.00050,    0.00043,    0.00039,    0.00035,    0.00030, 
    0.00027,    0.00023,    0.00021,    0.00018,    0.00016};

  static const int nbinsUV=91;
  static const double pUV[nbinsUV]={
    6.148,      6.005,      5.841,      5.636,      5.469,      5.311,      5.161,      4.982,      4.832,      4.657, 
    4.494,      4.343,      4.201,      4.068,      3.968,      3.815,      3.726,      3.591,      3.503,      3.410, 
    3.314,      3.214,      3.128,      3.054,      2.983,      2.902,      2.845,      2.783,      2.724,      2.668, 
    2.624,      2.572,      2.531,      2.487,      2.454,      2.431,      2.408,      2.381,      2.364,      2.342, 
    2.325,      2.305,      2.284,      2.260,      2.240,      2.221,      2.202,      2.180,      2.165,      2.147, 
    2.136,      2.119,      2.108,      2.091,      2.084,      2.071,      2.058,      2.048,      2.041,      2.032, 
    2.022,      2.013,      2.003,      1.997,      1.991,      1.982,      1.975,      1.969,      1.963,      1.957, 
    1.951,      1.946,      1.940,      1.934,      1.931,      1.922,      1.919,      1.914,      1.908,      1.902, 
    1.897,      1.894,      1.888,      1.883,      1.877,      1.872,      1.867,      1.861,      1.858,      1.853, 
    1.848}; 
  static const double qUV[nbinsUV]={
    0.09697,    0.10264,    0.11370,    0.12738,    0.13951,    0.15279,    0.16545,    0.18120,    0.19399,    0.20768, 
    0.21983,    0.23268,    0.24629,    0.25775,    0.26070,    0.27282,    0.27594,    0.27594,    0.27594,    0.27282, 
    0.26974,    0.26669,    0.26070,    0.25484,    0.24629,    0.23534,    0.22745,    0.21489,    0.20301,    0.18963, 
    0.17713,    0.16358,    0.15279,    0.13951,    0.12738,    0.11500,    0.10381,    0.09161,    0.08177,    0.07298, 
    0.06514,    0.05881,    0.05249,    0.04685,    0.04181,    0.03775,    0.03369,    0.03007,    0.02684,    0.02396, 
    0.02163,    0.01887,    0.01703,    0.01503,    0.01342,    0.01184,    0.01045,    0.00933,    0.00832,    0.00726, 
    0.00656,    0.00572,    0.00511,    0.00456,    0.00402,    0.00351,    0.00317,    0.00276,    0.00247,    0.00215, 
    0.00192,    0.00169,    0.00150,    0.00130,    0.00118,    0.00103,    0.00093,    0.00082,    0.00072,    0.00063, 
    0.00056,    0.00050,    0.00043,    0.00039,    0.00035,    0.00030,    0.00027,    0.00023,    0.00021,    0.00018, 
    0.00016};

//  static const int nPMTs=28;
  /*static const double effPMT[nPMTs]={
    1.26314, 1.14882,  0.897036, 1.01661,  1.21187, 0.974241, 0.82433, 0.919741, 1.05825, 1.13679, 
    1.06701, 1.04444,  1.23772,  1.19244,  1.19432, 1.05544,  1.09541, 1.00884,  1.18319, 0.870763,
    1.01581, 0.974612, 0.841163, 0.922537, 1.10571, 0.905572, 1.12075, 1.16352};*/
//  static const double effPMT[nPMTs]={
//    1.24404, 1.15731,  0.891041, 1.02305,  1.16932, 0.970796, 0.779371, 0.936425, 1.06004, 1.11025,
//    1.0508,  1.06841,  1.23461,  1.19373,  1.18116, 1.08872,  1.06598,  1.02615,  1.18402, 0.891628,
//    1.04522, 0.986448, 0.856585, 0.9269,   1.0875,  0.918675, 1.134,    1.21063};
	
  static const int    NIPXL   = 8 ;
  static const int    NJPXL   = 8 ;

  static const int    n_raws  = 8 ;
  static const int    n_cols  = 8 ;

  static const double PXDX[8] = { 6.26, 6.08, 6.08, 6.08, 6.08, 6.08, 6.08, 6.26 };
  static const double PXDY[8] = { 6.26, 6.08, 6.08, 6.08, 6.08, 6.08, 6.08, 6.26 };

  // std::vector<double> v_pixel_raw { 6.26, 6.08, 6.08, 6.08, 6.08, 6.08, 6.08, 6.26 }; 
  // std::vector<double> v_pixel_col { 6.26, 6.08, 6.08, 6.08, 6.08, 6.08, 6.08, 6.26 }; 

  static const double PhCsize_x( 6.26*2. + 6.08*6. );
  static const double PhCsize_y( 6.26*2. + 6.08*6. );

  static const double PXLDeadSpace_x(0.28);
  static const double PXLDeadSpace_y(0.28);

 // static const double PXLsize_x( ( PhCsize_x -7*PXLDeadSpace_x )/8 );
//  static const double PXLsize_y( ( PhCsize_y -7*PXLDeadSpace_y )/8 );
  
  //MC x axis is mirrored with respect the Lab reference system
  //shift the local reference center to up left corner (following configuration files)
  const double new_xpos = -xpos + PhCsize_x/2. ;
  const double new_ypos =  ypos + PhCsize_y/2. ;

  /* boolean to check if the photon hit the deadspace in y or x */
  bool _iDeadSpace = false ;

  //find the pixel row
  int _iraw = -1 ;

  double_t _x_progress = 0.0 ;
  for( int _ir = 0 ; _ir < n_raws; _ir++ ){

    if( new_xpos > _x_progress && new_xpos < _x_progress + PXDX[ _ir ] ){

      // cout << " x new pos is " << new_xpos << " and x_progress is " << _x_progress << " for ir = " << _ir << endl ;

      /* position of the hit in the pixel */
      double _x_pos_in_the_pixel = new_xpos - _x_progress ;

      _iraw = _ir + 1 ;

      /* verify it is not in the deadspace */
      if( ( _x_pos_in_the_pixel > ( PXDX[ _ir ] - PXLDeadSpace_x/2.) ) || ( _x_pos_in_the_pixel <  PXLDeadSpace_x/2. ) ) _iDeadSpace = true ;// _iraw = 0 ;

      break ;

    }

    _x_progress += PXDX[ _ir ] ;

  }


  //find the pixel column
  int _icol = -1 ;

  double_t _y_progress = 0.0 ;
  for( int _ic = 0 ; _ic < n_cols; _ic++ ){

    if( new_ypos > _y_progress && new_ypos < _y_progress + PXDY[ _ic ] ){

      // cout << " y new pos is " << new_ypos << " and y_progress is " << _y_progress << " for ir = " << _ic << endl ;

      /* position of the hit in the pixel */
      double _y_pos_in_the_pixel = new_ypos - _y_progress ;

      _icol = 8 - _ic ; 

      /* verify it is not in the deadspace */
      if( ( _y_pos_in_the_pixel > ( PXDY[ _ic ] - PXLDeadSpace_y/2.) ) || ( _y_pos_in_the_pixel <  PXLDeadSpace_y/2. ) ) _iDeadSpace = true ;// _iraw = 0 ;

      break ;

    }

    _y_progress += PXDY[ _ic ] ;

  }


  double ini=0.;
  //cout << " X --> " << new_xpos << endl;
  for (i=0; i<NIPXL; i++){
    if(new_xpos>ini && new_xpos<=ini+PXDX[i]){
      lx = new_xpos - ini - PXDX[i]/2.;
      if(abs(lx) < (PXDX[i]-PXLDeadSpace_x)/2.){
	iPMT=i+1;
	// cout << " 2 pixel x is " << iPMT << endl ;

      }else{
	iPMT=0;
      }
      //cout << "       --> " << i << " " << ini << " lx " << lx << " " << iPMT << endl;
      break;
    }
    ini+=PXDX[i];
  }

  //find the pixel column
  ini=0.;
  //cout << " Y --> " << new_ypos << endl;
  for (i=0; i<NJPXL; i++){
    if(new_ypos>ini && new_ypos<=ini+PXDY[i]){
      ly = new_ypos - ini - PXDY[i]/2.;
      if(abs(ly) < (PXDY[i]-PXLDeadSpace_y)/2.){
	jPMT=NJPXL-i;
	// cout << " 2 pixel y is " << jPMT << endl ;
      }else{
	jPMT=0;
      }
      //cout << "         --> " << i << " " << ini << " ly " << ly << " " << jPMT << endl;
      break;
    }
    ini+=PXDY[i];
  }

  // find the pixel number
//  int ijPMT2 = NIPXL*(jPMT-1) + iPMT;

  /* final pixel number */
  int ijPMT = NIPXL*( _icol - 1 ) + _iraw ;
  /* attribute a negative pixel number if the deadspace was hit */
  if( _iDeadSpace ) ijPMT = - ijPMT ;

  // if( ijPMT != ijPMT2 ) cout << "Pixel number is " << ijPMT<< "\t" << ijPMT2<< "\t" << iPMT<< "\t" << jPMT<< "\t" << _iraw<< "\t" << _icol << endl ;

  /* check now the quantum efficiency */
  if( ijPMT > 0 && ijPMT < 65 )
    {

      // cout << ijPMT << endl ;

      /* start checking the photon range - if outside the range, pixel number is set to 0 */
      if(Energy < pUV[0] && Energy >= pUV[ nbinsUV -1]){
	/* it was as below, but it seems to be wrong: && instead of || */
	// if(Energy < pUV[0] && Energy >= pUV[nbinsUV-1])
	double QE = -1;
	double QEUV = -1;
	/* standard - non-UV - quantum efficiency */
	for(int ie=0; ie<nbins-1; ie++)
	  {
	    if(Energy < p[ie] && Energy >= p[ie+1])
	      {
		QE = q[ie] + (Energy-p[ie])*(q[ie+1]-q[ie])/(p[ie+1]-p[ie]);
		break;
	      }
	  }
	/* UV quantum efficiency */
	for(int ie=0; ie<nbinsUV-1; ie++)
	  {
	    if(Energy < pUV[ie] && Energy >= pUV[ie+1])
	      {
		QEUV = qUV[ie] + (Energy-pUV[ie])*(qUV[ie+1]-qUV[ie])/(pUV[ie+1]-pUV[ie]);
		break;
	      }
	  }
	// pixel# += 100 if fails quauntum efficiency:
	// pixel# += 200 if fails extended quauntum efficiency:
	double QEtest = G4UniformRand(); 
	//if(G4UniformRand()>QE)
	//cout << " HitE " << Energy << " QE " << QE << " QEtest " << QEtest << endl;

	if(pad % 2 == 0){
	  //cout << " Hit Pad " << pad << " even QEUV " << endl;
	  //if(QEtest>QE*effPMT[pad-1])
	  if(QEtest>QE)
	    {
	      ijPMT += 100;
	      //if(QEtest>QEUV) ijPMT += 100;
	      //cout << "  --> Failed QE " << ijPMT << " E " << Energy << " QE " << QE << " vs " << QEtest << endl;
            }
	} else {
	  //cout << " Hit Pad " << pad << " odd QE" << endl;
	  //if(QEtest>QEUV*effPMT[pad-1])
	  if(QEtest>QE)
	    {
	      ijPMT += 100;
	      //if(QEtest>QEUV) ijPMT += 100;
	      //cout << "  --> Failed QE " << ijPMT << " E " << Energy << " QE " << QE << " vs " << QEtest << endl;
            }
	}
	/* if photon energy is outside the range */
      } else {
	//cout << "  --> Outside QE " << ijPMT << " E " << Energy << endl;
	ijPMT += 200;
      }
    }
  
    // G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    // if(pVVisManager)
    //   {

    // 	G4Circle circle(xyz);
    // 	circle.SetFillStyle(G4Circle::filled);
    // 	circle.SetScreenSize(7);

    // 	/* red hits are the good ones */
    // 	if (ijPMT > 0 && ijPMT < 65 ){
    // 	  G4Colour colour_touch (1.0, 0.0, 0.0);
    // 	  circle.SetVisAttributes(G4VisAttributes(colour_touch));
    // 	} else if ( ijPMT > 100 && ijPMT < 200 ){
    // 	  /* draw in green pixels corresponding to photons that didn't pass the quantum efficiency*/
    // 	  G4Colour colour_touch (0.0, 1.0, 0.0);
    // 	  circle.SetVisAttributes(G4VisAttributes(colour_touch));
    // 	  /* draw in gray pixels corresponding to photons hitting dead spaces */
    // 	} else if( ijPMT < 0 ){
    // 	  G4Colour colour_touch ( 0.5, 0.5, 0.5 );
    // 	  circle.SetVisAttributes(G4VisAttributes(colour_touch));
    // 	} else {

    // 	  G4Colour colour_touch ( 0.0, 0.0, 1.0 );
    // 	  circle.SetVisAttributes(G4VisAttributes(colour_touch));

    // 	}

    // 	pVVisManager->Draw(circle);

    //   }

    id2[2].id = ijPMT;
    //cout << " Pixalization end   " << id2[0].id << " " << id2[1].id << " " << id2[2].id << endl;
    return id2;

      }



  map< string, vector <int> >  rich_HitProcess :: multiDgt(MHit* aHit, int hitn)
  {
    map< string, vector <int> > MH;
	
    return MH;
  }


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> rich_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}









