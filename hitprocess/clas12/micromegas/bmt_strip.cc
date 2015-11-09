// gemc headers
#include "bmt_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>
#include <cstdlib>

void bmt_strip::fill_infos()
{
 // all dimensions are in mm
    Pi  = 3.14159265358;

    interlayer         = 15.0;
    Nsector            = 3;
    hDrift         = 3.0;
    hStrip2Det         = hDrift/2.;
    sigma_td_max   = 0.4;
    theta_L        = Pi*20./180.;
    w_i            = 25.0;


    interStripZ = 0.20;
    interStripC = 0.16;

 // pitch CR4Z -> CR4 = central tracker region 4
    pitchZ4 = 0.487; // and interpitch=0.201
 
 // pitch CR5Z -> CR5 = central tracker region 5
    pitchZ5 = 0.49; // none existing value. picked a rabdom value compatible with the geometry
 
 // pitch CR6Z -> CR6 = central tracker region 6
    double strip_widthZ = 0.328;
    pitchZ6 = interStripZ + strip_widthZ;
 
 // pitch CRC --> for the C-detectors, the strips are in bunches of equal pitches
    double CR4C_width[13]={0.345,0.28,0.225,0.175,0.17,0.21,0.26,0.31,0.37,0.44,0.515,0.605,0.7};     // width of the corresponding group of strips
    double CR4C_group[13]={32,32,32,32,624,32,32,32,32,32,32,32,32,896};                              // the number of strips with equal pitches
    // For CR5, no existing value. picked a random value compatible with the geometry
    double CR5C_width[1]={0.253};                                                                     // width of the corresponding group of strips
    double CR5C_group[1]={1024};                                                                      // the number of strips with equal pitches
    double CR6C_width[14]={0.38,0.32,0.27,0.23,0.17,0.18,0.22,0.25,0.29,0.33,0.37,0.41,0.46,0.51};    // width of the corresponding group of strips
    double CR6C_group[14]={32,32,32,32,704,64,32,32,32,32,32,32,32,32};                               // the number of strips with equal pitches
    int arraySize[] = {sizeof(CR4C_group) / sizeof(CR4C_group[0]),sizeof(CR5C_group) / sizeof(CR5C_group[0]),1sizeof(CR6C_group) / sizeof(CR6C_group[0])}; // size of arrays describing the C-detectors for a given MM region

    // pitch CR4C for first region of macromegas (CRC4)
    for(int i =0; i<sizeof(CR4C_group) / sizeof(CR4C_group[0]); i++)
        pitchC4.push_back(interStripC+CR4C_width[i]);
 
    // pitch CR5C for second region of macromegas (CRC5)
        for(int i =0; i<sizeof(CR5C_group) / sizeof(CR5C_group[0]); i++)
        pitchC5.push_back(interStripC+CR5C_width[i]);

    // pitch CR6C for third region of macromegas (CRC6)
    for(int i =0; i<sizeof(CR6C_group) / sizeof(CR6C_group[0]); i++)
        pitchC6.push_back(interStripC+CR6C_width[i]);

    // z of the upstream part of the layer
    Z0.push_back(-127.820); Z0.push_back(-127.820);
    Z0.push_back(-148.740); Z0.push_back(-148.740);
    Z0.push_back(-169.710); Z0.push_back(-169.710);
 
    // total z length of the layer of the layer
    DZ.push_back(372.75); DZ.push_back(372.75);
    DZ.push_back(423.99); DZ.push_back(423.99);
    DZ.push_back(444.96); DZ.push_back(444.96);

    // radii of layers
    R.push_back(145.731); R.push_back(R[0]+interlayer);
    R.push_back(175.731); R.push_back(R[2]+interlayer);
    R.push_back(205.731); R.push_back(R[4]+interlayer);

    // Number of strips (depends on radius, dead zones, and pitch!)
    // layers 0,2,4 are Z detectors, i.e. measure phi (strips in z direction)
    // layers 1,3,5 are C detectors, i.e. measure z (phi strips at z)
    Nstrips.push_back(14*64);         // CR4Z
    Nstrips.push_back(10*64);         // CR4C
    Nstrips.push_back(11*64);         // CR5Z
    Nstrips.push_back(16*64);        // CR5C
    Nstrips.push_back(768);          // CR6Z
    Nstrips.push_back(1152);         // CR6C

 
    // dead angles
    Inactivtheta.push_back((20/R[0])*(180./Pi));
    Inactivtheta.push_back((20/R[1])*(180./Pi));
    Inactivtheta.push_back((20/R[2])*(180./Pi));
    Inactivtheta.push_back((20/R[3])*(180./Pi));
    Inactivtheta.push_back((20/R[4])*(180./Pi));
    Inactivtheta.push_back((20/R[5])*(180./Pi));
 
    // dead zones
    DZ_inLength = 0; //for CRnZ
    DZ_inWidth  = 0; //for CRnC
    DZ4_inLength = ((R[1]*(120-Inactivtheta[1])*Pi/180)-(Nstrips[1]*pitchZ4+interStripZ))/2; //for CR4Z (((160.731*(120-7.129)*Pi/180)-(10*64*0.487)+0.201))/2 = 2.377
    DZ4_inWidth  = (DZ[0]-372.48-interStripC)/2; //for CR4C // (372.75 - 372.48 - 0.16)/2 = 0.055
    DZ5_inLength = ((R[2]*(120-Inactivtheta[2])*Pi/180)-(Nstrips[2]*pitchZ4+interStripZ))/2; //for CR5Z (((175.731*(120-6.521)*Pi/180)-(11*64*0.49)+0.201))/2 = 1.444
    DZ5_inWidth  = (DZ[3]-422.912-interStripC)/2; //for CR5C // (423.99 - 422.912 - 0.16)/2 = 0.459
    DZ6_inLength = ((R[4]*(120-Inactivtheta[4])*Pi/180)-(Nstrips[4]*pitchZ4+interStripZ))/2; //for CR6Z (((205.731*(120-5.57)*Pi/180)-(12*64*0.526)+0.201))/2 = 1.779
    DZ6_inWidth  = (DZ[5]-444.8-interStripC)/2; //for CR6C // (444.96 - 44.8 - 0.16)/2 = 0
 
    // mid angle of the sector
    MidTile.push_back(0);     MidTile.push_back(0);
    MidTile.push_back(0);     MidTile.push_back(0);
    MidTile.push_back(0);     MidTile.push_back(0);
}
/**
 * Method to get digi hits based on x,y,z position and E
 * layers from 0 to 5
 *   layers 0,2,4 are Z detectors (inner)
 *   layers 1,3,5 are C detectors (outer)
 */
vector<double>  bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
        // the return vector is always in pairs.
        // The first number is the ID,
        // the second number is the sharing percentage
        vector<double> strip_id;

        // dead zones
        if(layer == 0 || layer == 1){
                DZ_inLength = DZ4_inLength;
                DZ_inWidth  = DZ4_inWidth;
        }
        else if(layer == 2 || layer == 3){
                DZ_inLength = DZ5_inLength;
                DZ_inWidth  = DZ5_inWidth;
        }
        else if(layer == 4 || layer == 5){
                DZ_inLength = DZ6_inLength;
                DZ_inWidth  = DZ6_inWidth;
        }

        // number of electrons (Nt)
        Nel = (int) (1e6*Edep/w_i);

        // 1st define phi of the mean hit point
        double phi;
        if(x>0 && y>=0) phi = atan(y/x);
        else if(x>0 && y<0) phi = 2.*Pi+atan(y/x);
        else if(x<0) phi = Pi+atan(y/x);
        else if(x==0 && y>0) phi = Pi/2.;
        else if(x==0 && y<0) phi = 3.*Pi/2.;
        else phi = 0; // x = y = 0, phi not defined

        // now find the tile number (transverse diffusion will not change that)
        int ti=0;
        double theta_tmp=0;
        for(int t=0; t<Nsector; t++)
        {
                theta_tmp = MidTile[layer]+2.*t*Pi/Nsector;
                if(theta_tmp>2.*Pi) theta_tmp = theta_tmp - 2.*Pi;
                if(theta_tmp<0) theta_tmp = theta_tmp + 2.*Pi;
                if(fabs(phi-theta_tmp)<Pi/Nsector || fabs(2.*Pi-fabs(phi-theta_tmp))<Pi/Nsector) ti=t; // gives tile #
        }

        double phiij = MidTile[layer]+2.*ti*Pi/Nsector;
        if(phi-phiij<=-Pi/Nsector) phi = phi+2.*Pi;
        if(phi-phiij>Pi/Nsector) phi = phi-2.*Pi;
        if(phi-phiij<=-Pi/Nsector || phi-phiij>Pi/Nsector) cout << "WARNING: incorrect phi value in BMT: " << phi*180./Pi << " vs " << phiij*180./Pi << endl;
        int ClosestStrip=0;

        // now compute the sigma of the (transverse) dispersion for this interaction
        if(layer==1 || layer==3 || layer==5) sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/hDrift); // "C" det, transverse diffusion grows with square root of distance
        else sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/(cos(theta_L)*hDrift)); // same, but "Z" detectors, so Lorentz angle makes drift distance longer by 1./cos(theta_L) . Means sigma_td can be larger than sigma_td_max


        if(Nel>0)
        {
                for(int iel=0;iel<Nel;iel++)
                { // loop over (total) electrons
                // select the equal-size-pitch groups based on the layer

                        if(layer%2==1)
                        { //  for "C" layers, i.e. measuring z

                        Int_t CRC_group[arraySize[(layer-1)/2]];
                        Double_t CRC_width[arraySize[(layer-1)/2]];
                        if(layer==1) {
                                CRC_group[] = CR4C_group;
                                CRC_width[] = CR4C_width;
                        }
                        if(layer==3) {
                                CRC_group[] = CR5C_group;
                                CRC_width[] = CR5C_width;
                        }
                        if(layer==5) {
                                CRC_group[] = CR6C_group;
                                CRC_width[] = CR6C_width;
                        }

                        // the real z is generated using a gaussian number generator with sigma corresponding to the sigma of the dispersion
                        z_real = (double) (G4RandGauss::shoot(z,sigma_td));
                        // if this z is within the first group of equal-pitch strips then get its closest strip
                        if(z_real-Z0[layer]-DZ_inWidth>0 && z_real-Z0[layer]-DZ_inWidth<CRC_group[0]*(interStripC+CRC_width[0]))
                                ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/(interStripC+CRC_width[0])))+0.5);    // the closest strip
                        // is the effective z divided by the pitch +0.5
                        // if there is only one group -- stop here
                        // otherwise loop over the number of groups of equal-pitch strips and find the corresponding nearest strip from z
                        int nstripsPrevGrp =0;
                        if(arraySize[(layer-1)/2]>1)
                                for(int i =1; i<arraySize[(layer-1)/2]; i++) {
                                nstripsPrevGrp+=i*CRC_group[i];
                                if(z_real-Z0[layer]-DZ_inWidth>CRC_group[i-1]*(interStripC+CRC_width[i-1]) && z_real-Z0[layer]-DZ_inWidth<CRC_group[i]*(interStripC+CRC_width[i]))
                                        ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/(interStripC+CRC_width[i])))+0.5+nstripsPrevGrp);
                                }
                        }

                        //  for "Z" layers
                        int pitchZ;
                        // select the right pitch value based on the layer
                        if(layer==0)
                        pitchZ = pitchZ4;
                        if(layer==2)
                        pitchZ = pitchZ5;
                        if(layer==2)
                        pitchZ = pitchZ6;

                        if(layer%2==0)
                        { //  for "Z" layers, i.e. measuring phi
                                phi_real = phi + ((G4RandGauss::shoot(0,sigma_td))/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer]; // the sign of the 2nd term (Lorentz angle) should be a "-" as the B field is along +z and the MM are convex (and electrons are negatively charged)
                                ClosestStrip = (int) (floor(((R[layer]/pitchZ)*(phi_real-phiij+Pi/Nsector - (Inactivtheta[layer]/2.)*Pi/180. - DZ_inLength/R[layer]))+0.5));
                        }

                        // acceptance criteria
                        if(ClosestStrip>=0 && ClosestStrip<=Nstrips[layer] && z>=Z0[layer]+DZ_inWidth && z<=Z0[layer]+DZ[layer]-DZ_inWidth && phi>=phiij-Pi/Nsector+DZ_inLength/R[layer] && phi<=phiij+Pi/Nsector-DZ_inLength/R[layer])
                        { // strip is in the acceptance, check if new strip or not
                                for(int istrip=0;istrip< (int) (strip_id.size()/2);istrip++)
                                {
                                        if(strip_id[2*istrip]==ClosestStrip)
                                        {// already hit strip - add Edep
                                                strip_id[2*istrip+1]=strip_id[2*istrip+1]+1./((double) Nel); // no gain fluctuation yet
                                                ClosestStrip=-1; // not to use it anymore
                                        }
                                }
                                if(ClosestStrip>-1)
                                { // this is a new strip
                                        strip_id.push_back(ClosestStrip);
                                        strip_id.push_back(1./((double) Nel)); // no gain fluctuation yet
                                }
                        }
                        else
                        { // not in the acceptance
                                strip_id.push_back(-1);
                                strip_id.push_back(1);
                        }
                }
        }
        else
        { // Nel=0, consider the Edep is 0
                strip_id.push_back(-1);
                strip_id.push_back(1);
        }
        return strip_id;
}
