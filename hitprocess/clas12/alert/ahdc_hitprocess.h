#ifndef AHDC_HITPROCESS_H
#define AHDC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


class ahdcConstants
{
public:
	
	// Database parameters
	int    runNo;
	string date;
	string connection;
	char   database[80];
	
	// translation table
	TranslationTable TT;
};


// Class definition
/// \class ahdc_HitProcess
/// <b> Alert Drift Chamber Hit Process Routine</b>\n\n

class ahdc_HitProcess : public HitProcess
{
public:
	
	~ahdc_HitProcess(){;}
	
	// constants initialized with initWithRunNumber
	static ahdcConstants atc;
	
	void initWithRunNumber(int runno);
	
	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);
	
	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);
	
	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);
	
	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
	// creates the HitProcess
	static HitProcess *createHitClass() {return new ahdc_HitProcess;}
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
public:
	// AHDC geometry parameters
	float PAD_W, PAD_L, PAD_S, RTPC_L;
	float phi_per_pad;
	
	// parameters for drift and diffustion equations for drift time, 
	// drift angle, and diffusion in z
	float a_t, b_t, c_t, d_t;
	float a_phi, b_phi, c_phi, d_phi;
	float a_z, b_z;
	
	// variables for storing drift times and diffusion in time
	float t_2GEM2, t_2GEM3, t_2PAD, t_2END;
	float sigma_t_2GEM2, sigma_t_2GEM3, sigma_t_2PAD, sigma_t_gap;
	
	// variables for storing drift angle and diffusion in phi
	float phi_2GEM2, phi_2GEM3, phi_2PAD, phi_2END;
	float sigma_phi_2GEM2, sigma_phi_2GEM3, sigma_phi_2PAD, sigma_phi_gap;
	
	float z_cm;
	float TPC_TZERO;
	
	map<int, double> timeShift_map;
	double shift_t;
	
};

#include <string>
#include "CLHEP/GenericFunctions/Landau.hh"

/**
 * @class ahdcSignal
 * 
 * @brief ahdc signal simulation
 *
 * This class simulates the waveform of the ahdc signal and provide 
 * algorithms to extract relevant informations from this signal.
 *
 * @author Felix Touchte Codjo
 */
class ahdcSignal {
	// MHit or wires identifiers
	public : 
		int hitn; ///< n-th MHit of the event, also corresponds to the n-th activated wire
		int sector; ///< sector, first wire identifier
		int layer; ///< layer, second wire identifer
		int component; ///< component, third wire identifier
		int nsteps; ///< number of steps in this MHit, i.e number of Geant4 calculation steps in the sensitive area of the wire 
	// vectors
	private :
		std::vector<double> Edep; ///< array of deposited energy in each step [keV]
		std::vector<double> G4Time; ///< array of Geant4 time corresponding to each step [ns]
		std::vector<double> Doca; ///< array of distance of closest approach corresponding each step [mm]
		std::vector<double> DriftTime; ///< array of drift time corresponding each step [ns]
		
		/**
		 * @brief Fill the arrays Doca and DriftTime
		 * 
		 * Compute the doca corresponding to each step and
		 * deducte the driftime using a "time to distance"
		 * relation
		 *
		 * @param aHit an object derived from Geant4 "G4VHit" class
		 */
		void ComputeDocaAndTime(MHit * aHit);
		std::vector<short> Dgtz; ///< Array containing the samples of the simulated signal
		std::vector<short> Noise; ///< Array containing the samples of the simulated noise
	// setting parameters for digitization
	private : 
		const double tmin; ///< lower limit of the simulated time window
		const double tmax; ///< upper limit of the simulated time window
		const double timeOffset; ///< time offset for simulation purpose
		const double samplingTime; ///< sampling time [ns]
		const double Landau_width; ///< Width pararemeter of the Landau distribution
		double electronYield = 9500;   ///< ADC gain
		static const int ADC_LIMIT = 4095; ///< ADC limit, corresponds to 12 digits : 2^12-1
	// public methods
	public :
		/** @brief Default constructor */
		ahdcSignal() = default;
		
		/** @brief Constructor */
		ahdcSignal(MHit * aHit, int _hitn, double _tmin, double _tmax, double _timeOffset, double _samplingTime, double _Landau_width) 
		: tmin(_tmin), tmax(_tmax), timeOffset(_timeOffset), samplingTime(_samplingTime), Landau_width(_Landau_width) {
			// read identifiers
			hitn = _hitn;
			vector<identifier> identity = aHit->GetId();
			sector = 0;
			layer = 10 * identity[0].id + identity[1].id ; // 10*superlayer + layer
			component = identity[2].id;
			// fill vectors
			Edep = aHit->GetEdep();
			nsteps = Edep.size();
			for (int s=0;s<nsteps;s++){Edep.at(s) = Edep.at(s)*1000;} // convert MeV to keV
			G4Time = aHit->GetTime();
			this->ComputeDocaAndTime(aHit); // fills Doca and DriftTime
		}
		
		/** @brief Destructor */
		~ahdcSignal(){;}
		
		/** @brief Return the value of the attribut `electronYield` */
		double GetElectronYield() {return electronYield;}
		
		/** @brief Return the content of the attribut `Edep` */
		std::vector<double>                     GetEdep() 		{return Edep;}
		
		/** @brief Return the content of the attribut `G4Time` */
		std::vector<double>                     GetG4Time()		{return G4Time;}
		
		/** @brief Return the content of the attribut `Doca` */
		std::vector<double>                     GetDoca()		{return Doca;}
		
		/** @brief Return the content of the attribut `DriftTime` */
		std::vector<double>                     GetDriftTime()		{return DriftTime;}
		
		/** @brief Return the content of the attribut `Noise` */
		std::vector<short>                     GetNoise()              {return Noise;}
		
		/** @brief Return the content of the attribut `Dgtz` */
		std::vector<short> 			GetDgtz()		{return Dgtz;}
		
		/**
		 * @brief Set the electron yield. 
		 * 
		 * Only relevant before the use of the method `Digitize`
		 */
		void SetElectronYield(double electronYield_)		{electronYield = electronYield_;}
		
		/**
		 * @brief Overloaded `()` operator to get the value of the signal at a given time.
		 * 
		 * @param t Time at which to calculate the signal's value
		 *
		 * @return Value of the signal at the time `t`
		 */
		double operator()(double timePoint){
			using namespace Genfun;
			double signalValue = 0;
			for (int s=0; s<nsteps; s++){
				// setting Landau's parameters
				Landau L;
				L.peak() = Parameter("Peak",DriftTime.at(s),tmin,tmax); 
				L.width() = Parameter("Width",Landau_width,0,400); 
				signalValue += Edep.at(s)*L(timePoint-timeOffset);
			}
			return signalValue;
		}
		
		/**
		 * @brief Digitize the simulated signal
		 *
		 * This method perfoms several steps
		 * - step 1 : it produces samples from the simulated signal (using `samplingTime`)
		 * - step 2 : it converts keV/ns in ADC units (using `electronYield`)
		 * - step 3 : it adds noise
		 *
		 * @return Fill the attributs `Dgtz` and `Noise` 
		 */
		void Digitize();

		/** @brief Generate gaussian noise
		 *
		 * @return Fill the attribut `Noise`
		 */
		void GenerateNoise(double mean, double stdev);
		
		double GetMCTime(); // tmp
		double GetMCEtot(); // tmp
		
		/**
		 * @brief Extract various informations from the digitized signal 
		 *
		 * This method computes
		 * - `binMax`, `binOffset`, `adcMax`, `timeMax`, `integral` 
		 * - `timeRiseCFA`, `timeFallCFA`, `timeOverThresholdCFA`, `timeCFD` 
		 */
		std::map<std::string,double> Extract();
};


/**
 * @class ahdcExtractor
 * 
 * @author Felix Touchte Codjo
 */
class ahdcExtractor {
	public :
		float samplingTime; ///< time between two ADC bins
		int sparseSample = 0; ///< used to defined binOffset
		short adcOffset = 0; ///< pedestal or noise level
		long timeStamp = 0; ///< timeStamp timing informations (used to make fine corrections)
		float fineTimeStampResolution = 0; ///< precision of dream clock (usually 8) 
		static const short ADC_LIMIT = 4095; ///< Maximum value of ADC : 2^12-1
		float amplitudeFractionCFA; ///< amplitude fraction between 0 and 1
		int binDelayCFD; ///< CFD delay parameter
		float fractionCFD; ///< CFD fraction parameter between 0 and 1
	
	public :
		int binMax; ///< Bin of the max ADC over the pulse
                int binOffset; ///< Offset due to sparse sample
                float adcMax; ///< Max value of ADC over the pulse (fitted)
                float timeMax; ///< Time of the max ADC over the pulse (fitted)
                float integral; ///< Sum of ADCs over the pulse (not fitted)

		std::vector<short> samplesCorr; ///< Waveform after offset (pedestal) correction
                int binNumber; ///< Number of bins in one waveform

                float timeRiseCFA; ///< moment when the signal reaches a Constant Fraction of its Amplitude uphill (fitted)
                float timeFallCFA; ///< moment when the signal reaches a Constant Fraction of its Amplitude downhill (fitted)
                float timeOverThresholdCFA; ///< is equal to (timeFallCFA - timeRiseCFA)
                float timeCFD; ///< time extracted using the Constant Fraction Discriminator (CFD) algorithm (fitted)

		/** @brief Default constructor */
		ahdcExtractor() = default;

		/** @brief Constructor */
		ahdcExtractor(float _samplingTime,float _amplitudeFractionCFA, int _binDelayCFD, float _fractionCFD) :
			samplingTime(_samplingTime), amplitudeFractionCFA(_amplitudeFractionCFA), binDelayCFD(_binDelayCFD), fractionCFD(_fractionCFD) {}
		
		/** @brief Destructor */
		~ahdcExtractor(){;}
		
		/**
		 * This method extracts relevant informations from the digitized signal
		 * (the samples) and store them in a Pulse
		 *
		 * @param samples ADC samples
		 */
		std::map<std::string,double> extract(const std::vector<short> samples);

	private :
		/**
		* This method subtracts the pedestal (noise) from samples and stores it in : `samplesCorr`
		* It also computes a first value for : `adcMax`, `binMax`, `timeMax` and `integral`
		* This code is inspired by the one of coatjava/.../MVTFitter.java
		*/
		void waveformCorrection();

		/**
		 * This method gives a more precise value of the max of the waveform by computing the average of five points around the binMax
		 * It is an alternative to fitParabolic()
		 * The suitability of one of these fits can be the subject of a study
		 * Remark : This method updates adcMax but doesn't change timeMax
		 */
		void fitAverage();

		/**
		 * @brief Alternative to `fitAverage`
		 */
		void fitParabolic();
		
		/**
		 * From MVTFitter.java
		 * Make fine timestamp correction (using dream (=electronic chip) clock)
		 * 
		 * Parameter dependency :
		 * - `timeStamp` 
		 * - `fineTimeStampResolution`
		 */
		void fineTimeStampCorrection();
		//void fineTimeStampCorrection (long timeStamp, float fineTimeStampResolution);

		/**
		 * This method determines the moment when the signal reaches a Constant Fraction of its Amplitude (i.e amplitudeFraction*adcMax)
		 * It fills the attributs : `timeRiseCFA`, `timeFallCFA`, `timeOverThresholdCFA`
		 * 
		 * Parameter dependency :
		 * - `samplingTime` time between 2 ADC bins
		 * - `amplitudeFraction` amplitude fraction between 0 and 1
		 */
		void computeTimeAtConstantFractionAmplitude();

		/**
		 * This methods extracts a time using the Constant Fraction Discriminator (CFD) algorithm
		 * It fills the attribut : `timeCFD`
		 *
		 * Parameter dependency :
		 * - `samplingTime` time between 2 ADC bins
		 * - `fractionCFD` CFD fraction parameter between 0 and 1
		 * - `binDelayCFD` CFD delay parameter
		 */
		void computeTimeUsingConstantFractionDiscriminator();
	public:	
		std::vector<float> samplesCFD; ///< samples corresponding to the CFD signal
};

#endif




