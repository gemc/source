/// \file gbank.h
/// Defines the gemc bank class.\n
/// The output information, grouped by type, is dynamic
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef gbank_H
#define gbank_H 1

// gemc headers
#include "options.h"

// C++ headers
#include <sstream>
#include <set>
using namespace std;


// Banks. Let's put these in DB?

// simulation conditions NUM is 0 for the mother bank, 1 for the data bank
#define SIMULATION_CONDITIONS_BANK_TAG 5

// header bank
#define HEADER_BANK_TAG 10

// header bank
#define USER_HEADER_BANK_TAG 11

// RF bank
#define RF_BANK_TAG 30


// simulation conditions NUM is 0 for the mother bank, 1 for the data bank
#define GENERATED_PARTICLES_BANK_TAG 20

// simulation conditions NUM is 0 for the mother bank, 1 for the data bank
#define GENERATED_SUMMARY_BANK_TAG 21

// generated particle user info
#define GENERATED_USE_INFO_TAG 22

// FLUX Bank
#define FLUX_BANK_TAG 50

// Mirrors Bank
#define MIRRORS_BANK_TAG 60

// COUNTER Bank
#define COUNTER_BANK_TAG 70

// These Bank Types ID can be left hardcoded here
#define DETECTOR_BANK_ID 0

// true information, integrated over the hit
// this number adds to the TAG of the mother bank
// and it's added to the tag of the variables
// example:
// DC (100, 0)
//   - DC Raw (100 + RAWINT_ID, 0)
//     - DC pid (100 + RAWINT_ID, variable ID)
// identified by "R" in the variable type
// Variable ID must be different from zero
#define RAWINT_ID 1

// digitized information, integrated over the hit
// this number adds to the TAG of the mother bank
// and it's added to the tag of the variables
// example:
// DC (100, 0)
//   - DC Raw (100 + DGTINT_ID, 0)
//     - DC pid (100 + DGTINT_ID, variable ID)
// identified by "D" in the variable type
// Variable ID must be different from zero
#define DGTINT_ID 2

// true information, step by step
// identified by "D" in the variable type
#define RAWSTEP_ID 3

// digitized information, multi hit (step by step)
// identified by "M" in the variable type
#define DGTMULTI_ID 4

// charge time is a user-defined info for every step:
// it provides a (as seen by the PMT) charge and its timing
// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
#define CHARGE_TIME_ID 5

// quantum signal: it's the processed signal every # nanoseconds
// (bunch time parameter given by user)
// identified by "Q" in the variable type
#define QUANTUM_SIGNAL_ID 6


/// \class gBank
/// <b>gBank </b>\n\n
/// This class defines the general bank content.\n
/// The structure of the content is read
/// from the database or a text file.\n
class gBank
{
public:
	gBank(){;}
	gBank(int i, string bname, string d)
	{
		idtag        = i;
		bdescription = d;
		bankName = bname;
		name.clear();
		gid.clear();
		description.clear();
	}
	~gBank(){;}

public:
	int    idtag;                 ///< unique id for the bank
	string bdescription;          ///< bank description
	string bankName;              ///< name of the bank, it's also key in the map but we store it here as well

	vector<string>  name;         ///< Variable name.
	vector<int>     gid;           ///< Output variable identifier
	// variable type is 2 chars. The first char represent the type of bank:
	// N is for no level banks, the rest are defined above (N, R, D, S, M, V).
	// "i"nt , "d"ouble, "s"tring
	vector<string>  type;         ///< Variable type
	vector<string>  description;  ///< Variable description

	void load_variable(string, int, string, string);  // Load a variable in the bank definition

	// these two function return id and type
	// first occurance of string in the vector
	// there should be only one (not enforced right now)
	int getVarId(string);
	string getVarType(string);

	// returns the type of variable (rawInt, dgtInt, etc)
	int getVarBankType(string);

	// vector of names ordered by ID
	map<int, string> orderedNames;
	void orderNames();


	friend ostream &operator<<(ostream &stream, gBank);       ///< Overloaded "<<" for the class 'bank'

};


// get bank definitions (all)
gBank getBankFromMap(string, map<string, gBank>*);

// get dgt bank definitions
gBank getDgtBankFromMap(string, map<string, gBank>*);

// creates maps of banks based on sensitivity
map<string, gBank> read_banks(goptions gemcOpt, map<string, string> allSystems);


#endif





