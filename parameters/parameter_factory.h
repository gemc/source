#ifndef PARAMETER_FACTORY_H
#define PARAMETER_FACTORY_H

// C++ headers
#include <map>
#include <iostream>
using namespace std;

// gemc headers
#include "run_conditions.h"
#include "options.h"
#include "utils.h"


class parametersFactory
{
	public:
		virtual map<string, double> loadParameters(goptions, runConditions) = 0;    // Pure Virtual Method to initialize the parameters
		virtual ~parametersFactory(){}

		string factoryType;

		void initFactory(string ft)
		{
			cout << "  > gemc Init: " << ft << " Parameters Factory is Initialized "  << endl;
			factoryType = ft;
		}
};

typedef parametersFactory *(*parameterFactoryInMap)();                                 // Define parameterFactoryInMap as a pointer to a function that returns a pointer 

parametersFactory *getParameterFactory(map<string, parameterFactoryInMap> *, string);  // returns parameterFactory Function from Factory Map

map<string, parameterFactoryInMap> registerParameterFactories();                       // Registers parameterFactory in Factory Map

map<string, string> getParametersMap(map<string, double>);                             // Return Parameter Map

map<string, double> loadAllParameters(map<string, parameterFactoryInMap>, goptions, runConditions);


// from gtable to value
double get_par_value(gtable);

// log parameters 
void log_value(gtable, string);

#endif
