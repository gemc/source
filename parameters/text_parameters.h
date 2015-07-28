#ifndef TEXT_PARAMETERS_H
#define TEXT_PARAMETERS_H

#include "parameter_factory.h"

class text_parameters : public parametersFactory
{
	public:
		~text_parameters(){}
	
		// load all parameters that matches factorytype
		map<string, double> loadParameters(goptions, runConditions);  // Method to define the parameters
  
  static parametersFactory *createParametersFactory() 
  {
    return new text_parameters;
  }

};


#endif
