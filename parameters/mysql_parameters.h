#ifndef MYSQL_PARAMETERS_H
#define MYSQL_PARAMETERS_H

#include "parameter_factory.h"

class mysql_parameters : public parametersFactory
{
	public:
		~mysql_parameters(){}
	
		map<string, double> loadParameters(goptions, runConditions);  // Method to define the parameters
  
  static parametersFactory *createParametersFactory() 
  {
    return new mysql_parameters;
  }

};


#endif
