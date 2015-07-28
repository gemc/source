#ifndef MYSQL_MIRRORS_H
#define MYSQL_MIRRORS_H

#include "mirrors_factory.h"

class mysql_mirrors : public mirrors
{
	public:
		~mysql_mirrors(){}
	
		map<string, mirror*> initMirrors(runConditions, goptions);  // Method to define the mirrors
  
	  static mirrors *createMirrors()
	  {
		return new mysql_mirrors;
	  }

};


#endif
