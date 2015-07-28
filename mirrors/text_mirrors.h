#ifndef TEXT_MIRRORS_H
#define TEXT_MIRRORS_H

#include "mirrors_factory.h"

class text_mirrors : public mirrors
{
	public:
		~text_mirrors(){}
	
		map<string, mirror*> initMirrors(runConditions, goptions);  // Method to define the mirrors
  
	  static mirrors *createMirrors()
	  {
		return new text_mirrors;
	  }

};


#endif
