#ifndef CPP_MATERIALS_H
#define CPP_MATERIALS_H

#include "material_factory.h"

class cpp_materials : public materials
{
	public:
		~cpp_materials(){}
	
		map<string, G4Material*> initMaterials(runConditions, goptions);  // Method to define the G4 Materials
  
	  static materials *createMaterials() 
	  {
		return new cpp_materials;
	  }

};


#endif
