#ifndef GDML_MATERIALS_H
#define GDML_MATERIALS_H

#include "material_factory.h"

class gdml_materials : public materials
{
	public:
		~gdml_materials(){}
	
		map<string, G4Material*> initMaterials(runConditions, goptions);  // Method to define the G4 Materials
  
	  static materials *createMaterials() 
	  {
		return new gdml_materials;
	  }

};


#endif
