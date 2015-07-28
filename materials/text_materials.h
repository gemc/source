#ifndef TEXT_MATERIALS_H
#define TEXT_MATERIALS_H

#include "material_factory.h"

class text_materials : public materials
{
	public:
		~text_materials(){}
	
		map<string, G4Material*> initMaterials(runConditions, goptions);  // Method to define the G4 Materials
  
	  static materials *createMaterials() 
	  {
		return new text_materials;
	  }

};


#endif
