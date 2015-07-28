#ifndef MYSQL_MATERIALS_H
#define MYSQL_MATERIALS_H

#include "material_factory.h"

class mysql_materials : public materials
{
	public:
		~mysql_materials(){}
	
		map<string, G4Material*> initMaterials(runConditions, goptions);  // Method to define the G4 Materials
  
	  static materials *createMaterials() 
	  {
		return new mysql_materials;
	  }

};


#endif
