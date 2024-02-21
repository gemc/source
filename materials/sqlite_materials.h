#ifndef SQLITE_MATERIALS_H
#define SQLITE_MATERIALS_H

#include "material_factory.h"

class sqlitel_materials : public materials {
public:
    ~sqlitel_materials() {}

    map<string, G4Material *> initMaterials(runConditions, goptions);  // Method to define the G4 Materials

    static materials *createMaterials() {
        return new sqlitel_materials;
    }

};


#endif
