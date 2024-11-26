#ifndef SQLITE_MIRRORS_H
#define SQLITE_MIRRORS_H

#include "mirrors_factory.h"

class sqlite_mirrors : public mirrors {
public:
    ~sqlite_mirrors() {}

    map<string, mirror *> initMirrors(runConditions, goptions);  // Method to define the mirrors

    static mirrors *createMirrors() {
        return new sqlite_mirrors;
    }

};


#endif
