#ifndef SQLITECAD_DET_FACTORY_H
#define SQLITECAD_DET_FACTORY_H 1

// gemc headers
#include "detector_factory.h"

class sqlitecad_det_factory : public detectorFactory {
public:
    sqlitecad_det_factory() { ; }

    // load all detectors that matches factorytype
    map <string, detector> loadDetectors();

    // initialize factorytype, option and runcondition classes
    void initFactory(goptions, runConditions, string);

    static detectorFactory *createFactory() {
        return new sqlitecad_det_factory;
    }

    // check that file.(allowed extension) exist.
    // if it does, returns file.ext
    // otherwise returns "na"
    string checkFormat(string file);
};


#endif
