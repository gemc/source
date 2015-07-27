<b>GE</b>ant4 <b>M</b>onte-<b>C</b>arlo


<b>gemc</b> is an application based on Geant4 libraries to simulate the 
passage of particles through matter.



The simulation models are stored on databases (MYSQL for example). This makes <b>gemc</b> ideal for collaborations as any change in the models can be tested and used immediately by all users.


![](https://github.com/gemc/gemc.github.io/blob/master/img/gemcAbstract.png)

========


### Overview:

All the simulation parameters (geometry, fields, sensitivity, etc) are defined in external 
databases (MYSQL, TXT), At run time, options can be given to tilt objects, set conditions, etc. 

The API to build the model is based on simple scripts. No previous knowledge of C++ or geant4 is required, 
even for the most complex simulations:  the users can focus solely on the models geometry, materials, etc. 

Any change in the model is reflected immediately in the databses, and can be tested w/o having to re-compile code. In fact the same <b>gemc</b> executable is used for many different experiments.


========

Test

### Platform Supported:

* Linux (32, 64)
* Mac OS


========

### Documentation:
* <a href="gemc.jlab.org"> gemc home page </a>




