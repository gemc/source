<b>GE</b>ant4 <b>M</b>onte-<b>C</b>arlo


<b>gemc</b> is an application based on <a href="https://geant4.cern.ch"> Geant4 </a> to simulate the 
passage of particles through matter.

It provides:

* application independent geometry description
* easy interface to build / run experiments
* cad/gdml imports


<p align="center">
	<img src="https://github.com/gemc/docs/blob/master/webPage/source/beam.png"   height="170" width="170">
	<img src="https://github.com/gemc/docs/blob/master/webPage/source/eic_beam.png" height="170" width="170">
	<img src="https://github.com/gemc/docs/blob/master/webPage/source/eic.png"    height="170" width="170">
	<img src="https://github.com/gemc/docs/blob/master/webPage/source/bubble.png" height="170" width="170">
</p>
<i> examples of experiments using gemc. From left to right: a <a href="https://www.jlab.org/Hall-B/clas12-web/">CLAS12</a> beamline event; 
The electron and ion beamline and interation point detectors for the <a href="https://www.jlab.org/meic/">electron-ion collider</a>.
The gemc simulation of the Electron Ion Collider beamline and detectors.
10,000 electrons producing photons in the 6mm collimator in the bubble experiments at Jefferson Lab.</i>

========


### Overview:

gemc makes easy things trivial and hard things possible.

Users can build and run complex setups with minimal programming knowledge. 
See for example <a href="https://gemc.jlab.org/gemc/html/examples/paddles.html#simplepaddleexample">how to build a TOF with
few lines of code</a>.

Experiments can be loaded  using a combination of several available factories:

- MYSQL
- TEXT
- GDML
- CAD (STL, PLY, OBJ formats)
- C++ Plugin


<script src="https://embed.github.com/view/3d/gemc/detectors/master/humanBody/Upper_GI.stl?width=250"></script>
<script src="https://embed.github.com/view/3d/gemc/detectors/master/forFun/enterprise.stl?width=250"></script>





### Simulations are application independent

Once the user defined setup is loaded, it is translated in geant4. This includes:

- geometry
- materials
- mirrors
- physics list
- digitization
- electromagnetic fields

All particles are transported through matters and
produce radiation, hits, secondaries. The geant4 results are then collected and organized according to user preferences.



<p align="center">
	<img src="https://github.com/gemc/docs/blob/master/webPage/source/gemcArchitecture.png" width="90%">
</p>


========


### Platform Supported:

* Linux (32, 64)
* Mac OS


========

### Documentation:
* <a href="http://gemc.jlab.org"> gemc home page </a>

========


