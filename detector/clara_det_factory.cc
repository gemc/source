
// gemc headers
#include "clara_det_factory.h"
#include "utils.h"

#include <dlfcn.h>

typedef map<string, map<string,string> > volmap_t;

map<string, detector> clara_det_factory::loadDetectors()
{
	string hd_msg     = gemcOpt.optMap["LOG_MSG"].args + " CLARA Factory: >> ";
	// double verbosity  = gemcOpt.optMap["GEO_VERBOSITY"].arg;

	map<string, detector> dets;
	// first check if there's at least one detector with CLARA factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;
	
	// checking if the plugin directory exist 
	if(getenv ("GEMC_PLUGINS") == NULL)
	{
		cout << "  !!! Warning: the GEMC_PLUGINS env variable, needed for the CLARA plugin, is not set. " << endl;		
		cout << "  !!! Warning: for CLAS12, this is typically /group/clas12/lib " << endl;
		cout << "  !!! Warning: CLARA detectors won't be loaded. " << endl << endl;
		return dets;
	}

	string clasraPlugin = (string) getenv("GEMC_PLUGINS") + "/libclas12_geometry_gemc.so";
	
	cout << "  > Opening geometry plugin..." << endl;
	void* handle = dlopen(clasraPlugin.c_str(), RTLD_NOW);

	if (!handle)
	{
		cerr << "  !!! Error: Cannot open library: " << dlerror() << '\n';
		exit(0);
	}
	
	// reset errors
	dlerror();

	// extract function symbol from library
	typedef volmap_t (*get_volume_maps_t)(const map<string,string>&);
	get_volume_maps_t get_volume_maps = (get_volume_maps_t) dlsym(handle, "get_volume_maps");

	const char *dlsym_error = dlerror();
 	if (dlsym_error)
	{
		cerr << "  !!! Error: Cannot load symbol 'get_volume_maps': " << dlerror() << endl;
		dlclose(handle);
		exit(0);
	}

	
	// building detectors that are tagged with CLARA factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end(); it++)
	{
		if(it->second.get_factory() != factoryType )
			continue;
		
        
		// send request with a string with the detector name in it + "/volumes"
		// e.g. dc/volumes. 
		map<string,string> request;
		request["request"] = it->first + "/volumes";
 
		map<string, map<string,string> > volumes = get_volume_maps(request);
		
		cout << " >>>>  Loading: " << it->first << endl;
		
		for(map<string, map<string,string> >::iterator idet = volumes.begin(); idet != volumes.end(); idet++)
		{
			gtable gt;
			
			gt.add_data(idet->first);                   // 1 name
			gt.add_data(idet->second["mother"]);        // 2 mother volume 
			gt.add_data(idet->second["description"]);   // 3 description  
			gt.add_data(idet->second["pos"]);           // 4 position  
			gt.add_data(idet->second["rotation"]);      // 5 rotation  
			gt.add_data(idet->second["color"]);         // 6 color
			gt.add_data(idet->second["type"]);          // 7 type 
			gt.add_data(idet->second["dimensions"]);    // 8 dimensions 
			gt.add_data(idet->second["material"]);      // 9 material
			gt.add_data(idet->second["mfield"]);        // 10 magnetic field
			gt.add_data(idet->second["ncopy"]);         // 11 copy number
			gt.add_data(idet->second["pMany"]);         // 12 pmany
			gt.add_data(idet->second["exist"]);         // 13 activation flag
			gt.add_data(idet->second["visible"]);       // 14 visibility
			gt.add_data(idet->second["style"]);         // 15 style
			gt.add_data(idet->second["sensitivity"]);   // 16 sensitivity  
			gt.add_data(idet->second["hit_type"]);      // 17 hit_type
			gt.add_data(idet->second["identifiers"]);   // 18 identifiers
			gt.add_data(it->first);                     // 19 system
			gt.add_data((string) "CLARA");              // 20 factory
			gt.add_data(idet->second["variation"]);     // 21 variation
			gt.add_data(idet->second["run"]);           // 22 run number

			dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);


		}
	}
 
	// close the library
	dlclose(handle);

		
	return dets;
}





