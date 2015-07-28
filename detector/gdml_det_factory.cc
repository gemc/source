
// gemc headers
#include "gdml_det_factory.h"
#include "utils.h"
#include <cmath>
double verbosity;
map<string, detector> gdml_det_factory::loadDetectors()
{

	string hd_msg     = " >> GDML Factory: >> ";
	verbosity   = gemcOpt.optMap["GEO_VERBOSITY"].arg;
	map<string, detector> dets;
	
  	// first check if there's at least one detector with GDML factory
	if(!check_if_factory_is_needed(RC.detectorConditionsMap, factoryType))
		return dets;
	
  	// there is at least one build with this factory.
	// building all detectors that are tagged with GDML factory
	for(map<string, detectorCondition>::iterator it=RC.detectorConditionsMap.begin(); it != RC.detectorConditionsMap.end() ; it++)
	{
		if(it->second.get_factory() != factoryType )
			continue;
		
		string dname     = it->first;
		string fname     = dname + ".gdml";
		
		if(verbosity)
			cout <<  hd_msg << " Importing Detector: " <<  dname << " with " << factoryType << " factory. "  << endl;
		
		// If found, parse the gdml file
		QFile gdet(fname.c_str());
    
		if( !gdet.exists() )
		{
			cout << hd_msg << "  Failed to open geometry file " << fname << " for system: " << dname << ". Maybe the filename doesn't exist? Exiting." << endl;
			exit(0);
		}
		
		QDomDocument domDocument;
		// opening gcard and filling domDocument
		if(!domDocument.setContent(&gdet))
		{
			gdet.close();
			cout << hd_msg << "  Failed to open geometry file " << fname << " for system: " << dname << ". Wrong XML syntax. Exiting." << endl;
			exit(0);
		}
		gdet.close();
		
		
		// maps of various objects
		map<string, double> gconstantMap;
		map<string, gposition> gpositionsMap;
		map<string, grotation> grotationsMap;
		map<string, glogicV> glogicVMap;
		map<string, gsolid> gsolidsMap;
		map<string, gphysV> gphysVMap;
		map<string, goperation> operationMap;
		map<string, string> auxiliaryMap;
		map<string, gposition> widthOffsetMap;

		//add "unity" as default value for position and rotation.
		gpositionsMap["unity"] = gposition(0, 0, 0, "mm");
		grotationsMap["unity"] = grotation("0", "0", "0", "deg");
		//map<string, string> stringConstantMap;
		
		// fill with   <quantity name="wextent" type="length" value="10000.0" unit="mm"/>
		
		
		QDomElement docElem = domDocument.documentElement(); ///< reading gcard file
		QDomNode n = docElem.firstChild(); 
		while(!n.isNull())///< looping over options.
		{
			
			///< converts the node to an element.
			QDomElement e = n.toElement();
			///< the node really is an element.

			// selecting "define" nodes for definitions of parameters
			if(!e.isNull() && e.tagName().toStdString() == "define")
			{
				///<gets first child node of the element
				QDomNode nn= e.firstChild();
				
				while( !nn.isNull() )
				{
					QDomElement ee = nn.toElement();
					
					//parses constants and add to gconstantMap
					if(ee.tagName().toStdString() == "constant")
				 	{
						///<gets the string value of "name". "noname" is default value
						string constant_name= assignAttribute(ee, "name", "noname");
						string v = assignAttribute(ee, "value", "0");
				
						//call resolve_strings function to get double constant value
						gconstantMap[constant_name] = resolve_strings(constant_name, v, gconstantMap);
						if(verbosity>0) cout << "map value of " << constant_name << " " << gconstantMap[constant_name] << endl;

					}
					
					//calls resolve_string function for x, y and z respectively to get double values and add to gpositionsMap
					if(ee.tagName().toStdString() == "position")
				 	{
						
						string position_name= assignAttribute(ee, "name", "noname");
										
						if(verbosity>0)	cout << "!!!!!!!!!!!!!! position name : " << position_name << " !!!!!!!!!!!!!!!!! " << endl;
						
						string x = assignAttribute(ee, "x", "0");
						//resolve_string will calculate the string or change it into a double value
						double x_p = resolve_string(x, gconstantMap);
						if(verbosity>0) cout << "x_p is obtained " <<"x_p = " << x_p <<  endl;
						string y = assignAttribute(ee, "y", "0");
						double y_p = resolve_string(y, gconstantMap);
						if(verbosity>0) cout << "y_p is obtained " <<"y_p = " << y_p <<  endl;
						string z = assignAttribute(ee, "z", "0");
						double z_p = resolve_string(z, gconstantMap);
						if(verbosity>0) cout << "z_p is obtained " <<"z_p = " << z_p <<  endl;

						string lunit = assignAttribute(ee, "unit", "mm");
						
						gpositionsMap[position_name] = gposition(x_p, y_p, z_p, lunit);
						
					}
					//calls resolve_string function for x, y and z respectively to get double values and add to grotationsMap
					if(ee.tagName().toStdString() == "rotation")
					{
						
						string rotation_name= assignAttribute(ee, "name", "noname");
						
						if(verbosity>0) cout << "!!!!!!!!!! rotation name : " << rotation_name << "!!!!!!!!!!!!!!!!! " << endl;

						string a = rad_to_deg(assignAttribute(ee, "x", "0"));
						string a_r = stringify(resolve_string(a, gconstantMap));
						if(verbosity>0) cout << "a_r is obtained " << a_r << endl;
						string b = rad_to_deg(assignAttribute(ee, "y", "0"));
						string b_r = stringify(resolve_string(b, gconstantMap));
						if(verbosity>0) cout << "b_r is obtained " << b_r << endl;
						string c = rad_to_deg(assignAttribute(ee, "z", "0"));
						string c_r = stringify(resolve_string(c, gconstantMap));
						if(verbosity>0) cout << "c_r is obtained " << c_r << endl;
 						string aunit = assignAttribute(ee, "unit", "deg");
						
						//changes unit according to the original value before resolve_string is called
 						rad_to_deg_u (a, aunit);
 						rad_to_deg_u (b, aunit);
 						rad_to_deg_u (c, aunit);
						
 						grotationsMap[rotation_name] = grotation(a_r, b_r, c_r, aunit);
						
						
 					}
					
					
 					nn=nn.nextSibling();///<move to next sibling node after parsing each node
					
 				}
				
				
				
 			}
			
			
 			/// selecting "solids" node
 			if(!e.isNull() && e.tagName().toStdString() == "solids")	
 			{
	 			
				QDomNode nn= e.firstChild();
				vector<string> seconds;
	 			while( !nn.isNull() )
	 			{
	 				QDomElement ee = nn.toElement();
					//parses type and calls togType to get correct name of solids types
					string type = togType(ee.tagName().toStdString());
					
					//discriminates operations with type
	 				if(type == "union" || type == "subtraction")
					{
						string first, second;
						//define default value of position and rotation reference as unity
						string pos = "unity";
						string rot = "unity";
						string optype, renamed_second;
						
						string opname = assignAttribute(ee, "name", "noname");
						
						QDomNode nnn= ee.firstChild();
	 					while( !nnn.isNull() )
	 					{
	 						QDomElement eee = nnn.toElement();
							
							if(eee.tagName().toStdString() == "first")
								first = assignAttribute(eee, "ref", "noref");
							
							if(eee.tagName().toStdString() == "second")
							{
								second = assignAttribute(eee, "ref", "noref");
								//add every seconds to seconds vector to compare with other seconds
								seconds.push_back(second);
								//reused indicates the number of repetition as a seconds
								double reused = 0;
								for(unsigned int i=0; i<seconds.size(); i++)
								{
									
									if(seconds[i] == second)
										reused++;
								}
								//if a second is reused rename it with the number of repetition(reused)
								if(reused>1)
									renamed_second = second + "_" + stringify(reused);
								else 
									renamed_second = second;

							
							}
							
 							if(eee.tagName().toStdString() == "positionref")
								pos = assignAttribute(eee, "ref", "noref");
							
							if(eee.tagName().toStdString() == "rotationref")
								rot = assignAttribute(eee, "ref", "noref");
							
							nnn=nnn.nextSibling();
						}
						//defines types for union and subtraction with + and - respectively
						if(type == "union")
							optype = "Operation: " + first + "+" + second;
							
						if(type == "subtraction")
							optype = "Operation: " + first + "-" + second;
						
						//adds operations to gsolidsMap (type will be optype and dimension is 0)
						gsolidsMap[opname] = gsolid(optype, "0");
						//adds operation information to the map so that it can be used to make gtables for components.
						operationMap[opname] = goperation(first, second, renamed_second, pos, rot);
						
					}
					
					// not an operation
					else
					{
						string name = assignAttribute(ee, "name", "unnamed");
						
						cout << " " << endl;
						if(verbosity>0) cout << "-------solids name : " << name << " solids type : " <<type << "--------" <<  endl;
						//get_dimensions will return each solid dimension which is defined in gdml_solids.cc.
						gsolidsMap[name] = gsolid(type, get_dimensions(nn, gconstantMap));
						
					}
					
					nn=nn.nextSibling();
				}
			}
			cout << " " << endl;
			/// selecting "structure" nodes, they contain logical and physical volumes
			if(!e.isNull() && e.tagName().toStdString() == "structure")
			{
				QDomNode nn= e.firstChild();
				while( !nn.isNull() )
				{
					QDomElement ee = nn.toElement();
          
					/// selecting "volume" node
					if(!ee.isNull() && ee.tagName().toStdString() == "volume")
					{
						string lvol_name= assignAttribute(ee, "name", "unnamed");
						
						// childs of volumes to look for phys volumes, volumeref
						QDomNode nnn= ee.firstChild();
						string matref;
						string solidref;
						string auxtype, auxvalue;
						string replica_num;
						while( !nnn.isNull() )
						{
							QDomElement eee = nnn.toElement();
							
							// I suspect physappear volumes are child of volumes only if volume is "World"
							///< selecting "physical volume" nodes with World tag
							//lvol_name == "World" &&

							//selecting physical volume node
							if(eee.tagName().toStdString() == "physvol")
							{
								string posref = "unity";
								string rotref = "unity";
								string physvol_name;
								
								QDomNode nnnn= eee.firstChild();
								while( !nnnn.isNull() )
								{
									QDomElement eeee = nnnn.toElement();
									
									if(eeee.tagName().toStdString() == "volumeref")
										physvol_name=eeee.attributeNode("ref").value().toStdString();

									if(eeee.tagName().toStdString() == "positionref")
										posref=assignAttribute(eeee, "ref", "unity");

									if(eeee.tagName().toStdString() == "rotationref")
										rotref=assignAttribute(eeee, "ref", "unity");
									
									nnnn=nnnn.nextSibling();
								}
								//adds position and rotation reference to gphysVMap
								gphysVMap[physvol_name] = gphysV(posref, rotref);
							}
							
							
							if(eee.tagName().toStdString() == "materialref")
								matref=assignAttribute(eee, "ref", "");
							
							
							if(eee.tagName().toStdString() == "solidref")
								solidref=assignAttribute(eee, "ref", "");

							//selecting replica node
							if(eee.tagName().toStdString() == "replicavol")
							{
								//gets double replica value in gconstantMap with the string value
								replica_num=assignAttribute(eee, "number", "0");
								double j = gconstantMap[replica_num];
							
								
								QDomNode nnnn= eee.firstChild();
								string volume_ref;
								while( !nnnn.isNull() )
								{
									
									QDomElement eeee = nnnn.toElement();
									
									if(eeee.tagName().toStdString() == "volumeref")
										volume_ref=eeee.attributeNode("ref").value().toStdString();
										
									if(eeee.tagName().toStdString() == "replicate_along_axis")
									{
										QDomNode n_5= eeee.firstChild();
										
										string x_s, y_s, z_s, direction;
										double width = 0;
										double offset = 0;
										string width_ref, width_unit, offset_ref, offset_unit;
										
										while( !n_5.isNull() )
										{
											
											QDomElement e_5 = n_5.toElement();
											if(e_5.tagName().toStdString() == "direction")
											{
												x_s = assignAttribute(e_5, "x", "0");
												y_s = assignAttribute(e_5, "y", "0");
												z_s = assignAttribute(e_5, "z", "0");
											}
											if(e_5.tagName().toStdString() == "width")
											{
												width_ref = assignAttribute(e_5, "value", "no value");
												width_unit = assignAttribute(e_5, "unit", "mm");
												width = gconstantMap[width_ref];
											}
											if(e_5.tagName().toStdString() == "offset")
											{
												offset_ref = assignAttribute(e_5, "value", "no value");
												offset_unit = assignAttribute(e_5, "unit", "mm");
												offset = gconstantMap[offset_ref];
											}
											n_5 = n_5.nextSibling();
										}

									
										if(x_s == "1") direction = "1";
										if(y_s == "1") direction = "2";
										if(z_s == "1") direction = "3";
									
										//obtains replica type with all the replica information
										string replica_type = "Replica "+ volume_ref + "," +"root" +","+ direction+ "," + stringify(j) + "," +stringify(width)+"*"+width_unit + "," + stringify(offset)+"*"+offset_unit;
										//defines replica type with the string above
										gsolidsMap[solidref].type = replica_type;
								
										
									}
									nnnn=nnnn.nextSibling();
								}
							}
							//parese auxiliary node and fills the map with type and value
							if(eee.tagName().toStdString() == "auxiliary")
							{
								auxtype=assignAttribute(eee, "auxtype", "");
								auxvalue=assignAttribute(eee, "auxvalue", "");
							}
							
							nnn=nnn.nextSibling();
							
						}
						//fills glogicVMap with material and solid reference
						glogicVMap[lvol_name] = glogicV(matref, solidref);
						auxiliaryMap[auxtype] = auxvalue;
					}
					nn=nn.nextSibling();
				}
			}
			
 			// selecting "define" nodes for definitions of parameters
			
			
			n=n.nextSibling();
			
		}
		
		//makes gtable for every logical volume
		for(map<string, glogicV>::iterator itt=glogicVMap.begin(); itt != glogicVMap.end(); itt++)
		{
			gtable gt;
			
			string a ="root";
			string b = "no description";
			string c = "9999ff";
			string d = "no";
			string e = "1";
			string f = "0";
			
			gt.add_data(itt->first); // first item is name
			gt.add_data(a); // second element is mother volume (not found for now)
			gt.add_data(b);
			
			//position
			//If an element which is same as logical volume name(itt->first) doesn't exist, add default value.
			if(gphysVMap.find(itt->first) == gphysVMap.end())
			{
				if(verbosity>0)cout << " !! Attention: position ref doesn't exist. Defaulting to (0,0,0)" << endl;
				gt.add_data((string) "0, 0, 0");
			}
			
			else
			{
				string posref_name = gphysVMap[itt->first].positionRef;
				//If there is no position value for a position reference, add a default value
			 	if(gpositionsMap.find(posref_name) == gpositionsMap.end())
				{
					if(verbosity>0)cout << " !! Attention: position ref for " << posref_name << " not found. Defaulting to (0,0,0)" << endl;
					gt.add_data((string) "0, 0, 0");
				}
		
				else
					gt.add_data(gpositionsMap[posref_name].get_dimensions()); // position is obtained by get_dimensions() like " x*unit y*unit z*unit"
				if(verbosity>0)cout << "position get_dimensions() : " << gpositionsMap[posref_name].get_dimensions()<< endl;
			}
			
			//rotation
			//If an element which is same as logical volume name(itt->first) doesn't exist, add default value.
			if(gphysVMap.find(itt->first) == gphysVMap.end())
			{
				if(verbosity>0)cout << " !! Attention: rotation ref doesn't exist. Defaulting to (0,0,0)" << endl;
				gt.add_data((string) "0, 0, 0");
			}
			else
			{
				string rotref_name = gphysVMap[itt->first].rotationRef;
				//If there is no rotation value for a rotation reference, add a default value
				if(grotationsMap.find(rotref_name) == grotationsMap.end())
				{
					if(verbosity>0)cout << " !! Attention: rotation ref for " << rotref_name << " not found. Defaulting to (0,0,0)" << endl;
					gt.add_data((string) "0, 0, 0");
				}
				
				else
					gt.add_data(grotationsMap[rotref_name].get_dimensions());// rotation is obtained by get_dimensions() like " x*unit y*unit z*unit"
			
				if(verbosity>0)cout << "rotation get_dimensions() : " << grotationsMap[rotref_name].get_dimensions() << endl;
			}
			gt.add_data(c); // color
			
			string solid_type = gsolidsMap[itt->second.solidRef].type;
			gt.add_data(solid_type); // type
			
			string solid_dimensions = gsolidsMap[itt->second.solidRef].dimension;
			gt.add_data(solid_dimensions); // dimensions of solid
			
			gt.add_data(itt->second.materialRef); // material
			gt.add_data(d);  // magnetic field
			gt.add_data(e);  // copy number
			gt.add_data(e);  // pmany
			gt.add_data(e);  // 1 = active
			if (itt->first == "World")
				gt.add_data(f);  // 0 = invisible (World)
			else
 				gt.add_data(e);
			
			//gt.add_data(e);  // 1 = solid style ,  0 wirefeame
	 		gt.add_data(string("0"));  // 1 = solid style ,  0 wirefeame	
 			gt.add_data(d); // sensitivity
 			gt.add_data(d); // hitprocess
 			gt.add_data(f);   // identity
			gt.add_data(dname); //detector name
			gt.add_data((string) "GDML"); //GDML
			
			dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
			
			if(verbosity>0)cout << dets[gt.data[0]] << endl;
			
			
 		}
		
		//makes gtable for every component of operations
		for(map<string, goperation>::iterator itt=operationMap.begin(); itt != operationMap.end(); itt++)
		{
			//makes two gtables for every operation
			//One is for first component, the other is for second component
			for(unsigned int i=0; i<2; i++)
			{
				gtable gt;
			
				string a ="root";
				string b = "no description";
				string c = "9999ff";
				string d = "no";
				string e = "1";
				string f = "0";
				string g = "Component";
			
				string fs;
				
				// first item is name
				//makes a gtable for first component and then second component
				if(i == 0)
				{
					gt.add_data(itt->second.first);
					fs = "first";
					if(verbosity>0)cout << "-------------------first--------------------"<< endl;
				}
				else if(i == 1)
				{
					gt.add_data(itt->second.renamed_second); 
					fs = "second";
					if(verbosity>0)cout << "--------------------second----------------------" << endl;
				}
				gt.add_data(a); // second element is mother volume (not found for now)
				gt.add_data(b);
			
				//For position and rotation,
				//If the gtable if for first component, position is "0"
				//If it is for second component, get position from gpositionsMap and inline function get_dimensions()

				//position
				string posref_name;
				if(fs == "first")
					gt.add_data((string) "0, 0, 0");
				else if(fs == "second")
				{
					posref_name = itt->second.pos;
					//If there is no position value for a position reference, add a default value
					if(gpositionsMap.find(posref_name) == gpositionsMap.end())
					{
						if(verbosity>0)cout << " !! Attention: position ref for " << posref_name << " not found. Defaulting to (0,0,0)" << endl;
						gt.add_data((string) "0, 0, 0");
					}
					else
						gt.add_data(gpositionsMap[posref_name].get_dimensions());
					
					if(verbosity>0)cout << "position ref " << posref_name << " position " << gpositionsMap[posref_name].get_dimensions() << endl;
				}
			
				//rotation
				string rotref_name;
				if(fs == "first")
					gt.add_data((string) "0, 0, 0");
				else if(fs == "second")
				{
					rotref_name = itt->second.rot;
					//If there is no rotation value for a rotation reference, add a default value
				 	if(grotationsMap.find(rotref_name) == grotationsMap.end())
					{
						if(verbosity>0)cout << " !! Attention: position ref for " << rotref_name << " not found. Defaulting to (0,0,0)" << endl;
						gt.add_data((string) "0, 0, 0");
					}
		
					else
						gt.add_data(gpositionsMap[posref_name].get_dimensions()); 
			

					if(verbosity>0)cout << "rotation ref " << rotref_name << " rotation " << grotationsMap[rotref_name].get_dimensions() << endl;
				}
				
			
				gt.add_data(c); // color
			
				//type
				string lvol_name;
				//finds solid reference for first and second from operationMap
				//This solid reference(lvol_name) will be used to find type and dimension from gsolidsMap
				if(fs == "first")
					lvol_name = itt->second.first;
				else if(fs == "second")
					lvol_name = itt->second.second;

				string solid_type = gsolidsMap[lvol_name].type;
				gt.add_data(solid_type);
				
				//dimensions of solid
				string solid_dimensions = gsolidsMap[lvol_name].dimension;
				gt.add_data(solid_dimensions);
			
				gt.add_data(g); // material
				gt.add_data(d);  // magnetic field
				gt.add_data(e);  // copy number
				gt.add_data(e);  // pmany
				gt.add_data(e);  // 1 = active
				gt.add_data(e);  // 1 = visible

				//gt.add_data(e);  // 1 = solid style ,  0 wirefeame
	 			gt.add_data(string("0"));  // 1 = solid style ,  0 wirefeame
	 			gt.add_data(d); // sensitivity
	 			gt.add_data(d); // hitprocess
	 			gt.add_data(f);   // identity
				gt.add_data(dname); //detector name
				gt.add_data((string) "GDML"); //GDML
			
				dets[gt.data[0]] = get_detector(gt, gemcOpt, RC);
				if(verbosity>0)cout << dets[gt.data[0]] << endl;
				
			}
			
 		}


		
		
 	}
	
 	return dets;
	
}


ostream &operator<<(ostream &stream, map<string, glogicV> glogicVMap)
{
	for(map<string, glogicV>::iterator it = glogicVMap.begin(); it != glogicVMap.end(); it++)
		cout << "   Logic Volume :" << it->first << "     solidRef: " << it->second.solidRef << ",    materialRef: " << it->second.materialRef << endl;
	
	return stream;
}


ostream &operator<<(ostream &stream, map<string, gposition> gpositionsMap)
{
	for(map<string, gposition>::iterator it = gpositionsMap.begin(); it != gpositionsMap.end(); it++)
		cout << "   Position " << it->first << ":     x: " << it->second.x << ",    y: " << it->second.y << ",     z: " << it->second.z << " lunit: " << it->second.unit << endl;
	
	return stream;
}


ostream &operator<<(ostream &stream,map<string, grotation> grotationsMap )
{
	for(map<string,grotation>::iterator it = grotationsMap.begin(); it != grotationsMap.end(); it++)
		cout << "   Rotation " << it->first << ":     x: " << it->second.x << ",    y: " << it->second.y <<  ",     z: " << it->second.z << " aunit: " << it->second.unit << endl;
	
	return stream;
}


ostream &operator<<(ostream &stream,map<string, gsolid> gsolidsMap )
{
	for(map<string, gsolid>::iterator it = gsolidsMap.begin(); it != gsolidsMap.end(); it++)
		cout << "   Solid :" << it->first << "     type: " << it->second.type << ",    dimension: " << it->second.dimension << endl;
	
	return stream;
}


ostream &operator<<(ostream &stream,map<string, gphysV> gphysVMap )
{
	for(map<string, gphysV>::iterator it = gphysVMap.begin(); it != gphysVMap.end(); it++)
		cout << "   Physical Volume :" << it->first << "     positionRef: " << it->second.positionRef << ",    rotationRef: " << it->second.rotationRef << endl;
	
	return stream;
}



string togType(string type)
{
	if(type == "para")       return "Parallelepiped";
	if(type == "cone")       return "Cons";
	if(type == "box")        return "Box";
	if(type == "sphere")     return "Sphere";
	if(type == "ellipsoid")  return "Ellipsoid";
	if(type == "paraboloid") return "Paroboloid";
	if(type == "hype")       return "Hype";
	if(type == "tube")       return "Tube";
	if(type == "eltube")     return "EllipticalTube";
	if(type == "torus")      return "Torus";
	if(type == "trd")        return "Trd";
	if(type == "orb")        return "Orb";
	if(type == "tet")        return "Tet";
	
	return type;
}

double resolve_strings(string constant_name, string v, map<string, double> gconstantMap)
{
	
	double constantResult;

	// starts solving operations, and iteration will become smaller and smaller
	
	// If a string has "+","-","/" or "*", replace these to " + "," - "," / " and " * " so that it can be tokenized
	string w("+");	
	string x("-");
	string y("*");
	string z("/");
	v = replaceCharWithChars(v, w, " + ");
	v = replaceCharWithChars(v, x, " - ");
	v = replaceCharWithChars(v, y, " * ");
	v = replaceCharWithChars(v, z, " / ");

			vector<string> operation;
			vector<double> constant;
			vector<string> constantName;
		
			//puts constant names that exist in the resultMap by now to compare with operand tokens and change them with double values

			for(map<string, double>::iterator it=gconstantMap.begin(); it != gconstantMap.end(); it++)
				constantName.push_back(it->first);

			if (v == "0" || v == "0.0")
				constantResult = 0.0;
			
			else if ( atof(v.c_str())== 0 && atoi(v.c_str()) == 0)
			{
				if(verbosity>0)cout << "----------constant calculation---------- "<< endl;
				//tokenizes each constant value that needs to be substituted or calculated
				int size = v.length();
				char constantChar[size];
				strcpy(constantChar, v.c_str());
				char * pch;
				//tokenizes a string by space and "\n"
				pch = strtok (constantChar," \n");
				string token;
			
				//adds operators in operation vector and operands in constant vector after tokenizing
				while (pch != NULL)
				{
					 token = pch;
					 
					 if(token == "+" || token == "-" || token == "*" || token == "/")
					 	operation.push_back(token);
					 
					 //puts number tokens into the constant vector.
					 else if (atof(token.c_str()) != 0)
					 	constant.push_back(atof(token.c_str()));
					 
					 else
						constant.push_back(gconstantMap[token]);
					
					 
					if(verbosity>0) cout << "constant " <<  token << endl;
					 
					 pch = strtok (NULL, " \n");
				 }
				 
				 
				 //puts double constant value in the map if it doesn't include any operators
				 if (operation.size() == 0)
					 constantResult = constant[0];
				 //If calculation is needed.
				 else 
				 {
					 //gets a string with changed values and operators.
					 
					 unsigned int calSize = constant.size()*2-1;
					
					 vector<string> calculation;
					 calculation.resize(calSize);
					
					//make calculation vector in the same order of a formula (operand operator operand...)
					 for(unsigned int a=0;a<constant.size();a++)
					 	calculation[2*a] = stringify(constant[a]);
					 for(unsigned int a=0;a<constant.size()-1;a++)
					 	calculation[2*a+1] = operation[a];

					 if(verbosity>0) 
					 {
						 cout << " " << endl;
						 cout << "calculation" << endl;
						 cout << " " << endl;
					
						 for(unsigned int a=0;a<calSize;a++)
						 	cout << calculation[a] << endl;
					}

					//calls p_calculation to get the result of calculation.
					 calculation[0] = p_calculation(calculation);
					 
					if(verbosity>0) cout << "!!!!!calculated value!!!!!       " << calculation[0] << endl;
					 constantResult = atof(calculation[0].c_str());
					 
				 }
			
			}
			else
				constantResult = atof(v.c_str());

	return constantResult;
	
	
}


double resolve_string(string constString, map<string, double> gconstantMap)
{
	
	double result;
	
	//When there are these characters in a string, add space before and after the string
	string w("+");
	string x("-");
	string y("*");
	string z("/");
	string p1("(");
	string p2(")");
	string sine("sin");
	string cosine("cos");

	//For one character, use replaceCharWithChars which is defined in string_utilities.h
	constString = replaceCharWithChars(constString, w, " + ");
	constString = replaceCharWithChars(constString, x, " - ");
	constString = replaceCharWithChars(constString, y, " * ");
	constString = replaceCharWithChars(constString, z, " / ");
	constString = replaceCharWithChars(constString, p1, " ( ");
	constString = replaceCharWithChars(constString, p2, " ) ");
	
	//If it has more than one character, use replaceCharsWithChars also defined in string_utilities.h
	if(constString.size()>=sine.size())
		constString = replaceCharsWithChars(constString, sine, " sin ");
	if(constString.size()>=cosine.size())
		constString = replaceCharsWithChars(constString, cosine, " cos ");

	if(verbosity>0)cout << "constString after replacing : " << constString << endl;


	vector<string> calculation;
	
	vector<string> operation;
	vector<double> constant;
	
	//If constString is "0", return double value 0
	if (constString == "0" || constString == "0.0")
		result = 0.0;
	//If constString is not numbers and is a string,
	else if ( atof(constString.c_str())== 0 && atoi(constString.c_str()) == 0)
	{
		int size = constString.length();
		char constantChar[size];
		strcpy(constantChar, constString.c_str());
		char * pch;
		pch = strtok (constantChar," \n");
		string token;
		
		//tokenize the string and put values in gconstantMap	
		while (pch != NULL)
		{
			token = pch;

			//put each token to calculation vector
			if(token == "+" || token == "-" ||token == "*" ||token == "/" || token == "(" || token == ")" || token == "sin" || token == "cos" )
				calculation.push_back(token);
			//If a token is number, add double value of the number
			else if(atof(token.c_str()) != 0)
				calculation.push_back(token);

			//If the token is "pi", put "3.14" instead of it
			else if(token == "pi")
				calculation.push_back("3.14");
			//If the token is string, and gconstantMap has its value, put the value to calculation
			else 
				calculation.push_back(stringify(gconstantMap[token]));
	 		
	 		pch = strtok (NULL, " \n");
		}

		for(unsigned int i=0; i<calculation.size(); i++)
			if(verbosity>0)cout << "calculation[" << i << "] = " << calculation[i] << endl;		
		//If calculation has only one token, return it as a result
		if(calculation.size() == 1)
			result = atof(calculation[0].c_str());
		//If calculation vector has more than one element and needs calculation,
		else
		{
				
			vector<unsigned int> a;
			vector<unsigned int> b;
			unsigned int element_a, element_b;
			//makes vector for parenthesis and find the number of parenthesis
			for(unsigned int i=0; i<calculation.size();i++)
			{
				if(calculation[i] == "(")
					a.push_back(i);
				if(calculation[i] == ")")
					b.push_back(i);
			}
			unsigned int size_a = a.size();

			//If parenthesis exists, check if there is double parenthesis and calculate inside of every parenthesis pair
			if(size_a>0)
			{
				int doublep_count=0;
				vector<unsigned int> c;
				//finds double parenthesis by comparing the index of "(" and ")"
				for(unsigned int i=0; i<a.size()-1; i++)
				{
					//If a ")[i]" which is a pair with "([i]" is located in ahead of "([i+1]", there is double parenthesis  
					if(b[i]>a[i+1])
					{
						c.push_back(i+1);
						doublep_count++;
					}
				}
				unsigned int size_c = c.size();
	
				//If the string has one pair of parenthesis, calculate inside the parenthesis once using p_calculation.
				if(size_a>0 && doublep_count == 0)
				{
					element_a=a[0];
					element_b=b[0];
					
					while(size_a>0)//iterates until there is no parenthesis left
					{
						vector<string> paren_calculation;	
						//puts calculation elements which are between a pair of parenthesis to paren_calculation
						for(unsigned int i=element_a+1; i<element_b; i++)
							paren_calculation.push_back(calculation[i]);
						//calls p_calculation to calculate paren_calculation and add the result to the location of the corresponding "("
						calculation[element_a] = p_calculation(paren_calculation);	
						
						//reorganize calculation by shifting the next elements located after the parenthesis part and popping out as many times as the paren_calculation size+1
						for(unsigned int i = element_b+1; i<calculation.size(); i++)
							calculation[i-element_b+element_a] = calculation[i];
						for(unsigned int i = 0; i<element_b-element_a; i++)	
							calculation.pop_back();
						size_a--;
						//find the location(index) of each parenthesis again in the modified calculation vector
						for(unsigned int i=0; i<calculation.size();i++)
						{
							if(calculation[i] == "(")
								element_a = i;
							if(calculation[i] == ")")
								element_b = i;
						}
						for(unsigned int j=0; j<calculation.size(); j++)
							if(verbosity>0)cout << "calculation[" << j << "] = " << calculation[j] << endl;
					}
					
				}
				//If the string has double parenthesis, calculate twice using p_calculation.
				else if(doublep_count>0)
				{
					unsigned int m = c[0];
					//iterates calculation for double parenthesis until there is no double parenthesis. 
					//Same as single parenthesis calculation.
					while(size_c>0)
					{
						vector<string> paren_calculation;
						
						for(unsigned int n=a[m]+1; n<b[m-1]; n++)
						{
							paren_calculation.push_back(calculation[n]);
						}	
						for(unsigned int h=0; h<paren_calculation.size(); h++)
							if(verbosity>0)cout << "double_paren_calculation[i] = " << paren_calculation[h] << endl;

						calculation[a[m]] = p_calculation(paren_calculation);
						
						for(unsigned int n=b[m-1]+1; n<calculation.size(); n++)
							calculation[n-b[m-1]+a[m]] = calculation[n];
						for(unsigned int n=0; n<b[m-1]-a[m]; n++)
							calculation.pop_back();
						size_c--;
						size_a--;
						if(size_c>0)
						{
							vector<unsigned int>aa;
							vector<unsigned int>bb;
	
							for(unsigned int n=0;n<calculation.size();n++)
							{
								if(calculation[n] == "(")
									aa.push_back(n);
								if(calculation[n] == ")")
									bb.push_back(n);
								}
							for(unsigned int n=0;n<calculation.size()-1;n++)
							{
								if(bb[n]>aa[n+1])
									m=n;
							}
						}
							
					}

					//After removing double parenthesis, iterates again to remove single parenthesis
					vector<unsigned int> d;
					vector<unsigned int> e;
					unsigned int element_d, element_e;
					for(unsigned int i=0; i<calculation.size(); i++)
					{
						if(calculation[i] == "(")
								d.push_back(i);
						if(calculation[i] == ")")
							e.push_back(i);
					}
					element_d = d[0];
					element_e = e[0];
					unsigned int size_d = d.size();
					while(size_d>0)
					{
						vector<string>paren_calculation;
						for(unsigned int i=element_d+1; i<element_e; i++)	
							paren_calculation.push_back(calculation[i]);
						for(unsigned int h=0; h<paren_calculation.size(); h++)
							if(verbosity>0)cout << "__paren_calculation[i] = " << paren_calculation[h] << endl;
						calculation[element_d] = p_calculation(paren_calculation);
						
						for(unsigned int i=element_e+1; i<calculation.size(); i++)
							calculation[i-element_e+element_d] = calculation[i];
						for(unsigned int i=0; i<element_e-element_d; i++)
							calculation.pop_back();
						size_d--;
						for(unsigned int i=0; i<calculation.size();i++)
						{
							if(calculation[i] == "(")
								element_d = i;
							if(calculation[i] == ")")
								element_e = i;
						}
					}
				}
			}
			//If there is no parenthesis, go to the next step
			else
				if(verbosity>0)cout << "There is no parenthesis in calculation[]" << endl;
			//looks for sin and cos
			for(unsigned int i=0; i<calculation.size();i++)
			{
				//If there is sine or cosine in the calculation, calculate this trigonometric function and reorganize calculation by moving elements and popping out
				if(calculation[i] == "sin") 
				{
					calculation[i] = stringify(sin(atof(calculation[i+1].c_str())));
					for(unsigned int j=i+2; j<calculation.size();j++)
						calculation[j]=calculation[j-1];
					calculation.pop_back();
				}
	
				if(calculation[i] == "cos")
				{
					calculation[i] = stringify(cos(atof(calculation[i+1].c_str())));
					for(unsigned int j=i+2; j<calculation.size();j++)
						calculation[j]=calculation[j-1];
					calculation.pop_back();
				}
			}
			//If there is no parenthesis, sin and cos, calculate +-/* with p_calculation function
			result = atof(p_calculation(calculation).c_str());
			if(verbosity>0)cout << "return value result = " << result << endl;

			
		}
	}
	//If constString is composed of numbers, not "0"
	else 
		result = atof(constString.c_str());

return result;


}

//a function that calculates +, -, * and /
string p_calculation(vector<string> paren_calculation)
{	
	string paren_result = paren_calculation[0];

	for(unsigned int m=0; m<paren_calculation.size(); m++)
		if(verbosity>0)cout << "calculation in p_calculation : " << paren_calculation[m] << endl;

	//If "-" is the first element, (If there is a negative number)
	if(paren_calculation[0] == "-") 
	{
		//multiply -1 to the second element and reorganize calculation vector.
		paren_calculation[1] = stringify(-1*atof(paren_calculation[1].c_str()));
		if(verbosity>0)cout << "calculation[1] after multiplying -1 " << paren_calculation[1] << endl;
		for(unsigned int i =1;i<paren_calculation.size();i++)
			paren_calculation[i-1]=paren_calculation[i];
		paren_calculation.pop_back();
	}

	int count=0;
	unsigned int paren_calSize = paren_calculation.size();	
	
	//counts the number of "*" and "/"  
	for(unsigned int j=0;j<(paren_calSize-1)/2;j++)
	{
		if(paren_calculation[2*j+1] == "*" || paren_calculation[2*j+1] == "/")
			count++;
	}
	//If there is "*" or "/", calculate this operation first and reorganize calculation vector				 
	while(count>0)
	{
		for(unsigned int j=0;j<(paren_calculation.size()-1)/2;j++)
		{
				
			if(paren_calculation[2*j+1] == "*")
			{
							 
				paren_calculation[2*j] = stringify(atof(paren_calculation[2*j].c_str())*atof(paren_calculation[2*j+2].c_str()));
				//remove paren_calculation[2*j+1], [2*j+2] and shift other values
				for(unsigned int l=2*j+3;l<paren_calculation.size();l++)
					paren_calculation[l-2]=paren_calculation[l];
				//remove paren_calculation[end], paren_calculation[end-1]
				paren_calculation.pop_back();
				paren_calculation.pop_back();
				count=count-1;
							 
			}
			else if(paren_calculation[2*j+1] == "/")
			{
				paren_calculation[2*j] = stringify(atof(paren_calculation[2*j].c_str())/atof(paren_calculation[2*j+2].c_str()));
				//remove paren_calculation[2*j+1], [2*j+2] and shift other values
				for(unsigned int l=2*j+3;l<paren_calculation.size();l++)
					paren_calculation[l-2]=paren_calculation[l];
				//remove paren_calculation[end], paren_calculation[end-1]
				paren_calculation.pop_back();
				paren_calculation.pop_back();
				count=count-1;
								 
			}
		}
	}

	//After "*" and "/" calculation, calculate "+" and "-", the final result will be stored in paren_calculation[0]				 
	while(paren_calculation.size()>1)
	{
		int m = 0;
		if(paren_calculation[2*m+1] == "+")
		{
			paren_calculation[2*m] = stringify(atof(paren_calculation[2*m].c_str())+atof(paren_calculation[2*m+2].c_str()));
			for(unsigned int l=2*m+3;l<paren_calculation.size();l++)
				paren_calculation[l-2]=paren_calculation[l];
			//remove paren_calculation[end], paren_calculation[end-1]
			paren_calculation.pop_back();
			paren_calculation.pop_back();
		}
		if(paren_calculation[2*m+1] == "-")
		{
			paren_calculation[2*m] = stringify(atof(paren_calculation[2*m].c_str())-atof(paren_calculation[2*m+2].c_str()));
			for(unsigned int l=2*m+3;l<paren_calculation.size();l++)
				paren_calculation[l-2]=paren_calculation[l];
			//remove paren_calculation[end], paren_calculation[end-1]
			paren_calculation.pop_back();
			paren_calculation.pop_back();

		}	
						 

	}
	paren_result = paren_calculation[0];	

		return paren_result;

}

