#ifdef WINDOWS
#include <direct.h>
#define get_cwd _getcwd
#else
#include <unistd.h>
#define get_cwd getcwd
#endif


#include "ReadConfigFile.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

// Function for reading the configuration file
ReadConfigFile::ReadConfigFile(std::string& ifName, probParams& probParamsObject) {

	// obtaining the directory where the executable is stored
	getDirectory(probParamsObject.directory);

	// string variable for each line of the configuration file
	std::string line;
	// Configuration filename
	std::ifstream file(probParamsObject.directory + "\\" + ifName);

	// variables for finding first and last index of input variables
	int first, last;

	// Reading the lines of the configuration file and assigning variables accordingly
	// Check if file is opened
	if (file.is_open()){
		//Obtain line
		while (getline(file, line)){

			if(line.rfind("MESH_FILENAME")==0){
				findStringBounds(first,last,line);
				probParamsObject.mesh_ifName = line.substr(first,last-first);
			}
			else if(line.rfind("MESH_OUT_FILENAME")==0){
				findStringBounds(first,last,line);
				probParamsObject.mesh_ofName = line.substr(first,last-first);
			}
			else if(line.rfind("MARKER_BDRY")==0){
				getTags(line,probParamsObject.bdryTags);
			}
			else if(line.rfind("MARKER_MOVING")==0){
				getTags(line,probParamsObject.mTags);
			}
			else if(line.rfind("MARKER_PERIODIC")==0){
				getTags(line,probParamsObject.pTags);
			}
			else if(line.rfind("MARKER_INTERNAL")==0){
				getTags(line,probParamsObject.iTags);
			}
			else if(line.rfind("DEFORMATION_FILENAME")==0){
				findStringBounds(first,last,line);
				probParamsObject.dispFile =  line.substr(first,last-first);
			}
			else if(line.rfind("SLIDING_MODE")==0){
				findStringBounds(first,last,line);
				probParamsObject.smode =  line.substr(first,last-first);
			}
			else if(line.rfind("PERIODIC_MODE")==0){
				findStringBounds(first,last,line);
				probParamsObject.pmode =  line.substr(first,last-first);
			}
			else if(line.rfind("PERIODIC_TYPE") == 0){
				findStringBounds(first,last,line);
				probParamsObject.ptype = stoi(line.substr(first,last-first));
			}
			else if(line.rfind("PERIODIC_DIRECTION")==0){
				findStringBounds(first,last,line);

				bool invalid = false;
				std::string pDir = line.substr(first,last-first);
				for(int i = 0; i< last-first; i++){
					if(!isdigit(pDir[i])){
						invalid = true;
					}
				}

				probParamsObject.pDir = stoi(line.substr(first,last-first));

				try{
					if(probParamsObject.pDir > 2 || invalid){
						throw(pDir);
					}
				}
				catch(std::string& pDir){
					std::cout << "Invalid periodic direction parameter provided: " << pDir << std::endl;
					std::exit(0);
				}

			}
			else if(line.rfind("STEPS")==0){
				findStringBounds(first,last,line);
				probParamsObject.steps =  stoi(line.substr(first,last-first));
			}
			else if(line.rfind("CURVED")==0){
				std::string param = "CURVED";
				getBool(param, probParamsObject.curved, line);
			}

			else if(line.rfind("DATA_REDUCTION")==0){
				std::string param = "DATA_REDUCTION";
				getBool(param, probParamsObject.dataRed, line);
			}

			else if(line.rfind("DATA_RED_TOLERANCE")==0){
				findStringBounds(first,last,line);
				probParamsObject.tol =  stod(line.substr(first,last-first));
			}
			else if(line.rfind("INFLUENCE_FACTOR")==0){
				findStringBounds(first,last,line);
				probParamsObject.rFac =  stod(line.substr(first,last-first));
			}
			else if(line.rfind("GAMMA")==0){
				findStringBounds(first,last,line);
				probParamsObject.gamma =  stod(line.substr(first,last-first));
			}
			else if(line.rfind("MULTILEVEL")==0){
				std::string param = "MULTILEVEL";
				getBool(param, probParamsObject.multiLvl, line);
			}
			else if(line.rfind("DOUBLE_EDGE")==0){
				std::string param = "DOUBLE_EDGE";
				getBool(param, probParamsObject.doubleEdge, line);
			}
			else if(line.rfind("LVL_TOL_CRIT")==0){
				findStringBounds(first,last,line);
				probParamsObject.tolCrit = stod(line.substr(first,last-first));
			}
			else if(line.rfind("LVL_SIZE")==0){
				findStringBounds(first,last,line);
				probParamsObject.lvlSize =  stoi(line.substr(first,last-first));
			}
			else if(line.rfind("MULTI_CRIT")==0){
				findStringBounds(first,last,line);
				probParamsObject.mCrit =  line.substr(first,last-first);
			}
			else if(line.rfind("GENERATE_QUALITY")==0){
				std::string param = "GENERATE_QUALITY";
				getBool(param, probParamsObject.generateQuality, line);
			}
		}
		file.close();
	}else{std::cout << "Unable to open the configuration file: " << ifName << std::endl;
		std::exit(0);}

	// the periodic displacement with moving vertices requires a sliding method to be used.
	// check is performed whether that is the case
	try{
		if(probParamsObject.smode =="none" && (probParamsObject.pmode =="moving")){
			throw(probParamsObject.pmode);
		}
	}
	catch(std::string& pmode){
		std::cout << "Non-compatible configuration file settings. Periodicity setting: '" << pmode << "' is only applicable when sliding is allowed.";
		std::exit(0);
	}

	// if multi level greedy is used then dataRed must be set to YES also
	try{
		if(probParamsObject.dataRed == false && probParamsObject.multiLvl == true){
			throw(probParamsObject.multiLvl);
		}
	}
	catch(bool& mLvl){
		std::cout << "Multi-level was set to: " << mLvl << ", this is only possible when data reduction is enabled. \nMulti-level is reverted to FALSE."<< std::endl;
		probParamsObject.multiLvl = false;
	}

	std::cout << "Configuration file was read succesfully\n" << std::endl;
}


// function for removing whitespace of strings
void ReadConfigFile::findStringBounds(int& first, int& last, std::string& line){
	first = line.find("=")+1;
	last = line.size()-1;

	while(line[first] == ' '){
		first++;
	}

	while(line[last] == ' '){
		last = last -1;
	}
	last++;

	if(first == last){
		last++;
	}
}

// function for reading the provided tags in the configuration file
void ReadConfigFile::getTags(std::string& line, std::vector<std::string>& tagVec){
	int first = line.find('(')+1;
	int last = line.find(')');
	std::stringstream ss(line.substr(first,last-first));

	while(ss.good()){
		std::string marker;
		getline(ss, marker, ',');
		first = 0;
		last = marker.size()-1;

		while(marker[first] == ' '){
			first++;
		}
		while(marker[last-1] == ' '){
			last = last-1;
		}
		tagVec.push_back(marker.substr(first,last-first+1));
	}
}

// function for setting the bools according to the input in the configuration file
void ReadConfigFile::getBool(std::string& param, bool& boolean, std::string& line){
	int first, last;
	findStringBounds(first,last, line);
	std::string sstr = line.substr(first,last-first);

	try{
		if(sstr == "YES" || sstr == "yes"){
			boolean = true;
		}else if( sstr == "NO" || sstr == "no"){
			boolean = false;
		}else{
			throw(param);
		}
	}
	catch(const std::string& x){
		std::cout << "Invalid entry for parameter: " << param << std::endl;
		std::exit(0);
	}
}

// function for obtaining the directory where the executable or source code is located
void ReadConfigFile::getDirectory(std::string& dir){
	char fname_buff[FILENAME_MAX];
	get_cwd(fname_buff, FILENAME_MAX );
	std::string cwd(fname_buff);
	dir = cwd;
}
