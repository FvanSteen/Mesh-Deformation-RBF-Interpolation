#include "ReadConfigFile.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
ReadConfigFile::ReadConfigFile(std::string& ifName, probParams& probParamsObject) {
	std::cout << "Reading the configuration file" << std::endl;
//	std::string ifName = "config_file.txt";
	std::string line;							// string containing line obtained by getline() function
	std::ifstream file("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\ConfigFiles\\" + ifName);

	int first, last;

	// Check if file is opened
	if (file.is_open()){
		//Obtain line
		while (getline(file, line)){

			if(line.rfind("MESH_FILENAME")==0){
				findStringBounds(first,last,line);
				mesh_ifName = line.substr(first,last-first);
//				probParamsObject.convHistFile = mesh_ifName.substr(0,mesh_ifName.find(".")) + "_convHist.txt";
			}
			else if(line.rfind("MESH_OUT_FILENAME")==0){
				findStringBounds(first,last,line);
				mesh_ofName = line.substr(first,last-first);
			}
			else if(line.rfind("MARKER_BDRY")==0){
				getTags(line,bdryTags);
			}
			else if(line.rfind("MARKER_MOVING")==0){
				getTags(line,mTags);
			}
			else if(line.rfind("MARKER_PERIODIC")==0){
				getTags(line,pTags);
			}
			else if(line.rfind("DEFORMATION_FILENAME")==0){
				findStringBounds(first,last,line);
				probParamsObject.dispFile =  line.substr(first,last-first);
			}

			else if(line.rfind("SLIDING_MODE")==0){
				findStringBounds(first,last,line);
				smode =  line.substr(first,last-first);
				probParamsObject.smode = smode;
			}
			else if(line.rfind("PERIODIC_MODE")==0){
				findStringBounds(first,last,line);
				pmode =  line.substr(first,last-first);
				probParamsObject.pmode = pmode;
			}
			else if(line.rfind("PERIODIC_DIRECTION")==0){
				findStringBounds(first,last,line);
				pDir = line.substr(first,last-first);
				probParamsObject.pDir = pDir;

			}
			else if(line.rfind("STEPS")==0){
				findStringBounds(first,last,line);
				probParamsObject.steps =  stoi(line.substr(first,last-first));
			}
			else if(line.rfind("CURVED")==0){
				findStringBounds(first,last,line);
				if(line.substr(first,last-first) == "YES"){
					probParamsObject.curved = true;
				}else{
					probParamsObject.curved = false;
				}
			}
			else if(line.rfind("DATA_REDUCTION")==0){
				findStringBounds(first,last,line);
				if(line.substr(first,last-first) == "YES" || line.substr(first,last-first) == "yes"){
					dataRed = true;
				}else{
					dataRed = false;
				}
				probParamsObject.dataRed = dataRed;
			}
			else if(line.rfind("DATA_RED_TOLERANCE")==0){
				findStringBounds(first,last,line);
				tol =  stod(line.substr(first,last-first));
				probParamsObject.tol = tol;
			}
			else if(line.rfind("INFLUENCE_FACTOR")==0){
				findStringBounds(first,last,line);
				rFac =  stod(line.substr(first,last-first));
			}
			else if(line.rfind("GAMMA")==0){
				findStringBounds(first,last,line);
				gamma =  stod(line.substr(first,last-first));
				probParamsObject.gamma = gamma;
			}
			else if(line.rfind("MULTILEVEL")==0){
				findStringBounds(first,last,line);
				if(line.substr(first,last-first) == "YES"){
					probParamsObject.multiLvl = true;
				}else{
					probParamsObject.multiLvl = false;
				}
			}
			else if(line.rfind("LEVELSIZE")==0){
				findStringBounds(first,last,line);
				if(probParamsObject.multiLvl){
					probParamsObject.lvlSize = stoi(line.substr(first,last-first));
				}else{
					probParamsObject.lvlSize = 1;
				}
				probParamsObject.lvlSizeInit = probParamsObject.lvlSize;
				probParamsObject.convHistFile = mesh_ifName.substr(0,mesh_ifName.find(".")) + "_m" + std::to_string(probParamsObject.lvlSize) + ".txt";
			}
		}
		file.close();
	}else{std::cout << "Unable to open the configuration file: " << ifName << std::endl;}

	std::cout << "Done reading the configuration file" << std::endl;
}

void ReadConfigFile::findStringBounds(int& first, int& last, std::string& line){
	first = line.find("=")+1;
	last = line.size()-1;

	while(line[first] == ' '){
		first++;
	}

	while(line[last-1] == ' '){
		last = last -1;
	}
	last++;

	if(first == last){
		last++;
	}
}

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


