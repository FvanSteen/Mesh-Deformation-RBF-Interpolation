#include "ReadConfigFile.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
ReadConfigFile::ReadConfigFile(std::string& ifName) {
	std::cout << "Reading the displacement file" << std::endl;
//	std::string ifName = "config_file.txt";
	std::string line;							// string containing line obtained by getline() function
	std::ifstream file("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\ConfigFiles\\" + ifName); 	//opening file name stored in mFile object

	std::string dispFile,mesh_ifName,mesh_ofName,pDir,curved,dataRed;
	std::vector<std::string> bdryTags, mTags, pTags;
	int steps;
	double tol;

	int first, last;

	// Check if file is opened
	if (file.is_open()){
		//Obtain line
		while (getline(file, line)){

			if(line.rfind("MESH_FILENAME")==0){
				findStringBounds(first,last,line);
				mesh_ifName = line.substr(first,last-first);
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
				dispFile =  line.substr(first,last-first);
			}
			else if(line.rfind("PERIODIC_DIRECTION")==0){
				std::cout << line << std::endl;
				findStringBounds(first,last,line);
				pDir = line.substr(first,last-first);
				std::cout << "periodic Direction is " << pDir << std::endl;
			}
			else if(line.rfind("STEPS")==0){
				findStringBounds(first,last,line);
				steps =  stoi(line.substr(first,last-first));
			}
			else if(line.rfind("CURVED")==0){
				findStringBounds(first,last,line);
				curved =  line.substr(first,last-first);
			}
			else if(line.rfind("DATA_REDUCTION")==0){
				findStringBounds(first,last,line);
				dataRed =  line.substr(first,last-first);
			}
			else if(line.rfind("DATA_RED_TOLERANCE")==0){
				findStringBounds(first,last,line);
				tol =  stod(line.substr(first,last-first));
			}


		}
		std::cout << dispFile << std::endl;
		std::cout << mesh_ifName << std::endl;
		std::cout << mesh_ofName << std::endl;
		for(int x = 0; x<5; x++){
			std::cout << bdryTags[x] << std::endl;
		}

		std::cout << mTags[0] << std::endl;
		std::cout << pTags[0] << std::endl;
		std::cout << pTags[1] << std::endl;
		std::cout << std::endl;
		std::cout << pDir << std::endl;
		std::cout << steps << std::endl;
		std::cout << curved << std::endl;
		std::cout << dataRed << std::endl;
		std::cout << tol << std::endl;

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




