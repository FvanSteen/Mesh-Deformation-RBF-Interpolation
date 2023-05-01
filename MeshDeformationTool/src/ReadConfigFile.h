
#ifndef READCONFIGFILE_H_
#define READCONFIGFILE_H_
#include "probParams.h"
#include <string>
#include <vector>
class ReadConfigFile {
public:
	std::string mesh_ifName,mesh_ofName,pDir,curved,smode,pmode;
	std::vector<std::string> bdryTags, mTags, pTags;
	int steps;
	double tol, rFac, gamma;
	bool dataRed;

	ReadConfigFile(std::string& ifName, probParams& probParamsObject);
	void findStringBounds(int& first, int& last, std::string& line);
	void getTags(std::string& line, std::vector<std::string>& tagVec);
	void getBool(std::string& param, bool& boolean, std::string& line);
	void getDirectory(std::string& dir);
};

#endif /* READCONFIGFILE_H_ */
