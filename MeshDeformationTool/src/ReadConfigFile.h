#ifndef READCONFIGFILE_H_
#define READCONFIGFILE_H_
#include <string>
#include <vector>

#include "ProbParams.h"
class ReadConfigFile {

public:
	ReadConfigFile(std::string& ifName, probParams& probParamsObject);

private:
	void findStringBounds(int& first, int& last, std::string& line);
	void getTags(std::string& line, std::vector<std::string>& tagVec);
	void getBool(std::string& param, bool& boolean, std::string& line);
	void getDirectory(std::string& dir);
};

#endif /* READCONFIGFILE_H_ */
