
#ifndef READCONFIGFILE_H_
#define READCONFIGFILE_H_
#include <string>
#include <vector>
class ReadConfigFile {
public:
	ReadConfigFile(std::string& ifName);
	void findStringBounds(int& first, int& last, std::string& line);
	void getTags(std::string& line, std::vector<std::string>& tagVec);
};

#endif /* READCONFIGFILE_H_ */
