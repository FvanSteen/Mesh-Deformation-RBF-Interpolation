
#ifndef WRITERESULTS_H_
#define WRITERESULTS_H_
#include <string>

class WriteResults {
public:
	WriteResults();
	void setIntResults(int step, int lvl, double maxErr, double meanErr, double time, std::string& fName, int N);
	void createConvHistFile(std::string& fName);
};

#endif /* WRITERESULTS_H_ */
