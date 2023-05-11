
#ifndef WRITERESULTS_H_
#define WRITERESULTS_H_
#include <string>
#include "probParams.h"
class WriteResults {
public:
	WriteResults();
	void setIntResults(int step, int lvl, double maxErr, double meanErr, double time, std::string& fName, int N);
	void createConvHistFile(std::string& fName);
	void finalResult(probParams& p, long double& time, Eigen::Array3d& quals);
	bool existTest(std::string& fName);
};

#endif /* WRITERESULTS_H_ */
