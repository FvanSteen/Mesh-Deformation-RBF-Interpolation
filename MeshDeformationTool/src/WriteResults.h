
#ifndef WRITERESULTS_H_
#define WRITERESULTS_H_
#include <string>

#include "ProbParams.h"
class WriteResults {
public:
	WriteResults();
	void setIntResults(std::string& dir, int step, int lvl, double& maxErr, long double& time, int N);
	void createConvHistFile(std::string& dir);
	void finalResult(probParams& p, long double& time, Eigen::Array3d& quals);
	bool existTest(std::string& fName);
};

#endif /* WRITERESULTS_H_ */
