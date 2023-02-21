#ifndef PROBPARAMS_H_
#define PROBPARAMS_H_
#include <Eigen/Dense>
#include <string>
struct probParams{
//	Eigen::VectorXd dVec;
	int steps;
//	Eigen::VectorXd rotVec;
//	Eigen::RowVectorXd rotPnt;
	std::string smode;
	std::string pmode;
	bool curved;
	std::string pDir;
	bool dataRed;
	double tol;
	std::string dispFile;
	double gamma;
	bool multiLvl;
	int lvlSize;
	int lvlSizeInit;
	std::string convHistFile;
};


#endif /* PROBPARAMS_H_ */
