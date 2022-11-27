/*
 * probParams.h
 *
 *  Created on: 22 nov. 2022
 *      Author: floyd
 */

#ifndef PROBPARAMS_H_
#define PROBPARAMS_H_
#include <Eigen/Dense>
#include <string>
struct probParams{
	Eigen::VectorXd dVec;
	int steps;
	Eigen::VectorXd rotVec;
	Eigen::RowVectorXd rotPnt;
	std::string sMode;
	std::string pMode;
	bool curved;
	std::string pDir;
	bool dataRed;
	double tolerance;
};


#endif /* PROBPARAMS_H_ */
