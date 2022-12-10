/*
 * greedy.h
 *
 *  Created on: 17 nov. 2022
 *      Author: floyd
 */

#ifndef GREEDY_H_
#define GREEDY_H_
#include "rbfps.h"
#include "getNodeType.h"
#include "Mesh.h"
class greedy {
public:
	Eigen::ArrayXXd error;
	greedy();
	void getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, double& e, int& idxMax, std::string sMode, Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement);
	void correction(Mesh& m, getNodeType& n);
	void getNearestNode(Mesh& m, getNodeType& n, int& node, int& idxMin, double& dist);
	double rbfEval(double distance, double radius);
	double getMaxError();
};

#endif /* GREEDY_H_ */
