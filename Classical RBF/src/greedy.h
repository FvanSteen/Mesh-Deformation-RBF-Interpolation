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
	void getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, double& e, Eigen::ArrayXi& maxErrorNodes, std::string sMode, Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement,Eigen::VectorXd& pnVec);
	void correction(Mesh& m, getNodeType& n, double& gamma);
	void getNearestNode(Mesh& m, getNodeType& n, int& node, int& idxMin, double& dist);
	double rbfEval(double distance, double radius);
	double getMaxError();
	void project(Mesh& m, int& node, int& index, Eigen::ArrayXXd& disp,Eigen::ArrayXXd& error,Eigen::VectorXd& pnVec);
	void getDoubleEdgeError(Eigen::ArrayXXd& error, Eigen::ArrayXd& errorAngle, int& idxMax, int& N_i, int& doubleEdgeMax);
};

#endif /* GREEDY_H_ */
