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
#include "projection.h"
class greedy {
public:
//	Eigen::ArrayXXd error;
	Eigen::ArrayXXd delta;
	Eigen::ArrayXXi mNodesHist;
	Eigen::ArrayXXd alphaHist;
	Eigen::ArrayXXd error;


	greedy();
	void getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes,  Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp,Eigen::VectorXd& pnVec, projection* projPtr);
	void getErrorMultiLvl( getNodeType& n,  Eigen::ArrayXXd& d, double& e, Eigen::ArrayXi& maxErrorNodes,  Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement,Eigen::VectorXd& pnVec);
	void correction(Mesh& m, getNodeType& n, double& gamma);
	void getNearestNode(Mesh& m, getNodeType& n,  int& node, int& idxMin, double& dist);
	double rbfEval(double distance, double radius);

	void project(Mesh& m, int& node, int& index, Eigen::ArrayXXd& disp, Eigen::VectorXd& pnVec, projection* projPtr);
	int getDoubleEdgeError(Eigen::ArrayXd& errorAngle, int idxMax, int N_i);
	void setLevelParams(Mesh& m, getNodeType& n, int& lvl, int& lvlSize, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha);

	void getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXd& errorAngle, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection* projPtr);

};

#endif /* GREEDY_H_ */
