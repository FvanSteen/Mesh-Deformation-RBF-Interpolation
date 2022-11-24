/*
 * rbfGenFunc.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef RBFGENFUNC_H_
#define RBFGENFUNC_H_
#include "Mesh.h"
#include "probParams.h"
#include "getNodeType.h"
#include <Eigen/Dense>
class rbfGenFunc {
public:
	Mesh& m;
	probParams& params;

	Eigen::VectorXd pVec, pnVec;
	Eigen::Matrix3d rotMatX, rotMatY, rotMatZ;//todo differentiate between 2 and 3d!
	Eigen::Matrix2d rotMat;



	rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2);
	void getDefVec(Eigen::VectorXd& defVec, int& N, Eigen::ArrayXi& ibNodes);
	void getRotationalMat();
	void getNodeTypes();
	void getPeriodicParams();
	double rbfEval(double distance);

};

#endif /* RBFGENFUNC_H_ */
