/*
 * rbfds.h
 *
 *  Created on: 17 nov. 2022
 *      Author: floyd
 */

#ifndef RBFDS_H_
#define RBFDS_H_
#include "Mesh.h"
#include "rbfGenFunc.h"

#include "getNodeType.h"
#include "greedy.h"


#include <string>
#include <Eigen/Dense>

class rbf_ds : public rbfGenFunc{
public:

//	projection* p;


	rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi, PhiStruct* PhiPtr, Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b);
	void getPhiDS(Eigen::MatrixXd& Phi,PhiStruct* PhiPtr ,getNodeType& n);
	void getIdxSlidingNodes(Eigen::ArrayXi* sPtr, Eigen::ArrayXi& idx, int type);
	void setDefVec_b(Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, getNodeType& n, PhiStruct* PhisPtr);
private:

};

#endif /* RBFDS_H_ */
