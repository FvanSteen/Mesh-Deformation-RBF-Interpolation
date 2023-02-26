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
#include "projection.h"
#include "getNodeType.h"
#include "greedy.h"


#include <string>
#include <Eigen/Dense>

class rbf_ds : public rbfGenFunc{
public:

//	projection* p;


	rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_cc,Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, projection& p, Eigen::ArrayXi& iNodes,Eigen::ArrayXi& iNodesGrdy, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s);
	void getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_cc,Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& nVec, Eigen::ArrayXXd& tVec ,getNodeType& n);
	void getIdxSlidingNodes(Eigen::ArrayXi* sPtr, Eigen::ArrayXi& idx, int type);
	void getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, getNodeType& n);
	void performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, getNodeType& n);
	void setDefVec_b(Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, getNodeType& n, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_ss);
private:
//	Eigen::ArrayXi iNodes,mNodes,sNodes,mNodesStd;
//	int N_i,N_m,N_s,N_mStd;
};

#endif /* RBFDS_H_ */
