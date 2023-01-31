/*
 * rbfds.h
 *
 *  Created on: 17 nov. 2022
 *      Author: floyd
 */

#ifndef RBFDS_H_
#define RBFDS_H_
#include "Mesh.h"
#include "rbfstd.h"
#include "projection.h"
#include "getNodeType.h"
#include "greedy.h"


#include <string>
#include <Eigen/Dense>

class rbf_ds : public rbf_std{
public:

//	projection* p;


	rbf_ds(Mesh& meshObject, struct probParams& probParamsObject);
	void perform_rbf(getNodeType& n);
	void performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, projection* proPnt, Eigen::ArrayXi& iNodes,Eigen::ArrayXi& iNodesGrdy, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s);
	void getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t,int& N_m, int& N_s, Eigen::ArrayXi& sNodes);
	void getIdxSlidingNodes(Eigen::ArrayXi& sNodes, Eigen::ArrayXi& idx, int type);
	void getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, getNodeType& n);
	void performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, getNodeType& n);
private:
//	Eigen::ArrayXi iNodes,mNodes,sNodes,mNodesStd;
//	int N_i,N_m,N_s,N_mStd;
};

#endif /* RBFDS_H_ */
