/*
 * rbfstd.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef RBFSTD_H_
#define RBFSTD_H_
#include "Mesh.h"
#include "probParams.h"
#include "rbf.h"
#include "rbfGenFunc.h"
#include "getNodeType.h"



class rbf_std : public rbfGenFunc
//class rbf_std
{
public:

//	Eigen::ArrayXi iNodes,mNodes;
//	int N_i, N_m;
	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;
//	rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir);
	rbf_std(Mesh& meshObject, struct probParams& probParamsObject);
	virtual void perform_rbf(getNodeType& n);
	void performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes, int& N);
	void updateNodes(Eigen::MatrixXd& Phi_imGreedy, getNodeType& n, Eigen::VectorXd& defVec);
	void getExactDef(getNodeType& n, Eigen::VectorXd& exactDeformation);
	virtual ~rbf_std(){};
};

#endif /* RBFSTD_H_ */
