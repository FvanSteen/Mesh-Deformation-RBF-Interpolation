
#ifndef RBFSTD_H_
#define RBFSTD_H_
#include "Mesh.h"
#include "probParams.h"
#include "rbf.h"
#include "rbfGenFunc.h"
#include "getNodeType.h"
#include "projection.h"





class rbf_std : public rbfGenFunc
//class rbf_std
{
public:


	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;

	rbf_std(Mesh& meshObject, struct probParams& probParamsObject);
	virtual void perform_rbf(getNodeType& n);
	void performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes, int& N);
	void updateNodes(Eigen::MatrixXd& Phi_imGreedy, getNodeType& n, Eigen::VectorXd& defVec , Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);
	virtual ~rbf_std(){};
};

#endif /* RBFSTD_H_ */
