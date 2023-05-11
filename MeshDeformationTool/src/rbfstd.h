
#ifndef RBFSTD_H_
#define RBFSTD_H_
#include "Mesh.h"
#include "probParams.h"
#include "rbfGenFunc.h"
#include "getNodeType.h"
#include "greedy.h"



class rbf_std : public rbfGenFunc
//class rbf_std
{
public:



	CoordTransform transform;
	rbf_std(struct probParams& probParamsObject, Mesh& meshObject,  getNodeType& n);
	void perform_rbf(getNodeType& n);
	void perform_rbf(getNodeType& n, greedy& g);

private:
	Eigen::VectorXd defVec;
//	void performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* movingNodes, Eigen::ArrayXi* internalNodes, int& N);
//	void updateNodes(Eigen::MatrixXd& Phi_icGreedy, getNodeType& n, Eigen::VectorXd& defVec , Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);
//	virtual ~rbf_std(){};
};

#endif /* RBFSTD_H_ */
