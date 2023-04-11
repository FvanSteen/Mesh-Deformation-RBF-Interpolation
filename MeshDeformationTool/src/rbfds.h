#ifndef RBFDS_H_
#define RBFDS_H_
#include "Mesh.h"
#include "rbfGenFunc.h"

#include "getNodeType.h"
#include "greedy.h"
#include <string>
#include <Eigen/Dense>
#include "SPDS.h"

class rbf_ds : public rbfGenFunc{
public:


	SPDS p;

	rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void perform_rbf(getNodeType& n, greedy& g);
	void performRBF_DS(getNodeType& n, PhiStruct* PhiPtr);
	void getPhiDS(getNodeType& n, PhiStruct* PhiPtr);
	void getIdxSlidingNodes(Eigen::ArrayXi* sPtr, Eigen::ArrayXi& idx, Eigen::ArrayXi& sNodesInit);
	void setDefVec_all(getNodeType& n, PhiStruct* PhisPtr);
private:
	Eigen::VectorXd defVec_ds, defVec_all;
	Eigen::MatrixXd Phi;
};

#endif /* RBFDS_H_ */
