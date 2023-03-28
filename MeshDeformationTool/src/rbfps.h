
#ifndef RBFPS_H_
#define RBFPS_H_
#include "Mesh.h"

#include "getNodeType.h"
#include "greedy.h"
#include "rbfGenFunc.h"
#include <string>
#include <Eigen/Dense>
#include "SPDS.h"



//class rbf_ps : public rbf_std, public projection
class rbf_ps : public rbfGenFunc
{
public:
	// spatial partial data structere class
	SPDS p;

	rbf_ps(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void performRBF_PS(PhiStruct* PhiPtr, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec_b, getNodeType& n);


};


#endif /* RBFPS_H_ */
