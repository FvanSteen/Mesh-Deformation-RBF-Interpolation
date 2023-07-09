#ifndef RBFSTD_H_
#define RBFSTD_H_
#include "Mesh.h"
#include "rbfGenFunc.h"
#include "getNodeType.h"
#include "greedy.h"
#include "ProbParams.h"
#include "WriteResults.h"


class rbf_std : public rbfGenFunc
{
public:
	CoordTransform transform;
	WriteResults w;
	rbf_std(struct probParams& probParamsObject, Mesh& meshObject,  getNodeType& n);
	void perform_rbf(getNodeType& n);
	void perform_rbf(getNodeType& n, greedy& g);

private:
	Eigen::VectorXd defVec;
};

#endif /* RBFSTD_H_ */
