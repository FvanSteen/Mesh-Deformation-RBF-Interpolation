
#ifndef RBFPS_H_
#define RBFPS_H_
#include "Mesh.h"

#include "getNodeType.h"
#include "greedy.h"
#include "rbfGenFunc.h"
#include <string>
#include <Eigen/Dense>
#include "SPDS.h"
#include "CoordTransform.h"
#include "WriteResults.h"


class rbf_ps : public rbfGenFunc
{
public:

	SPDS p;
	CoordTransform transform;
	WriteResults w;
	rbf_ps(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void perform_rbf(getNodeType& n, greedy& g);
	void pseudo_sliding_edge(PhiStruct* PhiPtr, getNodeType& n);
	void pseudo_sliding_surf(PhiStruct* PhiPtr, getNodeType& n);

private:

	Eigen::VectorXd defVec_m, defVec_me, defVec_all;
	Eigen::ArrayXXd free_disp_edge, proj_disp_edge, free_disp_all, proj_disp_all;
};


#endif /* RBFPS_H_ */
