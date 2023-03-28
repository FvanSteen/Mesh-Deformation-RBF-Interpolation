
#ifndef SPDS_H_
#define SPDS_H_
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <vector>
#include "Mesh.h"
#include "getNodeType.h"
#include "nanoflann.hpp"


class SPDS {
public:

	SPDS();

	void kdt_NNSearch(Eigen::ArrayXi& bdryIndex, Eigen::ArrayXi& intIndex, Eigen::ArrayXXd& coords, const size_t dim, double& gamma, double&maxError, Eigen::ArrayXXd* ePtr);
	double rbfEval(double distance, double radius);

	void project(Mesh& m, getNodeType& n, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec);
	void projectEdge(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec, size_t startIdx, size_t endIdx, int project);
	void projectSurf(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec, size_t startIdx, size_t endIdx, int project);

};

#endif /* SPDS_H_ */


