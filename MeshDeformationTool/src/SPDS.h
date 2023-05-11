
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

	void project(Mesh& m, getNodeType& n, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, int ptype);
	void projectEdge(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, size_t startIdx, size_t endIdx, int project, int ptype);
	void projectSurf(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, size_t startIdx, size_t endIdx, int project, int ptype);

	void getDistNearestNeighbour(Eigen::ArrayXXd& midPnts,  Eigen::VectorXd& d, std::vector<size_t> idx, Eigen::ArrayXd& query, int ptype);
	void transformProjection(Eigen::RowVectorXd& project_i, Eigen::ArrayXd& query);
	void projectSlidingVertex(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int project, int ptype);
	void projectSlidingEdge(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int idxPerEdge, int project, int ptype);
};

#endif /* SPDS_H_ */


