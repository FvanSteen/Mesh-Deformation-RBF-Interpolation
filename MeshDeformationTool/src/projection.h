/*
 * projection.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "Mesh.h"
#include "getNodeType.h"
class projection{
public:

	projection(Eigen::VectorXd& pVec);
	void project(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec);
	void projectIter(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, int N_se);
	void projectFun(Mesh& m,  Eigen::RowVectorXd& projection, Eigen::ArrayXXd& dist, int edge);

	void kd_tree(Mesh& m, getNodeType& n, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec);
private:
	Eigen::VectorXd pVec;
};

#endif /* PROJECTION_H_ */
