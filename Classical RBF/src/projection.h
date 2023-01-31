/*
 * projection.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "Mesh.h"
class projection{
public:
	projection();
	void project(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec);
	void projectIter(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, int N_se);
	void projectFun(Mesh& m,  Eigen::RowVectorXd& projection, Eigen::ArrayXXd& dist, int edge);
};

#endif /* PROJECTION_H_ */
