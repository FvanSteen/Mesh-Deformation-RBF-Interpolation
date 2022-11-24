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
};

#endif /* PROJECTION_H_ */
