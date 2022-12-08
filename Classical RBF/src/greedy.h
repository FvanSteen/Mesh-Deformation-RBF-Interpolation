/*
 * greedy.h
 *
 *  Created on: 17 nov. 2022
 *      Author: floyd
 */

#ifndef GREEDY_H_
#define GREEDY_H_
#include "rbfps.h"
#include "getNodeType.h"
#include "Mesh.h"
class greedy {
public:
	greedy();
	void getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, double& e, int& idxMax, std::string sMode, Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement);

};

#endif /* GREEDY_H_ */
