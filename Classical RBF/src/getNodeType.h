/*
 * getNodeType.h
 *
 *  Created on: 19 nov. 2022
 *      Author: floyd
 */

#ifndef GETNODETYPE_H_
#define GETNODETYPE_H_
#include "Mesh.h"
#include <Eigen/Dense>

class getNodeType {
public:
	Mesh& m;
	Eigen::ArrayXi iNodes,mNodes,sNodes,mNodesStd,ibNodes,esNodes;
	int N_i,N_m,N_s,N_mStd,N_ib,N_es;
	getNodeType(Mesh& meshOb);
	void assignNodeTypes();
	void greedyNodes(int node, std::string smode);
	void GreedyInit();

};

#endif /* GETNODETYPE_H_ */
