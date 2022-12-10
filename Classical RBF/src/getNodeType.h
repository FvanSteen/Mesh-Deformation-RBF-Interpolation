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
	Eigen::ArrayXi iNodes, mNodes, sNodes,mNodesStd,ibNodes,esNodes;

	Eigen::ArrayXi* mPtr;
	Eigen::ArrayXi* iPtr;
	Eigen::ArrayXi* iPtrGrdy;


	int N_i,N_m,N_s,N_mStd,N_ib,N_es,N_i_grdy;
	getNodeType(Mesh& meshOb, bool& dataRed);
	void assignNodeTypes();
	void addControlNode(int& node, std::string& smode);
//	void greedyNodes(int node, std::string smode);
	void assignNodeTypesGreedy();
};

#endif /* GETNODETYPE_H_ */
