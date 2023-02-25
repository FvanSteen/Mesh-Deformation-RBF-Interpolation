#ifndef GETNODETYPE_H_
#define GETNODETYPE_H_
#include "Mesh.h"
#include "probParams.h"
#include <Eigen/Dense>

class getNodeType {
public:
	Eigen::ArrayXi iNodes, mNodes, seNodes,mNodesStd,ibNodes,esNodes, sNodes,ssNodes;

	Eigen::ArrayXi* mPtr;
	Eigen::ArrayXi* iPtr;
	Eigen::ArrayXi* iPtrGrdy;
	Eigen::ArrayXi* sePtr;
	Eigen::ArrayXi* mStdPtr;
	Eigen::ArrayXi* sPtr;
	Eigen::ArrayXi* ssPtr;



	int N_i,N_m,N_se,N_mStd,N_ib,N_es,N_i_grdy, N_s, N_ss;



	getNodeType(probParams& params, Mesh& m);
	void assignNodeTypes(Mesh& m);
	void addControlNode(int node, std::string& smode, Mesh& m);

	void assignNodeTypesGrdy(Mesh& m);
};

#endif /* GETNODETYPE_H_ */
