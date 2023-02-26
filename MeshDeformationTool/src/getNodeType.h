#ifndef GETNODETYPE_H_
#define GETNODETYPE_H_
#include "Mesh.h"
#include "probParams.h"
#include <Eigen/Dense>

class getNodeType {
public:


	Eigen::ArrayXi* cPtr;
	Eigen::ArrayXi* bPtr;
	Eigen::ArrayXi* iPtr;
	Eigen::ArrayXi* iPtrGrdy;
	Eigen::ArrayXi* sePtr;
	Eigen::ArrayXi* mStdPtr;
	Eigen::ArrayXi* sPtr;
	Eigen::ArrayXi* ssPtr;

	Eigen::ArrayXi iNodesIdx, cNodesIdx;

	int N_i,N_m,N_se,N_mStd,N_ib,N_es,N_iGrdy, N_s, N_ss;
	int N_c, N_b;


	getNodeType(probParams& params, Mesh& m);
	void assignNodeTypes(Mesh& m);
	void addControlNode(int node, std::string& smode, Mesh& m);

	void assignNodeTypesGrdy(Mesh& m);
private:
	Eigen::ArrayXi iNodes, bNodes, cNodes, mNodes, seNodes,mNodesStd,ibNodes,esNodes, sNodes,ssNodes;

};

#endif /* GETNODETYPE_H_ */
