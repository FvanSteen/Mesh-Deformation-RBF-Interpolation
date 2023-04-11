#ifndef GETNODETYPE_H_
#define GETNODETYPE_H_
#include "Mesh.h"
#include "probParams.h"
#include <Eigen/Dense>

class getNodeType {
public:


	Eigen::ArrayXi* mPtr;
	Eigen::ArrayXi* cPtr;
	Eigen::ArrayXi* iPtr;
	Eigen::ArrayXi* iPtrGrdy;
	Eigen::ArrayXi* sePtr;
	Eigen::ArrayXi* ssPtr;
	Eigen::ArrayXi* iPtr_reduced;

	Eigen::ArrayXi iNodesIdx, cNodesIdx;


	int N_i,N_se,N_iGrdy, N_ss;
//	int N_c, N_b;
	int N_m, N_c;

	struct addedNodesData{
		Eigen::ArrayXi idx, idx_i, type;
	};
	addedNodesData addedNodes;
	getNodeType(probParams& params, Mesh& m);
	void assignNodeTypes(Mesh& m);

	void addControlNodes(Eigen::ArrayXi& nodes, std::string& smode, Mesh& m);
	void addControlNode(int node, std::string& smode, Mesh& m, int i);

	void assignNodeTypesGrdy(Mesh& m);
	void assignNodeTypesGrdy(Mesh& m, std::string& smode);
private:
	bool pseudo;
	Eigen::ArrayXi iNodes, cNodes, mNodes, seNodes,ssNodes, iNodesReduced;


};

#endif /* GETNODETYPE_H_ */
