#ifndef MESH_H_
#define MESH_H_

#include <string>
#include <Eigen/Dense>
#include <vector>

#include "ProbParams.h"
class Mesh
{
public:
	Mesh(probParams& params);

	int nNodes, nDims;

	Eigen::ArrayXXd coords, coords_polar_cylindrical;
	Eigen::ArrayXXd surfMidPnts, edgeMidPnts,edgeMidPnts_polar_cylindrical, surfMidPnts_polar_cylindrical;

	Eigen::ArrayXXd* ptrCoords;
	Eigen::ArrayXXd* edgeMidPntPtr;
	Eigen::ArrayXXd* surfMidPntPtr;

	Eigen::ArrayXXd t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss;
	Eigen::ArrayXXd surfMidPntNormals, edgeMidPntNormals1,edgeMidPntNormals2, periodicEdgeNormals;

	Eigen::ArrayXi mNodes,seNodes,iNodes, ssNodes, verticesNodes, periodicVerticesNodes, periodicEdgeNodes,internalEdgeNodes;

	Eigen::MatrixXd periodicVecs;


	double r,periodic_length;
	int N_i, N_nonzeroDisp, N_se, N_ss, N_m, N_periodic_vertices, N_pe;

	void getVecs(probParams& params);
	void getMidPnts(probParams& params);

	void writeMeshFile(std::string& directory, std::string& ifName, std::string& ofName);


private:



	Eigen::ArrayXXi edgeConnectivity, edgeConnectivityPeriodic, surfConnectivity, bdryNodesMat, extBdryEdgeSegments;
	Eigen::ArrayXi nrElemsBdry;
	std::vector<std::string> srtdTags;

	void readMeshFile(probParams& params);

	void getNodeTypes(probParams& params);
	void getIntNodes();
	void removeDuplicates(Eigen::ArrayXi& arr);

	void getEdgeConnectivity(std::string& pmode, Eigen::ArrayXXi& edgeConnectivity,Eigen::ArrayXi& seNodes, int size);
	void getSurfConnectivity();
	void findStringBounds(int& first, int& last, std::string& line);

	void getPeriodicVector(probParams& params);

	void removeMutualNodes(Eigen::ArrayXi& array_in, int& size, Eigen::ArrayXi& to_remove_nodes, int end_idx);
	void getCharPerLength(probParams& params);
	double getCharDomLength();

	void getPerpVecs(Eigen::ArrayXXd& vecs, Eigen::ArrayXXd& p1, Eigen::ArrayXXd& p2);
	void getEdgeTan(Eigen::ArrayXXd& t,Eigen::ArrayXXi& edgeConnectivity, Eigen::ArrayXi& seNodes);
	void getSurfNormal();
	void getSurfNormalPeriodic(int directSliding);


	void getExtBdryEdgeSegments();

};

#endif /* MESH_H_ */
