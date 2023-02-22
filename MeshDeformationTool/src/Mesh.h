#ifndef MESH_H_
#define MESH_H_
//#include "ReadConfigFile.h"
#include <string>
#include <Eigen/Dense>
#include <vector>
#include "probParams.h"
class Mesh
{
public:
//	Mesh(const std::string& inputFileName,const std::string& outputFileName, const std::vector<std::string>& Tags,const double& rFac,const int& debugLvl, const std::string& slidingMode, const std::string& periodicMode,const std::vector<std::string>& movingTags, const std::vector<std::string>& periodicTags);
	Mesh(probParams& params, const int& debugLvl);
//	const std::string ifName, ofName;
//	const std::vector<std::string>& Tags;
	int nNodes, nDims,nElem, lvl;
	Eigen::ArrayXXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXXi extBdryNodesMat,bdryNodesMat;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXXi intBdryNodesMat;
	Eigen::ArrayXi intNodes;
	Eigen::ArrayXXd nVecs;
	Eigen::ArrayXXd tVecs;
	Eigen::ArrayXXd midPnts, edgeMidPnts;
	Eigen::ArrayXi extStaticNodes, slidingSurfNodes, slidingEdgeNodes, periodicNodes;
	Eigen::ArrayXi nrElemsBdry;
	Eigen::ArrayXXi edgeConnectivity, surfConnectivity;
	Eigen::ArrayXXd n, t;
	Eigen::ArrayXXd t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss;
	Eigen::ArrayXi extBdryEdgeNodes;
	Eigen::ArrayXXi extBdryEdgeSegments;
	Eigen::ArrayXXd midPntNormals, edgeMidPntNormals1,edgeMidPntNormals2;
	Eigen::ArrayXi mNodes,seNodes,iNodes, staticNodes, ssNodes;
	Eigen::ArrayXi ibIndices, ebIndices;
	Eigen::ArrayXi mNodesStd;

	Eigen::ArrayXi intCorNodes;

//	struct Domains{
//		Eigen::ArrayXi d1;
//		Eigen::ArrayXi d2;
//		Eigen::ArrayXi d3;
//		Eigen::ArrayXi d4;
//	};
//
//	struct Domains subDoms;


	double r; // support radius
	int N_i, N_ib, N_eb, N_se, N_es, N_ss, N_p,N_m, N_mStd;
//	const std::string smode, pmode;
//	const std::vector<std::string> perBdry;
//	const std::vector<std::string> mTags,pTags;
	std::vector<std::string> srtdTags;
	double lambda;

	double charLength();
	void removeDuplicates(Eigen::ArrayXi& arr);
	void getIntNodes();
	void readMeshFile(probParams& params);
//	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile(std::string& ifName, std::string& ofName);
	void getExtBdryData();
	void getBdryNodes(Eigen::ArrayXXi& bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int& nBdryNodes, int& nBdryElems);
	void getNodeVecs(Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);
	void getNormals(Eigen::ArrayXXd& n);
	void getNodeTypes(probParams& params);
	void getNodeType(Eigen::ArrayXi& nrElemsExtBdry, const std::vector<std::string>& ebTags);
	void getEdgeConnectivity();
	void getSurfConnectivity();
	void getVecs3D(Eigen::ArrayXXd& t_se, Eigen::ArrayXXd& n1_se, Eigen::ArrayXXd& n2_se, Eigen::ArrayXXd& n_ss, Eigen::ArrayXXd& t1_ss, Eigen::ArrayXXd& t2_ss);
	void getPerpVecs(Eigen::ArrayXXd& vecs, Eigen::ArrayXXd& p1, Eigen::ArrayXXd& p2);


	void getVecs();
	void getEdgeTan(Eigen::ArrayXXd& t);
	void getSurfNormal();
	void getExtBdryEdgeSegments();
	void getMidPnts(probParams& params);

	void getSubDomains(Eigen::ArrayXi& subDomains, Eigen::ArrayXi& subDomLen, Eigen::ArrayXi& subDomBdry, Eigen::ArrayXi& subDomBdryLen);
	void getIntCorNodes(double& gamma, double& tol);
	void getCharLength(std::string& pDir);

	void findStringBounds(int& first, int& last, std::string& line);
//	void getInternalCorrectionNodes(Eigen::ArrayXi& subDomains,  Eigen::ArrayXi& subDomsLen, Eigen::ArrayXXd& bdryCoord);

};

#endif /* MESH_H_ */
