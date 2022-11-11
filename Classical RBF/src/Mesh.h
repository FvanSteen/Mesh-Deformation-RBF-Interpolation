#include <string>
#include <Eigen/Dense>
#include <vector>
#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(const std::string& inputFileName,const std::string& outputFileName, const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags,const double& rFac,const int& debugLvl, const std::string& slidingMode, const std::vector<std::string>& periodicBdry, const std::string& periodicMode);
	const std::string ifName, ofName;

	int nNodes, nDims,nElem, lvl;
	Eigen::ArrayXXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXXi extBdryNodesMat;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXXi intBdryNodesMat;
	Eigen::ArrayXi intNodes;
	Eigen::ArrayXXd nVecs;
	Eigen::ArrayXXd tVecs;
	Eigen::ArrayXXd midPnts;
	Eigen::ArrayXi extStaticNodes, slidingSurfNodes, slidingEdgeNodes, periodicNodes;
	Eigen::ArrayXi nrElemsExtBdry;
	Eigen::ArrayXXi edgeConnectivity, surfConnectivity;
	Eigen::ArrayXXd n, t;
	Eigen::ArrayXXd t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss;
	Eigen::ArrayXi extBdryEdgeNodes;
	Eigen::ArrayXXi extBdryEdgeSegments;
	Eigen::ArrayXXd midPntNormals;
	double r; // support radius
	int N_i, N_ib, N_eb, N_se, N_es, N_ss, N_p;
	const std::string smode, pmode;
	const std::vector<std::string> perBdry;


	double charLength();
	void removeDuplicates(Eigen::ArrayXi& arr);
	void getIntNodes();
	void readMeshFile(const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags);
//	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile();
	void getExtBdryData();
	void getBdryNodes(Eigen::ArrayXXi& bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int& nBdryNodes, int& nBdryElems);
	void getNodeVecs(Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);
	void getNormals(Eigen::ArrayXXd& n);

	void getNodeType(Eigen::ArrayXi& nrElemsExtBdry, const std::vector<std::string>& ebTags);
	void getEdgeConnectivity(Eigen::ArrayXi& nrElemsExtBdry);
	void getSurfConnectivity();
	void getVecs3D(Eigen::ArrayXXd& t_se, Eigen::ArrayXXd& n1_se, Eigen::ArrayXXd& n2_se, Eigen::ArrayXXd& n_ss, Eigen::ArrayXXd& t1_ss, Eigen::ArrayXXd& t2_ss);
	void getPerpVecs(std::string& type);


	void getVecs();
	void getEdgeTan(Eigen::ArrayXXd& t);
	void getSurfNormal();
	void getExtBdryEdgeSegments(Eigen::ArrayXi& nrElemsExtBdry);
	void getMidPnts();

};

#endif /* MESH_H_ */
