#include <string>
#include <Eigen/Dense>
#include <vector>
#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(const std::string& inputFileName,const std::string& outputFileName, const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags,const double& rFac,const int& debugLvl);
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
	Eigen::ArrayXi movingNodes, slidingNodes;
	double r; // support radius


	double charLength();
	void removeDuplicates(Eigen::ArrayXi& arr);
	void getIntNodes();
	void readMeshFile(const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags,const double& rFac);
//	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile();
	void getSlidingNodes();
	void getExtBdryData();
	void getMovingNodes(const std::vector<std::string>& ebTags, Eigen::ArrayXi& extBdryEndsIdx);
	void getBdryNodes(Eigen::ArrayXXi& bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int& nBdryNodes, int& nBdryElems);
	void getNodeVecs(Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);
};

#endif /* MESH_H_ */
