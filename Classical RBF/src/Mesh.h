#include <string>
#include <Eigen/Dense>
#include <vector>
#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(std::string inputFileName, std::string outputFileName, std::vector<std::string> ibTags,std::vector<std::string> ebTags, double rFac,int debugLvl);
	std::string ifName, ofName;
	int nNodes, nDims,nElem, lvl;
	Eigen::ArrayXXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXi intNodes;
	Eigen::ArrayXXd nVecs;
	Eigen::ArrayXXd tVecs;
	Eigen::ArrayXXd nVecsT;
	Eigen::ArrayXXd tVecsT;
	Eigen::ArrayXXd midPnts;
	Eigen::ArrayXXd midPntsT;
	Eigen::ArrayXi movingNodes, slidingNodes;
	double r; // support radius
	double charLength();

	Eigen::ArrayXi removeDuplicates(Eigen::ArrayXi& arr);
	void getIntNodes();
	void readMeshFile(std::vector<std::string> ibTags,std::vector<std::string> ebTags,double rFac);
//	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile(Eigen::MatrixXd& newCoords);
	void getNormals(Eigen::VectorXi nodes, int& cntExtElems);
	void getNodeVecs(Eigen::ArrayXi& idxs, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);
	void getSlidingNodes();
	void getNormalsTest(Eigen::ArrayXXi nodes);
	void getMovingNodes(std::vector<std::string> ebTags, Eigen::ArrayXXi extBdryNodesMat, Eigen::ArrayXi extBdryEndsIdx);
	void getBdryNodes(Eigen::ArrayXXi bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int nBdryNodes, int nBdryElems);
};

#endif /* MESH_H_ */
