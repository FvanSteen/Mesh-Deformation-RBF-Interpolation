#include <string>
#include <Eigen/Dense>
#include <vector>
#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(std::string fileName, double supportRadius);
	std::string fName;
	int pntsIdx, nPnts, nDims,nElem, bcNodesIdx, nBdryNodes, intNodesIdx;
	Eigen::MatrixXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXi intNodes;
	Eigen::MatrixXd newCoords;
	double r; // support radius
	void findProblemChars(std::vector<std::string> ibTags,std::vector<std::string> ebTags);
	void obtainCoords();
	Eigen::MatrixXd interpMat(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2);
	double rbfEval(double distance);
	Eigen::VectorXi UniqueElems(Eigen::ArrayXi arr);
	Eigen::ArrayXi obtainIntNodes();
	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void WriteMeshFile(std::string fName);
	void wmf(std::string ifName,std::string ofName);

};

#endif /* MESH_H_ */
