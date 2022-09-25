#include <string>
#include <Eigen/Dense>
#include <vector>
#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(std::string fileName, double supportRadius,int debugLvl);
	std::string fName;
	int nPnts, nDims,nElem, nBdryNodes,lvl;
	Eigen::MatrixXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXi intNodes;
	Eigen::MatrixXd newCoords;
	double r; // support radius


	Eigen::MatrixXd interpMat(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2);
	double rbfEval(double distance);
	Eigen::VectorXi UniqueElems(Eigen::ArrayXi arr);
	Eigen::ArrayXi obtainIntNodes();
	void readMeshFile(std::vector<std::string> ibTags,std::vector<std::string> ebTags);
	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile(std::string ifName,std::string ofName);

};

#endif /* MESH_H_ */
