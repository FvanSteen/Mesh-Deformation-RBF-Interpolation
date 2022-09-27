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
	Eigen::MatrixXd coords;
	Eigen::ArrayXi extBdryNodes;
	Eigen::ArrayXi bdryNodes; // Remove this one later
	Eigen::ArrayXi intBdryNodes;
	Eigen::ArrayXi intNodes;
	double r; // support radius
	double charLength();

	Eigen::VectorXi UniqueElems(Eigen::ArrayXi& arr);
	void obtainIntNodes();
	void readMeshFile(std::vector<std::string> ibTags,std::vector<std::string> ebTags,double rFac);
	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);
	void writeMeshFile(Eigen::MatrixXd& newCoords);

};

#endif /* MESH_H_ */
