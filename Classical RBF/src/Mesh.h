#include <string>
#include <Eigen/Dense>

#ifndef MESH_H_
#define MESH_H_

class Mesh{
public:
	Mesh(std::string fileName, double supportRadius);
	std::string fName;
	int pntsIdx, nPnts, nDims,nElem, bcNodesIdx, nBdryNodes;
	Eigen::MatrixXd coords;
	Eigen::ArrayXi bdryNodes;
	Eigen::ArrayXi intNodes;
	Eigen::MatrixXd newCoords;
	double r; // support radius
	void findProblemChars();
	void obtainCoords();
	Eigen::MatrixXd interpMat(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2);
	double rbfEval(double distance);
	Eigen::VectorXi UniqueElems();
	Eigen::ArrayXi obtainIntNodes();
	void updateNodes(Eigen::VectorXd dVec, Eigen::VectorXd disp);
	void WriteMeshFile(std::string fName);

};

#endif /* MESH_H_ */
