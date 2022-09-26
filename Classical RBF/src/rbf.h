#include "Mesh.h"
#include <string>
#include <Eigen/Dense>
#include <vector>

#ifndef RBF_H_
#define RBF_H_

class rbf {
public:

	rbf(Mesh meshOb);
	Eigen::MatrixXd newCoords;
	void performRbfInterpolation();
	Eigen::MatrixXd getPhi(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2);
	double rbfEval(double distance);
	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);

private:
	Mesh m;

};

#endif /* RBF_H_ */
