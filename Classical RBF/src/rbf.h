#include "Mesh.h"
#include <string>
#include <Eigen/Dense>
#include <vector>

#ifndef RBF_H_
#define RBF_H_

class rbf {
public:
	Mesh m;
	rbf(Mesh meshOb);
	Eigen::MatrixXd newCoords;

	void performRbfInterpolation(const double& xDef, const double& yDef, const double& rotDefDeg);
	Eigen::MatrixXd getPhi(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2);
	Eigen::MatrixXd getDefVec(double xDef, double yDef, double rotDefDeg);
	Eigen::MatrixXd getRotDef(double rotDef);
	double rbfEval(double distance);
	void updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp);

private:


};

#endif /* RBF_H_ */
