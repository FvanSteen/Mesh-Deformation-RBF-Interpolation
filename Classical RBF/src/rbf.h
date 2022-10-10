#include "Mesh.h"
#include <string>
#include <Eigen/Dense>
#include <vector>

#ifndef RBF_H_
#define RBF_H_

class rbf {
public:
	Mesh m;
	const double dx, dy;
	Eigen::RowVectorXd rotPnt;
	const int steps;
	Eigen::Matrix2d rotMat;
	Eigen::ArrayXi mNodes;

	rbf(Mesh& meshOb,const double xDef, const double yDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd rotationPnt);
//	Eigen::MatrixXd newCoords;

	void performRbfDS();

	Eigen::MatrixXd getDefVecDS(double xDef, double yDef, double rotDefDeg, Eigen::VectorXd rotPnt,Eigen::ArrayXi intN);
	void performRbfInterpolation();
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2);

	void getDefVec(Eigen::VectorXd& defVec);
	void getRotDef();
	double rbfEval(double distance);
	void getDisplacement(Eigen::MatrixXd& Phi, Eigen::VectorXd& a_x, Eigen::VectorXd& a_y, Eigen::VectorXd& defVec);

private:


};

#endif /* RBF_H_ */
