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
	int N_m, N_se;
	Eigen::Matrix2d rotMat;
	Eigen::ArrayXi mNodes;
	const std::string mode;
	// todo can rotationPnt be passed by reference?
	rbf(Mesh& meshOb,const double xDef, const double yDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd rotationPnt, const std::string& mode);
//	Eigen::MatrixXd newCoords;

	void performRbfDS();
	void performRbfInterpolation();

	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2);
	void getDefVec(Eigen::VectorXd& defVec);
	void getRotDef();
	// todo check the rbf eval function
	double rbfEval(double distance);
	void getDisplacement(Eigen::MatrixXd& Phi, Eigen::VectorXd& a_x, Eigen::VectorXd& a_y, Eigen::VectorXd& defVec);
	void getPhiDS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);
	void getDisplacementDS(Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is,Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& alpha, Eigen::VectorXd& defVec);

	void performRbfPS();

private:


};

#endif /* RBF_H_ */
