#include "Mesh.h"
#include <string>
#include <Eigen/Dense>
#include <vector>

#ifndef RBF_H_
#define RBF_H_

class rbf {
public:
	Mesh m;
	const double dx, dy, dz;
	Eigen::RowVectorXd rotPnt;
	const int steps;
	int N_m;
	int N_mPro; // subset used in the first projection step
	Eigen::ArrayXi mNodesPro;
	Eigen::Matrix2d rotMat;
	Eigen::ArrayXi mNodes;
	const std::string mode;
	Eigen::VectorXd dVec;
	rbf(Mesh& meshOb,const double xDef, const double yDef, const double zDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd& rotationPnt, const std::string& slidingMode);

	void RBFMain();

	void RBF_standard();
	void performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2);
	void getDefVec(Eigen::VectorXd& defVec, int& N);

	void RBF_DS();
	void performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha);
	void getPhiDS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);

	void RBF_PS();
	void performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& n, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec);
	//	double rbfEval(double distance);
private:


};

#endif /* RBF_H_ */
