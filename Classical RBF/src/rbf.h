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
	int N_m,N_s;
	int N_mPro; // subset used in the first projection step
	Eigen::ArrayXi mNodesPro;
	Eigen::Matrix3d rotMatX, rotMatY, rotMatZ;//todo differentiate between 2 and 3d!
	Eigen::Matrix2d rotMat;
	Eigen::ArrayXi mNodes, iNodes,sNodes;
	const std::string smode, pmode;
	Eigen::VectorXd dVec;
	Eigen::VectorXd rotVec;
	rbf(Mesh& meshOb,const double xDef, const double yDef, const double zDef, Eigen::VectorXd& rVec, const int steps, Eigen::RowVectorXd& rotationPnt, const std::string& slidingMode);

	void RBFMain();

	void RBF_standard();
	void performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2);
	void getDefVec(Eigen::VectorXd& defVec, int& N);
	void getRotationalMat();

	void RBF_DS();
	void RBF_DS_3D();
	void getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss);
	void performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha);
	void performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha);
	void getPhiDS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t);

	void RBF_PS();
	void performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec);
	//	double rbfEval(double distance);
private:


};

#endif /* RBF_H_ */
