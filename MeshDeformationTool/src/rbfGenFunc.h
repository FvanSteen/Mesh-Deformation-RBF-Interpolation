
#ifndef RBFGENFUNC_H_
#define RBFGENFUNC_H_
#include "Mesh.h"
#include "getNodeType.h"
#include <Eigen/Dense>
#include "ProbParams.h"
#include "CoordTransform.h"

class rbfGenFunc {
public:
	Mesh& m;
	probParams& params;
	CoordTransform transform;


	Eigen::ArrayXi movingIndices;
	Eigen::ArrayXXd exactDisp;
	Eigen::ArrayXXd exactDisp_polar_cylindrical;
	Eigen::ArrayXXd* exactDispPtr;

	Eigen::ArrayXi* dispIdx;
	Eigen::ArrayXXd* disp;


	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;

	rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject);

	void getPhis(getNodeType& n, int iter);

	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2);


	void getDefVec(Eigen::VectorXd& defVec, int N_m, Eigen::ArrayXi* mPtr);
	void getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef, int N, int N_init);
	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_m);

	void readDisplacementFile();

	double getDistance(int node1, int node2);
	double rbfEval(double distance);

	void performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int N);
	void updateNodes(getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);


	void getPhisReduced(getNodeType& n);
	void getPhisFull(getNodeType& n);
	void getReducedPhi(Eigen::MatrixXd& Phi, getNodeType& n);
	void adjustPhi(Eigen::MatrixXd& Phi, getNodeType& n, int type);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int idx, int type);
private:



protected:
	struct PhiStruct{
		Eigen::MatrixXd Phi_mm, Phi_em;
		Eigen::MatrixXd Phi_cc, Phi_ic;
		Eigen::MatrixXd Phi_mc, Phi_ec, Phi_sc;
		Eigen::MatrixXd Phi_meme, Phi_sme, Phi_ic_reduced;
	};
	PhiStruct Phis;
	PhiStruct* PhiPtr;

};

#endif /* RBFGENFUNC_H_ */
