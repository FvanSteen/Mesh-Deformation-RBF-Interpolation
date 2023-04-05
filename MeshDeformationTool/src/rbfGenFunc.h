
#ifndef RBFGENFUNC_H_
#define RBFGENFUNC_H_
#include "Mesh.h"
#include "probParams.h"
#include "getNodeType.h"
#include <Eigen/Dense>

class rbfGenFunc {
public:
	Mesh& m;
	probParams& params;


	Eigen::ArrayXi movingIndices;
	Eigen::ArrayXXd exactDisp;

	Eigen::ArrayXi* dispIdx;
	Eigen::ArrayXXd* disp;

	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;

	rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject);

	void getPhis(getNodeType& n, int iter);

	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2);

//	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, int lvl, Eigen::ArrayXXd& errorPrevLvl);
	void getDefVec(Eigen::VectorXd& defVec, int N_m, Eigen::ArrayXi* mPtr);
	void getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef);

	// for the multi levels
	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_m);
	void readDisplacementFile();

	void getNodeTypes();

	double rbfEval(double distance);
	void getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N, Eigen::ArrayXi* mPtr);
	void performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int& N);
	void updateNodes(getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);




	void getPhisReduced(getNodeType& n);
	void getPhisFull(getNodeType& n);
	double getDistance(int node1, int node2);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int newNode, int type);
	void getReducedPhi(Eigen::MatrixXd& Phi, getNodeType& n);

	void adjustPhi(Eigen::MatrixXd& Phi, getNodeType& n, int type);
	void getPhi2(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int idx, int type);
private:



protected:
	struct PhiStruct{
		Eigen::MatrixXd Phi_mm, Phi_me, Phi_em, Phi_ee, Phi_im;
		Eigen::MatrixXd Phi_ms, Phi_sm, Phi_ss, Phi_es,Phi_se;
		Eigen::MatrixXd Phi_cc, Phi_ic;
		Eigen::MatrixXd Phi_mc, Phi_ec, Phi_sc;
	};
	PhiStruct Phis;
	PhiStruct* PhiPtr;

};

#endif /* RBFGENFUNC_H_ */
