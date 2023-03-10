/*
 * rbfGenFunc.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

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

	Eigen::VectorXd periodicVec, periodicNormalVec1,periodicNormalVec2;
	Eigen::ArrayXi movingIndices;
	Eigen::ArrayXXd exactDisp;

	Eigen::ArrayXi* dispIdx;
	Eigen::ArrayXXd* disp;

	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;

	rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject);
	void getPhis(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::ArrayXi* mPtr, Eigen::ArrayXi* iPtr);
	void getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_bb, Eigen::MatrixXd& Phi_ib, Eigen::ArrayXi* cPtr, Eigen::ArrayXi* sPtr, Eigen::ArrayXi* bPtr, Eigen::ArrayXi* iPtr);
	void getPhis(getNodeType& n);
	void getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ce, Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_ec, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::MatrixXd& Phi_ic, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is , getNodeType& n);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2);

//	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, int lvl, Eigen::ArrayXXd& errorPrevLvl);
	void getDefVec(Eigen::VectorXd& defVec, int N_c, Eigen::ArrayXi* cPtr);
	void getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef);

	// for the multi levels
	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_c);
	void readDisplacementFile();

	void getNodeTypes();
	void getPeriodicParams();
	double rbfEval(double distance);
	void getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N, Eigen::ArrayXi* mPtr);
	void performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int& N);
	void updateNodes(getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);
private:



protected:
	struct PhiStruct{
		Eigen::MatrixXd Phi_cc, Phi_cs, Phi_sc, Phi_ss, Phi_ic, Phi_is;
		Eigen::MatrixXd Phi_ce, Phi_ec, Phi_es,Phi_ee,Phi_se,Phi_ie;
		Eigen::MatrixXd Phi_bb, Phi_ib;
	};
	PhiStruct Phis;
	PhiStruct* PhiPtr;

};

#endif /* RBFGENFUNC_H_ */
