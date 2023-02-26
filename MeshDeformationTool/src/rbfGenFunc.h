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

	Eigen::VectorXd pVec, pnVec;
	Eigen::ArrayXi movingIndices;
	Eigen::ArrayXXd exactDisp;

	Eigen::ArrayXi* dispIdx;
	Eigen::ArrayXXd* disp;

	Eigen::ArrayXXd d;
	Eigen::VectorXd alpha;

	rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject);
	void getPhis(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::ArrayXi* mPtr, Eigen::ArrayXi* iPtr);
	void getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_bb, Eigen::MatrixXd& Phi_ib, Eigen::ArrayXi* cPtr, Eigen::ArrayXi* sPtr, Eigen::ArrayXi* bPtr, Eigen::ArrayXi* iPtr);
	void getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_sc,Eigen::MatrixXd&  Phi_ss, Eigen::MatrixXd& Phi_ic, Eigen::MatrixXd& Phi_is, getNodeType& n);
	void getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2);

//	void getDefVec(Eigen::VectorXd& defVec, getNodeType& n, int lvl, Eigen::ArrayXXd& errorPrevLvl);
	void getDefVec(Eigen::VectorXd& defVec, int N_c, Eigen::ArrayXi* cPtr);
	void getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef);
	void readDisplacementFile();

	void getNodeTypes();
	void getPeriodicParams();
	double rbfEval(double distance);
	void getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N, Eigen::ArrayXi*& mPtr);
	void performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int& N);
	void updateNodes(Eigen::MatrixXd& Phi_icGrdy, getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr);

};

#endif /* RBFGENFUNC_H_ */
