/*
 * greedy.h
 *
 *  Created on: 17 nov. 2022
 *      Author: floyd
 */

#ifndef GREEDY_H_
#define GREEDY_H_
#include "rbfps.h"
#include "getNodeType.h"
#include "Mesh.h"
#include "projection.h"
class greedy {
public:

	Eigen::ArrayXXd delta;
	Eigen::ArrayXXd deltaInternal;

	Eigen::ArrayXXd error;
	Eigen::ArrayXXd errorPrevLvl;

	Eigen::VectorXd alphaGrdy;
	Eigen::ArrayXi ctrlNodesAll;
	Eigen::MatrixXd alphaTotal;

	double maxErrorPrevLvl;


	Eigen::VectorXd* alpha_step;
	Eigen::ArrayXXd* d_step;
	Eigen::ArrayXi* ctrlPtr;


	greedy(probParams& params, Eigen::VectorXd& alpha, Eigen::ArrayXXd& d);
	void getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes,  Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp,Eigen::VectorXd& pVec, projection& p, bool multiLvl, int lvl, bool doubleEdge);
	void getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXXd& errorAngle, Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pVec, projection& p);
	void getErrorMultiLvl( getNodeType& n,  Eigen::ArrayXXd& errorAngle,Eigen::ArrayXXd& d, Mesh& m, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& pVec, projection& p);

	void correction(Mesh& m, getNodeType& n, double& gamma, bool& multiLvl);
	void getNearestNode(Mesh& m, getNodeType& n,  int& node, int& idxMin, double& dist);
	double rbfEval(double distance, double radius);
	void project(Mesh& m, int& node, int& index, Eigen::ArrayXXd& disp, Eigen::VectorXd& pnVec, projection& p, int edge);
	int getDoubleEdgeError(Eigen::ArrayXXd& errorAngle, int idxMax, int N_i, Eigen::ArrayXXd& error);
	void setLevelParams(Mesh& m, getNodeType& n, int lvl, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha, double maxError, Eigen::VectorXd& defVec, Eigen::ArrayXi* cPtr, int N_c);
	void setInitMaxErrorNodes(Mesh& m, Eigen::ArrayXXd& coords, Eigen::ArrayXXd& disp, Eigen::ArrayXi& mIdx,Eigen::ArrayXi& maxErrorNodes, bool doubleEdge);
	void setMaxErrorNodes(Mesh& m, Eigen::ArrayXi& maxErrorNodes);
	void getAlphaIdx(Eigen::ArrayXi& mNodes, Eigen::ArrayXi* mNodesGrdy, int N, Eigen::ArrayXi& idxAlpha);
	void getAlphaVector();


	void getErrorMovingNodes(Eigen::ArrayXi* nodes,Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, size_t N);
	void getErrorSlidingEdgeNodes(Mesh& m, Eigen::ArrayXi* nodes, Eigen::ArrayXXd& d, Eigen::VectorXd& pVec, size_t N, size_t startIdx);
	void getErrorSlidingSurfNodes(Mesh& m, Eigen::ArrayXi* nodes, Eigen::ArrayXXd& d, size_t N, size_t startIdx);
	void getErrorAngle(Eigen::ArrayXXd& errorAngle, size_t dims, size_t N);


	//todo can be deleted
	void getErrorSingleLvlOld(Mesh& m, getNodeType& n, Eigen::ArrayXXd& errorAngle, Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection& p);
};

#endif /* GREEDY_H_ */
