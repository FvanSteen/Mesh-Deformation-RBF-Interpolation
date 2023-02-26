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

	greedy();
	void getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes,  Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp,Eigen::VectorXd& pnVec, projection& p, bool multiLvl, int lvl, bool doubleEdge);
	void getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXd& errorAngle, Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection& p);
	void getErrorMultiLvl( getNodeType& n,  Eigen::ArrayXd& errorAngle,Eigen::ArrayXXd& d, Mesh& m, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& pnVec, projection& p);

	void correction(Mesh& m, getNodeType& n, double& gamma);
	void getNearestNode(Mesh& m, getNodeType& n,  int& node, int& idxMin, double& dist);
	double rbfEval(double distance, double radius);
	void project(Mesh& m, int& node, int& index, Eigen::ArrayXXd& disp, Eigen::VectorXd& pnVec, projection& p, int edge);
	int getDoubleEdgeError(Eigen::ArrayXd& errorAngle, int idxMax, int N_i, Eigen::ArrayXXd& error);
	void setLevelParams(Mesh& m, getNodeType& n, int& lvl, int& lvlSize, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha, double maxError);
	void setInitMaxErrorNodes(Mesh& m, Eigen::ArrayXXd& coords, Eigen::ArrayXXd& disp, Eigen::ArrayXi& mIdx,Eigen::ArrayXi& maxErrorNodes, bool doubleEdge);
	void setMaxErrorNodes(Mesh& m, Eigen::ArrayXi& maxErrorNodes);
	void getAlphaIdx(Eigen::ArrayXi& mNodes, Eigen::ArrayXi* mNodesGrdy, int N, Eigen::ArrayXi& idxAlpha);
	void getAlphaVector();

};

#endif /* GREEDY_H_ */
