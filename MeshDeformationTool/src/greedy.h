#ifndef GREEDY_H_
#define GREEDY_H_
#include "getNodeType.h"
#include "Mesh.h"
#include "CoordTransform.h"
#include "SPDS.h"
class greedy {
public:


	double maxError;
	Eigen::ArrayXi maxErrorNodes;
	Eigen::ArrayXi* errIdxPtr;

	Eigen::ArrayXXd delta;
	Eigen::ArrayXXd deltaInternal;

	Eigen::ArrayXXd error, errorPolarCylindrical;
	Eigen::ArrayXXd errorPrevLvl;

	Eigen::VectorXd alphaGrdy;
	Eigen::ArrayXi ctrlNodesAll;
	Eigen::MatrixXd alphaTotal;

	double maxErrorPrevLvl;

	Eigen::VectorXd* alpha_step;
	Eigen::ArrayXXd* d_step;
	Eigen::ArrayXi* ctrlPtr;

	Eigen::ArrayXXd* exctDispPtr;
	Eigen::ArrayXi* mIdxPtr;

	probParams* p;

	Mesh* mPtr;
	CoordTransform transform;
	SPDS SPDSobj;

	greedy(Mesh& m, probParams& params, Eigen::ArrayXXd* disp, Eigen::ArrayXi& movingIndices,  Eigen::VectorXd& alpha, Eigen::ArrayXXd& d);


	void getError(getNodeType& n, Eigen::ArrayXXd& d, int lvl);
	void getErrorSingleLvl(getNodeType& n,  Eigen::ArrayXXd& d);
	void getErrorMultiLvl( getNodeType& n,  Eigen::ArrayXXd& d);

	void correction(Mesh& m, getNodeType& n, double& gamma, bool& multiLvl);


	int getDoubleEdgeError(int idxMax, int N_i, Eigen::ArrayXXd& error);

	void setLevelParams(getNodeType& n, int lvl, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha, Eigen::VectorXd& defVec, Eigen::ArrayXi* cPtr, int N_c);
	void setInitMaxErrorNodes();

	void getAlphaVector();


	void getErrorMovingNodes(Eigen::ArrayXi* nodes,Eigen::ArrayXXd& d,  size_t N);
	void getErrorAngle();
private:
	Eigen::ArrayXXd errorAngle;
};

#endif /* GREEDY_H_ */
