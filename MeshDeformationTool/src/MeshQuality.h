
#ifndef MESHQUALITY_H_
#define MESHQUALITY_H_
#include "probParams.h"
#include "Mesh.h"

class MeshQuality {
public:
	MeshQuality(probParams& p, Eigen::ArrayXXd& coords);
	void getDeformedMeshQual(Eigen::ArrayXXd& coords);
	Eigen::Array3d defQuals;
private:
	Eigen::ArrayXXd* alpha_ptr;
	Eigen::ArrayXXd alpha_init, alpha;
	Eigen::ArrayXXd lambda_11, lambda_22, lambda_33, lambda_12, lambda_23, lambda_13;
	Eigen::ArrayXd qual;
	Eigen::ArrayXXi k;
	int nElem, nDim;
	Eigen::ArrayXXi elems;
	Eigen::ArrayXi elemType;

	std::string fName;

	bool setQuadDir = false;
	bool quadClockWise;
	void getInitialMeshQualParams(probParams& p, Eigen::ArrayXXd& coords);

	void getMeshQual();
	void getElemConnectivity(probParams& p);
	void getQualParams(Eigen::ArrayXXd& coords);
	void writeQualFile();
	bool existTest(std::string& fName);
	int setElemTypeParams(int type);
	double get_f_skew(int type, int i);
	void getQuadRotation(int i, Eigen::ArrayXXd& coords);
};

#endif /* MESHQUALITY_H_ */
