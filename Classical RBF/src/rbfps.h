/*;
 * rbfps.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef RBFPS_H_
#define RBFPS_H_
#include "Mesh.h"
#include "rbfstd.h"
#include "projection.h"
#include "getNodeType.h"
#include "greedy.h"
#include <string>
#include <Eigen/Dense>




//class rbf_ps : public rbf_std, public projection
class rbf_ps : public rbf_std
{
public:
//	Eigen::ArrayXi mNodes,iNodes,sNodes,mNodesPro;
//	int N_m,N_i,N_s, N_mPro;
	projection* p;

	rbf_ps(Mesh& meshObject, struct probParams& probParamsObject);
	void perform_rbf(getNodeType& n);
	void performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec, getNodeType& n);
	void getExactDef(getNodeType& n, Eigen::VectorXd& exactDeformation);
};


#endif /* RBFPS_H_ */
