/*;
 * rbfps.h
 *
 *  Created on: 16 nov. 2022
 *      Author: floyd
 */

#ifndef RBFPS_H_
#define RBFPS_H_
#include "Mesh.h"
#include "projection.h"
#include "getNodeType.h"
#include "greedy.h"
#include "rbfGenFunc.h"
#include <string>
#include <Eigen/Dense>




//class rbf_ps : public rbf_std, public projection
class rbf_ps : public rbfGenFunc
{
public:

	rbf_ps(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n);
	void perform_rbf(getNodeType& n);
	void performRBF_PS(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_bb, Eigen::MatrixXd& Phi_ib, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec_b, getNodeType& n, projection& p);

};


#endif /* RBFPS_H_ */
