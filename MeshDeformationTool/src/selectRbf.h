/*
 * selectRbf.h
 *
 *  Created on: 19 nov. 2022
 *      Author: floyd
 */

#ifndef SELECTRBF_H_
#define SELECTRBF_H_
#include <iostream>
template<typename base>
class selectRbf : public base
{
public:
	selectRbf(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir);
//	rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
};

template <typename base>
selectRbf<base>::selectRbf(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
:rbf_std(meshObject, dVec, rotPnt, rotVec, steps, smode, curved, pDir)
{
	std::cout << "testing the rbf selector " << std::endl;
}




#endif /* SELECTRBF_H_ */
