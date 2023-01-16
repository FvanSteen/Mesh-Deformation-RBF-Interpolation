

#include "projection.h"

#include <Eigen/Dense>
#include <iostream>
#include <chrono>

projection::projection() {

}

void projection::project(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec){
//	std::cout << "Doing Projection" << std::endl;

	Eigen::RowVectorXd d;
	Eigen::ArrayXd dist, projection;
	Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(m.midPnts.rows(),0,m.midPnts.rows()-1);

	for(int i=0; i<sNodes.size(); i++){
		if(std::find(std::begin(m.staticNodes),std::end(m.staticNodes),sNodes(i))  != std::end(m.staticNodes)){
			finalDef.row(i) = delta.row(i)*pVec.transpose().array();
		}else{
			// distance to all midpoints
			dist = (m.midPnts.rowwise()-(m.coords.row(sNodes(i)) + delta.row(i))).rowwise().norm();

			// finding closest midpoint by sorting them in ascending order
			std::sort(index.begin(), index.end(),[&](const int& a, const int& b) {
				return (dist[a] < dist[b]);
				}
			);

			// vector with distance from closest midpoint to displaced sliding node.
			d = m.midPnts.row(index(0)) - (m.coords.row(sNodes(i)) + delta.row(i));

			// Projection required to bring node back on the external boundary
			projection = d.dot(m.midPntNormals.row(index(0)).matrix())*m.midPntNormals.row(index(0));

			// Final deformation is the initial displacement plus the projection to bring the node back on the boundary
			finalDef.row(i) = delta.row(i) + projection.transpose();
		}

	}
//	std::cout << "projection is done " << std::endl;
}

void projection::projectIter(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef){
//	std::cout << "in the projection iteration function" << std::endl;

	Eigen::ArrayXXd dist;
	Eigen::RowVectorXd finalProjection;

	for(int i = 0; i<sNodes.size(); i++){

		dist = m.midPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i));
		projectFun(m,delta,finalProjection, dist);

		finalDef.row(i) = delta.row(i) + finalProjection.array();
	}
}

void projection::projectFun(Mesh& m, Eigen::ArrayXXd& delta, Eigen::RowVectorXd& projection, Eigen::ArrayXXd& dist){

	int idxMin;
	Eigen::RowVectorXd project_i;

	projection = Eigen::RowVectorXd::Zero(m.nDims);
	dist.rowwise().norm().minCoeff(&idxMin);

	double residual = 1;
	double tol = 1e-12;
	while(fabs(residual) > tol){
		project_i = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix()) * m.midPntNormals.row(idxMin);

		projection += project_i;

		dist.rowwise() -= project_i.array();
		dist.rowwise().norm().minCoeff(&idxMin);

		residual = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix());
	}

}


