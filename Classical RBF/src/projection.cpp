

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
	Eigen::ArrayXd projection,finalProjection;
	int idxMin;
	for(int i = 0; i<sNodes.size(); i++){

//		std::cout << "Node in question: " << sNodes(i) << std::endl;



		finalProjection = Eigen::ArrayXd::Zero(m.nDims);
//		std::cout << "initial position: " << m.coords.row(sNodes(i)) << std::endl;
//		std::cout << "displaced position: " << m.coords.row(sNodes(i)) + delta.row(i) << std::endl;

		dist = m.midPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i));
		projectFun(m,delta,finalProjection, dist);
		std::cout << delta.row(i) + finalProjection.transpose() << std::endl;
		std::cout << "done" << std::endl;
		std::exit(0);


		dist.rowwise().norm().minCoeff(&idxMin);
//		std::cout << dist << std::endl;
//		std::cout << "nearest midPoint: " << m.midPnts.row(idxMin) << std::endl;

		double residual = 1;
		double tol = 1e-12;

		while(fabs(residual) > tol){
			projection = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix()) * m.midPntNormals.row(idxMin);

			finalProjection += projection;

			dist.rowwise() -= projection.transpose();

			dist.rowwise().norm().minCoeff(&idxMin);

			residual = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix());
//			std::cout <<"residual: " << residual << std::endl;

		}

		finalDef.row(i) = delta.row(i) + finalProjection.transpose();
		std::cout << finalDef.row(i) << std::endl;

		std::cout << "done" << std::endl;
		std::exit(0);
//		m.coords.row(sNodes(i)) += finalDef.row(i);
	}
}

void projection::projectFun(Mesh& m, Eigen::ArrayXXd& delta, Eigen::ArrayXd& projection, Eigen::ArrayXXd& dist){
	int idxMin;
	Eigen::ArrayXd project_i;

	projection = Eigen::ArrayXd::Zero(m.nDims);
	dist.rowwise().norm().minCoeff(&idxMin);

	double residual = 1;
	double tol = 1e-12;
	while(fabs(residual) > tol){
		project_i = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix()) * m.midPntNormals.row(idxMin);

		projection += project_i;

		dist.rowwise() -= project_i.transpose();
		dist.rowwise().norm().minCoeff(&idxMin);

		residual = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix());
	}
}


