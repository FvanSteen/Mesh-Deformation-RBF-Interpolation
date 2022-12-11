

#include "projection.h"

#include <Eigen/Dense>
#include <iostream>

projection::projection() {

}

void projection::project(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec){
	std::cout << "Doing Projection" << std::endl;

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
	std::cout << "projection is done " << std::endl;
}


