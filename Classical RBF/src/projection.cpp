

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

void projection::projectIter(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef){
	std::cout << "in the projection iteration function" << std::endl;

	Eigen::ArrayXXd dist;
	Eigen::ArrayXd projection;
	int idxMin,idxMinNew;
	for(int i = 0; i<sNodes.size(); i++){
//		idxMin = -1;
		idxMinNew = -1;
//
//		projection << 0,0;

		dist = m.midPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i));
		dist.rowwise().norm().minCoeff(&idxMin);

		projection = dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix()) * m.midPntNormals.row(idxMin);

//		int iter = 0;
		while(idxMin != idxMinNew){
			dist = m.midPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i) + projection.transpose());
			dist.rowwise().norm().minCoeff(&idxMinNew);

			projection += dist.row(idxMinNew).matrix().dot(m.midPntNormals.row(idxMinNew).matrix()) * m.midPntNormals.row(idxMinNew);
			idxMin = idxMinNew;
//			iter++;
//			std::cout << iter << std::endl;
		}
		finalDef.row(i) = delta.row(i) + projection.transpose();
//		std::cout << "Node: " << sNodes(i) << std::endl;
//		std::cout << (m.coords.row(sNodes(i)) + delta.row(i)) << std::endl;

//		std::cout << dist << std::endl;







//		std::cout << dist.row(idxMin) << std::endl;
//		std::cout << m.midPntNormals.row(idxMin) << std::endl;
//		std::cout << "midpoint: "  <<  m.midPnts.row(idxMin) << std::endl;
//		std::cout << std::endl;
//		std::cout << projection << std::endl;

//		finalDef.row(i) = delta.row(i);
//		std::cout << idxMin << std::endl;

//		std::cout << dist(idxMin) << std::endl;

//		dist = m.midPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i) + projection.transpose());
//		int idxMin2;
//		dist.rowwise().norm().minCoeff(&idxMin2);
//		if (idxMin2 != idxMin){
////			std::cout << "these nodes are not the same" << std::endl;
//			std::cout << sNodes(i) << ", ";
//		}


//		m.coords.row(sNodes(i)) += finalDef.row(i);

	}
//	std::cout << std::endl;

//	m.writeMeshFile();
//	std::exit(0);
}


