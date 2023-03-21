

#include "projection.h"

#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include "nanoflann.hpp"


projection::projection(Eigen::VectorXd& pVec)
:pVec(pVec)
{}

void projection::project(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec){
	std::cout << "Doing Projection CHECK WERE THIS IS CALLED" << std::endl;
	std::exit(0);

	Eigen::RowVectorXd d;
	Eigen::ArrayXd dist, projection;
	Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(m.edgeMidPnts.rows(),0,m.edgeMidPnts.rows()-1);

	for(int i=0; i<sNodes.size(); i++){
		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),sNodes(i))  != std::end(m.periodicVerticesNodes)){
			finalDef.row(i) = delta.row(i)*pVec.transpose().array();
		}else{
			// distance to all midpoints
			dist = (m.edgeMidPnts.rowwise()-(m.coords.row(sNodes(i)) + delta.row(i))).rowwise().norm();

			// finding closest midpoint by sorting them in ascending order
			std::sort(index.begin(), index.end(),[&](const int& a, const int& b) {
				return (dist[a] < dist[b]);
				}
			);

			// vector with distance from closest midpoint to displaced sliding node.
			d = m.edgeMidPnts.row(index(0)) - (m.coords.row(sNodes(i)) + delta.row(i));

			// Projection required to bring node back on the external boundary
			projection = d.dot(m.edgeMidPntNormals1.row(index(0)).matrix())*m.edgeMidPntNormals1.row(index(0));

			// Final deformation is the initial displacement plus the projection to bring the node back on the boundary
			finalDef.row(i) = delta.row(i) + projection.transpose();
		}

	}
//	std::cout << "projection is done " << std::endl;
}

void projection::projectIter(Mesh& m, Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, int N_se){
//	std::cout << "in the projection iteration function" << std::endl;


	Eigen::ArrayXXd dist;
	Eigen::RowVectorXd finalProjection;

	int edge;

	for(int i = 0; i < sNodes.size(); i++){
		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),sNodes(i))  != std::end(m.periodicVerticesNodes)){
			finalDef.row(i) = delta.row(i)*pVec.transpose().array();
		}
//		else if(std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),sNodes(i))  != std::end(m.periodicEdgeNodes)){
//
//		}
		else{
			if(i < N_se){
				edge = 1;
				dist = m.edgeMidPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i));
			}else{
				edge = 0;
				dist = m.edgeMidPnts.rowwise() - (m.coords.row(sNodes(i)) + delta.row(i));
			}
			projectFun(m,finalProjection, dist, edge);

			finalDef.row(i) = delta.row(i) + finalProjection.array();
		}
	}
}


void projection::projectFun(Mesh& m,  Eigen::RowVectorXd& projection, Eigen::ArrayXXd& dist, int edge){

	int idxMin;
	Eigen::RowVectorXd project_i;

	projection = Eigen::RowVectorXd::Zero(m.nDims);
	dist.rowwise().norm().minCoeff(&idxMin);

	double residual = 1;
	double tol = 1e-12;


	if(edge){
		while(fabs(residual) > tol){
			project_i = dist.row(idxMin).matrix().dot(m.edgeMidPntNormals1.row(idxMin).matrix()) * m.edgeMidPntNormals1.row(idxMin) + dist.row(idxMin).matrix().dot(m.edgeMidPntNormals2.row(idxMin).matrix()) * m.edgeMidPntNormals2.row(idxMin);

			projection += project_i;

			dist.rowwise() -= project_i.array();
			dist.rowwise().norm().minCoeff(&idxMin);

			residual = dist.row(idxMin).matrix().dot(m.edgeMidPntNormals1.row(idxMin).matrix()) + dist.row(idxMin).matrix().dot(m.edgeMidPntNormals2.row(idxMin).matrix());
		}
	}else{
		while(fabs(residual) > tol){
			project_i = dist.row(idxMin).matrix().dot(m.surfMidPntNormals.row(idxMin).matrix()) * m.surfMidPntNormals.row(idxMin);

			projection += project_i;

			dist.rowwise() -= project_i.array();
			dist.rowwise().norm().minCoeff(&idxMin);

			residual = dist.row(idxMin).matrix().dot(m.surfMidPntNormals.row(idxMin).matrix());
		}
	}



}

void projection::kd_tree(Mesh& m, getNodeType& n, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec){
	std::cout << "kd tree projection" << std::endl;

	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;

	kdt mat_index(m.nDims, std::cref(m.edgeMidPnts), 10 /* max leaf */);


	Eigen::ArrayXd query(m.nDims);
	Eigen::VectorXd distNN;
	Eigen::RowVectorXd project_i;


	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);
//	double dist

	Eigen::RowVectorXd projection;

	nanoflann::KNNResultSet<double> resultSet(N);

	double tol = 1e-6;

	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	std::vector<Eigen::ArrayXXd*> ptr(m.nDims-1);
	ptr[0] = &m.edgeMidPntNormals1;
	ptr[1] = &m.edgeMidPntNormals2;

	for(size_t i = 0; i < size_t(n.N_se); i++){
		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),(*n.sePtr)(i)) != std::end(m.periodicVerticesNodes)){
//			std::cout << delta.row(i) << std::endl;
//			std::cout << pVec << std::endl;
//
			finalDef.row(i) = delta.row(i)*pVec.transpose().array();
//			std::cout << finalDef.row(i) << std::endl;
//			std::exit(0);
		}
		else{

			projection = Eigen::RowVectorXd::Zero(m.nDims);
	//		std::cout << "node in question: " << sNodes(i) << std::endl;

			query = m.coords.row((*n.sePtr)(i)) + delta.row(i);
	//		std::cout << "init query: \n" << query << std::endl;

			resultSet.init(&idx[0], &distSqrd[0]);

			mat_index.index_->findNeighbors(resultSet, &query(0));

			distNN = m.edgeMidPnts.row(idx[0]) - query.transpose();
	//    	projectionMagnitude =  distNN.dot(m.edgeMidPntNormals1.row(idx[0]).matrix());

			for(size_t j = 0; j < size_t(m.nDims-1); j++){
				projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
			}
			while(abs(projectionMagnitude).sum() > tol){

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					project_i =  projectionMagnitude(j) *(*ptr[j]).row(idx[0]);
					projection += project_i;
					query += project_i.array();
				}


	//	    	std::cout << "query\n" << query << std::endl;
				mat_index.index_->findNeighbors(resultSet, &query(0));

				distNN = m.edgeMidPnts.row(idx[0]) - query.transpose();

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
				}
	//	    	projectionMagnitude =  distNN.dot(m.edgeMidPntNormals1.row(idx[0]).matrix());
	//	    	std::cout << "Projection Magnitude: " << projectionMagnitude << std::endl;
			}
	//	    std::cout << projection << std::endl;
	//	    std::cout << "DONE" << std::endl;


			finalDef.row(i) = delta.row(i) + projection.array();
	//	    std::cout << delta.row(i) << std::endl;
	//	    std::cout << projection << std::endl;
	//	    std::cout << finalDef.row(i) << std::endl;
	//	    std::exit(0);
		}
	}


	if(m.N_ss >0){
		kdt mat_index2(m.nDims, std::cref(m.surfMidPnts), 10 /* max leaf */);

		projectionMagnitude.resize(1);

		for(size_t  i = n.N_se; i < size_t(n.N_s); i++){

			int idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*n.ssPtr)(i-n.N_se)));
			if(idxNode != m.N_pe){
	//			std::cout << "node: " << (*n.ssPtr)(i-n.N_se) << std::endl;
	//			std::cout << delta.row(i) << std::endl;
	//			std::cout << m.n_ss.row(m.N_ss-m.N_pe+idxNode) << std::endl;

				finalDef.row(i)  = delta.row(i) + delta.row(i)*m.n_ss.row(m.N_ss-m.N_pe+idxNode);
	//			std::cout << finalDef.row(i) << std::endl;

			}
			else{
				projection = Eigen::RowVectorXd::Zero(m.nDims);

				query = m.coords.row((*n.ssPtr)(i-n.N_se)) + delta.row(i);


				resultSet.init(&idx[0], &distSqrd[0]);

				mat_index2.index_->findNeighbors(resultSet, &query(0));

				distNN = m.surfMidPnts.row(idx[0]) - query.transpose();

				projectionMagnitude(0) =  distNN.dot( m.surfMidPntNormals.row(idx[0]).matrix());

				while(abs(projectionMagnitude).sum() > tol){

					project_i =  projectionMagnitude(0) * m.surfMidPntNormals.row(idx[0]);
					projection += project_i;
					query += project_i.array();

					mat_index2.index_->findNeighbors(resultSet, &query(0));

					distNN = m.surfMidPnts.row(idx[0]) - query.transpose();


					projectionMagnitude(0) =  distNN.dot(m.surfMidPntNormals.row(idx[0]).matrix());

				}
				finalDef.row(i) = delta.row(i) + projection.array();
			}

		}
	}

	std::cout << "DONE KD tree pojection\n";
}


