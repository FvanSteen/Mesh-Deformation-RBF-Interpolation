#include "SPDS.h"

SPDS::SPDS() {}


void SPDS::kdt_NNSearch(Eigen::ArrayXi& bdryIndex, Eigen::ArrayXi& intIndex, Eigen::ArrayXXd& coords, const size_t dim, double& gamma, double&maxError, Eigen::ArrayXXd* ePtr){
	std::cout << "NN search  function" << std::endl;

	Eigen::ArrayXXd bdryPnts(bdryIndex.size(),dim);
	bdryPnts << coords(bdryIndex, Eigen::all);

	using kdt_bdry = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;
//	std::clock_t c_start = std::clock();
	kdt_bdry mat_index(dim, std::cref(bdryPnts), 10 /* max leaf */);
//	std::clock_t c_end = std::clock();
//	long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
//	std::cout << "CPU time used building kd tree: " << time_elapsed_ms << " ms\n";

	Eigen::ArrayXd query(dim);

	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);
	double dist;



	for(size_t i = 0; i < intIndex.size(); i++){

		query = coords.row(intIndex(i));

	    nanoflann::KNNResultSet<double> resultSet(N);

	    resultSet.init(&idx[0], &distSqrd[0]);

	    mat_index.index_->findNeighbors(resultSet, &query(0));

	    dist = sqrt(distSqrd[0]);

	    coords.row(intIndex[i]) -= (*ePtr).row(idx[0])*rbfEval(dist,gamma*maxError);
	}
}



double SPDS::rbfEval(double distance, double radius){
	double xi = distance/radius;	// distance scaled by support radius
	double f_xi;
	if(xi > 1){
		f_xi = 0;
	}else{
		f_xi = pow((1-(xi)),4)*(4*(xi)+1);
	}
	return f_xi;
}

//todo remove pVec as argument since m is also included
void SPDS::project(Mesh& m, getNodeType&n, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, int ptype){
//	std::cout << "SPDS project\n";
	if(n.N_se > 0){
		projectEdge(m, n.sePtr, array_in, array_out, 0, n.N_se, 1, ptype);
//		std::cout << "edge projection is done\n";
	}

//	if(n.N_ss > 0){
	if(array_in.rows() > n.N_se){
		projectSurf(m, n.ssPtr, array_in, array_out, n.N_se, n.N_se + n.N_ss, 1, ptype);
//		std::cout << "surf projection is done\n";
	}

}

void SPDS::projectSurf(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, size_t startIdx, size_t endIdx, int project, int ptype){
	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;

	// midpoints in cartesian coords
	kdt mat_index2(m.nDims, std::cref(m.surfMidPnts), 10 /* max leaf */);

	Eigen::ArrayXd query(m.nDims);
	Eigen::VectorXd distNN(m.nDims);
	Eigen::RowVectorXd project_i;


	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);

	Eigen::RowVectorXd projection;

	nanoflann::KNNResultSet<double> resultSet(N);

	double tol = 1e-12;
	//todo does not have to be an array probably
	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	projectionMagnitude.resize(1);

	for(size_t  i = startIdx; i < endIdx; i++){
		int idxNode;

		//todo check why this difference is being made for not projection purposes
		if(project){
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i-startIdx)));
		}else{
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i)));
		}
		if(idxNode != m.N_pe){
			if(project){
				projectSlidingEdge(m, array_out, array_in, (*nodesPtr)(i-startIdx), i, idxNode, project, ptype);
			}
			else{
				projectSlidingEdge(m, array_out, array_in, (*nodesPtr)(i), i, idxNode, project, ptype);
			}

		}
		else{

			projection = Eigen::RowVectorXd::Zero(m.nDims);

			if(project){
				query = m.coords.row((*nodesPtr)(i-startIdx)) + array_in.row(i);
			}else{
				query = m.coords.row((*nodesPtr)(i)) + array_in.row(i);
			}

			resultSet.init(&idx[0], &distSqrd[0]);

			mat_index2.index_->findNeighbors(resultSet, &query(0));

			getDistNearestNeighbour(m.surfMidPnts, distNN, idx, query, ptype);

			projectionMagnitude(0) =  distNN.dot( m.surfMidPntNormals.row(idx[0]).matrix());

			while(abs(projectionMagnitude).sum() > tol){

				project_i =  projectionMagnitude(0) * m.surfMidPntNormals.row(idx[0]);

				if(ptype){
					transformProjection(project_i, query);
				}

				projection += project_i;
				query += project_i.array();

				mat_index2.index_->findNeighbors(resultSet, &query(0));

				getDistNearestNeighbour(m.surfMidPnts, distNN, idx, query, ptype);


				projectionMagnitude(0) =  distNN.dot(m.surfMidPntNormals.row(idx[0]).matrix());

			}
			if(project){
				array_out.row(i) = array_in.row(i) + projection.array();
			}
			else{
				array_out.row(i) = -projection;
			}
		}

	}
}

void SPDS::projectEdge(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, size_t startIdx, size_t endIdx, int project, int ptype){

	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;
	kdt mat_index(m.nDims, std::cref(m.edgeMidPnts), 10 /* max leaf */);


	Eigen::ArrayXd query(m.nDims);
	Eigen::VectorXd distNN(m.nDims);
	Eigen::RowVectorXd project_i;


	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);

	Eigen::RowVectorXd projection;

	nanoflann::KNNResultSet<double> resultSet(N);

	double tol = 1e-12;

	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	std::vector<Eigen::ArrayXXd*> ptr(m.nDims-1);
	ptr[0] = &m.edgeMidPntNormals1;
	ptr[1] = &m.edgeMidPntNormals2;

	for(size_t i = startIdx; i < endIdx; i++){

		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),(*nodesPtr)(i)) != std::end(m.periodicVerticesNodes)){
//			std::cout << "SLINDING A PERIODIC VERTEX\n";
			projectSlidingVertex(m, array_out, array_in, (*nodesPtr)(i), i, project, ptype);

		}
		else{
			projection = Eigen::RowVectorXd::Zero(m.nDims);
//			std::cout << "node: " << (*nodesPtr)(i) << std::endl;
			query = m.coords.row((*nodesPtr)(i)) + array_in.row(i);
//			std::cout << query << std::endl;
			resultSet.init(&idx[0], &distSqrd[0]);

			mat_index.index_->findNeighbors(resultSet, &query(0));

			getDistNearestNeighbour(m.edgeMidPnts, distNN, idx, query, ptype);
//			std::cout << distNN << std::endl;
//			std::cout << "nearest midpoint: \n" << m.edgeMidPnts.row(idx[0]) << std::endl;

			for(size_t j = 0; j < size_t(m.nDims-1); j++){
				projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
			}


			while(abs(projectionMagnitude).sum() > tol){

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					project_i =  projectionMagnitude(j) *(*ptr[j]).row(idx[0]);
					if(ptype){
						transformProjection(project_i, query);
					}
//					std::cout <<"\n" <<  project_i << std::endl;
					projection += project_i;
					query += project_i.array();
//					std::cout << "\nquery:\n" << query << std::endl;
				}


				mat_index.index_->findNeighbors(resultSet, &query(0));

				getDistNearestNeighbour(m.edgeMidPnts, distNN, idx, query, ptype);

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
				}
			}
			if(project)
				array_out.row(i) = array_in.row(i) + projection.array();
			else
				array_out.row(i) = -projection;
		}
//		std::cout << i << '\t' << array_out.row(i) << std::endl;
	}
}

void SPDS::projectSlidingEdge(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int idxPerEdge, int project, int ptype){
	Eigen::VectorXd delta(m.nDims);
//	std::cout << "node: " << node << std::endl;

	if(ptype){
		double dr  = sqrt( pow(m.coords(node,0)+array_in(idx,0),2) + pow(m.coords(node,1)+array_in(idx,1),2)) - sqrt( pow(m.coords(node,0),2) + pow(m.coords(node,1),2));
		double dtheta = atan2(m.coords(node,1)+array_in(idx,1),m.coords(node,0)+array_in(idx,0)) - atan2(m.coords(node,1),m.coords(node,0));

		if(m.nDims == 2)
			delta << dr, dtheta;
		else
			delta << dr, dtheta, array_in(idx,2);
	}else{
		delta = array_in.row(idx);
	}


	Eigen::ArrayXd projection;
	if(project){
		projection = delta.array().transpose() - delta.dot(m.periodicEdgeNormals.row(idxPerEdge).matrix())*m.periodicEdgeNormals.row(idxPerEdge);
	}else{
		projection = delta.dot(m.periodicEdgeNormals.row(idxPerEdge).matrix())*m.periodicEdgeNormals.row(idxPerEdge);
	}

//	if(project == 1){
//		std::exit(0);
//	}

//	std::cout << delta.dot(m.periodicEdgeNormals.row(idxPerEdge).matrix())*m.periodicEdgeNormals.row(idxPerEdge) << std::endl;
//	std::cout << projection << "\n\n" << std::endl;
//
//
//	std::cout << "\nProjection:\n";
//	std::cout << projection << std::endl;
//	std::cout << "\nedge normal:\n";
//	std::cout << m.periodicEdgeNormals.row(idxPerEdge) << std::endl;
//
//	std::cout << "\n corresponding node: " << m.periodicEdgeNodes(idxPerEdge) << std::endl;
//	std::cout << "oldposition:\n" << m.coords_polar_cylindrical.row(node)<< std::endl;
//	std::cout << "newposition:\n" << m.coords_polar_cylindrical.row(node) + delta.array().transpose() - projection.transpose() << std::endl;
	array_out(idx,0) = (m.coords_polar_cylindrical(node,0) + projection(0))*cos(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*cos(m.coords_polar_cylindrical(node,1));
	array_out(idx,1) = (m.coords_polar_cylindrical(node,0) + projection(0))*sin(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*sin(m.coords_polar_cylindrical(node,1));
	if(m.nDims == 3)
		array_out(idx,2) = projection(2);

//	std::cout << "array out:\n" << array_out.row(idx) << std::endl;


//	if(project == 0){
//		std::exit(0);
//	}

}

void SPDS::projectSlidingVertex(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int project, int ptype){
	Eigen::ArrayXd delta(m.nDims);

	if(ptype){
		double dr  = sqrt( pow(m.coords(node,0)+array_in(idx,0),2) + pow(m.coords(node,1)+array_in(idx,1),2)) - sqrt( pow(m.coords(node,0),2) + pow(m.coords(node,1),2));
		double dtheta = atan2(m.coords(node,1)+array_in(idx,1),m.coords(node,0)+array_in(idx,0)) - atan2(m.coords(node,1),m.coords(node,0));

		if(m.nDims == 2)
			delta << dr, dtheta;
		else
			delta << dr, dtheta, array_in(idx,2);
	}else{
		delta = array_in.row(idx);
	}

	Eigen::ArrayXd projection;
	if(project)
		projection = delta*m.periodicVecs.col(0).array();
	else
		projection = delta - delta*m.periodicVecs.col(0).array();



	array_out(idx,0) = (m.coords_polar_cylindrical(node,0) + projection(0))*cos(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*cos(m.coords_polar_cylindrical(node,1));
	array_out(idx,1) = (m.coords_polar_cylindrical(node,0) + projection(0))*sin(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*sin(m.coords_polar_cylindrical(node,1));
	if(m.nDims == 3)
		array_out(idx,2) = projection(2);
}

void SPDS::transformProjection(Eigen::RowVectorXd& project_i, Eigen::ArrayXd& query){
	Eigen::RowVectorXd transformedProjection(project_i.size());

	double r = sqrt( pow(query(0),2) + pow(query(1),2));
	double theta = atan2(query(1),query(0));

	if(project_i.size() == 3){
		transformedProjection << (project_i(0)+r)*cos(project_i(1)+theta) - query(0), (project_i(0)+r)*sin(project_i(1)+theta) - query(1), project_i(2);
	}else{
		transformedProjection << (project_i(0)+r)*cos(project_i(1)+theta) - query(0), (project_i(0)+r)*sin(project_i(1)+theta) - query(1);
	}

	project_i = transformedProjection;

}

void SPDS::getDistNearestNeighbour(Eigen::ArrayXXd& midPnts,  Eigen::VectorXd& d, std::vector<size_t> idx, Eigen::ArrayXd& query, int ptype){
	//todo input argument should be the midpoints either edge or surf
	// todo m can be removed and replaced by query.size()
	if(ptype){
		double dr = sqrt(pow(midPnts(idx[0],0),2) + pow(midPnts(idx[0],1),2)) - sqrt(pow(query(0),2) + pow(query(1),2));
		double dtheta = atan2(midPnts(idx[0],1), midPnts(idx[0],0)) - atan2(query(1), query(0));

		if(query.size() == 3){
			double dz = midPnts(idx[0],2) - query(2);
			d << dr, dtheta, dz;
		}else{
			d << dr, dtheta;
		}

	}else{
		d = midPnts.row(idx[0]) - query.transpose();
	}
}
