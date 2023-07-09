#include "SPDS.h"

SPDS::SPDS() {}


/* kdt_NNSearch
 *
 * This function uses a KD-tree to find the nearest neighbour boundary node for all internal nodes
 * This nearest neighbour is used to perform the local boundary correction
 */

void SPDS::kdt_NNSearch(Eigen::ArrayXi& bdryIndex, Eigen::ArrayXi& intIndex, Eigen::ArrayXXd& coords, const size_t dim, double& gamma, double&maxError, Eigen::ArrayXXd* ePtr){

	// Set of the boundary points
	Eigen::ArrayXXd bdryPnts(bdryIndex.size(),dim);
	bdryPnts << coords(bdryIndex, Eigen::all);

	// setting up the KD tree, see https://github.com/jlblancoc/nanoflann for more information about the nanoflann header
	using kdt_bdry = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;
	kdt_bdry mat_index(dim, std::cref(bdryPnts), 10 /* max leaf */);

	// the query points
	Eigen::ArrayXd query(dim);

	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);
	double dist;

	// looping through the internal nodes
	for(int i = 0; i < intIndex.size(); i++){
		//setting query
		query = coords.row(intIndex(i));

		// finding nearest neighbour
	    nanoflann::KNNResultSet<double> resultSet(N);
	    resultSet.init(&idx[0], &distSqrd[0]);
	    mat_index.index_->findNeighbors(resultSet, &query(0));

	    // determine distance
	    dist = sqrt(distSqrd[0]);

	    // apply the correction. See eq. 2.30 of the manuscript
	    coords.row(intIndex[i]) -= (*ePtr).row(idx[0])*rbfEval(dist,gamma*maxError);
	}
}

/* rbfEval
 *
 * Function that evaluates the RBF
 *
 */

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


/* project
 *
 * performs the projection of the sliding nodes
 *
 */

void SPDS::project(Mesh& m, getNodeType&n, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, int ptype){


	// projection of the edge nodes
	if(n.N_se > 0){
		projectEdge(m, n.sePtr, array_in, array_out, 0, n.N_se, 1, ptype);
	}

	// projection of the surface nodes
	if(array_in.rows() > n.N_se){
		projectSurf(m, n.ssPtr, array_in, array_out, n.N_se, n.N_se + n.N_ss, 1, ptype);
	}

}

/* projectSurf
 *
 * Performs the projection of the surface nodes
 */

void SPDS::projectSurf(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, size_t startIdx, size_t endIdx, int project, int ptype){
	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;

	// KD tree with the surface midpoints
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
	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	projectionMagnitude.resize(1);

	// looping through the nodes
	for(size_t  i = startIdx; i < endIdx; i++){
		int idxNode;

		// check if the considered node is a periodic edge node.
		// in case of projection (project) the moving nodes at the start of the array_in array should be skipped.
		if(project){
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i-startIdx)));
		}else{
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i)));
		}

		// in case the node is among the periodic edge nodes then its projected as a sliding edge node
		if(idxNode != m.N_pe){
			if(project){
				projectSlidingEdge(m, array_out, array_in, (*nodesPtr)(i-startIdx), i, idxNode, project, ptype);
			}
			else{
				projectSlidingEdge(m, array_out, array_in, (*nodesPtr)(i), i, idxNode, project, ptype);
			}

		}

		// projection of the surface node
		else{

			projection = Eigen::RowVectorXd::Zero(m.nDims);

			// establish the query point for the nearest neighbour search
			if(project){
				query = m.coords.row((*nodesPtr)(i-startIdx)) + array_in.row(i);
			}else{
				query = m.coords.row((*nodesPtr)(i)) + array_in.row(i);
			}

			resultSet.init(&idx[0], &distSqrd[0]);

			mat_index2.index_->findNeighbors(resultSet, &query(0));

			// determine the nearest neighbour
			getDistNearestNeighbour(m.surfMidPntPtr, distNN, idx, query, ptype);

			// projection magnitude
			projectionMagnitude(0) =  distNN.dot( m.surfMidPntNormals.row(idx[0]).matrix());

			// iteratively perform the projection until the tolerance is reached
			while(abs(projectionMagnitude).sum() > tol){

				// the projection
				project_i =  projectionMagnitude(0) * m.surfMidPntNormals.row(idx[0]);

				// transformation to cartesian coordinates if required
				if(ptype){
					transformProjection(project_i, query);
				}

				// add to total projection
				projection += project_i;

				// update query points
				query += project_i.array();

				// find new nearest neighbour
				mat_index2.index_->findNeighbors(resultSet, &query(0));

				// determine distance nearest neighbour
				getDistNearestNeighbour(m.surfMidPntPtr, distNN, idx, query, ptype);


				projectionMagnitude(0) =  distNN.dot(m.surfMidPntNormals.row(idx[0]).matrix());

			}

			// determine the projected displacement
			if(project){
				array_out.row(i) = array_in.row(i) + projection.array();
			}
			// or the error
			else{
				array_out.row(i) = -projection;
			}
		}

	}
}

/* projectEdge
 *
 * performs the projection of the edge nodes
 */

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

		// if the considered node is a periodic vertex then its projected as such
		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),(*nodesPtr)(i)) != std::end(m.periodicVerticesNodes)){

			projectSlidingVertex(m, array_out, array_in, (*nodesPtr)(i), i, project, ptype);

		}
		else{
			projection = Eigen::RowVectorXd::Zero(m.nDims);

			// query points
			query = m.coords.row((*nodesPtr)(i)) + array_in.row(i);

			resultSet.init(&idx[0], &distSqrd[0]);

			// find new nearest neighbour index
			mat_index.index_->findNeighbors(resultSet, &query(0));

			// determine distance nearest neighbour
			getDistNearestNeighbour(m.edgeMidPntPtr, distNN, idx, query, ptype);


			// get magnitude of the projection
			for(size_t j = 0; j < size_t(m.nDims-1); j++){
				projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
			}


			while(abs(projectionMagnitude).sum() > tol){

				for(size_t j = 0; j < size_t(m.nDims-1); j++){

					// find the projection of the current iteration
					project_i =  projectionMagnitude(j) *(*ptr[j]).row(idx[0]);

					if(ptype){
						transformProjection(project_i, query);
					}
					// update the total projection
					projection += project_i;
					// update query point
					query += project_i.array();

				}

				// find new nearest neighbour index
				mat_index.index_->findNeighbors(resultSet, &query(0));

				// determine distance nearest neighbour
				getDistNearestNeighbour(m.edgeMidPntPtr, distNN, idx, query, ptype);

				// update the projection magnitude
				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
				}
			}

			// determine the projected displacement
			if(project){
				array_out.row(i) = array_in.row(i) + projection.array();
			}
			// or the error
			else{
				array_out.row(i) = -projection;
			}
		}
	}
}

/* projectSlidingEdge
 *
 * Performs the projection of the sliding edges
 */

void SPDS::projectSlidingEdge(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int idxPerEdge, int project, int ptype){
	Eigen::VectorXd delta(m.nDims);


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


	if(ptype){
		array_out(idx,0) = (m.coords_polar_cylindrical(node,0) + projection(0))*cos(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*cos(m.coords_polar_cylindrical(node,1));
		array_out(idx,1) = (m.coords_polar_cylindrical(node,0) + projection(0))*sin(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*sin(m.coords_polar_cylindrical(node,1));
		if(m.nDims == 3)
			array_out(idx,2) = projection(2);
	}else{
		array_out.row(idx) = projection;
	}
}

/* projectSlidingVertex
 *
 * Performs the projection of the vertices
 */

void SPDS::projectSlidingVertex(Mesh& m, Eigen::ArrayXXd& array_out,Eigen::ArrayXXd& array_in, int node, int idx, int project, int ptype){
	Eigen::ArrayXd delta(m.nDims);
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
	if(project)
		projection = delta*m.periodicVecs.col(0).array();
	else
		projection = delta - delta*m.periodicVecs.col(0).array();

	if(ptype){
		array_out(idx,0) = (m.coords_polar_cylindrical(node,0) + projection(0))*cos(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*cos(m.coords_polar_cylindrical(node,1));
		array_out(idx,1) = (m.coords_polar_cylindrical(node,0) + projection(0))*sin(m.coords_polar_cylindrical(node,1) + projection(1)) - m.coords_polar_cylindrical(node,0)*sin(m.coords_polar_cylindrical(node,1));
		if(m.nDims == 3)
			array_out(idx,2) = projection(2);
	}else{
		array_out.row(idx) = projection;
	}

}


/* transformProjection
 *
 * Function that transforms the projection from polar/ cylindrical coordinates to Cartesian coordinates
 *
 */

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


/* getDistNearestNeighbour
 *
 *  function that determines the distance to the nearest neighbour
 */

void SPDS::getDistNearestNeighbour(Eigen::ArrayXXd* midPnts,  Eigen::VectorXd& d, std::vector<size_t> idx, Eigen::ArrayXd& query, int ptype){
	if(ptype){

		d = (*midPnts).row(idx[0]);
		d(0) -= sqrt(pow(query(0),2) + pow(query(1),2));
		d(1) -= atan2(query(1),query(0));
		if(query.size() == 3){
			d(2) -= query(2);
		}

	}else{
		d = (*midPnts).row(idx[0]) - query.transpose();
	}

}
