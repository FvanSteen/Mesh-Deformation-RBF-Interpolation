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
void SPDS::project(Mesh& m, getNodeType&n, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec){
//	std::cout << "SPDS project\n";
	if(n.N_se > 0){
		projectEdge(m, n.sePtr, array_in, array_out, pVec, 0, n.N_se, 1);
//		std::cout << "edge projection is done\n";
	}

//	if(n.N_ss > 0){
	if(array_in.rows() > n.N_se){
		projectSurf(m, n.ssPtr, array_in, array_out, pVec, n.N_se, n.N_se + n.N_ss, 1);
//		std::cout << "surf projection is done\n";
	}

}

void SPDS::projectSurf(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec, size_t startIdx, size_t endIdx, int project){
	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;
	kdt mat_index2(m.nDims, std::cref(m.surfMidPnts), 10 /* max leaf */);


	Eigen::ArrayXd query(m.nDims);
	Eigen::VectorXd distNN;
	Eigen::RowVectorXd project_i;


	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);

	Eigen::RowVectorXd projection;

	nanoflann::KNNResultSet<double> resultSet(N);

	double tol = 1e-6;

	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	projectionMagnitude.resize(1);

	for(size_t  i = startIdx; i < endIdx; i++){
		int idxNode;
		if(project){
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i-startIdx)));
		}else{
			idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),(*nodesPtr)(i)));
		}
		if(idxNode != m.N_pe){
			if(project){
				array_out.row(i)  = array_in.row(i) + array_in.row(i)*m.n_ss.row(m.N_ss-m.N_pe+idxNode);
			}else{
				array_out.row(i)  = -array_in.row(i)*m.n_ss.row(m.N_ss-m.N_pe+idxNode);
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
			if(project){
				array_out.row(i) = array_in.row(i) + projection.array();
			}
			else{
				array_out.row(i) = -projection;
			}
		}
	}
}

void SPDS::projectEdge(Mesh& m, Eigen::ArrayXi* nodesPtr, Eigen::ArrayXXd& array_in, Eigen::ArrayXXd& array_out, Eigen::VectorXd& pVec, size_t startIdx, size_t endIdx, int project){

	using kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;

	kdt mat_index(m.nDims, std::cref(m.edgeMidPnts), 10 /* max leaf */);


	Eigen::ArrayXd query(m.nDims);
	Eigen::VectorXd distNN;
	Eigen::RowVectorXd project_i;


	const size_t N = 1;
	std::vector<size_t> idx(N);
	std::vector<double> distSqrd(N);

	Eigen::RowVectorXd projection;

	nanoflann::KNNResultSet<double> resultSet(N);

	double tol = 1e-6;

	Eigen::ArrayXd projectionMagnitude(m.nDims-1);

	std::vector<Eigen::ArrayXXd*> ptr(m.nDims-1);
	ptr[0] = &m.edgeMidPntNormals1;
	ptr[1] = &m.edgeMidPntNormals2;

	for(size_t i = startIdx; i < endIdx; i++){
		if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),(*nodesPtr)(i)) != std::end(m.periodicVerticesNodes)){
			if(project)
				array_out.row(i) = array_in.row(i)*pVec.transpose().array();
			else
				array_out.row(i) =   array_in.row(i) - array_in.row(i)*pVec.transpose().array();

		}
		else{
			projection = Eigen::RowVectorXd::Zero(m.nDims);


			query = m.coords.row((*nodesPtr)(i)) + array_in.row(i);


			resultSet.init(&idx[0], &distSqrd[0]);

			mat_index.index_->findNeighbors(resultSet, &query(0));

			distNN = m.edgeMidPnts.row(idx[0]) - query.transpose();


			for(size_t j = 0; j < size_t(m.nDims-1); j++){
				projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
			}
			while(abs(projectionMagnitude).sum() > tol){

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					project_i =  projectionMagnitude(j) *(*ptr[j]).row(idx[0]);
					projection += project_i;
					query += project_i.array();
				}

				mat_index.index_->findNeighbors(resultSet, &query(0));

				distNN = m.edgeMidPnts.row(idx[0]) - query.transpose();

				for(size_t j = 0; j < size_t(m.nDims-1); j++){
					projectionMagnitude(j) =  distNN.dot((*ptr[j]).row(idx[0]).matrix());
				}
			}
			if(project)
				array_out.row(i) = array_in.row(i) + projection.array();
			else
				array_out.row(i) = -projection;
		}
	}
}
