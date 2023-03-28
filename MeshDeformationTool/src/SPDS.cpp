#include "SPDS.h"

SPDS::SPDS() {
	std::cout << "in the space partitioning data structure class" << std::endl;

}


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


//	std::clock_t c_start2 = std::clock();
	for(size_t i = 0; i < intIndex.size(); i++){

		query = coords.row(intIndex(i));

	    nanoflann::KNNResultSet<double> resultSet(N);

	    resultSet.init(&idx[0], &distSqrd[0]);

	    mat_index.index_->findNeighbors(resultSet, &query(0));

	    dist = sqrt(distSqrd[0]);

	    coords.row(intIndex[i]) -= (*ePtr).row(idx[0])*rbfEval(dist,gamma*maxError);

	}
//	std::clock_t c_end2 = std::clock();
//	long double time_elapsed_ms2 = 1000.0 * (c_end2-c_start2) / CLOCKS_PER_SEC;
//	std::cout << "CPU time used single qu: " << time_elapsed_ms2 << " ms\n";
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

