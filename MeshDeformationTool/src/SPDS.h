
#ifndef SPDS_H_
#define SPDS_H_
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <vector>
#include "nanoflann.hpp"
//#include "utils.h"


class SPDS {
public:

	SPDS();
	template <typename num_t>
	void kdtree(const size_t N, const size_t dim);

	template <typename Der>
	void generateRandomPointCloud(Eigen::ArrayBase<Der>& mat, const size_t N, const size_t dim, const typename Der::Scalar max_range);


	void kdt_NNSearch(Eigen::ArrayXi& bdryIndex, Eigen::ArrayXi& intIndex, Eigen::ArrayXXd& coords, const size_t dim, double& gamma, double&maxError, Eigen::ArrayXXd* ePtr);
	double rbfEval(double distance, double radius);
//	void generateKDtree(Eigen::ArrayXXd& bdryIndex, const size_t dim);





//	using kdt_bdry = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::ArrayXXd>;
//	kdt_bdry bdry_index;
};





template<typename num_t>
void SPDS::kdtree(const size_t N, const size_t dim){
	std::cout << "IN THE FUNCTION" << std::endl;
	using matrix_t = Eigen::Array<num_t, Eigen::Dynamic, Eigen::Dynamic>;

	matrix_t mat(N, dim);

	const num_t max_range = 10;

	generateRandomPointCloud(mat, N, dim, max_range);

	std::cout << "random point cloud matrix: \n" << mat << std::endl;

    // Query point:
//    std::vector<num_t> query_pt(dim);

    using vec = Eigen::Array<num_t, Eigen::Dynamic,1>;
    vec query(dim);


    for (size_t d = 0; d < dim; d++){
//        query_pt[d] = max_range * (rand() % 1000) / num_t(1000);
//    	query(d) = query_pt[d];
    	query(d) = max_range * (rand() % 1000) / num_t(1000);
    }
//    std::cout << "Query point:\n" << query_pt[0] << '\t' << query_pt[1] << std::endl;
    std::cout << "QUERY: " << query << std::endl;

    using my_kd_tree_t = nanoflann::KDTreeEigenMatrixAdaptor<matrix_t>;

    my_kd_tree_t mat_index(dim, std::cref(mat), 10 /* max leaf */);

    // do a knn search
    const size_t        num_results = 3;
    std::vector<size_t> ret_indexes(num_results);
    std::vector<num_t>  out_dists_sqr(num_results);

    nanoflann::KNNResultSet<num_t> resultSet(num_results);

    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
//    mat_index.index_->findNeighbors(resultSet, &query_pt[0]);
    mat_index.index_->findNeighbors(resultSet, &query(0));

    std::cout << "knnSearch(nn=" << num_results << "): \n";

    for (size_t i = 0; i < resultSet.size(); i++){
    	std::cout << "ret_index[" << i << "]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
    }

	std::cout << "DONE" << std::endl;

}

template <typename Der>
void SPDS::generateRandomPointCloud(Eigen::ArrayBase<Der>& mat, const size_t N, const size_t dim, const typename Der::Scalar max_range){
//	error is in the last argument.
	std::cout << "generating " << N << " random points...\n";




	for (size_t i = 0; i < N; i++)
		for (size_t d = 0; d < dim; d++)
			mat(i, d) = max_range * (rand() % 1000) / typename Der::Scalar(1000);


}


#endif /* SPDS_H_ */


