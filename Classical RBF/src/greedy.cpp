#include "greedy.h"
#include <iostream>
#include <Eigen/Dense>
greedy::greedy()  {
//
////	std::cout << ptr->mNodes << std::endl;
////	std::cout << ptr->iNodes <<std::endl;
////	std::cout << ptr->sNodes(2) << std::endl;
//	std::cout << "In the greedy class" << std::endl;
////	int idx = 2;
//
////	std::cout << n.sNodes << std::endl;
//
//	n.GreedyInit();
//
//	n.greedyNodes(rbf.m.intBdryNodes(0));
//
//
////	rbf.params.steps = 1;
//
////	std::cout <<rbf.params.steps << std::endl;
//	rbf.m.getVecs();
//	rbf.perform_rbf(n);
////	std::cout << rbf.m.coords << std::endl;
//	int idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//	n.greedyNodes(n.iNodes(idxMax));
//
//
////	std::cout << n.iNodes << std::endl;
////	rbf.params.dataRed = false;
//
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	std::cout << n.iNodes << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	std::cout << n.iNodes << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
//
//	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//
//	std::cout << n.mNodesStd << std::endl;
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;

//	std::cout << n.iNodes << std::endl;

//	std::cout << rbf.d.rows() << std::endl;
//	Eigen::ArrayXi freeNodesIdx;
//	freeNodesIdx <<

//	ptr->mNodes.resize(ptr->m.N_ib + ptr->m.N_es + 1);
//	ptr->mNodes << ptr->m.intBdryNodes, ptr->m.extStaticNodes, ptr->m.slidingEdgeNodes(idx);
//
//	ptr->sNodes.resize(1);
//	ptr->sNodes << ptr->m.slidingEdgeNodes(idx);
//
//	ptr->iNodes.resize(ptr->m.N_se-1);
//	ptr->iNodes << ptr->m.slidingEdgeNodes(Eigen::seqN(0,idx)), ptr->m.slidingEdgeNodes(Eigen::seqN(idx+1,ptr->m.N_se-idx-1));
//	std::cout << ptr->iNodes << std::endl;
//	ptr->N_m = ptr->mNodes.size();
//	ptr->N_s = ptr->sNodes.size();
//	ptr->N_i = ptr->iNodes.size();
//	std::cout << ptr->mNodes << std::endl;
//	std::cout << ptr->sNodes << std::endl;
//
//	std::cout << '\n' << ptr->iNodes << std::endl;

//	ptr->iNodes.conservativeResize(iNodes.size()+)
//	ptr->mNodes.conservativeResize

//	Eigen::ArrayXXd coordsInit;
//	coordsInit = ptr->m.coords(ptr->iNodes,Eigen::all);
//	ptr->m.getVecs();
//
//	ptr->perform_rbf();
//
//	Eigen::ArrayXXd diff;
//	diff = coordsInit-ptr->m.coords(ptr->iNodes,Eigen::all);
//	std::cout << diff << std::endl;
//	getError(diff,ptr);

}


void greedy::getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, Eigen::VectorXd& exactDeformation, double& e, int& idxMax){
	std::cout << "obtaining error" << std::endl;
//	Eigen::VectorXd defVec;
//	int N = rbf.m.N_ib-n.N_ib;
//	defVec = Eigen::VectorXd::Zero(rbf.m.nDims*N);
////	std::cout << n.iNodes(Eigen::seq(rbf.m.N_i+rbf.m.N_p, rbf.m.N_i+rbf.m.N_p+rbf.m.N_ib-n.N_ib-1)) << std::endl;
//
//
//	Eigen::ArrayXi ibNodes = n.iNodes(Eigen::seq(rbf.m.N_i+rbf.m.N_p, rbf.m.N_i+rbf.m.N_p+rbf.m.N_ib-n.N_ib-1));
//	rbf.getDefVec(defVec, N, ibNodes);

	Eigen::ArrayXd error;
	error = Eigen::ArrayXd::Zero(n.N_i - meshOb.N_i-meshOb.N_p);
	int i,j,k;
	for(i=0; i<meshOb.N_ib-n.N_ib;i++){

		for(int dim=0; dim<meshOb.nDims;dim++){
//			std::cout << rbf.d(rbf.m.N_i+rbf.m.N_p + i,dim)-defVec(N*dim+i) << std::endl;
			error(i) += pow(d(meshOb.N_i+meshOb.N_p + i,dim)-exactDeformation((meshOb.N_ib-n.N_ib)*dim+i),2);

		}
		error(i) = sqrt(error(i));
//		std::cout << n.iNodes(rbf.m.N_i+rbf.m.N_p + i) << std::endl;

	}
	std::cout << "check" << std::endl;
	for(j = 0; j < meshOb.N_es - n.N_es ; j++){
//		std::cout << i + j << std::endl;
		for(int dim=0; dim<meshOb.nDims;dim++){
			error(i+j) += pow(d(meshOb.N_i+meshOb.N_p + i + j,dim),2);
		}
		error(i+j) = sqrt(error(i+j));
	}
	std::cout << "check" << std::endl;

//	std::cout << rbf.m.n << std::endl;
	int idx;
	for(k=0; k< meshOb.N_se - n.N_s; k++){
//		std::cout << n.iNodes(rbf.m.N_i+rbf.m.N_p+i+j+k) << std::endl;
		idx = std::distance(std::begin(meshOb.slidingEdgeNodes),std::find(std::begin(meshOb.slidingEdgeNodes),std::end(meshOb.slidingEdgeNodes),n.iNodes(meshOb.N_i+meshOb.N_p+i+j+k)));

//		std::cout << idx << std::endl;
//		std::cout << i+j+k << std::endl;
//		std::cout << rbf.m.n.row(k) << std::endl;
//		std::cout << rbf.d.row(rbf.m.N_i+rbf.m.N_p + i + j+k) << std::endl;
//		std::cout << std::abs(rbf.m.n.row(k).matrix().dot(rbf.d.row(rbf.m.N_i+rbf.m.N_p + i + j+k).matrix())) << std::endl;
		error(i+j+k) = std::abs(meshOb.n.row(idx).matrix().dot(d.row(meshOb.N_i+meshOb.N_p + i + j + k).matrix()));

	}

	std::cout << "check" << std::endl;

	std::cout << "error Array: \n" << error << '\n' << std::endl;
//	std::exit(0);

	error.maxCoeff(&idxMax);

//	std::cout << n.iNodes << std::endl;

//	std::cout << error(idxMax) << std::endl;
//	std::cout << meshOb.N_i+meshOb.N_p+idxMax << std::endl;
	e = error(idxMax);
//	std::cout << rbf.d << std::endl;
//	std::cout << rbf.d << std::endl;
//	std::cout << defVec << std::endl;
//	Eigen::ArrayXXd error;
//
//	std::cout << ptr->m.n << std::endl;

//	error = ptr->m.n

	idxMax = n.iNodes(meshOb.N_i+meshOb.N_p+idxMax);


}

