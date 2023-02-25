
#include "getNodeType.h"
#include <iostream>
#include <iterator>
getNodeType::getNodeType(probParams& params, Mesh& m)
{
	if(params.dataRed){
		assignNodeTypesGrdy(m);
	}else{
		assignNodeTypes(m);
	}
}

void getNodeType::assignNodeTypes(Mesh& m){
	N_c = m.N_m;
	N_i = m.N_i;

	cPtr = &m.mNodes;
	iPtr = &m.iNodes;

	N_se = m.N_se;
	sePtr = &m.seNodes;

	N_ss = m.N_ss;
	ssPtr = &m.ssNodes;

	N_mStd = m.N_m + m.N_se + m.N_ss;
	mNodesStd.resize(N_mStd);
	mNodesStd << m.mNodes, m.seNodes, m.ssNodes;
	bPtr = &mNodesStd;


	N_s = m.N_se+m.N_ss;
	sNodes.resize(N_s);
	sNodes << m.seNodes, m.ssNodes;
	sPtr = &sNodes;
}



void getNodeType::assignNodeTypesGrdy(Mesh& m){

	N_i = m.N_m + m.N_se + m.N_ss;
	iNodes.resize(N_i);
	iNodes << m.mNodes, m.seNodes, m.ssNodes;

	iPtr = &iNodes;

	N_c = 0;
	cNodes.resize(N_c);
	cPtr = &cNodes;

	N_se = 0;
	seNodes.resize(N_se);
	sePtr = &seNodes;

	N_ss = 0;
	ssNodes.resize(N_ss);
	ssPtr = &ssNodes;


	iPtrGrdy = &m.iNodes;
	N_iGrdy = m.N_i;

	N_b = 0;
	bNodes.resize(N_b);
	bPtr = &bNodes;

	N_s = 0;
	sNodes.resize(N_s);
	sPtr = &sNodes;
}


void getNodeType::addControlNode(int node, std::string& smode, Mesh& m){

	if (std::find(std::begin(m.mNodes),std::end(m.mNodes),node) != std::end(m.mNodes)){
		N_c++;
		cNodes.conservativeResize(N_c);
		cNodes(N_c-1) = node;



	}else if(std::find(std::begin(m.seNodes),std::end(m.seNodes),node) != std::end(m.seNodes)){

		N_se++;
		seNodes.conservativeResize(N_se);
		seNodes(N_se-1) = node;
//		std::cout << seNodes << std::endl;
		N_s++;
		sNodes.resize(N_s);
		sNodes << seNodes, ssNodes;
	}else if(std::find(std::begin(m.ssNodes), std::end(m.ssNodes),node) != std::end(m.ssNodes)){

		N_ss++;
		ssNodes.conservativeResize(N_ss);
		ssNodes(N_ss-1) = node;

		N_s++;
		sNodes.resize(N_s);
		sNodes << seNodes, ssNodes;

	}



	if(smode != "none"){
		N_mStd++;
		mNodesStd.resize(N_mStd);
		mNodesStd << mNodes,seNodes, ssNodes;

	}

	// reducing the iNodes array:
	int idx;
	idx = std::distance(std::begin(iNodes), std::find(std::begin(iNodes), std::end(iNodes),node));
	N_i --;
	iNodes(Eigen::seqN(0,N_i)) << iNodes(Eigen::seqN(0,idx)), iNodes(Eigen::seq(idx+1,N_i));
	iNodes.conservativeResize(N_i);

//	std::cout << "control node: " << node << " is added"<< std::endl;
//	std::cout << "control node is added " << std::endl;

}





//void getNodeType::greedyNodes(int node,std::string smode){
////	std::cout << "in the greedyNodes function" << std::endl;
////	std::cout << "Node added to greedy control nodes: \t" << node << std::endl;
////	std::distance(std::begin(extBdryNodesMat.col(col)), std::find( std::begin(extBdryNodesMat.col(col)), std::end(extBdryNodesMat.col(col)), slidingSurfNodes(i)));
//	auto iterI = std::find(std::begin(m.intBdryNodes),std::end(m.intBdryNodes),node);
//	auto iterS = std::find(std::begin(m.slidingEdgeNodes),std::end(m.slidingEdgeNodes),node);
//	auto iterES = std::find(std::begin(m.extStaticNodes),std::end(m.extStaticNodes),node);
//	int idx;
//	if( iterI != std::end(m.intBdryNodes)){
//		// in case the added node is an internal boundary node
////		std::cout << "yes int bdry" << std::endl;
//
//		// update the included internal boundary nodes
//		N_ib++;
//		ibNodes.conservativeResize(N_ib);
//		ibNodes(N_ib-1) = node;
//
//
//		// find index of the node in the iNodes array
//		idx = std::distance(std::begin(iNodes),std::find(std::begin(iNodes),std::end(iNodes),node));
//
//		// shifting the nodes to remove the node from iNodes
//		for(int i = idx; i<N_i-1;i++ ){
//			iNodes(i) = iNodes(i+1);
//		}
//
//		// updating size of the iNodes array
//		N_i -= 1;
//		iNodes.conservativeResize(N_i);
////		std::cout << iNodes << std::endl;
//
//		// updating the moving nodes array
//		N_m++;
//		mNodes.resize(N_m);
//		if(smode == "none"){
//			mNodes << ibNodes, esNodes, sNodes;
//		}else{
//			mNodes << ibNodes, esNodes;
//		}
//
//		// updating the nodes that are used in the std rbf interpolation
//		N_mStd++;
//		mNodesStd.resize(N_mStd);
//		mNodesStd << ibNodes, esNodes, sNodes;
//
//	}else if( iterS != std::end(m.slidingEdgeNodes)){
////		std::cout << "yes slide" << std::endl;
//
//		N_s++;
//		sNodes.conservativeResize(N_s);
//		sNodes(N_s-1) = node;
//
//		idx = std::distance(std::begin(iNodes),std::find(std::begin(iNodes),std::end(iNodes),node));
//		// shifting the nodes to remove the node from iNodes
//		for(int i = idx; i<N_i-1;i++ ){
//			iNodes(i) = iNodes(i+1);
//		}
//		// updating size of the iNodes array
//		N_i -= 1;
//		iNodes.conservativeResize(N_i);
//
//		// updating the moving nodes array
//		if(smode == "none"){
//
//			N_m++;
//			mNodes.resize(N_m);
//			mNodes << ibNodes,esNodes, sNodes;
//		}
////		N_m++;
////		mNodes.resize(N_m);
////		mNodes << ibNodes, esNodes;
//
//		// updating the nodes that are used in the std rbf interpolation
//		//todo this part might not be used in case of standard rbf.
//
//		N_mStd++;
//		mNodesStd.resize(N_mStd);
//		mNodesStd << ibNodes, esNodes, sNodes;
//
////		std::cout << iNodes << std::endl;
////		std::cout << std::endl;
////		std::cout << sNodes << std::endl;
////		std::cout << std::endl;
////		std::cout << mNodes << std::endl;
////		std::cout << std::endl;
////		std::cout << mNodesStd << std::endl;
//
//	}else if( iterES != std::end(m.extStaticNodes)){
//		// in case the added node is an internal boundary node
////		std::cout << "yes ext stat" << std::endl;
//
//		// update the included internal boundary nodes
//		N_es++;
//		esNodes.conservativeResize(N_es);
//		esNodes(N_es-1) = node;
//
//		// find index of the node in the iNodes array
//		idx = std::distance(std::begin(iNodes),std::find(std::begin(iNodes),std::end(iNodes),node));
//
//		// shifting the nodes to remove the node from iNodes
//		for(int i = idx; i<N_i-1;i++ ){
//			iNodes(i) = iNodes(i+1);
//		}
//
//		// updating size of the iNodes array
//		N_i -= 1;
//		iNodes.conservativeResize(N_i);
////		std::cout << iNodes << std::endl;
//
//		// updating the moving nodes array
//		N_m++;
//		mNodes.resize(N_m);
//		if(smode == "none"){
//			mNodes << ibNodes, esNodes , sNodes;
//		}else{
//			mNodes << ibNodes, esNodes;
//		}
//
//		// updating the nodes that are used in the std rbf interpolation
//		N_mStd++;
//		mNodesStd.resize(N_mStd);
//		mNodesStd << ibNodes, esNodes, sNodes;
//
////				std::cout << iNodes << std::endl;
////				std::cout << std::endl;
////				std::cout << sNodes << std::endl;
////				std::cout << std::endl;
////				std::cout << mNodes << std::endl;
////				std::cout << std::endl;
////				std::cout << mNodesStd << std::endl;
//	}
//
//
//
//}
//
