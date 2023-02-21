/*
 * getNodeType.cpp
 *
 *  Created on: 19 nov. 2022
 *      Author: floyd
 */

#include "getNodeType.h"
#include <iostream>
#include <iterator>
getNodeType::getNodeType(Mesh& meshOb, bool& dataRed)
:m(meshOb)
{
	std::cout << "done assigning node types" << std::endl;
	if(dataRed){
		assignNodeTypesGreedy();
	}else{
		assignNodeTypes();
	}
}

void getNodeType::assignNodeTypes(){
//
//	mNodes = m.mNodes;
//	N_m = m.N_m;
//	iNodes = m.iNodes;
//	N_i = m.N_i;
//	std::cout << *ptr << std::endl;


	N_m = m.N_m;
	N_i = m.N_i;
	mPtr = &m.mNodes;
	iPtr = &m.iNodes;

	N_se = m.N_se;
	sePtr = &m.seNodes;

	N_ss = m.N_ss;
	ssPtr = &m.ssNodes;

	N_mStd = m.N_m + m.N_se + m.N_ss;
	mNodesStd.resize(N_mStd);
	mNodesStd << m.mNodes, m.seNodes, m.ssNodes;
	mStdPtr = &mNodesStd;


	N_s = m.N_se+m.N_ss;
	sNodes.resize(N_s);
	sNodes << m.seNodes, m.ssNodes;
	sPtr = &sNodes;


//todo move this to a seperate function
//	N_ib = m.N_ib;
//	ibNodes.resize(N_ib);
//	ibNodes = m.intBdryNodes;
//
//	N_i = m.N_i;
//	iNodes.resize(N_i);
//	iNodes << m.iNodes;
//
////	N_es = m.N_es;
////	esNodes.resize(N_es);
////	esNodes << m.extStaticNodes;
//
//	if(m.smode == "none"){
//		N_m = m.N_m;
//		mNodes.resize(N_m);
//		mNodes << m.mNodes;
//
//	}else{
//		// only in case pmode != moving
//		N_m = m.N_m;
//		mNodes.resize(N_m);
//		mNodes << m.mNodes;
//
//		N_s = m.N_se;
//		sNodes.resize(N_s);
//		sNodes << m.seNodes;
//
//		N_mStd = m.N_mStd;
//		mNodesStd.resize(N_mStd);
//		mNodesStd << m.mNodesStd;
//	}
//		else if(m.smode == "ds"){ // only in case of non moving periodic boundaries.
//		N_m = m.N_ib + m.N_es;
//		mNodes.resize(N_m);
//		mNodes << m.intBdryNodes, m.extStaticNodes;
//
//		N_s = m.N_se;
//		sNodes.resize(N_s);
//		sNodes << m.slidingEdgeNodes;
//
//		// todo rename to make clearer
//		N_mStd = m.N_ib + m.N_es + m.N_se;
//		mNodesStd.resize(N_mStd);
//		mNodesStd << m.intBdryNodes,m.extStaticNodes, m.slidingEdgeNodes;
//
//	}

	std::cout << "done" << std::endl;


}



void getNodeType::assignNodeTypesGreedy(){

	// for the non sliding rbf

	N_i = m.N_m + m.N_se + m.N_ss;
	iNodes.resize(N_i);
	iNodes << m.mNodes, m.seNodes, m.ssNodes;

	iPtr = &iNodes;

	N_m = 0;
	mNodes.resize(N_m);
	mPtr = &mNodes;

	N_se = 0;
	seNodes.resize(N_se);
	sePtr = &seNodes;

	N_ss = 0;
	ssNodes.resize(N_ss);
	ssPtr = &ssNodes;


	iPtrGrdy = &m.iNodes;
	N_i_grdy = m.N_i;

	N_mStd = 0;
	mNodesStd.resize(N_mStd);
	mStdPtr = &mNodesStd;

	N_s = 0;
	sNodes.resize(N_s);
	sPtr = &sNodes;

//	N_m = 0;
//	N_s = 0;
//	N_mStd = 0;
//	N_ib = 0;
//	N_es = 0;
//	N_i = m.N_i+m.N_p+m.N_ib+ m.N_es + m.N_se;
//
//	iNodes.resize(N_i);
//	iNodes <<  m.intNodes, m.periodicNodes, m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;
//
//	esNodes.resize(N_es);
//	mNodes.resize(N_m);
//	sNodes.resize(N_s);
//	mNodesStd.resize(N_mStd);
//	ibNodes.resize(N_ib);

}


void getNodeType::addControlNode(int node){



	if (std::find(std::begin(m.mNodes),std::end(m.mNodes),node) != std::end(m.mNodes)){
		N_m++;
		mNodes.conservativeResize(N_m);
		mNodes(N_m-1) = node;



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



	if(m.smode != "none"){


		N_mStd++;
		mNodesStd.resize(N_mStd);
//		mNodesStd(N_mStd-1) = node;
		mNodesStd << mNodes,seNodes, ssNodes;
//		std::cout << mNodesStd << std::endl;

	}


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
