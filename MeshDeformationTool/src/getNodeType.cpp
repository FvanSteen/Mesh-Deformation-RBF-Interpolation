
#include "getNodeType.h"
#include <iostream>
#include <iterator>
getNodeType::getNodeType(probParams& params, Mesh& m)
{
	if(params.dataRed){
//		if(params.smode == "ps" && m.nDims ==3){
//			assignNodeTypesGrdy(m, params.smode);
//		}else{
		assignNodeTypesGrdy(m);
//		}
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

	N_b = m.N_m + m.N_se + m.N_ss;
	bNodes.resize(N_b);
	bNodes << m.mNodes, m.seNodes, m.ssNodes;
	bPtr = &bNodes;


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

	iNodesIdx = Eigen::ArrayXi::LinSpaced(N_i, 0, N_i-1);
	cNodesIdx.resize(0);

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
	// check if the node is a moving node with known displacement
	if (std::find(std::begin(m.mNodes),std::end(m.mNodes),node) != std::end(m.mNodes)){
		N_c++;
		cNodes.conservativeResize(N_c);
		cNodes(N_c-1) = node;
	// check if the node is among the sliding edge nodes
	}else if(std::find(std::begin(m.seNodes),std::end(m.seNodes),node) != std::end(m.seNodes)){
		N_se++;
		seNodes.conservativeResize(N_se);
		seNodes(N_se-1) = node;

		N_s++;
		sNodes.resize(N_s);
		sNodes << seNodes, ssNodes;
	// check is the node is a silding surface node
	}else if(std::find(std::begin(m.ssNodes), std::end(m.ssNodes),node) != std::end(m.ssNodes)){

		N_ss++;
		ssNodes.conservativeResize(N_ss);
		ssNodes(N_ss-1) = node;

		N_s++;
		sNodes.resize(N_s);
		sNodes << seNodes, ssNodes;

	}



	if(smode != "none"){
		N_b++;
		bNodes.resize(N_b);
		bNodes << cNodes, seNodes, ssNodes;
	}


	// removing the node from the iNodes array
	int idx;
	idx = std::distance(std::begin(iNodes), std::find(std::begin(iNodes), std::end(iNodes),node));

	N_i --;
	iNodes(Eigen::seqN(0,N_i)) << iNodes(Eigen::seqN(0,idx)), iNodes(Eigen::seq(idx+1,N_i));
	iNodes.conservativeResize(N_i);

	// cNodesIdx contains the order of the control node set (ordered as moving nodes, sliding edge nodes, sliding surf nodes)
	cNodesIdx.conservativeResize(N_c+N_s);

	if(idx <= (m.N_m-N_c)){
		cNodesIdx(Eigen::seqN(N_c,N_s)) = cNodesIdx(Eigen::seqN(N_c-1,N_s)).eval();
		cNodesIdx(N_c-1) = iNodesIdx(idx);
	}else if(idx <= (m.N_m-N_c) + (m.N_se-N_se)){
		cNodesIdx(Eigen::seqN(N_c+N_se,N_ss)) = cNodesIdx(Eigen::seqN(N_c+N_se-1,N_ss)).eval();
		cNodesIdx(N_c+N_se-1) = iNodesIdx(idx);
	}else{
		cNodesIdx(N_c+N_se+N_ss-1) = iNodesIdx(idx);
	}

	// remaining indices of the unselected boundary nodes
	iNodesIdx(Eigen::seqN(0,N_i)) << iNodesIdx(Eigen::seqN(0,idx)), iNodesIdx(Eigen::seq(idx+1,N_i));
	iNodesIdx.conservativeResize(N_i);

}




