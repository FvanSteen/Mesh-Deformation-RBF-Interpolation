#include "getNodeType.h"
#include <iostream>
#include <iterator>


/* getNodeType
 *
 * This class assigns the various node types based on the problem parameters
 *
 */

getNodeType::getNodeType(probParams& params, Mesh& m)
{
	// For the pseudo sliding in 3D a reduced set of internal nodes will be defined
	// this boolean triggers that statement
	if(params.smode == "ps" && m.nDims == 3){
		pseudo = true;
	}

	// assign node types based on whether data reduction is applied
	if(params.dataRed){
		assignNodeTypesGrdy(m);
	}else{
		assignNodeTypes(m);
	}

}

/* assignNodeTypes
 *
 * set various points to the moving, internal, sliding edge, sliding surface and control nodes
 *
 */
void getNodeType::assignNodeTypes(Mesh& m){
	N_m = m.N_m;
	N_i = m.N_i;

	mPtr = &m.mNodes;
	iPtr = &m.iNodes;

	N_se = m.N_se;
	sePtr = &m.seNodes;

	N_ss = m.N_ss;
	ssPtr = &m.ssNodes;


	N_c = m.N_m + m.N_se + m.N_ss;
	cNodes.resize(N_c);
	cNodes << m.mNodes, m.seNodes, m.ssNodes;
	cPtr = &cNodes;
}

/* assignNodeTypes
 *
 * set various points to the moving, internal, sliding edge, sliding surface and control nodes
 * In this case the list of control nodes is empty and all nodes are included in the internal node array
 * The greedy algorithm will add control nodes, which will be removed from the internal node array
 */

void getNodeType::assignNodeTypesGrdy(Mesh& m){


	N_i = m.N_m + m.N_se + m.N_ss;
	iNodes.resize(N_i);
	iNodes << m.mNodes, m.seNodes, m.ssNodes;
	iPtr = &iNodes;

	if(pseudo){
		iNodesReduced = iNodes(Eigen::seqN(0,m.N_m+m.N_se));
		iPtr_reduced = &iNodesReduced;
	}

	iNodesIdx = Eigen::ArrayXi::LinSpaced(N_i, 0, N_i-1);
	cNodesIdx.resize(0);

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
	N_iGrdy = m.N_i;

	N_c = 0;
	cNodes.resize(N_c);
	cPtr = &cNodes;
}

/* addControlNodes
 *
 * Function for adding an array of control nodes
 */

void getNodeType::addControlNodes(Eigen::ArrayXi& nodes, std::string& smode, Mesh& m){

	addedNodes.idx.resize(nodes.size());
	addedNodes.idx_i.resize(nodes.size());
	addedNodes.type.resize(nodes.size());
	for(int i = 0; i < nodes.size(); i++){
		addControlNode(nodes[i], smode, m, i);
	}
}


/* addControlNodes
 *
 * Function for the addition of a single control node
 * The control nodes are always ordered by the sequence moving nodes, sliding edge nodes and sliding surface nodes
 *
 */

void getNodeType::addControlNode(int node, std::string& smode, Mesh& m, int i){
	// check if the node is a moving node with known displacement
	if (std::find(std::begin(m.mNodes),std::end(m.mNodes),node) != std::end(m.mNodes)){
		N_m++;
		mNodes.conservativeResize(N_m);
		mNodes(N_m-1) = node;

		// saving the index and type of added control node required for the update of the interpolation matrices
		addedNodes.idx[i] = N_m-1;
		addedNodes.type[i] = 0;

	// check if the node is among the sliding edge nodes
	}else if(std::find(std::begin(m.seNodes),std::end(m.seNodes),node) != std::end(m.seNodes)){
		N_se++;
		seNodes.conservativeResize(N_se);
		seNodes(N_se-1) = node;

		// saving the index and type of added control node required for the update of the interpolation matrices
		addedNodes.idx[i] = N_se-1;
		addedNodes.type[i] = 1;


	// check is the node is a silding surface node
	}else if(std::find(std::begin(m.ssNodes), std::end(m.ssNodes),node) != std::end(m.ssNodes)){
		N_ss++;
		ssNodes.conservativeResize(N_ss);
		ssNodes(N_ss-1) = node;

		// saving the index and type of added control node required for the update of the interpolation matrices
		addedNodes.idx[i] = N_ss-1;
		addedNodes.type[i] = 2;

	}


	// update the control node size and array
	N_c++;
	cNodes.resize(N_c);
	cNodes << mNodes, seNodes, ssNodes;


	// removing the node from the iNodes array
	// saving the index of the removed internal node for updating the interpolation matrices
	int idx;
	idx = std::distance(std::begin(iNodes), std::find(std::begin(iNodes), std::end(iNodes),node));
	addedNodes.idx_i[i] = idx;

	N_i --;
	iNodes(Eigen::seqN(0,N_i)) << iNodes(Eigen::seqN(0,idx)), iNodes(Eigen::seq(idx+1,N_i));
	iNodes.conservativeResize(N_i);


	// Saving the indices of the control nodes in the initial iNodes array for the multi level greedy algorithm
	// cNodesIdx contains the order of the control node set (ordered as moving nodes, sliding edge nodes, sliding surf nodes)
	cNodesIdx.conservativeResize(N_m+N_se+N_ss);

	switch(addedNodes.type[i]){
		case 0:
			cNodesIdx(Eigen::seqN(N_m,N_se+N_ss)) = cNodesIdx(Eigen::seqN(N_m-1,N_se+N_ss)).eval();
			cNodesIdx(N_m-1) = iNodesIdx(idx);
			break;

		case 1:
			cNodesIdx(Eigen::seqN(N_m+N_se,N_ss)) = cNodesIdx(Eigen::seqN(N_m+N_se-1,N_ss)).eval();
			cNodesIdx(N_m+N_se-1) = iNodesIdx(idx);
			break;

		case 2:
			cNodesIdx(N_m+N_se+N_ss-1) = iNodesIdx(idx);
			break;

	}

	// remaining indices of the unselected boundary nodes
	iNodesIdx(Eigen::seqN(0,N_i)) << iNodesIdx(Eigen::seqN(0,idx)), iNodesIdx(Eigen::seq(idx+1,N_i));
	iNodesIdx.conservativeResize(N_i);

	// In case of pseudo sliding a reduced set of internal nodes is used to seperately perform the sliding on edges and surfaces
	if(pseudo){
		iNodesReduced = iNodes(Eigen::seqN(0,N_i- m.N_ss + N_ss));
	}
}




