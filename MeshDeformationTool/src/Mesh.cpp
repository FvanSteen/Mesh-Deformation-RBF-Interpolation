#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <chrono>
#include <Eigen/Dense>
#include "SPDS.h"
#include "CoordTransform.h"



Mesh::Mesh(probParams& params, const int& debugLvl)
:lvl(debugLvl)
{
readMeshFile(params);
r = params.rFac*getCharDomLength();
getPeriodicParams(params);
if(params.ptype){
	coords_polar_cylindrical.resize(nNodes,nDims);
	ptrCoords = &coords_polar_cylindrical;
}else{
	ptrCoords = &coords;
}
}

// Main function for reading the .su2 mesh files

void Mesh::readMeshFile(probParams& params){
	if(lvl>=1){
		std::cout << "Reading mesh file: " << params.mesh_ifName << std::endl;
	}

	int lineNo = 0;								// line number counter
	int nBdryElemsTotal = 0;					// stores total number of boundary elements

	int markerIdx = -2;							// contain line index containing the latest marker found

//	int nBdryNodes = 0;							// counter of the total number of boundary nodes

	int nodeCnt = 0;							// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is

	int bdryElemCnt = 0; 						// counting the elements of the boundaries
	int markerElems;							// locally stores how many elements are in that boundary
	int nMarker = 0; 							// Counts the number of external boundary markers
	nNodes = -1; // set to a default so some if statements are not triggered //todo check which ones
	int idxDelta;
	int firstIdx, lastIdx;
	bdryNodesMat.resize(0,3);//, intBdryNodesMat.resize(0,3);	// the int/ ext boundary node arrays have a minimum of 3 columns. One for the node type and at least two node indices.
																// The array will be adjusted to appropriate size depending on the boundary element type.
	srtdTags.resize(params.bdryTags.size());
	nrElemsBdry.resize(params.bdryTags.size());		// Array containing the sizes of each ext boundary
	std::string line;							// string containing line obtained by getline() function
	std::ifstream mFile("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes\\" + params.mesh_ifName); 	//opening file name stored in mFile object
	// Check if file is opened
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){
			if (line.rfind("NDIME= ",0)==0){							// save number of dimensions
				nDims = stoi(line.substr(7));
			}
			else if (line.rfind("NPOIN= ",0)==0){// save nr of points
				nNodes = stoi(line.substr(7));
				pntsIdx = lineNo;										// save line nr.
				coords.resize(nNodes, nDims);							// resizing the array containing the coordinates
				std::cout << "Saving node coordinates" << std::endl;
			}
//			else if (line.rfind("NELEM= ",0)==0){						// save nr of elements
//				nElem = stoi(line.substr(7));
//
//			}

			// Checking whether provided boundary tags equals the amount in the mesh file
			else if (line.rfind("NMARK=",0)==0){
				findStringBounds(firstIdx,lastIdx,line);
				try{
					if(stoi(line.substr(firstIdx, lastIdx-firstIdx)) != int(params.bdryTags.size())){
						throw(stoi(line.substr(firstIdx, lastIdx-firstIdx)));
					}
				}
				catch(int nTag){
					std::cout << "Number of tags provided in config file (" << params.bdryTags.size() << ") does not match number of tags found in mesh file (" << nTag << ")." << std::endl;
					std::exit(0);
				}
			}

			// Finding tags of the boundaries
			else if (line.rfind("MARKER_TAG=",0)==0){
				// start here
				findStringBounds(firstIdx, lastIdx,line);
				std::string tag =  line.substr(firstIdx, lastIdx-firstIdx);
				try{
					if(std::find(std::begin(params.bdryTags), std::end(params.bdryTags), tag) != std::end(params.bdryTags)){
						if(lvl >=2){
							std::cout << "Saving nodes of boundary: " << tag << std::endl;
						}
						markerIdx = lineNo;
						srtdTags[nMarker] = tag;
						if(tag == "SEND_RECEIVE"){
							idxDelta = 2;
						}else{
							idxDelta = 1;
						}

					}
					else throw(tag);
				}
				// If the provided tag does not match any tag in the meshing file an error is thrown
				catch(std::string& tag){
					std::cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << std::endl;
					std::exit(0);
				}
			}

			// Check whether line corresponds to ext. boundary
			else if(lineNo == markerIdx+1){

				findStringBounds(firstIdx,lastIdx,line);
				markerElems = stoi(line.substr(firstIdx, lastIdx-firstIdx));

					// Updating number of external boundary elements
				nBdryElemsTotal+= markerElems;
					// Saving number of elements of each boundary in an array
				nrElemsBdry(nMarker) = markerElems;
					// resizing the array containing the external boundary data
				bdryNodesMat.conservativeResize(nBdryElemsTotal,bdryNodesMat.cols());
				nMarker++;

				}

			else if(markerIdx > 0 && lineNo > markerIdx+idxDelta && lineNo <= markerIdx+idxDelta + markerElems){


				// split line by '\t' character
				std::istringstream is(line);
				// elements on the line are assigned to data in a while loop. lineElem counts the number of elements per line
				int data, lineElem = 0;
				while(is >> data){
					// First integer on line describes the type of element.
					// Triangle: 5,	Quadrilateral: 9

					// following if statements check the elementtype and if needed resize the intBdryNodesMat array.
					if(lineElem == 0 && data == 5 && bdryNodesMat.cols() < 4){
						bdryNodesMat.conservativeResize(nBdryElemsTotal,4);
					}
					else if(lineElem == 0 && data == 9 && bdryNodesMat.cols() < 5){
						bdryNodesMat.conservativeResize(nBdryElemsTotal,5);
					}

					// assigning the node information to the array
					bdryNodesMat(bdryElemCnt,lineElem) = data;

					lineElem++;
				}
				bdryElemCnt++;
			}
			// Check if line corresponds to line containing node coordinates
			if (lineNo > pntsIdx && nodeCnt < nNodes){

				// split line by '\t' character
				std::istringstream is(line);
				// based on the number of dimensions 2 or 3 coordinates are assigned to the coords array
				switch(nDims){
					case 2:
						is >> coords(nodeCnt,0) >> coords(nodeCnt,1);
						break;
					case 3:
						is >> coords(nodeCnt,0) >> coords(nodeCnt,1) >> coords(nodeCnt,2);
						break;
				}
				nodeCnt++;

			}
			lineNo++;

		}
		//closing meshing file
		mFile.close();
	}
	// If the file is not opened then the following error message will be displayed
	else std::cout << "Not able to open input mesh file";


	getNodeTypes(params, nMarker);
	std::cout << "established node types\n";


	int idx = srtdTags.size();
	int cnt = 0;
	std::vector<int> SEND_RECEIVE_idx;
	while(cnt < idx ){
		if(srtdTags[cnt] == "SEND_RECEIVE"){
			SEND_RECEIVE_idx.push_back(cnt);
		}
		cnt++;
	}




	if(SEND_RECEIVE_idx.size() > 0 ){

		Eigen::ArrayXi oldIdx(nrElemsBdry[SEND_RECEIVE_idx[0]]/2), newIdx(nrElemsBdry[SEND_RECEIVE_idx[1]]/2);

		int startIdxOld = nrElemsBdry(Eigen::seqN(0,SEND_RECEIVE_idx[0])).sum();
		int startIdxNew = nrElemsBdry(Eigen::seqN(0,SEND_RECEIVE_idx[1])).sum()+oldIdx.size();

		for(int i = 0; i < oldIdx.size(); i++){
			oldIdx(i) = bdryNodesMat(startIdxOld+i,1);
			newIdx(i) = bdryNodesMat(startIdxNew+i,1);
		}

		removeMutualNodes(mNodes, N_m, newIdx);
		removeMutualNodes(seNodes, N_se, newIdx);
		removeMutualNodes(ssNodes, N_ss, periodicEdgeNodes);
		/*
		Eigen::ArrayXi seNodesIdx(seNodes.size());
		int cnt = 0;
		for(int i = 0; i < seNodes.size(); i++ ){
			if(std::find(std::begin(newIdx), std::end(newIdx), seNodes(i)) == std::end(newIdx)){
				seNodesIdx(cnt) = i;
				cnt++;
			}
		}

		seNodes(Eigen::seqN(0,cnt)) = seNodes(seNodesIdx(Eigen::seqN(0,cnt))).eval();
		seNodes.conservativeResize(cnt);
		N_se = cnt;




		Eigen::ArrayXi mNodesIdx(mNodes.size());
		cnt = 0;
		for(int i = 0; i < mNodes.size(); i++ ){
			if(std::find(std::begin(newIdx), std::end(newIdx), mNodes(i)) == std::end(newIdx)){
				mNodesIdx(cnt) = i;
				cnt++;
			}
		}

		mNodes(Eigen::seqN(0,cnt)) = mNodes(mNodesIdx(Eigen::seqN(0,cnt))).eval();
		mNodes.conservativeResize(cnt);
		N_m = cnt;

		Eigen::ArrayXi ssNodesIdx(ssNodes.size());
		cnt = 0;
		for(int i = 0; i < ssNodes.size(); i++ ){
			if(std::find(std::begin(seNodes), std::end(seNodes), ssNodes(i)) == std::end(seNodes)){
				ssNodesIdx(cnt) = i;
				cnt++;
			}
		}

		ssNodes(Eigen::seqN(0,cnt)) = ssNodes(ssNodesIdx(Eigen::seqN(0,cnt))).eval();
		ssNodes.conservativeResize(cnt);
		N_ss = cnt;
		*/
	}







//	std::cout <<"moving Nodes: \n" <<  mNodes << std::endl;
//	std::cout <<"sliding edge Nodes: \n" <<  seNodes << std::endl;
//	std::cout <<"sliding surf Nodes: \n" <<  ssNodes << std::endl;


//	std::cout << "internal Nodes: \n" << iNodes << std::endl;

	if(nDims == 3 && (params.smode == "ds" || params.smode == "ps") && (params.pmode == "fixed" || params.pmode == "moving")){
		// adding the periodic edge nodes to the sliding surface nodes
		N_ss += N_pe;
		ssNodes.conservativeResize(N_ss);
		ssNodes(Eigen::lastN(N_pe)) = periodicEdgeNodes;


		// removing the periodic edge nodes from the sliding edge node selection
		removeMutualNodes(seNodes, N_se, periodicEdgeNodes);
	}

	removeMutualNodes(seNodes, N_se, mNodes);

	removeMutualNodes(ssNodes, N_ss, seNodes);

	std::cout <<  "MOVING\n";
	for(int i = 0; i < mNodes.size(); i++){
		std::cout << mNodes(i) << ", ";
	}
	std::cout << std::endl;
	std::cout <<  "surf\n";
	for(int i = 0; i < ssNodes.size(); i++){
		std::cout << ssNodes(i) << ", ";
	}
	std::cout << std::endl;
	std::cout <<  "per edge\n";
	for(int i = 0; i < periodicEdgeNodes.size(); i++){
		std::cout << periodicEdgeNodes(i) << ", ";
	}
	std::cout << std::endl;
	std::cout <<  "per vertices\n";
	for(int i = 0; i < periodicVerticesNodes.size(); i++){
		std::cout << periodicVerticesNodes(i) << ", ";
	}
	std::cout << std::endl;
	std::cout <<  "vertices\n";
		for(int i = 0; i < verticesNodes.size(); i++){
			std::cout << verticesNodes(i) << ", ";
		}
	std::cout << std::endl;

	std::cout <<  "edge\n";
	for(int i = 0; i < N_se; i++){
		std::cout << seNodes(i) << ", ";
	}
	std::cout << std::endl;

	std::cout << N_m << '\t' << N_se << '\t'<< N_ss << std::endl;

//	Eigen::ArrayXi allNodes(N_m+N_se+N_ss);
//	allNodes << mNodes, seNodes, ssNodes;
//	std::cout << "size: " << allNodes.size() << std::endl;
//	removeDuplicates(allNodes);
//	std::cout << "size: " << allNodes.size() << std::endl;

	getIntNodes();

	if(params.smode != "none"){
		getEdgeConnectivity(params.pmode, edgeConnectivity,seNodes,  N_se - N_periodic_vertices);
		if(nDims == 3){
			getExtBdryEdgeSegments(); // todo check if N_pe should be considered or not
			getSurfConnectivity();
			if(params.pmode == "fixed" || "moving"){
				getEdgeConnectivity(params.pmode, edgeConnectivityPeriodic, periodicEdgeNodes, N_pe);
			}
		}
	}



	if(params.pmode != "none"){
		getCharPerLength(params);
	}

//	if(params.dataRed){
//		std::cout << mNodes << std::endl;
//		std::cout << "HERE" << std::endl;
//		std::exit(0);
//	}



//todo next part should be replaced by a kd tree implementation
//	if(params.dataRed){
//		getIntCorNodes();
//	}
	std::cout << "Mesh file read successfully" << std::endl;
}


void Mesh::removeMutualNodes(Eigen::ArrayXi& array_in, int& size, Eigen::ArrayXi& to_remove_nodes){
	Eigen::ArrayXi idx_keep(array_in.size());
	int cnt = 0;
	for(int i = 0; i < array_in.size(); i++ ){
		if(std::find(std::begin(to_remove_nodes), std::end(to_remove_nodes), array_in(i)) == std::end(to_remove_nodes)){
			idx_keep(cnt) = i;
			cnt++;
		}
	}

	array_in(Eigen::seqN(0,cnt)) = array_in(idx_keep(Eigen::seqN(0,cnt))).eval();
	array_in.conservativeResize(cnt);
	size = cnt;
}

/* getNodeTypes function
 *
 * This function will go through each boundary specified in the mesh file individually. *
 * In case of a 3D mesh the sliding surface nodes are identified by the fact that these nodes will be part of 3 or more surfaces and are therefore
 * as many times present in the array containing all external boundary nodes.
 * The sliding edge nodes are identified by the fact that these nodes only appear twice in the external boundary nodes array.
 * The remaining elements only appear once in the array and are therefore static external nodes.  *
 * For 2D meshes the sliding surface nodes don't have to be found.
 */
void Mesh::getNodeTypes(probParams& params, int nMarker){

	if(lvl>=1){
		std::cout << "Obtaining node types" << std::endl;
	}

	Eigen::ArrayXi bdryNodesArr; 						// 1D array that will contain the all nodes for each respective boundary
	Eigen::ArrayXi idxMoving, idxSlidingEdge, idxSlidingSurf;		// Arrays containing specific type of nodes
	int cntMoving, cntSlidingEdge, cntSlidingSurf;								// counters for the number of sliding surface (SS) sliding edge (SE), static (Stat) and periodic (Per) nodes

	N_nonzeroDisp = 0;	//number of nodes with nonzero displacement
	N_periodic_vertices = 0;
	verticesNodes.resize(10);
	periodicVerticesNodes.resize(10);
	int verticesCnt = 0;
//	int edgeNodeCnt = 0;
	bool periodic, moving;										// boolean that is set based on whether its a periodic boundary element or not.


	// for each external boundary the sliding edge, sliding surface and static nodes are identified
	for(int elem = 0; elem < nMarker; elem++){
//		std::cout << "MARKER: " << srtdTags[elem] << std::endl;
		// resizing the array that will contain all the boundary nodes of that boundary

		if(srtdTags[elem] == "SEND_RECEIVE"){
			bdryNodesArr.resize(0);
		}else{
			bdryNodesArr.resize(nrElemsBdry(elem)*(bdryNodesMat.cols()-1));

			// next two loops ensure that all nodes are included in the 1D node array.
			for(int i =0; i<nrElemsBdry(elem);i++){
				for(int j=0; j< bdryNodesMat.cols()-1; j++ ){
					bdryNodesArr(j+i*(bdryNodesMat.cols()-1)) = bdryNodesMat(i+ nrElemsBdry(Eigen::seqN(0,elem)).sum(),j+1);
				}
			}
		}

		// sorting the array such that the nodes are ascending
		std::sort(std::begin(bdryNodesArr),std::end(bdryNodesArr));

		// Set size of the arrays of the various nodes equal to the size of the bdryNodesArr.
		if(idxMoving.size() != bdryNodesArr.size()){
			idxMoving.resize(bdryNodesArr.size());
			idxSlidingEdge.resize(bdryNodesArr.size());
			idxSlidingSurf.resize(bdryNodesArr.size());
		}

		// setting counters to zero
		cntMoving = 0, cntSlidingEdge = 0, cntSlidingSurf = 0;

		moving = false;
		periodic = false;

		if(std::find(std::begin(params.mTags),std::end(params.mTags),srtdTags[elem]) != std::end(params.mTags)){
			std::cout << srtdTags[elem] << " is a moving boundary" << std::endl;
			moving = true;

		}else if((params.pmode == "moving" || params.pmode == "fixed") && params.smode != "none" && std::find(std::begin(params.pTags),std::end(params.pTags),srtdTags[elem]) != std::end(params.pTags)){
			std::cout << srtdTags[elem] << " is periodic boundary" << std::endl;
			periodic = true;
		}

		// if the marker is not periodic
		if(periodic == false){
			// if the marker is part of a moving boundary then save it as such
			if(moving){
				idxMoving(Eigen::seqN(0,bdryNodesArr.size())) = bdryNodesArr;
				cntMoving += bdryNodesArr.size();
				removeDuplicates(bdryNodesArr);
				N_nonzeroDisp += bdryNodesArr.size();
			}
			// else its either a sliding or moving node
			else{
				for (int i= 0; i< bdryNodesArr.size();i++){

					// for 3D, in case of 4 subsuquent equal nodes its a surface node
					if(nDims == 3 && (i< bdryNodesArr.size()-2 && bdryNodesArr(i) == bdryNodesArr(i+2))){


						int n = 2;
						while(bdryNodesArr(i) == bdryNodesArr(i+n) && i+n < bdryNodesArr.size()){
							n++;
						}

						// if no sliding is allowed then its a moving node, else its a sliding surface node
						if(params.smode == "none"){
							idxMoving(cntMoving) = bdryNodesArr(i);
							cntMoving++;
						}else{
							idxSlidingSurf(cntSlidingSurf) = bdryNodesArr(i);
							cntSlidingSurf++;
						}
						i +=(n-1);

					}

					// else in case 2 subsequent nodes are equal its an edge node
					else if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){
						// if no sliding is allowed then its a moving node, else its a sliding edge node
						if(params.smode == "none"){
							idxMoving(cntMoving) = bdryNodesArr(i);
							cntMoving++;

						}else{
							idxSlidingEdge(cntSlidingEdge) = bdryNodesArr(i);
							cntSlidingEdge++;
						}
						i++;
					}
					// else in case there is a single occurence of the node in that bdry, its a corner node.
					else{

						if(params.pmode !=  "moving" || params.smode == "none"){
							idxMoving(cntMoving) = bdryNodesArr(i);
							cntMoving++;
						}
							verticesNodes(verticesCnt) = bdryNodesArr(i);
							verticesCnt++;
							if(verticesCnt == verticesNodes.size()){
								verticesNodes.conservativeResize(verticesNodes.size()+10);
							}
							if(params.pmode == "moving"){
								periodicVerticesNodes(N_periodic_vertices) = bdryNodesArr(i);
								N_periodic_vertices++;
								if(N_periodic_vertices == periodicVerticesNodes.size()){
									periodicVerticesNodes.conservativeResize(periodicVerticesNodes.size()+10);
								}
//								std::cout << "deze: " << bdryNodesArr(i) << std::endl;
//							}
						}
					}
				}
			}
		}else{
			Eigen::ArrayXi idxPeriodicEdgeNodes(bdryNodesArr.size());
			int cntPeriodicEdgeNodes = 0;
//			for(int x = 0; x < bdryNodesArr.size(); x++){
//				std::cout << bdryNodesArr(x) << ", ";
//			}

			for(int i = 0; i < bdryNodesArr.size(); i++){

				if ( (i< bdryNodesArr.size()-2 && bdryNodesArr(i) == bdryNodesArr(i+2))){

					int n = 2;
					while(bdryNodesArr(i) == bdryNodesArr(i+n) && i+n < bdryNodesArr.size()){
						n++;
					}

					i+= (n-1);
				}else if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){
					idxPeriodicEdgeNodes(cntPeriodicEdgeNodes) = bdryNodesArr(i);
					cntPeriodicEdgeNodes++;
					i++;
				}else{
					idxMoving(cntMoving) = bdryNodesArr(i);
					cntMoving++;
				}
			}
			periodicEdgeNodes.conservativeResize(periodicEdgeNodes.size()+cntPeriodicEdgeNodes);
			periodicEdgeNodes(Eigen::lastN(cntPeriodicEdgeNodes)) = idxPeriodicEdgeNodes(Eigen::seqN(0,cntPeriodicEdgeNodes));

		}

//		std::cout << "\n\n\n" << bdryNodesArr << std::endl;
//		std::cout << srtdTags[elem] << '\t' << cntMoving << '\t' << cntSlidingEdge << '\t' << cntSlidingSurf << std::endl;
//		for(int x = 0; x < bdryNodesArr.size(); x++){
//			std::cout << bdryNodesArr(x) << ", ";
//		}
//		std::cout << std::endl;


//		std::cout << "moving nodes: \n" << idxMoving(Eigen::seqN(0,cntMoving)) << std::endl;
//		std::cout << "slidingEdge: \n" << idxSlidingEdge(Eigen::seqN(0,cntSlidingEdge)) << std::endl;
//		std::cout << "sliding surf: \n" << idxSlidingSurf(Eigen::seqN(0,cntSlidingSurf)) << std::endl;
//		std::cout << "static: \n" << idxStatic(Eigen::seqN(0,cntStatic)) << std::endl;

		// Resize the arrays containing the different type of nodes according to how many are found in this boundary
		// And assigning last n-elements to the found nodes.

		mNodes.conservativeResize(mNodes.size()+cntMoving);
		mNodes(Eigen::lastN(cntMoving)) = idxMoving(Eigen::seqN(0,cntMoving));

		seNodes.conservativeResize(seNodes.size()+cntSlidingEdge);
		seNodes(Eigen::lastN(cntSlidingEdge)) = idxSlidingEdge(Eigen::seqN(0,cntSlidingEdge));

		ssNodes.conservativeResize(ssNodes.size()+cntSlidingSurf);
		ssNodes(Eigen::lastN(cntSlidingSurf)) = idxSlidingSurf(Eigen::seqN(0,cntSlidingSurf));

	}

	// Calling a function that removes any duplicate elements from the arrays.
	if(ssNodes.size() != 0){
		removeDuplicates(ssNodes);
	}

	if(seNodes.size() != 0 ){
		removeDuplicates(seNodes);
	}

	removeDuplicates(mNodes);



	verticesNodes.conservativeResize(verticesCnt);
	removeDuplicates(verticesNodes);

	periodicVerticesNodes.conservativeResize(N_periodic_vertices);
	if(N_periodic_vertices!= 0){
		removeDuplicates(periodicVerticesNodes);
		N_periodic_vertices = periodicVerticesNodes.size();
	}

	if(params.pmode == "moving"){
		removeMutualNodes(mNodes,N_m,periodicVerticesNodes);
		seNodes.conservativeResize(seNodes.size()+N_periodic_vertices);
		seNodes(Eigen::lastN(N_periodic_vertices)) = periodicVerticesNodes;
	}

	N_ss = ssNodes.size();
	N_se = seNodes.size();
	N_m = mNodes.size();
	N_pe = periodicEdgeNodes.size();
}




/* removeDuplicates function
 *
 * This function determines the unique elements in the array provided as argument.
 * This is done by initialising a temporary array that will save the unique elements.
 * The provided array is sorted and subsequent values are compared.
 * Since the for loop does not take care of final element, a separate comparison is done
 */


void Mesh::removeDuplicates(Eigen::ArrayXi& arr){
	// sorting the given array
	std::sort(std::begin(arr), std::end(arr));
	// Initialise array for the unique elements, could at most contain arr.size() unique elements
	Eigen::ArrayXi uniqueElems(arr.size());

	int cnt = 0;											// counter for the amount of unique elements
	for(int x = 0; x < arr.size()-1; x++){					// for loop over all but last element
		if(arr(x+1)!=arr(x)){								// if x+1-th is not equal to x-th element
			uniqueElems(cnt) = arr(x);						// include x-th element
			cnt++;											// update count of unique elements
		}
	}

	if(arr(arr.size()-1) != uniqueElems(cnt-1)){			// check if last element is already included
		uniqueElems(cnt) = arr(arr.size()-1);				// if not then include that element
	}

	// cut off all elements that were unassigned in the initialised array.
	arr = uniqueElems(Eigen::seq(0,cnt));
}

/* getIntNodes function
 *
 * This function takes the previously found boundary nodes and compares them to an integer array containing
 * all nodes from 0 to nNodes. If a node is not among the boundary nodes then the node is included in
 * the intNodes.
 */

void Mesh::getIntNodes(){
	if(lvl>=1){
		std::cout << "Determining internal nodes " << std::endl;
	}

	Eigen::ArrayXi bdryNodes(N_m+N_se+N_ss);
	bdryNodes << mNodes,seNodes,ssNodes;

	std::sort(std::begin(bdryNodes), std::end(bdryNodes));		// bdryNodes need to be sorted for this algorithm to work
	iNodes.resize(nNodes-N_m-N_se-N_ss);						// Resizing iNodes accordingly

	int cnt = 0, i=0, j=0;					// cnt keeps track of index of intNodes, i loops over all nodes, j over the boundary nodes


	while(i < nNodes){							// while loop over all nodes

		if(j < (N_m+N_se+N_ss)){				// check if j is within the size of the boundary nodes

			if(i == bdryNodes(j)){				// if node i is equal to the j-th boundary elements then i is not a internal node
				i++; j++;						// going to next node i and element j of the boundary nodes
			}
			else if(i < bdryNodes(j)){			// if i is smaller than the j-th element of the boundary nodes
				iNodes(cnt) = i;				// then i is an internal node
				i++; cnt++;						// update index of internal nodes and i

			}
			else if(i > bdryNodes(j)){			// if i is larger than the j-th boundary element
				j++;							// then next boundary element should be checked.
			}
		}
		else{									// if all j boundary elements are checked and i < nNodes
			iNodes(cnt) = i;					// include remaining points as internal nodes
			cnt++; i++;							// update index and i
		}

	}

	N_i = iNodes.size();
}

/*charLength function
 *
 * This function determines the maximum and minimal values in the coordinates.
 * Based on this the characteristic length of the domain is obtained. *
 */

double Mesh::getCharDomLength(){
	double charLength = 0;
	for(int i=0; i<nDims;i++){
		double length =  coords.col(i).maxCoeff()-coords.col(i).minCoeff();
		if (length > charLength){
			charLength = length;
		}
	}
	return charLength;
}


/* writeMeshFile function
 *
 * The writeMeshFile function writes the outputted mesh file containing the updated coordinates. This is done by reading the initial mesh file again
 * and copying its contents regarding the element connectivity and boundaries to the new mesh file. The part containing the coordinates is replaced by
 * the newly found coordinates.
 */

void Mesh::writeMeshFile(std::string& ifName, std::string& ofName){
	std::cout << "Writing output file " << std::endl;

	std::ofstream outputF;		// Making an output stream class to operate on the output file
	outputF.precision(15);		// sets precision of the floats in the file

	// opening existing or creating new output file. In its respective folder.
	outputF.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes\\" + ofName, std::ios::out); // ios::out allows output to file

	// Reopening the initial mesh file
	std::ifstream inputF("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes\\" + ifName);
	// string containing the contents of each line
	std::string line;
	// boolean that will be set to true whenever the new coordinates have to be specified.
	bool printFlag = false;
	// counter that keeps track of number of coordinates that have been outputted in the new file.
	int cnt = 0;

	// obtaining each line in a while loop
	while (getline(inputF, line)){
		// if the coordinates have to be outputted and the total number of nodes has not been reached.
		if(printFlag && cnt < nNodes){
			// set the coordinates based on the number of dimensions of the mesh.
			if(nDims == 2){
				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t'<< cnt << std::endl;
//				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t'<< cnt << std::endl;
			}
			else if(nDims == 3){
				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t' << coords(cnt,2) << '\t' << cnt << std::endl;
//				outputF << coords_polar_cylindrical(cnt,0)<< '\t' << coords_polar_spherical(cnt,1) << '\t' << coords_polar_spherical(cnt,2) << '\t' << cnt << std::endl;
			}
			cnt++;

		} else outputF << line << std::endl;	// if line does not correspond to any coordinates then simply copy that line to the new file.

		// if statement for finding the point in the file where the new coordinates have to be specified.
		if (line.rfind("NPOIN= ",0)==0){
			printFlag = true;
		}
	}

	// closing both files
	inputF.close();
	outputF.close();
	std::cout << "Done writing mesh file: " << ofName << std::endl;
}


/* getEdgeConnectivity function
 *
 * This function obtains for each sliding edge node, the other two nodes that make up the line segments eminating from the sliding node.
 * Each sliding node is part of two boundary elements. If one of those elements has an static nodes then that node is one of the connectivity nodes.
 * Otherwise, the sliding edge node will be connected to two other sliding edge nodes.
 *
 * For each sliding edge nodes it corresponding element is found by iterating through each column of the external boundary nodes array.
 * For the resulting row that corresponds to an element each node is first checked whether it is an external static node.
 * If not then they are checked to find which other node is a sliding edge node.
 */

void Mesh::getEdgeConnectivity(std::string& pmode, Eigen::ArrayXXi& edgeConnectivity, Eigen::ArrayXi& seNodes, int size){
	if(lvl>=1){
		std::cout << "Obtaining edge connectivity" << std::endl;
	}


//	int size = N_se - N_periodic_vertices;

	// resizing the array containing the connectivity information
	edgeConnectivity.resize(size,2);

	// row and col are the respective row and column index in the boundary node array.
	// idx is index in case either a static or sliding node is found.
	// cnt is a counter for the number of found connected nodes
	int row, col, idx,  cnt;

	// for-loop going through each sliding node
	for(int node = 0; node < size; node++){
		if(std::find(std::begin(verticesNodes), std::end(verticesNodes), seNodes(node)) == std::end(verticesNodes)){
	//		std::cout << node << '\t' << seNodes(node) << std::endl;

			col = 1;	// start from second column since first contains the element type
			cnt = 0;	// set count to zero, will go to 2 when both nodes connected by a line segment are found


	//		While all columns have not been searched and while 2 nodes are not found
			while(col<bdryNodesMat.cols() && cnt < 2){

				// find the row that contains the sliding edge node
				row = std::distance(std::begin(bdryNodesMat.col(col)), std::find(std::begin(bdryNodesMat.col(col)), std::end(bdryNodesMat.col(col)), seNodes(node)));

				// if the sliding edge node is found than the row index is smaller than the number of rows of the entire array.
				if (row!=bdryNodesMat.rows()){
	//				std::cout << "row: " << bdryNodesMat.row(row) << std::endl;
					// setting the index equal to the number of moving nodes, since this array includes the static nodes (zero movement).
					idx = verticesNodes.size();

					// starting from the second column, as first contains element type information
					int j = 1;

					// going through all columns until the idx changes or all columns have been done.
					while(j < bdryNodesMat.cols() && idx == verticesNodes.size()){
	//					std::cout << j << " / " << bdryNodesMat.cols() << std::endl;;
						idx = std::distance(std::begin(verticesNodes), std::find(std::begin(verticesNodes), std::end(verticesNodes),bdryNodesMat(row,j)));
						j++;
					}
					// if the index is not equal to the size of mNodes than a static node has been found.
					// this node should be included if either the count is zero or that node is not equal to the previous found node.
	//				if(idx!= mNodes.size() && ( cnt == 0 || mNodes(idx) != edgeConnectivity(node,cnt-1) ) ){
					if(idx!= verticesNodes.size()){
						if(cnt == 0 || verticesNodes(idx) != edgeConnectivity(node,cnt-1)){
							edgeConnectivity(node,cnt) = verticesNodes(idx);
//							std::cout << "found static node: "<< verticesNodes(idx) << std::endl;
							cnt++;
						}
					}
					// if no static node was found than, the function iterates through the columns to find the sliding nodes
					else{
						//starting again from the second column
						j=1;
						// setting index equal to the number of sliding edge nodes
						idx = seNodes.size();

						// going though all columns untill all columns are searched or the idx changes.
						while(j < bdryNodesMat.cols() && idx == seNodes.size()){
							// find index in case that row contains a sliding edge node
							idx = std::distance(std::begin(seNodes), std::find(std::begin(seNodes), std::end(seNodes),bdryNodesMat(row,j)));
							// if the found sliding edge node is the node that is considered then the idx is set back to the number of sliding edge nodes

							if(idx != seNodes.size() && seNodes(idx) == seNodes(node)){
								idx = seNodes.size();
							}
							j++;
						}

						// if an index is found and the count is either zero or that node has not been included than the node should be included.
						if(idx != seNodes.size() && ( cnt == 0 || seNodes(idx) != edgeConnectivity(node,cnt-1) ) ){
							edgeConnectivity(node,cnt) = seNodes(idx);
	//						std::cout << "found  sliding node: "<< seNodes(idx) << std::endl;
							cnt++;
						}
					}
				}
			// Moving to the next column
			col++;
			}
		}
	}
//	 debug message stating each sliding edge node and its found connecting nodes.
//	if(lvl>=3){
//		for(int i=0; i < edgeConnectivity.rows(); i++){
//			std::cout << seNodes(i) << '\t'  << edgeConnectivity.row(i)<<std::endl;
//		}
//	}
}

/* getSurfConnectivity function
 *
 * This function aims at finding the boundary elements that include the sliding surface nodes.
 *
 * These are found by finding the sliding surface node in each column of extBdryNodesMat, the corresponding row consists of the nodes making up the element.
 * The found indices that correspond to an element which contain the i-th sliding surface node are saved as the i-th row in the surface connectivity information array.
 */


void Mesh::getSurfConnectivity(){
//	if(lvl>=1){
//		std::cout << "Obtaining surface connectivity" << std::endl;
//	}

	// resizing of the surface connectivity information array. In case of an hexahedral mesh each sliding surface node is involved in four surfaces.
	surfConnectivity.resize(N_ss-N_pe, 4);

	// col is the variable used to iterate through the columns of the external boundary nodes array.
	// cnt keeps track of how many elements are found.
	// idx is index or row of the found element.
	int col,cnt,idx;

	// for-loop going through each sliding surface node.
	for(int i=0; i<N_ss-N_pe; i++){
		col = 1;	// First column contain information on the element type and is therefore skipped.
		cnt = 0;	// set count to zero.

		// iteratively going through the columns
		while(col < bdryNodesMat.cols()){

			// Finding the index of the row containing a sliding surface node.
			idx = std::distance(std::begin(bdryNodesMat.col(col)), std::find( std::begin(bdryNodesMat.col(col)), std::end(bdryNodesMat.col(col)), ssNodes(i)));

			// a sliding surface node is found when idx is not equal to the total amount of rows in extBdryNodesMat
			if(idx!=bdryNodesMat.rows()){

				// including the found index in the surface connectivity array and updating the count.
				surfConnectivity(i,cnt) = idx;
				cnt++;
			}
		col++;
		}
	}
	// debug message outputting the entire surfConnectivity array
//	if(lvl>=3){
//		std::cout << surfConnectivity << std::endl;
//	}
}

/* getVecs function
 *
 * The getVecs function finds the normal and tangential vectors required for the sliding algorithms.
 * In case of 2 dimensional problems only the tangential and normal vectors of the boundary line segments are found at the sliding edge nodes.
 *
 * In case of 3 dimensional problems the tangential vector is found for the sliding edge nodes on the boundary line segments and subsequently two perpendicular vectors are found.
 * For the surface sliding surface nodes the normal vector is found and subsequently two perpendicular vectors to it.
 *
 */

void Mesh::getVecs(probParams& params){
	if(lvl>=1){
//		std::cout << "Obtaining normal and tangential vectors " << std::endl;
	}

	// in case of 2D
	if(nDims == 2){
		// resizing of the normal and tangential vector arrays.
		n1_se.resize(edgeConnectivity.rows(), nDims); t_se.resize(edgeConnectivity.rows(), nDims);
		// Calling function to obtain the tangential vectors along the line segments at the sliding boundary node.
		getEdgeTan(t_se, edgeConnectivity, seNodes);


		// Finding normal of 2D problems.
		n1_se << t_se.col(1), -t_se.col(0);

		if(params.pmode == "moving"){
			n1_se.conservativeResize(N_se,nDims);
			t_se.conservativeResize(N_se,nDims);
			for(int i = N_se - periodicVerticesNodes.size(); i< N_se; i++){
				n1_se.row(i) = periodicVecs.col(1);
				t_se.row(i) = periodicVecs.col(0);
			}
		}
	}
	// else in case of 3D
	else if(nDims == 3){
		// resizing of the sliding edge normal and tangential vectors.
		t_se.resize(N_se,nDims); n1_se.resize(N_se,nDims); n2_se.resize(N_se,nDims);

		// calling function to obtain the tangential vectors along the boundary line segments at each sliding edge node.
		getEdgeTan(t_se, edgeConnectivity, seNodes);

		// set the type of element to edge and obtaining 2 vectors perpendicular to the tangential vectors.
		getPerpVecs(t_se, n1_se,n2_se);

		// resizing of the sliding surface normal and tangential vectors.
		n_ss.resize(N_ss-N_pe,nDims); t1_ss.resize(N_ss-N_pe,nDims); t2_ss.resize(N_ss-N_pe,nDims);

		// calling function to get the normal vectors of the surface sliding nodes.
		getSurfNormal();
		// set the type of element to surface and calling a function to obtain the two vector perpendicular to the normal vectors.
		getPerpVecs(n_ss, t1_ss,t2_ss);


		if(N_pe != 0){
			getSurfNormalPeriodic();
		}

		if(params.pmode == "moving"){
			n1_se.conservativeResize(N_se,nDims);
			n2_se.conservativeResize(N_se, nDims);
			t_se.conservativeResize(N_se, nDims);
			for(int i = N_se - periodicVerticesNodes.size(); i< N_se; i++){
				n1_se.row(i) = periodicVecs.col(1);
				n2_se.row(i) = periodicVecs.col(2);
				t_se.row(i) = periodicVecs.col(0);
			}
		}
	}

}

void Mesh::getSurfNormalPeriodic(){
	Eigen::ArrayXXd t(N_pe,nDims);
	getEdgeTan(t, edgeConnectivityPeriodic, periodicEdgeNodes);

	Eigen::ArrayXd periodicVec(nDims);
	periodicVec << 0,1,0; //todo remove this hardcoded periodic vector

	Eigen::ArrayXXd tp(N_pe,nDims), n(N_pe, nDims);
	for(int i = 0; i < N_pe; i++ ){
		tp.row(i) = periodicVec;
		t.row(i) =   t.row(i) - periodicVec.transpose() * t.row(i);
	}



	for(int i =0; i < nDims; i++){
		n.col(i) =  t.col((i+1)%3)*tp.col((i+2)%3) - t.col( (i+2)%3 )*tp.col((i+1)%3);
	}


	n_ss.conservativeResize(N_ss,nDims);
	t1_ss.conservativeResize(N_ss,nDims);
	t2_ss.conservativeResize(N_ss,nDims);
	n_ss(Eigen::lastN(N_pe),Eigen::all) = n;

	t1_ss(Eigen::lastN(N_pe), Eigen::all) = t;
	t2_ss(Eigen::lastN(N_pe), Eigen::all) = tp;
}

/* getEdgeTan function
 *
 * This function has as argument an empty array that will be filled with the tangential vectors.
 * The tangential vectors are found by using the edge connectivity information.
 * Two vectors are defined for each line segment connected to the sliding edge node.
 * The tangential vector at the node itself is found by taking a weighted average of the two vectors based on the length of the line segments.
 *
 * (this would result in the same vector if one would do it based on the distance from the sliding node to the midpoint of the line segment.
 * 		since the length of the vector is twice the distance to the midpoint.)
 */

void Mesh::getEdgeTan(Eigen::ArrayXXd& t, Eigen::ArrayXXi& edgeConnectivity, Eigen::ArrayXi& seNodes){

//	if(lvl>=2){
//		std::cout << "Obtaining edge tangential vectors" << std::endl;
//	}

	// initialising vectors to store intermediate results.
	Eigen::VectorXd v1(nDims), v2(nDims), tan(nDims);

	// looping through each sliding edge node
	for(int i = 0; i < edgeConnectivity.rows(); i++){
		// defining vectors along the line segments. Its important to note that these vectors need to have the same direction.
		// Otherwise, the vector will (partially) cancel each other out when taking the average if they are opposite to each other.
		// So one vector goes from connectivity node 1 to the sliding node and the second vector goes from the sliding node to connectivity node 2.

		v1 = (*ptrCoords).row(edgeConnectivity(i,0)) - (*ptrCoords).row(seNodes(i));
		v2 = (*ptrCoords).row(seNodes(i)) - (*ptrCoords).row(edgeConnectivity(i,1));

		tan = (v1/v1.norm() + v2/v2.norm());

		// Transforming the tangential vector in a unit vector and assigning it to its respective row in the tangential vector array.
		t.row(i) = tan/tan.norm();
	}

}

/* getSurfNormal function
 *
 * The getSurfNormal function is used to find a weighted average of the normal vector of the sliding surface node.
 *
 * For each surface that includes a sliding surface node the normal vector is found by defining two vectors to define its plane and then finding a vector normal to it.
 * The normal vector is found by the fact that the cross product of the two vectors is perpendicular to both vectors.
 * The resulting normal vectors are multiplied with the inverse relative distance from the sliding node to the midpoint of the element to obtain a weighted average.
 * Finally the vector is made of unit length and stored in the normal vector array.
 *
 */

// NOTE: FOR EFFICIENCY IT WILL BE BETER TO FIRST OBTAIN ALL NORMAL VECTORS FOR ALL BOUNDARY SURFACES AND AFTERWARDS DO THE WEIGHTED SUM.
// NOW THE NORMAL VECTOR OF A SINGLE ELEMENT IS CALCULATED MULTIPLE TIMES AS IT IS INVOLVED IN DIFFERENT SLIDING NODES.


void Mesh::getSurfNormal(){
	if(lvl>=2){
			std::cout << "Obtaining surface normal vectors" << std::endl;
		}

	// variables to store intermediate results
	Eigen::VectorXd n(nDims), vec1(nDims), vec2(nDims);
	// invDist is the inverse distance from the sliding surface node to the midpoint of the element.
	double invDist;

	// array dMidPnt will contain the relative distance from the sliding surface node to an element
	Eigen::ArrayXd dMidPnt(nDims);

	// looping through all sliding surface nodes
	for(int i = 0; i<N_ss-N_pe; i++){
		// initialise a zero normal vector
		n = Eigen::VectorXd::Zero(nDims);

		// for all j nodes corresponding to a boundary element index in the surface connectivity array.
		for(int j=0; j < surfConnectivity.cols(); j++){


			// two vector are set up to define a plane of the boundary element.
			// one vector from the first node to the second node and a second vector from the first node to the last node.
			vec1 = (*ptrCoords).row(bdryNodesMat(surfConnectivity(i,j),2)) - (*ptrCoords).row(bdryNodesMat(surfConnectivity(i,j),1));
			vec2 = (*ptrCoords).row(bdryNodesMat(surfConnectivity(i,j),4)) - (*ptrCoords).row(bdryNodesMat(surfConnectivity(i,j),1));

			// initialise a zero distance vector
			dMidPnt = Eigen::ArrayXd::Zero(nDims);

			// summing all coordinates of the nodes of the element and storing it in dMidPnt
			for(int l = 1; l< bdryNodesMat.cols(); l++){
				dMidPnt += (*ptrCoords).row(bdryNodesMat(surfConnectivity(i,j),l));
			}

			// taking the average and substracting the coordinates of the sliding surface node to obtain the relative distance
			// from the midpoint of the element to the sliding surface node
			dMidPnt = dMidPnt/(bdryNodesMat.cols()-1) - (*ptrCoords).row(ssNodes(i)).transpose();

			// calculating the inverse of the distance
			invDist = 1/dMidPnt.matrix().norm();

			// calculating the vector normal to the surface by using the fact that the cross product of vec1 and vec2 is perpendicular to both vec1 and vec2
			// the += operator is used since all surfaces that include the sliding surface node are summed.
			// The result is multiplied with the inverse of the distance in order to obtain a weighted average of the different elements.
			for(int k=0;k<nDims;k++){
				n(k) += (vec1((k+1)%3)*vec2((k+2)%3) - vec1((k+2)%3)*vec2((k+1)%3))*invDist;
			}
		}

		// Dividing the normal vector with its length to obtain a unit vector and assigning it in the normal vector array.
		// Usually one would divide by the inverse of the summed relative midpoint distances,
		// however this would simply scale the vector and since it is made of unit length anyway this operation can be left out.
		n_ss.row(i) =  n/n.norm();
	}
}

/* getPerpVecs function
 *
 * The getPerpVecs function aims at finding the vectors perpendicular to the one provided.
 * in 2D this involves swapping the x and y coordinate and introducing one minus sign. I.e. if v = [x, y] then its normal would be n = [y,-x]
 *
 * in 3D its slightly more complicating. A single vector is provided and two perpendicular vectors have to be found.
 * The first is found by using the fact that the dot product of perpendicular vectors is zero.
 * a second vector is found by the fact that the cross product of the first vector and the provided vector is perpendicular to both vectors.
 *
 *  A distinction is made for the type of element that is considered such that the correct vectors are calculated.
 */

void Mesh::getPerpVecs(Eigen::ArrayXXd& vecs, Eigen::ArrayXXd& p1, Eigen::ArrayXXd& p2){




		// Dot product for finding the first vector.
		for(int i = 0; i < nDims; i++){
			p1.col(i) = vecs.col((i+1)%3) - vecs.col((i+2)%3);
		}

		// Cross product for finding the second vector.
		for(int i =0; i < nDims; i++){
			p2.col(i) =  vecs.col((i+1)%3)*p1.col((i+2)%3) - vecs.col( (i+2)%3 )*p1.col((i+1)%3);
		}

		// Dividing the vectors by its length to ensure they are of unit length.
		for(int x = 0; x<p1.rows();x++){
			p1.row(x) = p1.row(x)/ p1.row(x).matrix().norm();
			p2.row(x) = p2.row(x)/ p2.row(x).matrix().norm();
		}

		// in case of an edge type element, the two tangentials to an normal vector should be found.
//		else if(type == "surface"){
//
//			// Dot product for finding the first vector
//			for(int i = 0; i < nDims; i++){
//				t1_ss.col(i) = n_ss.col((i+1)%3) - n_ss.col((i+2)%3);
//			}
//
//			// Cross product for finding the second vector
//			for(int i =0; i < nDims; i++){
//				t2_ss.col(i) =  n_ss.col((i+1)%3)*t1_ss.col((i+2)%3) - n_ss.col( (i+2)%3 )*t1_ss.col((i+1)%3);
//			}
//
//			// Dividing vectors by its length to ensure they are of unit length.
//			for(int x = 0; x<t1_ss.rows();x++){
//				t1_ss.row(x) = t1_ss.row(x)/ t1_ss.row(x).matrix().norm();
//				t2_ss.row(x) = t2_ss.row(x)/ t2_ss.row(x).matrix().norm();
//			}
//		}
//	}
}



/* getExtBdryEdgeSegment function
 *
 * This function takes an array with the external boundary edge nodes in the order as specified in the mesh file and
 * creates an array (extBdryEdgeSegments) where the rows will contain the nodes making up all the boundary line segments.
 * This array can then be further used to find the midpoints of these line segments and its normals.
 * Therefore, the midPnts and midPntNormals are already resized accordingly in this function
 */

void Mesh::getExtBdryEdgeSegments(){
	extBdryEdgeSegments.resize(N_se*2, 2);

	int idx,cnt = 0;
	Eigen::ArrayXi inclIdx;

	//todo check if this will still work for the new seNodes definitions

	// loop through sliding edge nodes
	for(int i = 0; i<N_se-N_periodic_vertices;i++){
		// nodes already included in the second column
		inclIdx = extBdryEdgeSegments(Eigen::seqN(0,cnt),1);

		//check if node seNodes(i) is already in the second column
		idx = std::distance(std::begin(inclIdx),std::find( std::begin(inclIdx) ,std::end(inclIdx),seNodes(i)));

		// if it is already in the second column
		if( idx != cnt){
			// if the value in the first column is equal to the first value of the edge connectivity then the second value should be added
			if(extBdryEdgeSegments(idx,0) == edgeConnectivity(i,0)){

				extBdryEdgeSegments.row(cnt) << seNodes(i), edgeConnectivity(i,1);

			}
			// else the first value should be added
			else{
				extBdryEdgeSegments.row(cnt) << seNodes(i), edgeConnectivity(i,0);
			}
			cnt++;
		}else{
			extBdryEdgeSegments(Eigen::seqN(cnt,2),Eigen::all)  << seNodes(i), edgeConnectivity(i,0) , seNodes(i), edgeConnectivity(i,1);
			cnt+= 2;
		}

	}
	extBdryEdgeSegments.conservativeResize(cnt,2);
}

/* getMidPnts function
 *
 * This function obtains the midpoints of the external boundary edge segments by using the array
 * containing the nodes making up the line segments extBdryEdgeSegments.
 * This array is also used to find the tangential of the line segment,
 * which allows for obtaining the normal of the line segments.
 *
 */

void Mesh::getMidPnts(probParams& params){

	int nrElems = 0, size = 0;

	Eigen::ArrayXi sTag(nrElemsBdry.size());
	for (int i = 0; i<nrElemsBdry.size();i++){
		if(std::find(std::begin(params.mTags),std::end(params.mTags),srtdTags[i]) == std::end(params.mTags) && (std::find(std::begin(params.pTags),std::end(params.pTags),srtdTags[i]) == std::end(params.pTags) || params.pmode == "none" || params.pmode == "periodic")){
			nrElems += nrElemsBdry(i);
			sTag(i) = 1;
		}else{
			sTag(i) = 0;
		}
	}

	if(nDims == 2){
		edgeMidPnts.resize(nrElems,nDims);
		edgeMidPntNormals1.resize(nrElems,nDims);
	}else if (nDims==3){
		surfMidPnts.resize(nrElems,nDims);
		surfMidPntNormals.resize(nrElems,nDims);

		size = extBdryEdgeSegments.rows();
		edgeMidPnts.resize(size,nDims);
		edgeMidPntNormals1.resize(size,nDims);
		edgeMidPntNormals2.resize(size,nDims);
	}


	int startIdx = 0;
	Eigen::ArrayXi indices(nrElems);
	for(int i= 0; i<sTag.size(); i++){
		if(sTag(i)){
			indices(Eigen::seqN(startIdx, nrElemsBdry(i))) = Eigen::ArrayXi::LinSpaced(nrElemsBdry(i),nrElemsBdry(Eigen::seqN(0,i)).sum(),nrElemsBdry(Eigen::seqN(0,i)).sum()+nrElemsBdry(i));
			startIdx += nrElemsBdry(i);
		}
	}
//	std::cout << indices << std::endl;
	Eigen::VectorXd tan(nDims);

	if(nDims == 2){
		for(int i = 0; i < nrElems; i++){

			edgeMidPnts.row(i) = ((*ptrCoords).row(bdryNodesMat(indices(i),1)) + (*ptrCoords).row(bdryNodesMat(indices(i),2)))/2;

			tan = (*ptrCoords).row(bdryNodesMat(indices(i),2)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));

			edgeMidPntNormals1.row(i) << tan(1),-tan(0);

			edgeMidPntNormals1.row(i) = edgeMidPntNormals1.row(i)/tan.norm();
		}
	}else{

		Eigen::ArrayXXd midPntTan(size,nDims);
		for(int i = 0 ; i < size; i++ ){
			edgeMidPnts.row(i) = ((*ptrCoords).row(extBdryEdgeSegments(i,0)) + (*ptrCoords).row(extBdryEdgeSegments(i,1)))/2;
			midPntTan.row(i) = (*ptrCoords).row(extBdryEdgeSegments(i,1)) - (*ptrCoords).row(extBdryEdgeSegments(i,0));
		}

		getPerpVecs(midPntTan,edgeMidPntNormals1,edgeMidPntNormals2);

		Eigen::VectorXd n(nDims), vec1(nDims), vec2(nDims);

		for(int i = 0; i < nrElems; i++){

			n = Eigen::VectorXd::Zero(nDims);
			// todo what if bdryNodesMat has less than 5 columns/?
			surfMidPnts.row(i) = ((*ptrCoords).row(bdryNodesMat(indices(i),1)) + (*ptrCoords).row(bdryNodesMat(indices(i),2)) + (*ptrCoords).row(bdryNodesMat(indices(i),3))+ (*ptrCoords).row(bdryNodesMat(indices(i),4)))/4;
			vec1 = (*ptrCoords).row(bdryNodesMat(indices(i),2)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));
			vec2 = (*ptrCoords).row(bdryNodesMat(indices(i),4)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));
//			std::cout << vec1 << std::endl;
//			std::cout << vec2 << std::endl;
			for(int k=0;k<nDims;k++){
				n(k) += (vec1((k+1)%3)*vec2((k+2)%3) - vec1((k+2)%3)*vec2((k+1)%3));
			}
//			std::cout << n << std::endl;
			surfMidPntNormals.row(i) = n/n.norm();

		}
	}
//	std::cout << midPnts << std::endl;
//	std::cout << midPntNormals << std::endl;


}





void Mesh::getCharPerLength(probParams& params){

	Eigen::ArrayXXd vertices(verticesNodes.size(), nDims);
	if(params.ptype){
		Eigen::ArrayXXd verticesCart(verticesNodes.size(),nDims);
		verticesCart = coords(verticesNodes,Eigen::all);

		CoordTransform transform;
		transform.cart_to_polar_cylindrical(verticesCart, vertices);
	}else{
		vertices = coords(verticesNodes,Eigen::all);
	}


	std::vector<double>  minVals;
	std::vector<int>  cols;
	for(int i = 0; i < nDims; i++){
		if(i != params.pDir){
			cols.push_back(i);
			minVals.push_back(vertices.col(i).minCoeff());
		}
	}




	Eigen::ArrayXd periodicVals(vertices.rows());
	int cnt = 0;

	for(int i = 0; i < vertices.rows(); i++){
		if(vertices(i,cols[0]) == minVals[0] && (nDims == 2 || vertices(i,cols[1]) == minVals[1])){

			periodicVals(cnt) = vertices(i,params.pDir);
			cnt++;
		}
	}



	periodicVals.conservativeResize(cnt);

	periodic_length = periodicVals.maxCoeff()-periodicVals.minCoeff();

}

void Mesh::findStringBounds(int& first, int& last, std::string& line){
	first = line.find("=")+1;
	last = line.size();

	while(line[first] == ' '){
		first++;
	}

	while(line[last-1] == ' '){
		last = last -1;
	}

	if(first == last){
		last++;
	}
}


void Mesh::getPeriodicParams(probParams& params){
	// each column of periodicVecs contains the periodic vectors for either the cartesian/ polar/ spherical coordinates
	// the first column is the normal vector in periodic direction
	// second and third column contains the two vectors normal to the periodic direction
	periodicVecs = Eigen::MatrixXd::Zero(nDims,nDims);

	int cnt = 1;
	for(int row = 0; row < nDims; row++){
		if(row == params.pDir){
			periodicVecs(row,0) = 1;
		}else{
			periodicVecs(row, cnt) = 1;
			cnt++;
		}
	}
}
