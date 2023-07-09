#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
//#include <math.h>
#include <Eigen/Dense>



Mesh::Mesh(probParams& params)
{
	// reading the input mesh file
	readMeshFile(params);

	// establish the nodetype for each of the boundary nodes
	getNodeTypes(params);

	// for sliding determine the connectivity of the boundary elements required for finding the local normal and tangential vector
	if(params.smode != "none"){
		// establish edge connectivity
		getEdgeConnectivity(params.pmode, edgeConnectivity,seNodes,  N_se - N_periodic_vertices);

		if(nDims == 3){
			// establish the 3d edge connectivity
			getExtBdryEdgeSegments();
			//establish surface connectivity
			getSurfConnectivity();
			// in case of periodic displacement, determine connectivity of the periodic edge nodes
			if(params.pmode == "fixed" || "moving"){
				getEdgeConnectivity(params.pmode, edgeConnectivityPeriodic, periodicEdgeNodes, N_pe);
			}
		}
	}

	// determing the support radius of the RBF, factor times the max domain length
	r = params.rFac*getCharDomLength();

	// determining the periodic vectors
	getPeriodicVector(params);


	// set pointer to either Cartesian or polar/cylindrical coordinates
	if(params.ptype){
		// resize the polar/ cylindrical coordinates
		coords_polar_cylindrical.resize(nNodes,nDims);
		ptrCoords = &coords_polar_cylindrical;
	}else{
		ptrCoords = &coords;
	}

	// determining the periodic length of the domain
	if(params.pmode != "none"){
		getCharPerLength(params);
	}


}

// The readMeshFile function reads the input mesh file and stores the required information
void Mesh::readMeshFile(probParams& params){


	int lineNo = 0;								// line number counter
	int nBdryElemsTotal = 0;					// stores total number of boundary elements

	int markerIdx = -1;							// contain line index containing the latest marker found

	int nodeCnt = 0;							// counter for number of points
	int pntsIdx = 0;								// int that stores the line where "NPOIN= " is

	int bdryElemCnt = 0; 						// counting the elements of the boundaries
	int markerElems = 0;						// locally stores how many elements are in that boundary
	int markerCnt = 0; 							// Counts the number of external boundary markers
	nNodes = 0; 								// Number of nodes

	int firstIdx, lastIdx;						// used to remove whitespace from input variables
	bdryNodesMat.resize(0,3);					// thet boundary node array saves information on the boundary elements. First column contains element type,
												// remaining columns contain nodes of the element


	srtdTags.resize(params.bdryTags.size());		// Tags sorted by appearance in the input mesh file
	nrElemsBdry.resize(params.bdryTags.size());		// Array containing the sizes of each boundary

	std::string line;														// string containing line obtained by getline() function
	std::ifstream mFile(params.directory + "\\" +  params.mesh_ifName); 	//opening file name and storing in mFile object

	// Check if file is opened
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){

			// Save number of dimensions
			if (line.rfind("NDIME= ",0)==0){
				findStringBounds(firstIdx, lastIdx,line);
				nDims = stoi(line.substr(firstIdx, lastIdx-firstIdx));
			}
			// Save number of nodes, index where the coordinates start and resize the coordinates array
			else if (line.rfind("NPOIN= ",0)==0){
				findStringBounds(firstIdx, lastIdx,line);
				nNodes = stoi(line.substr(firstIdx, lastIdx-firstIdx));
				pntsIdx = lineNo;
				coords.resize(nNodes, nDims);
			}

			// Checking whether provided boundary tags equal the amount in the mesh file
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
				findStringBounds(firstIdx, lastIdx,line);
				std::string tag =  line.substr(firstIdx, lastIdx-firstIdx);

				// if found tag is among those provided by user then store them in the order as they are found. Else throw error.
				try{
					if(std::find(std::begin(params.bdryTags), std::end(params.bdryTags), tag) != std::end(params.bdryTags)){
						// saving the index where the tag is found
						markerIdx = lineNo;
						// saving tags in the order found in the mesh file
						srtdTags[markerCnt] = tag;
					}
					else throw(tag);
				}
				// If the provided tag does not match any tag in the meshing file an error is thrown
				catch(std::string& tag){
					std::cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << std::endl;
					std::exit(0);
				}
			}

			// Line after the tag contains the number of boundary elements
			else if(lineNo == markerIdx+1){

				// store number of elements of the boundary
				findStringBounds(firstIdx,lastIdx,line);
				markerElems = stoi(line.substr(firstIdx, lastIdx-firstIdx));

				// Updating total number of boundary elements
				nBdryElemsTotal+= markerElems;

				// Saving number of elements of each boundary in an array
				nrElemsBdry(markerCnt) = markerElems;

				// resizing the array containing the boundary node data
				bdryNodesMat.conservativeResize(nBdryElemsTotal,bdryNodesMat.cols());

				//updating the marker count
				markerCnt++;

				}

			// Lines describing boundary elements
			else if(lineNo > markerIdx+1 && lineNo <= markerIdx+1+ markerElems){

				// split line by '\t' character
				std::istringstream is(line);

				// elements on the line are assigned to data in a while loop. lineElem counts the number of elements per line
				int data, elemCnt = 0;
				while(is >> data){
					// First integer on line describes the type of element according to VTK format. See SU2 documentation for more details on the mesh file

					// resizing of the bdryNodesMat array to accomodate all the boundary nodes if required
					if(elemCnt == 0 && data == 5 && bdryNodesMat.cols() < 4){
						bdryNodesMat.conservativeResize(nBdryElemsTotal,4);
					}
					else if(elemCnt == 0 && data == 9 && bdryNodesMat.cols() < 5){
						bdryNodesMat.conservativeResize(nBdryElemsTotal,5);
					}

					// assigning the node information to the array
					bdryNodesMat(bdryElemCnt,elemCnt) = data;
					// updating the line element count
					elemCnt++;
				}
				// updating the count of boundary elements
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
				// update node count
				nodeCnt++;
			}
			// update line number
			lineNo++;

		}
		//closing meshing file
		mFile.close();
	}
	// If the file is not opened then the following error message will be displayed
	else{ std::cout << "Not able to open input mesh file\n";
		std::exit(0);
	}
}



/*
 * removeMutualNodes function
 *
 * This function removes any nodes in array_in that are shared with the array of to_remove_nodes
 */

void Mesh::removeMutualNodes(Eigen::ArrayXi& array_in, int& size, Eigen::ArrayXi& to_remove_nodes, int end_idx){
	Eigen::ArrayXi idx_keep(array_in.size());
	int cnt = 0;
	for(int i = 0; i < array_in.size(); i++ ){
		if(std::find(std::begin(to_remove_nodes(Eigen::seqN(0,end_idx))), std::end(to_remove_nodes(Eigen::seqN(0,end_idx))), array_in(i)) == std::end(to_remove_nodes(Eigen::seqN(0,end_idx)))){
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
 * This function will go through each boundary specified in the mesh file individually and establish based on the problem parameters
 * and occurences of the node what type of node it is.
 * Based on the parameters it will then assign the various node types
 * Lastly, the internal nodes are found
 */

void Mesh::getNodeTypes(probParams& params){

	// array containing all nodes of the boundary
	Eigen::ArrayXi bdryNodesArr;

	// arrays containing the indices of different types of nodes
	Eigen::ArrayXi idxMoving, idxEdge, idxSurf;

	// counters for the node types
	int cntMoving, cntEdge, cntSurf;

	// number of nodes with nonzero displacement
	N_nonzeroDisp = 0;

	// number of periodic vertices in the domain
	N_periodic_vertices = 0;

	// number of periodic edge ndoes
	N_pe = 0;

	// array with nodes corresponding to vertices
	verticesNodes.resize(10);

	// array with nodes that are periodic vertices
	periodicVerticesNodes.resize(10);

	// vertices counter
	int verticesCnt = 0;

	// booleans set when a boundary is periodic moving or internal
	bool periodic, moving, internal;


	// for loop going through the boundaries
	for(int elem = 0; elem < int(params.bdryTags.size()); elem++){

		// resizing the array that will contain all the boundary nodes of that boundary
		bdryNodesArr.resize(nrElemsBdry(elem)*(bdryNodesMat.cols()-1));

		// next two loops ensure that all nodes are included in the 1D node array.
		for(int i =0; i<nrElemsBdry(elem);i++){
			for(int j=0; j< bdryNodesMat.cols()-1; j++ ){
				bdryNodesArr(j+i*(bdryNodesMat.cols()-1)) = bdryNodesMat(i+ nrElemsBdry(Eigen::seqN(0,elem)).sum(),j+1);
			}
		}


		// sorting the array
		std::sort(std::begin(bdryNodesArr),std::end(bdryNodesArr));


		// Set size of the arrays of the various nodes equal to the size of the bdryNodesArr.
		if(idxMoving.size() != bdryNodesArr.size()){
			idxMoving.resize(bdryNodesArr.size());
			idxEdge.resize(bdryNodesArr.size());
			idxSurf.resize(bdryNodesArr.size());
		}

		// setting counters to zero
		cntMoving = 0, cntEdge = 0, cntSurf = 0;

		// setting all booleans to false
		moving = false;
		periodic = false;
		internal = false;

		// check if the boundary is a moving boundary
		if(std::find(std::begin(params.mTags),std::end(params.mTags),srtdTags[elem]) != std::end(params.mTags)){
//			std::cout << srtdTags[elem] << " is a moving boundary" << std::endl;
			moving = true;
		}
		// check if boundary is periodic, in case the RBF interpolation is periodic
		else if((params.pmode == "fixed" || params.pmode == "moving") && std::find(std::begin(params.pTags),std::end(params.pTags),srtdTags[elem]) != std::end(params.pTags)){
//			std::cout << srtdTags[elem] << " is periodic boundary" << std::endl;
			periodic = true;
		}
		// check if the boundary should be treated as internal nodes
		else if(std::find(std::begin(params.iTags),std::end(params.iTags),srtdTags[elem]) != std::end(params.iTags)){
//			std::cout << srtdTags[elem] << " is a boundary treated as internal nodes" << std::endl;
			internal = true;
		}


		// if moving boundary then save all nodes as moving nodes
		if(moving){
			idxMoving(Eigen::seqN(0,bdryNodesArr.size())) = bdryNodesArr;
			cntMoving += bdryNodesArr.size();

			// establish  the number of moving nodes that have a nonero displacement
			removeDuplicates(bdryNodesArr);
			N_nonzeroDisp += bdryNodesArr.size();
		}

		// if not moving node then its a sliding node
		else{
			// looping through the nodes
			for (int i= 0; i< bdryNodesArr.size();i++){

				// for 3D, in case of 4 or more subsuquent equal nodes its a surface node
				if(nDims == 3 && (i< bdryNodesArr.size()-4 && bdryNodesArr(i) == bdryNodesArr(i+3))){

					// check the amount of equal subsequent nodes
					int n = 4;
					while(i+n < bdryNodesArr.size() && bdryNodesArr(i) == bdryNodesArr(i+n)){
						n++;
					}

					// save the node in question
					idxSurf(cntSurf) = bdryNodesArr(i);

					//update surface node count
					cntSurf++;

					// skip subsequent equal nodes
					i +=(n-1);
				}

				// else in case of 2 or 3 subsequent equal nodes then its an edge node
				else if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){

					// check the number of equal subsequent nodes
					int n = 2;
					while(i+n < bdryNodesArr.size() && bdryNodesArr(i) == bdryNodesArr(i+n)){
						n++;
					}

					// save as sliding edge node
					idxEdge(cntEdge) = bdryNodesArr(i);

					// update the count
					cntEdge++;

					// skip equal nodes in the loop
					i +=(n-1);
				}
				// else in case there is a single occurence of the node in that bdry, its a vertex of the domain
				else{

					// save as moving node
					idxMoving(cntMoving) = bdryNodesArr(i);

					// update moving node count
					cntMoving++;

					// save as vertex and update count
					verticesNodes(verticesCnt) = bdryNodesArr(i);
					verticesCnt++;

					// if array size is reached then increase the size with 10
					if(verticesCnt == verticesNodes.size()){
						verticesNodes.conservativeResize(verticesNodes.size()+10);
					}

					// in case of periodic displacement with moving vertices, the vertex is also saved as periodic vertex
					if(params.pmode == "moving"){
						periodicVerticesNodes(N_periodic_vertices) = bdryNodesArr(i);
						N_periodic_vertices++;
						if(N_periodic_vertices == periodicVerticesNodes.size()){
							periodicVerticesNodes.conservativeResize(periodicVerticesNodes.size()+10);
						}
					}
				}
			}
		}


		// if the boundary should be considerd as internal nodes, then save edge nodes as internal edge nodes
		if(internal){
			internalEdgeNodes.conservativeResize(internalEdgeNodes.size() + cntEdge);
			internalEdgeNodes(Eigen::lastN(cntEdge)) = idxEdge(Eigen::seqN(0,cntEdge));

		// if periodic boundary then save edge nodes as periodic edge nodes
		}else if(periodic){
			periodicEdgeNodes.conservativeResize(periodicEdgeNodes.size()+cntEdge);
			periodicEdgeNodes(Eigen::lastN(cntEdge)) = idxEdge(Eigen::seqN(0,cntEdge));


		}else{
			// save moving nodes in the mNodes matrix
			mNodes.conservativeResize(mNodes.size()+cntMoving);
			mNodes(Eigen::lastN(cntMoving)) = idxMoving(Eigen::seqN(0,cntMoving));

			//  in case of no sliding then also the edge and surface nodes are saved as moving nodes
			if(params.smode == "none"){
				mNodes.conservativeResize(mNodes.size()+cntEdge);
				mNodes(Eigen::lastN(cntEdge)) = idxEdge(Eigen::seqN(0,cntEdge));

				if(cntSurf != 0){
					mNodes.conservativeResize(mNodes.size()+cntSurf);
					mNodes(Eigen::lastN(cntSurf)) = idxSurf(Eigen::seqN(0,cntSurf));
				}

			// saving sliding edge nodes
			}else{
				seNodes.conservativeResize(seNodes.size()+cntEdge);
				seNodes(Eigen::lastN(cntEdge)) = idxEdge(Eigen::seqN(0,cntEdge));

				// saving sliding surface nodes
				if(cntSurf != 0){
					ssNodes.conservativeResize(ssNodes.size()+cntSurf);
					ssNodes(Eigen::lastN(cntSurf)) = idxSurf(Eigen::seqN(0,cntSurf));
				}
			}
		}
	}

	// Calling a function that removes any duplicate elements from the arrays.
	if(ssNodes.size() != 0){
		removeDuplicates(ssNodes);
	}

	if(seNodes.size() != 0 ){
		removeDuplicates(seNodes);
	}

	removeDuplicates(mNodes);


	if(internalEdgeNodes.size() !=0){
		removeDuplicates(internalEdgeNodes);
	}

	verticesNodes.conservativeResize(verticesCnt);
	if(verticesCnt != 0){
		removeDuplicates(verticesNodes);
	}


	periodicVerticesNodes.conservativeResize(N_periodic_vertices);


	if(N_periodic_vertices!= 0){
		removeDuplicates(periodicVerticesNodes);
		N_periodic_vertices = periodicVerticesNodes.size();
	}

	if(periodicEdgeNodes.size()!=0){
		removeDuplicates(periodicEdgeNodes);
	}


	// establish the amount of the various node types
	N_ss = ssNodes.size();
	N_se = seNodes.size();
	N_m = mNodes.size();
	N_pe = periodicEdgeNodes.size();


	// in case of 3D and periodic displacement of the boundaries, the periodic edge nodes will be treated as sliding surface nodes
	if(nDims == 3 && (params.smode == "ds" || params.smode == "ps") && (params.pmode == "fixed" || params.pmode == "moving")){

		// adding periodic edge nodes to the sliding surface nodes
		N_ss += N_pe;
		ssNodes.conservativeResize(N_ss);
		ssNodes(Eigen::lastN(N_pe)) = periodicEdgeNodes;

		// removing the periodic edge nodes from the sliding edge node selection
		removeMutualNodes(seNodes, N_se, periodicEdgeNodes, N_pe);
	}

	// in case of periodic displacement with moving vertices, the periodic vertices are treated as sliding edge nodes
	if(params.pmode == "moving"){
		// remove periodic vertices from the moving nodes
		removeMutualNodes(mNodes,N_m,periodicVerticesNodes, N_periodic_vertices);

		// add periodic vertices to the sliding edge nodes
		N_se += N_periodic_vertices;
		seNodes.conservativeResize(N_se);
		seNodes(Eigen::lastN(N_periodic_vertices)) = periodicVerticesNodes;
	}


	// in case the considered object intersects with an external boundary,
	// then the duplicates of moving nodes have to be removed from the sliding edge and surface nodes
	removeMutualNodes(seNodes, N_se, mNodes, N_m);
	removeMutualNodes(ssNodes, N_ss, mNodes, N_m);

	// and the sliding edge nodes from the sliding surface ndoes
	removeMutualNodes(ssNodes, N_ss, seNodes, N_se);

	// call function to establish the internal nodes
	getIntNodes();

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
 * the newly computed coordinates.
 */

void Mesh::writeMeshFile(std::string& directory, std::string& ifName, std::string& ofName){

	std::ofstream outputF;		// Making an output stream class to operate on the output file
	outputF.precision(15);		// sets precision of the floats in the file

	// opening existing or creating new output file.
	outputF.open(directory + "\\" + ofName, std::ios::out); // ios::out allows output to file

	// Reopening the initial mesh file
	std::ifstream inputF(directory + "\\" + ifName);
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
			}
			else if(nDims == 3){
				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t' << coords(cnt,2) << '\t' << cnt << std::endl;
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
	std::cout << "Output mesh file: " << ofName << " has been generated." <<  std::endl;
}


/* getEdgeConnectivity function
 *
 * This function obtains for each sliding edge node, the other two nodes that make up the line segments eminating from the sliding node.
 * Each sliding node is part of two boundary elements. If one of those elements has a static node then that node is one of the connectivity nodes.
 * Otherwise, the sliding edge node will be connected to two other sliding edge nodes.
 *
 * For each sliding edge nodes it corresponding element is found by iterating through each column of the boundary nodes array.
 * For the resulting row that corresponds to an element each node is first checked whether it is an static node.
 * If not then they are checked to find which other node is a sliding edge node.
 */

void Mesh::getEdgeConnectivity(std::string& pmode, Eigen::ArrayXXi& edgeConnectivity, Eigen::ArrayXi& seNodes, int size){


	// resizing the array containing the connectivity information
	edgeConnectivity.resize(size,2);

	// row and col are the respective row and column index in the boundary node array.
	// idx is index in case either a static or sliding node is found.
	// cnt is a counter for the number of found connected nodes
	int row, col, idx,  cnt;

	// for-loop going through each sliding node
	for(int node = 0; node < size; node++){
		if(std::find(std::begin(verticesNodes), std::end(verticesNodes), seNodes(node)) == std::end(verticesNodes)){

			col = 1;	// start from second column since first contains the element type
			cnt = 0;	// set count to zero, will go to 2 when both nodes connected by a line segment are found

	//		While all columns have not been searched and while 2 nodes are not found
			while(col<bdryNodesMat.cols() && cnt < 2){

				// find the row that contains the sliding edge node
				row = std::distance(std::begin(bdryNodesMat.col(col)), std::find(std::begin(bdryNodesMat.col(col)), std::end(bdryNodesMat.col(col)), seNodes(node)));

				// if the sliding edge node is found than the row index is smaller than the number of rows of the entire array.
				if (row!=bdryNodesMat.rows()){

					// setting the index equal to the number of vertex nodes, since this array includes the static nodes (zero movement).
					idx = verticesNodes.size();

					// starting from the second column, as first contains element type information
					int j = 1;

					// going through all columns until the idx changes or all columns have been done.
					while(j < bdryNodesMat.cols() && idx == verticesNodes.size()){
						idx = std::distance(std::begin(verticesNodes), std::find(std::begin(verticesNodes), std::end(verticesNodes),bdryNodesMat(row,j)));
						j++;
					}

					// if the index is not equal to the number of vertices than a static node has been found.
					// this node should be included if either the count is zero or that node is not equal to the previous found node.
					if(idx!= verticesNodes.size()){
						if(cnt == 0 || verticesNodes(idx) != edgeConnectivity(node,cnt-1)){
							edgeConnectivity(node,cnt) = verticesNodes(idx);
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
							cnt++;
						}
					}
				}
			col++;
			}
		}
	}

}

/* getSurfConnectivity function
 *
 * This function aims at finding the boundary elements that include the sliding surface nodes.
 *
 * These are found by finding the sliding surface node in each column of bdryNodesMat, the corresponding row consists of the nodes making up the element.
 * The found indices that correspond to a surface element which contain the i-th sliding surface node are saved as the i-th row in the surface connectivity information array.
 */


void Mesh::getSurfConnectivity(){

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

			// Finding the index of the row containing the sliding surface node.
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

	// in case of 2D
	if(nDims == 2){
		// resizing of the normal and tangential vector arrays.
		n1_se.resize(edgeConnectivity.rows(), nDims); t_se.resize(edgeConnectivity.rows(), nDims);
		// Calling function to obtain the tangential vectors along the line segments at the sliding boundary node.
		getEdgeTan(t_se, edgeConnectivity, seNodes);

		// Finding normal of 2D problems. Eq. 2.27 from the manuscript
		n1_se << t_se.col(1), -t_se.col(0);

		// in case of periodic displacement with moving vertices, the tangential vector is the one tangential to the periodic direction
		// and the normal is a vector normal to the periodic direction
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


		// In case of periodic displacement the periodic edge nodes are treated as sliding surface nodes.
		// therefore, their vectors are added to the sliding surface vector arrays in case of direct sliding
		if(N_pe != 0){
			int ds = 0;
			if(params.smode == "ds"){
				ds = 1; 	 				// set to 1 in order to include the vectors to the sliding surface vectors
			}

			// determining normal and tangential vectors of the periodic edge nodes
			getSurfNormalPeriodic(ds);
		}

		// in case of periodic displacement with moving vertices the normal and tangential vectors of the vertices are included in the sliding edge nodes
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


/* getSurfNormalPeriodic
 *
 * Obtains the tangential vectors and normal vector of the periodic edge nodes, which are treated as sliding surface nodes
 */

void Mesh::getSurfNormalPeriodic(int directSliding){

	// resizing the tangential vector array
	Eigen::ArrayXXd t(N_pe,nDims);

	// obtaining tangential vectors
	getEdgeTan(t, edgeConnectivityPeriodic, periodicEdgeNodes);

	// vector in the periodic direction
	Eigen::ArrayXd periodicVec(nDims);
	periodicVec = periodicVecs.col(0);


	Eigen::ArrayXXd tp(N_pe,nDims);

	// resizing of the periodic edge normals
	if(periodicEdgeNormals.rows() != N_pe){
		periodicEdgeNormals.resize(N_pe,nDims);
	}

	// First tangential vector is in the periodic direction
	// second is found by determing a vector perpendicular to it
	for(int i = 0; i < N_pe; i++ ){
		tp.row(i) = periodicVec;
		t.row(i) =   t.row(i) - periodicVec.transpose() * t.row(i);
	}

	// find the vector normal to the 2 tangential vectors
	for(int i =0; i < nDims; i++){
		periodicEdgeNormals.col(i) =  t.col((i+1)%3)*tp.col((i+2)%3) - t.col( (i+2)%3 )*tp.col((i+1)%3);
	}

	// normalising to obtain unit vector
	for(int row = 0; row < periodicEdgeNormals.rows(); row++){
		periodicEdgeNormals.row(row) = periodicEdgeNormals.row(row)/(periodicEdgeNormals.row(row)).matrix().norm();
	}


	// in case of direct sliding adding them to the sliding surface vector arrays
	if(directSliding){
		n_ss.conservativeResize(N_ss,nDims);
		t1_ss.conservativeResize(N_ss,nDims);
		t2_ss.conservativeResize(N_ss,nDims);
		n_ss(Eigen::lastN(N_pe),Eigen::all) = periodicEdgeNormals;

		t1_ss(Eigen::lastN(N_pe), Eigen::all) = t;
		t2_ss(Eigen::lastN(N_pe), Eigen::all) = tp;
	}
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


	// initialising vectors to store intermediate results.
	Eigen::VectorXd v1(nDims), v2(nDims), tan(nDims);

	// looping through each sliding edge node
	for(int i = 0; i < edgeConnectivity.rows(); i++){

		// defining vectors along the line segments. Its important to note that these vectors need to have the same direction.
		// Otherwise, the vector will (partially) cancel each other out when taking the average if they are opposite to each other.
		// So one vector goes from connectivity node 1 to the sliding node and the second vector goes from the sliding node to connectivity node 2.
//
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


void Mesh::getSurfNormal(){

	// variables to store intermediate results
	Eigen::VectorXd n(nDims), vec1(nDims), vec2(nDims);
	// invDist is the inverse distance from the sliding surface node to the midpoint of the element.
	double invDist = 0.;

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

}



/* getExtBdryEdgeSegment function
 *
 * This function takes an array with the external boundary edge nodes in the order as specified in the mesh file and
 * creates an array (extBdryEdgeSegments) where the rows will contain the nodes making up all the boundary line segments.
 * This array can then be further used to find the midpoints of these line segments and its normals.
 */

void Mesh::getExtBdryEdgeSegments(){
	extBdryEdgeSegments.resize(N_se*2, 2);

	int idx,cnt = 0;
	Eigen::ArrayXi inclIdx;


	// loop through sliding edge nodes
	for(int i = 0; i<N_se-N_periodic_vertices;i++){

		// nodes already included in the second column
		inclIdx = extBdryEdgeSegments(Eigen::seqN(0,cnt),1);

		//check if edge node is already in the second column
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

	// Determining the the total number of sliding elements considered and saving which boundaries are sliding boundaries
	Eigen::ArrayXi sTag(nrElemsBdry.size());
	for (int i = 0; i<nrElemsBdry.size();i++){
		if(std::find(std::begin(params.mTags),std::end(params.mTags),srtdTags[i]) == std::end(params.mTags) && (std::find(std::begin(params.pTags),std::end(params.pTags),srtdTags[i]) == std::end(params.pTags) || params.pmode == "none" || params.pmode == "periodic")){
			nrElems += nrElemsBdry(i);
			sTag(i) = 1;
		}else{
			sTag(i) = 0;
		}
	}

	// resizing of midpoint arrays and setting pointers
	if(nDims == 2){
		edgeMidPnts.resize(nrElems,nDims);
		if(params.ptype){
			edgeMidPnts_polar_cylindrical.resize(nrElems,nDims);
			edgeMidPntPtr = &edgeMidPnts_polar_cylindrical;
		}else{
			edgeMidPntPtr = &edgeMidPnts;
		}
		edgeMidPntNormals1.resize(nrElems,nDims);
	}else if (nDims==3){
		surfMidPnts.resize(nrElems,nDims);
		surfMidPntNormals.resize(nrElems,nDims);

		// setting the size of the edge midpoints
		size = extBdryEdgeSegments.rows();
		edgeMidPnts.resize(size,nDims);

		edgeMidPntNormals1.resize(size,nDims);
		edgeMidPntNormals2.resize(size,nDims);

		if(params.ptype){
			surfMidPnts_polar_cylindrical.resize(nrElems,nDims);
			surfMidPntPtr = &surfMidPnts_polar_cylindrical;

			edgeMidPnts_polar_cylindrical.resize(size,nDims);
			edgeMidPntPtr = &edgeMidPnts_polar_cylindrical;
		}else{
			edgeMidPntPtr = &edgeMidPnts;
			surfMidPntPtr = &surfMidPnts;
		}
	}


	// setting the indices of the elements for which the midpoints have to be determined
	int startIdx = 0;
	Eigen::ArrayXi indices(nrElems);
	for(int i= 0; i<sTag.size(); i++){
		if(sTag(i)){
			indices(Eigen::seqN(startIdx, nrElemsBdry(i))) = Eigen::ArrayXi::LinSpaced(nrElemsBdry(i),nrElemsBdry(Eigen::seqN(0,i)).sum(),nrElemsBdry(Eigen::seqN(0,i)).sum()+nrElemsBdry(i));
			startIdx += nrElemsBdry(i);
		}
	}

	// tangential vector
	Eigen::VectorXd tan(nDims);

	if(nDims == 2){
		// looping through the sliding elements
		for(int i = 0; i < nrElems; i++){

			// determining edge midpoint
			if(params.ptype){
				edgeMidPnts_polar_cylindrical.row(i) = (coords_polar_cylindrical.row(bdryNodesMat(indices(i),1)) + coords_polar_cylindrical.row(bdryNodesMat(indices(i),2)))/2;
			}
			edgeMidPnts.row(i) = (coords.row(bdryNodesMat(indices(i),1)) + coords.row(bdryNodesMat(indices(i),2)))/2;

			// finding tangential vector
			tan = (*ptrCoords).row(bdryNodesMat(indices(i),2)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));

			// establish normal vector
			edgeMidPntNormals1.row(i) << tan(1),-tan(0);

			// normalise normal vector
			edgeMidPntNormals1.row(i) = edgeMidPntNormals1.row(i)/tan.norm();
		}
	}else{


		// determining the midpoints and tangential vector of the edge segments
		Eigen::ArrayXXd midPntTan(size,nDims);
		for(int i = 0 ; i < size; i++ ){
			if(params.ptype){
				edgeMidPnts_polar_cylindrical.row(i) = (coords_polar_cylindrical.row(extBdryEdgeSegments(i,0)) + coords_polar_cylindrical.row(extBdryEdgeSegments(i,1)))/2;
			}
			edgeMidPnts.row(i) = (coords.row(extBdryEdgeSegments(i,0)) + coords.row(extBdryEdgeSegments(i,1)))/2;
			midPntTan.row(i) = (*ptrCoords).row(extBdryEdgeSegments(i,1)) - (*ptrCoords).row(extBdryEdgeSegments(i,0));
		}
		// determining the normals of the edge elements
		getPerpVecs(midPntTan,edgeMidPntNormals1,edgeMidPntNormals2);


		// finding normals and tangentials of the surface midpoints
		Eigen::VectorXd n(nDims), vec1(nDims), vec2(nDims);
		for(int i = 0; i < nrElems; i++){
			// normal vector
			n = Eigen::VectorXd::Zero(nDims);

			if(params.ptype){
				surfMidPnts_polar_cylindrical.row(i) = (coords_polar_cylindrical.row(bdryNodesMat(indices(i),1)) + coords_polar_cylindrical.row(bdryNodesMat(indices(i),2)) + coords_polar_cylindrical.row(bdryNodesMat(indices(i),3))+ coords_polar_cylindrical.row(bdryNodesMat(indices(i),4)))/4;
			}

			surfMidPnts.row(i) = (coords.row(bdryNodesMat(indices(i),1)) + coords.row(bdryNodesMat(indices(i),2)) + coords.row(bdryNodesMat(indices(i),3))+ coords.row(bdryNodesMat(indices(i),4)))/4;

			// finding two vector that make up the plane of the surface element
			vec1 = (*ptrCoords).row(bdryNodesMat(indices(i),2)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));
			vec2 = (*ptrCoords).row(bdryNodesMat(indices(i),4)) - (*ptrCoords).row(bdryNodesMat(indices(i),1));

			// determining the vector normal to both vec1 and vec2
			for(int k=0;k<nDims;k++){
				n(k) += (vec1((k+1)%3)*vec2((k+2)%3) - vec1((k+2)%3)*vec2((k+1)%3));
			}

			// normalising
			surfMidPntNormals.row(i) = n/n.norm();

		}

		// in case of periodic displacement the normal and tangential vector at the periodic edge nodes are required for the projection
		if(params.smode == "ps" && N_pe != 0){
			int ds = 0;					// Found vectors are not included in the sliding surface vector arrays
			getSurfNormalPeriodic(ds);
		}
	}
}


/* getCharPerLength function
 *
 * Due to periodicity a periodic vertex pair will have the same coordinate in the non-periodic directions and
 * the difference in coordinates is just the periodic length. This function aims at finding a coordinate pair by determining
 * the minimum coordinate values of the vertices in the non-periodic directions and then finding the vertices that share these values.
 * That will be a periodic pair and the difference in the periodic direction will be equal to the periodic length of the domain
 *
 */

void Mesh::getCharPerLength(probParams& params){

	// array containing the vertices of the domain
	Eigen::ArrayXXd vertices(verticesNodes.size(), nDims);

	// in case of rotational periodic domain
	if(params.ptype){
		// vertices in cartesian coordinates
		Eigen::ArrayXXd verticesCart(verticesNodes.size(),nDims);
		verticesCart = coords(verticesNodes,Eigen::all);

		// transformation to polar/cylindrical coordinates
		vertices.col(0) = verticesCart.leftCols(2).rowwise().norm();
		for(int row = 0; row < verticesCart.rows(); row++){
			vertices(row,1) = atan2(verticesCart(row,1), verticesCart(row,0));
		}
		if(verticesCart.cols() == 3){
			vertices.col(2) = verticesCart.col(2);
		}

	}
	// vertices in cartesian coordinates
	else{
		vertices = coords(verticesNodes,Eigen::all);
	}

	// finding the minimum values of the coordinates in the non-periodic directions
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

	// finding vertices that share the same minimum values in the non-periodic direction and saving their coordinate in periodic direction
	for(int i = 0; i < vertices.rows(); i++){
		if( fabs(vertices(i,cols[0]) - minVals[0]) < 1e-6 && (nDims == 2 || fabs(vertices(i,cols[1]) - minVals[1]) < 1e-6)){
			periodicVals(cnt) = vertices(i,params.pDir);
			cnt++;
		}
	}

	// establishing periodic length by taking difference between the max and min value of the periodic coordinate
	periodicVals.conservativeResize(cnt);
	periodic_length = periodicVals.maxCoeff()-periodicVals.minCoeff();
}


/* findStringBounds function
 *
 * This function removes any whitespaces that could be present in the input meshfile
 */

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

/* getPeriodicVector function
 *
 * This function establishes the vector in the periodic direction and two perpendicular vectors normal to this direction
 */

void Mesh::getPeriodicVector(probParams& params){
	// each column of periodicVecs contains the periodic vectors for either the cartesian/ polar/ spherical coordinates
	// the first column is the vector tangential to the periodic direction
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
