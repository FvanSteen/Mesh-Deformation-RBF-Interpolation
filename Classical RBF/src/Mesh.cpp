#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <chrono>
using Eigen::MatrixXd;

Mesh::Mesh(ReadConfigFile& cfg, const int& debugLvl)
:ReadConfigFile(cfg),lvl(debugLvl)
{
readMeshFile();
r = rFac*charLength();
std::cout << "RADIUS: " << r << std::endl;
}

// Main function for reading the .su2 mesh files
void Mesh::readMeshFile(){

	if(lvl>=1){
		std::cout << "Reading mesh file: " << mesh_ifName << std::endl;
	}

	int lineNo = 0;								// line number counter
	int nBdryElems = 0;	// variables that sum the amount of marker elements
//	bool bdry = true;		// booleans to acknowledge if on int/ext boundary
	int markerIdx = -2;
	int nBdryNodes = 0;	// counters for int/ext boundary nodes
	int nPnts = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is

	int bdryElemCnt = 0; // counting the elements of the boundaries
	int MarkerElems;							// locally stores how many elements are in that boundary
	int nMarker = 0; 						// Counts the number of external boundary markers

	bdryNodesMat.resize(0,3), intBdryNodesMat.resize(0,3);	// the int/ ext boundary node arrays have a minimum of 3 columns. One for the node type and at least two node indices.
																// The array will be adjusted to appropriate size depending on the boundary element type.
	srtdTags.resize(bdryTags.size());
	nrElemsBdry.resize(bdryTags.size());		// Array containing the sizes of each ext boundary
	std::string line;							// string containing line obtained by getline() function
	std::ifstream mFile("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" + mesh_ifName); 	//opening file name stored in mFile object
	// Check if file is opened
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){

			if (line.rfind("NDIME= ",0)==0){							// save number of dimensions
				nDims = stoi(line.substr(7));

			}
			else if (line.rfind("NPOIN= ",0)==0){						// save nr of points

				nNodes = stoi(line.substr(7));
				pntsIdx = lineNo;										// save line nr.
				coords.resize(nNodes, nDims);							// resizing the array containing the coordinates
				std::cout << "Saving node coordinates" << std::endl;
			}
			else if (line.rfind("NELEM= ",0)==0){						// save nr of elements

				nElem = stoi(line.substr(7));
			}

			// Checking whether provided boundary tags equals the amount in the mesh file
			else if (line.rfind("NMARK= ",0)==0 && lvl >= 2){

				try{
					if(stoi(line.substr(7)) == int(bdryTags.size())){
						std::cout << "Number of boundary tags are matched" << std::endl;
					} else{
						throw(stoi(line.substr(7)));
					}
				}
				catch(int nTag){
					std::cout << "Number of tags provided (" << bdryTags.size() << ") does not match number of tags found in file (" << nTag << ")." << std::endl;
					std::exit(0);
				}
			}

			// Finding tags of the boundaries
			else if (line.rfind("MARKER_TAG= ",0)==0){
				int cnt = 12;
				while(line[cnt] == ' '){
					cnt++;
				}

				std::string tag =  line.substr(cnt);

				try{
					if(std::find(std::begin(bdryTags), std::end(bdryTags), tag) != std::end(bdryTags)){
						if(lvl >=2){
							std::cout << "Saving nodes of boundary: " << tag << std::endl;
						}
						markerIdx = lineNo;
						srtdTags[nMarker] = tag;
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

//				if(line.rfind("MARKER_ELEMS= ",0)==0){
					// Obtaining number of elements in the boundary
				MarkerElems = stoi(line.substr(13));

					// Updating number of external boundary elements
				nBdryElems+= MarkerElems;
					// Saving number of elements of each boundary in an array
				nrElemsBdry(nMarker) = MarkerElems;
					// resizing the array containing the external boundary data
				bdryNodesMat.conservativeResize(nBdryElems,bdryNodesMat.cols());
				nMarker++;

				}

			else if(markerIdx > 0 && lineNo > markerIdx+1 && lineNo <= markerIdx+1 + MarkerElems){


				// split line by '\t' character
				std::istringstream is(line);
				// elements on the line are assigned to data in a while loop. lineElem counts the number of elements per line
				int data, lineElem = 0;
				while(is >> data){
					// First integer on line describes the type of element.
					// Triangle: 5,	Quadrilateral: 9

					// following if statements check the elementtype and if needed resize the intBdryNodesMat array.
					if(lineElem == 0 && data == 5 && bdryNodesMat.cols() < 4){
						bdryNodesMat.conservativeResize(nBdryElems,4);
					}
					else if(lineElem == 0 && data == 9 && bdryNodesMat.cols() < 5){
						bdryNodesMat.conservativeResize(nBdryElems,5);
					}

					// assigning the node information to the array
					bdryNodesMat(bdryElemCnt,lineElem) = data;
					// counting of the number of internal boundary nodes°
					if(lineElem>0){
						nBdryNodes++;
					}
					lineElem++;
				}
				bdryElemCnt++;
			}

			// Check if line corresponds to line containing node coordinates
			if (lineNo > pntsIdx && nPnts < nNodes){
				// split line by '\t' character
				std::istringstream is(line);
				// based on the number of dimensions 2 or 3 coordinates are assigned to the coords array
				switch(nDims){
					case 2:
						is >> coords(nPnts,0) >> coords(nPnts,1);
						break;
					case 3:
						is >> coords(nPnts,0) >> coords(nPnts,1) >> coords(nPnts,2);
						break;
				}
				nPnts++;
			}
			lineNo++;
		}
		//closing meshing file
		mFile.close();
	}
	// If the file is not opened then the following error message will be displayed
	else std::cout << "Not able to open input mesh file";

	getNodeTypes();
	N_se = seNodes.size();
	if(pmode == "moving"){
		std::cout << "adjusting the sliding edge nodes by adding the static nodes" << std::endl;
//		std::cout << N_se << std::endl;
		N_se += staticNodes.size();
//		std::cout << N_se << std::endl;
		seNodes.conservativeResize(N_se);
//		std::cout << seNodes << std::endl;
		seNodes(Eigen::lastN(staticNodes.size())) = staticNodes;
	}

//	std::cout <<"moving Nodes: \n" <<  mNodes << std::endl;
//	std::cout <<"sliding Nodes: \n" <<  seNodes << std::endl;
//	std::cout << "static Nodes: \n " << staticNodes << std::endl;
	N_m = mNodes.size();

//	std::cout << intBdryNodes << std::endl;


	ibIndices.resize(intBdryNodes.size());

	for(int x=0; x< intBdryNodes.size(); x++){
		ibIndices(x) = std::distance(std::begin(mNodes), std::find(std::begin(mNodes), std::end(mNodes),intBdryNodes(x)));
	}


//	std::cout << '\n' << mNodes << std::endl;


	getIntNodes();
//	std::cout << "internal Nodes: \n" << iNodes << std::endl;




//	N_ss = slidingSurfNodes.size();
//	N_es = extStaticNodes.size();
//	N_p = periodicNodes.size();


	std::cout << "node types are obtained" << std::endl;
//	std::cout << N_se << '\t' << N_ss  << '\t' << N_es  << '\t' << N_p << std::endl;
	//	std::cout << extBdryNodesMat << std::endl;


	getEdgeConnectivity();
//	std::cout << edgeConnectivity << std::endl;
	if(nDims == 3){
		getSurfConnectivity();
	}
//	}

	std::cout << "getEdgeConnectivity is done " << std::endl;



	std::cout << "Boundary nodes are stored " << std::endl;
//	std::cout << intBdryNodes << std::endl;
//	std::cout << "\n" << extBdryNodes << std::endl;

//	midPnts.resize(nExtBdryElems,nDims);
//	nVecs.resize(nExtBdryElems,nDims);
//	tVecs.resize(nExtBdryElems,nDims);

//	getIntNodes();
//	std::cout << "internal nodes are obtained" << std::endl;
	N_i = iNodes.size();
	N_ib = intBdryNodes.size();

//	N_eb = extBdryNodes.size();


//	std::cout << "Internal boundary nodes: \n" << intBdryNodes <<std::endl;
//	std::cout << "External boundary nodes: \n" << extBdryNodes <<std::endl;
//	std::cout << "Static external boundary nodes: \n" << extStaticNodes <<std::endl;
//	std::cout << "Sliding boundary nodes: \n" << slidingEdgeNodes << std::endl;
//	std::cout << "Internal nodes: \n " << intNodes << std::endl;
//	std::cout << "Periodic Nodes: \n " << periodicNodes << std::endl;


	// obtaining the nodes that make up the line segments of the external boundary
	getExtBdryEdgeSegments();
//	std::cout << extBdryEdgeSegments << std::endl;

	if(smode != "none"){
		N_mStd = N_m+N_se;
		mNodesStd.resize(N_mStd);
		mNodesStd << mNodes, seNodes;
	}
//	std::cout << mNodesStd << std::endl;

	if(pmode != "none"){
//		std::cout << "obtaining characteristic length" <<std::endl;
		getCharLength();
//		std::cout << "done" << std::endl;
	}


	if(dataRed){
		getIntCorNodes();
	}
	std::cout << "Mesh file read successfully" << std::endl;


}

/* getNodeType function
 *
 * This function will go through each boundary specified in the mesh file individually. *
 * In case of a 3D mesh the sliding surface nodes are identified by the fact that these nodes will be part of 3 or more surfaces and are therefore
 * as many times present in the array containing all external boundary nodes.
 * The sliding edge nodes are identified by the fact that these nodes only appear twice in the external boundary nodes array.
 * The remaining elements only appear once in the array and are therefore static external nodes.  *
 * For 2D meshes the sliding surface nodes don't have to be found.
 */
void Mesh::getNodeTypes(){
	if(lvl>=1){
		std::cout << "Obtaining node types" << std::endl;
	}

	Eigen::ArrayXi bdryNodesArr; 						// 1D array that will contain the all nodes for each respective boundary
	Eigen::ArrayXi idxMoving, idxSliding, idxStatic;		// Arrays containing specific type of nodes
	int cntMoving, cntSliding, cntStatic;								// counters for the number of sliding surface (SS) sliding edge (SE), static (Stat) and periodic (Per) nodes


	// array that will contain the external edge boundary nodes in the order of the meshing file.
	// This will be used to establish the line segments of the boundary and the corresponding midpoints.
	int nrSegments = 0;
	for(int elem = 0; elem < nrElemsBdry.size(); elem++){
		if(std::find(std::begin(mTags),std::end(mTags),srtdTags[elem]) == std::end(mTags) && (std::find(std::begin(pTags),std::end(pTags),srtdTags[elem]) == std::end(pTags) || pmode == "none" || pmode == "periodic" )){
			nrSegments += nrElemsBdry(elem) + 1;
		}
	}

	extBdryEdgeNodes.resize(nrSegments);

	int edgeNodeCnt = 0;
	bool periodic, moving;										// boolean that is set based on whether its a periodic boundary element or not.


	// for each external boundary the sliding edge, sliding surface and static nodes are identified
	for(int elem = 0; elem < nrElemsBdry.size(); elem++){

		// resizing the array that will contain all the boundary nodes of that boundary
		bdryNodesArr.resize(nrElemsBdry(elem)*(bdryNodesMat.cols()-1));

		// next two loops ensure that all nodes are included in the 1D node array.
		for(int i =0; i<nrElemsBdry(elem);i++){
			for(int j=0; j< bdryNodesMat.cols()-1; j++ ){
				bdryNodesArr(j+i*(bdryNodesMat.cols()-1)) = bdryNodesMat(i+ nrElemsBdry(Eigen::seqN(0,elem)).sum(),j+1);
			}
		}

		// sorting the array such that the nodes are ascending
		std::sort(std::begin(bdryNodesArr),std::end(bdryNodesArr));



		// Set size of the arrays of the various nodes equal to the size of the bdryNodesArr.

		if(idxMoving.size() != bdryNodesArr.size()){
			idxMoving.resize(bdryNodesArr.size());
			idxSliding.resize(bdryNodesArr.size());
			idxStatic.resize(bdryNodesArr.size());
		}

		// setting counters to zero
		cntMoving = 0, cntSliding = 0, cntStatic = 0;



		// In 2D only sliding edge and static nodes have to be considered for the sliding algorithms
		if(nDims ==2){
			moving = false;
			periodic = false;

			if(std::find(std::begin(mTags),std::end(mTags),srtdTags[elem]) != std::end(mTags)){
				std::cout << srtdTags[elem] << " is a moving boundary" << std::endl;
				moving = true;

			}else if((pmode == "moving" || pmode == "fixed") && std::find(std::begin(pTags),std::end(pTags),srtdTags[elem]) != std::end(pTags)){
				std::cout << srtdTags[elem] << " is periodic boundary" << std::endl;
				periodic = true;
			}

//			std::cout << "array: \n"<< bdryNodesArr << std::endl;
			if (periodic == false){
				for (int i= 0; i< bdryNodesArr.size();i++){

					if(moving){
						idxMoving(cntMoving) = bdryNodesArr(i);
						cntMoving++;
					}
	//				std::cout << bdryNodesArr(i) << std::endl;
					else{
					// checking if 2 subsequent nodes are equal. If so then its a sliding edge node
						if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){

							if(smode == "none"){
								idxMoving(cntMoving) = bdryNodesArr(i);
								cntMoving++;
							}else{

								idxSliding(cntSliding) = bdryNodesArr(i);
								cntSliding++;
							}

							i++;
						// else its a static node
						}else{
							if(pmode == "moving"){
//								idxSliding(cntSliding) = bdryNodesArr(i);
//								cntSliding++;
								// todo find a better name for this as these do slide in the periodic vector direction!
								idxStatic(cntStatic) = bdryNodesArr(i);
								cntStatic++;

//								std::cout << "elem: " << elem<< std::endl;
//								std::cout << "node: " << bdryNodesArr(i) << std::endl;
//								std::cout << "index: " << cntSliding-1 << std::endl;
							}else{
								idxMoving(cntMoving) = bdryNodesArr(i);
								cntMoving++;


							}
						}
						extBdryEdgeNodes(edgeNodeCnt) = bdryNodesArr(i);
						edgeNodeCnt++;
					}

					// saving the nodes in the order they were found to be able to establish the line segments of the external boundary
					// and their corresponding midpoints.

				}
			}
		}
		// Resize the arrays containing the different type of nodes according to how many are found in this boundary
		// And assigning last n-elements to the found nodes.


		if(moving){
			intBdryNodes.conservativeResize(intBdryNodes.size()+cntMoving);
			intBdryNodes(Eigen::lastN(cntMoving)) = idxMoving(Eigen::seqN(0,cntMoving));
		}

		if(cntStatic != 0){
			staticNodes.conservativeResize(staticNodes.size()+cntStatic);
			staticNodes(Eigen::lastN(cntStatic)) = idxStatic(Eigen::seqN(0,cntStatic));

		}

		mNodes.conservativeResize(mNodes.size()+cntMoving);
		mNodes(Eigen::lastN(cntMoving)) = idxMoving(Eigen::seqN(0,cntMoving));

		seNodes.conservativeResize(seNodes.size()+cntSliding);
		seNodes(Eigen::lastN(cntSliding)) = idxSliding(Eigen::seqN(0,cntSliding));

	}


	// Calling a function that removes any duplicate elements from the arrays.
	if(nDims==3){
		removeDuplicates(slidingSurfNodes);
	}

	removeDuplicates(mNodes);
	removeDuplicates(intBdryNodes);
	if(seNodes.size() != 0 ){
		removeDuplicates(seNodes);
	}


	// In case there is just a single external boundary then the final element is equal to the first element
	// This ensures that the line segments making up the external boundary is closed.

	if(int(bdryTags.size()-mTags.size())==1){
		extBdryEdgeNodes.conservativeResize(edgeNodeCnt+1);
		extBdryEdgeNodes(edgeNodeCnt) = extBdryEdgeNodes(0);
	}






}



/* getBdryNodes function
 *
 * This function is given an array that contains the information of the internal/external boundary.
 * Each row of this array starts with a integer that specifies the type of boundary element. Followed by the nodes that correspond to the element itself.
 * This function translates all these nodes to a 1D array.
 *
 */

//void Mesh::getBdryNodes(Eigen::ArrayXXi& bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int& nBdryNodes, int& nBdryElems){
//	// Resize the array containing the boundary nodes based on how many inputs where read in the mesh file.
//	bdryNodesArr.resize(nBdryNodes);
//
//	// set counter to zero
//	int count = 0;
//	// for each row corresponding to an element in the array.
//	for(int i = 0; i < nBdryElems; i++){
//		int n;
//		// Determine how many elements have to be added based on the integer determining the element type.
//		switch (bdryNodesMat(i,0)){
//			case 3:
//				n = 2;
//				break;
//			case 5:
//				n = 3;
//				break;
//			case 9:
//				n = 4;
//				break;
//		}
//		// adding the nodes to the array.
//		int j = 0;
//		while(j<n){
//			bdryNodesArr(count) = bdryNodesMat(i,j+1);
//			j++;count++;
//		}
//	}
//	// calling a function that removes any duplicates in the array.
//	removeDuplicates(bdryNodesArr);
//}

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

	// resize the array containing all boundary nodes and assign the int + ext boundary nodes to it.
	bdryNodes.resize(mNodes.size()+seNodes.size());
	bdryNodes << mNodes, seNodes;

	std::sort(std::begin(bdryNodes), std::end(bdryNodes));		// bdryNodes need to be sorted for this algorithm to work
	iNodes.resize(nNodes-bdryNodes.size());						// Resizing iNodes accordingly

	int cnt = 0, i=0, j=0;					// cnt keeps track of index of intNodes, i loops over all nodes, j over the boundary nodes


	while(i < nNodes){							// while loop over all nodes

		if(j < (bdryNodes.size())){				// check if j is within the size of the boundary nodes

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
//		std::cout << cnt << std::endl;

	}
	bdryNodes << mNodes,seNodes;	// Elsewhere, the code requires bdryNodes to be unsorted so intBdryNodes and extBdryNodes are assigned again

}

/*charLength function
 *
 * This function determines the maximum and minimal values in the coordinates.
 * Based on this the characteristic length of the domain is obtained. *
 */

double Mesh::charLength(){
	Eigen::Index row[2],col[2];
	Eigen::VectorXd maxVal(2);
	Eigen::VectorXd minVal(2);
	for(int i=0; i<2;i++){
		maxVal(i) = coords.col(i).maxCoeff(&row[i], &col[i]);
		minVal(i) = coords.col(i).minCoeff(&row[i], &col[i]);
	}
	double charLength = std::max(maxVal(0) - minVal(0),maxVal(1) - minVal(1));
	return charLength;
}

/* writeMeshFile function
 *
 * The writeMeshFile function writes the outputted mesh file containing the updated coordinates. This is done by reading the initial mesh file again
 * and copying its contents regarding the element connectivity and boundaries to the new mesh file. The part containing the coordinates is replaced by
 * the newly found coordinates.
 */

void Mesh::writeMeshFile(){
	std::cout << "Writing output file " << std::endl;

	std::ofstream outputF;		// Making an output stream class to operate on the output file
	outputF.precision(15);		// sets precision of the floats in the file

	// opening existing or creating new output file. In its respective folder.
	outputF.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" + mesh_ofName, std::ios::out); // ios::out allows output to file

	// Reopening the initial mesh file
	std::ifstream inputF("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" + mesh_ifName);
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
	std::cout << "Done writing mesh file: " << mesh_ofName << std::endl;
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

void Mesh::getEdgeConnectivity(){
	if(lvl>=1){
		std::cout << "Obtaining edge connectivity" << std::endl;
	}

	// resizing the array containing the connectivity information

	//todo make some counter for the number of elements
	edgeConnectivity.resize(N_se-staticNodes.size(),2);


	// row and col are the respective row and column index in the boundary node array.
	// idx is index in case either a static or sliding node is found.
	// cnt is a counter for the number of found connected nodes
	int row, col, idx,  cnt;

	// for-loop going through each sliding node
	for(int node = 0; node < N_se-staticNodes.size(); node++){


		col = 1;	// start from second column since first contains the element type
		cnt = 0;	// set count to zero, will go to 2 when both nodes connected by a line segment are found


//		While all columns have not been searched and while 2 nodes are not found
		while(col<bdryNodesMat.cols() && cnt < 2){

			// find the row that contains the sliding edge node
			row = std::distance(std::begin(bdryNodesMat.col(col)), std::find(std::begin(bdryNodesMat.col(col)), std::end(bdryNodesMat.col(col)), seNodes(node)));

			// if the sliding edge node is found than the row index is smaller than the number of rows of the entire array.
			if (row!=bdryNodesMat.rows()){

				// setting the index equal to the number of moving nodes, since this array includes the static nodes (zero movement).
				idx = mNodes.size();

				// starting from the second column, as first contains element type information
				int j = 1;

				// going through all columns until the idx changes or all columns have been done.
				while(j < bdryNodesMat.cols() && idx == mNodes.size()){
					idx = std::distance(std::begin(mNodes), std::find(std::begin(mNodes), std::end(mNodes),bdryNodesMat(row,j)));
					j++;
				}

				// if the index is not equal to the size of mNodes than a static node has been found.
				// this node should be included if either the count is zero or that node is not equal to the previous found node.
				if(idx!= mNodes.size() && ( cnt == 0 || mNodes(idx) != edgeConnectivity(node,cnt-1) ) ){
					edgeConnectivity(node,cnt) = mNodes(idx);
//					std::cout << "found node: "<<mNodes(idx) << std::endl;
					cnt++;
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
						if(seNodes(idx) == seNodes(node)){
							idx = seNodes.size();
						}
						j++;
					}

					// if an index is found and the count is either zero or that node has not been included than the node should be included.
					if(idx != seNodes.size() && ( cnt == 0 || seNodes(idx) != edgeConnectivity(node,cnt-1) ) ){
						edgeConnectivity(node,cnt) = seNodes(idx);
//						std::cout << "found node: "<< seNodes(idx) << std::endl;
						cnt++;
					}
				}
			}
		// Moving to the next column
		col++;
		}
	}
//	 debug message stating each sliding edge node and its found connecting nodes.
	if(lvl>=3){
		for(int i=0; i < edgeConnectivity.rows(); i++){
			std::cout << seNodes(i) << '\t'  << edgeConnectivity.row(i)<<std::endl;
		}
	}




}

/* getSurfConnectivity function
 *
 * This function aims at finding the boundary elements that include the sliding surface nodes.
 *
 * These are found by finding the sliding surface node in each column of extBdryNodesMat, the corresponding row consists of the nodes making up the element.
 * The found indices that correspond to an element which contain the i-th sliding surface node are saved as the i-th row in the surface connectivity information array.
 */


void Mesh::getSurfConnectivity(){
	if(lvl>=1){
			std::cout << "Obtaining surface connectivity" << std::endl;
		}

	// resizing of the surface connectivity information array. In case of an hexahedral mesh each sliding surface node is involved in four surfaces.
	surfConnectivity.resize(N_ss, 4);

	// col is the variable used to iterate through the columns of the external boundary nodes array.
	// cnt keeps track of how many elements are found.
	// idx is index or row of the found element.
	int col,cnt,idx;

	// for-loop going through each sliding surface node.
	for(int i=0; i<N_ss; i++){
		col = 1;	// First column contain information on the element type and is therefore skipped.
		cnt = 0;	// set count to zero.

		// iteratively going through the columns
		while(col < extBdryNodesMat.cols()){

			// Finding the index of the row containing a sliding surface node.
			idx = std::distance(std::begin(extBdryNodesMat.col(col)), std::find( std::begin(extBdryNodesMat.col(col)), std::end(extBdryNodesMat.col(col)), slidingSurfNodes(i)));

			// a sliding surface node is found when idx is not equal to the total amount of rows in extBdryNodesMat
			if(idx!=extBdryNodesMat.rows()){

				// debug message stating the found row.
				if(lvl>=4){
					std::cout << idx << '\t' << extBdryNodesMat.row(idx) << std::endl;
				}

				// including the found index in the surface connectivity array and updating the count.
				surfConnectivity(i,cnt) = idx;
				cnt++;
			}
		col++;
		}
	}
	// debug message outputting the entire surfConnectivity array
	if(lvl>=3){
		std::cout << surfConnectivity << std::endl;
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

void Mesh::getVecs(){
	if(lvl>=1){
//		std::cout << "Obtaining normal and tangential vectors " << std::endl;
	}

	// string that contains the type of element that is considered
	std::string type;

	// in case of 2D
	if(nDims == 2){
		// resizing of the normal and tangential vector arrays.
		n.resize(N_se-staticNodes.size(), nDims); t.resize(N_se-staticNodes.size(), nDims);
		// Calling function to obtain the tangential vectors along the line segments at the sliding boundary node.
		getEdgeTan(t);

		// set the type of element to edge and calling a function to obtain the normal vector perpendicular to the tangential vectors.
		type = "edge";
		getPerpVecs(type);

	}
	// else in case of 3D
	else if(nDims == 3){
		// resizing of the sliding edge normal and tangential vectors.
		t_se.resize(N_se,nDims); n1_se.resize(N_se,nDims); n2_se.resize(N_se,nDims);

		// calling function to obtain the tangential vectors along the boundary line segments at each sliding edge node.
		getEdgeTan(t_se);

		// set the type of element to edge and obtaining 2 vectors perpendicular to the tangential vectors.
		type = "edge";
		getPerpVecs(type);

		// resizing of the sliding surface normal and tangential vectors.
		n_ss.resize(N_ss,nDims); t1_ss.resize(N_ss,nDims); t2_ss.resize(N_ss,nDims);

		// calling function to get the normal vectors of the surface sliding nodes.
		getSurfNormal();

		// set the type of element to surface and calling a function to obtain the two vector perpendicular to the normal vectors.
		type = "surface";
		getPerpVecs(type);
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

void Mesh::getEdgeTan(Eigen::ArrayXXd& t){
	if(lvl>=2){
//		std::cout << "Obtaining edge tangential vectors" << std::endl;
	}

	// initialising vectors to store intermediate results.
	Eigen::VectorXd v1(nDims), v2(nDims), tan(nDims);

	// looping through each sliding edge node
	for(int i = 0; i < N_se-staticNodes.size(); i++){
		// defining vectors along the line segments. Its important to note that these vectors need to have the same direction.
		// Otherwise, the vector will (partially) cancel each other out when taking the average if they are opposite to each other.
		// So one vector goes from connectivity node 1 to the sliding node and the second vector goes from the sliding node to connectivity node 2.

		v1 = coords.row(edgeConnectivity(i,0)) - coords.row(seNodes(i));
		v2 = coords.row(seNodes(i)) - coords.row(edgeConnectivity(i,1));

		// Calculating the tangential vector based on a weighted average.
		tan = (v1/v1.norm() + v2/v2.norm())/(1/v1.norm()+1/v2.norm());
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
	for(int i = 0; i<N_ss; i++){
		// initialise a zero normal vector
		n = Eigen::VectorXd::Zero(nDims);

		// for all j nodes corresponding to a boundary element index in the surface connectivity array.
		for(int j=0; j < surfConnectivity.cols(); j++){

			// two vector are set up to define a plane of the boundary element.
			// one vector from the first node to the second node and a second vector from the first node to the last node.
			vec1 = coords.row(extBdryNodesMat(surfConnectivity(i,j),2)) - coords.row(extBdryNodesMat(surfConnectivity(i,j),1));
			vec2 = coords.row(extBdryNodesMat(surfConnectivity(i,j),4)) - coords.row(extBdryNodesMat(surfConnectivity(i,j),1));

			// initialise a zero distance vector
			dMidPnt = Eigen::ArrayXd::Zero(nDims);

			// summing all coordinates of the nodes of the element and storing it in dMidPnt
			for(int l = 1; l< extBdryNodesMat.cols(); l++){
				dMidPnt += coords.row(extBdryNodesMat(surfConnectivity(i,j),l));
			}

			// taking the average and substracting the coordinates of the sliding surface node to obtain the relative distance
			// from the midpoint of the element to the sliding surface node
			dMidPnt = dMidPnt/(extBdryNodesMat.cols()-1) - coords.row(slidingSurfNodes(i)).transpose();

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

void Mesh::getPerpVecs(std::string& type){
	// Finding normal of 2D problems.
	if(nDims == 2){
		n << t.col(1), -t.col(0);
	}

	// finding perpendicular vectors 3D problems
	if(nDims ==3){

		// in case of an edge type element, the two normals to an tangential vector should be found.
		if(type == "edge"){

			// Dot product for finding the first vector.
			for(int i = 0; i < nDims; i++){
				n1_se.col(i) = t_se.col((i+1)%3) - t_se.col((i+2)%3);
			}

			// Cross product for finding the second vector.
			for(int i =0; i < nDims; i++){
				n2_se.col(i) =  t_se.col((i+1)%3)*n1_se.col((i+2)%3) - t_se.col( (i+2)%3 )*n1_se.col((i+1)%3);
			}

			// Dividing the vectors by its length to ensure they are of unit length.
			for(int x = 0; x<n1_se.rows();x++){
				n1_se.row(x) = n1_se.row(x)/ n1_se.row(x).matrix().norm();
				n2_se.row(x) = n2_se.row(x)/ n2_se.row(x).matrix().norm();
			}
		}

		// in case of an edge type element, the two tangentials to an normal vector should be found.
		else if(type == "surface"){

			// Dot product for finding the first vector
			for(int i = 0; i < nDims; i++){
					t1_ss.col(i) = n_ss.col((i+1)%3) - n_ss.col((i+2)%3);
				}

			// Cross product for finding the second vector
			for(int i =0; i < nDims; i++){
				t2_ss.col(i) =  n_ss.col((i+1)%3)*t1_ss.col((i+2)%3) - n_ss.col( (i+2)%3 )*t1_ss.col((i+1)%3);
			}

			// Dividing vectors by its length to ensure they are of unit length.
			for(int x = 0; x<t1_ss.rows();x++){
				t1_ss.row(x) = t1_ss.row(x)/ t1_ss.row(x).matrix().norm();
				t2_ss.row(x) = t2_ss.row(x)/ t2_ss.row(x).matrix().norm();
			}
		}
	}
}



/* getExtBdryEdgeSegment function
 *
 * This function takes an array with the external boundary edge nodes in the order as specified in the mesh file and
 * creates an array (extBdryEdgeSegments) where the rows will contain the nodes making up all the boundary line segments.
 * This array can then be further used to find the midpoints of these line segments and its normals.
 * Therefore, the midPnts and midPntNormals are already resized accordingly in this function
 */

void Mesh::getExtBdryEdgeSegments(){
	if(lvl>=2){
//		std::cout << "obtaining the external boundary edge segments" << std::endl;
	}

	int nrSegments = 0;
	Eigen::ArrayXi sTag(nrElemsBdry.size());

	for(int i = 0; i<nrElemsBdry.size(); i++){
		if(std::find(std::begin(mTags),std::end(mTags),srtdTags[i]) == std::end(mTags) && (std::find(std::begin(pTags),std::end(pTags),srtdTags[i]) == std::end(pTags) || pmode == "none" || pmode == "periodic")){
			nrSegments += nrElemsBdry(i);

			sTag(i) = 1;

		}
		else{
			sTag(i) = 0;
		}
	}

	extBdryEdgeSegments.resize(nrSegments,2);


//	std::exit(0);
//
//	if(pmode == "fixed" || pmode == "moving"){
//		extBdryEdgeSegments.resize(extBdryEdgeNodes.size()-bdryTags.size() + mTags.size() + pTags.size(),2);
//	}else{
//		extBdryEdgeSegments.resize(extBdryEdgeNodes.size()-bdryTags.size() + mTags.size(),2);
//	}

	// starting index and segment counter
	int startIdx = 0,segment=0;
	// going through each external boundary seperately
	for(int i=0;i<nrElemsBdry.size();i++){
		if(sTag(i)){
			// for all boundary elements specified for that boundary
			// for a boundary having n elements there are n+1 nodes that make up n line segments.
			for(int j=0;j<nrElemsBdry(i);j++){

				// each row consists of the two nodes that make up a line segment
				extBdryEdgeSegments.row(segment) << bdryNodesMat(startIdx+j,1), bdryNodesMat(startIdx+j,2);

				segment++;
			}
		}
		startIdx += (nrElemsBdry(i));
	}

	// resizing the arrays containing the midpoints of the line segments and their normals
	midPnts.resize(extBdryEdgeSegments.rows(),nDims);
	midPntNormals.resize(extBdryEdgeSegments.rows(),nDims);

}

/* getMidPnts function
 *
 * This function obtains the midpoints of the external boundary edge segments by using the array
 * containing the nodes making up the line segments extBdryEdgeSegments.
 * This array is also used to find the tangential of the line segment,
 * which allows for obtaining the normal of the line segments.
 *
 */

void Mesh::getMidPnts(){
	if(lvl>=2){
//		std::cout << "obtaining the midpoints and normals of the external boundary edge segments" << std::endl;
	}



	// vector containing the tangential of a line segment
	Eigen::VectorXd tan(nDims);

	// for each line segment on the external boundary
	for(int i=0;i<extBdryEdgeSegments.rows();i++){

		// obtaining midpoint by averaging the coordinates of both nodes
		midPnts.row(i) = (coords.row(extBdryEdgeSegments(i,0)) + coords.row(extBdryEdgeSegments(i,1)))/2;

		// defining the tangential vector by moving from one node to the other.
		tan = coords.row(extBdryEdgeSegments(i,1)) - coords.row(extBdryEdgeSegments(i,0));

		// a normal of a tangential [x,y] is [y,-x]
		midPntNormals.row(i) << tan(1),-tan(0);

		// making the normal unit length
		midPntNormals.row(i) = midPntNormals.row(i)/tan.norm();

	}

//	std::cout << "midpoints and normals are obtained" << std::endl;
}




void Mesh::getSubDomains(Eigen::ArrayXi& subDomains, Eigen::ArrayXi& subDomLen, Eigen::ArrayXi& subDomBdry, Eigen::ArrayXi& subDomBdryLen){
//	std::cout << "dividing the subdomains" << std::endl;

	double max_x = coords.col(0).maxCoeff();
	double min_x = coords.col(0).minCoeff();

	double max_y = coords.col(1).maxCoeff();
	double min_y = coords.col(1).minCoeff();
//	std::cout << "max x " << max_x << " min x " << min_x << std::endl;
//	std::cout << "max y " << max_y << " min y " << min_y << std::endl;
	Eigen::ArrayXi domain1(nNodes), domain2(nNodes), domain3(nNodes), domain4(nNodes);
	Eigen::ArrayXi bdrydomain1(N_m+N_se), bdrydomain2(N_m+N_se), bdrydomain3(N_m+N_se), bdrydomain4(N_m+N_se);


	for(int i = 0; i < nNodes; i++){
//		std::cout << i << std::endl;
		if(std::find(std::begin(bdryNodes),std::end(bdryNodes),i) == std::end(bdryNodes)){
			if(coords(i,0) <= (max_x+min_x)/2){
				if(coords(i,1) <= (max_y+min_y)/2){
					domain1(subDomLen(0)) = i;
					subDomLen(0)++;
				}else{
					domain3(subDomLen(2)) = i;
					subDomLen(2)++;
				}
			}else{
				if(coords(i,1) <= (max_y+min_y)/2){
					domain2(subDomLen(1)) = i;

					subDomLen(1)++;
				}else{
					domain4(subDomLen(3)) = i;
					subDomLen(3)++;
				}
			}
		}else{
			if(coords(i,0) <= (max_x+min_x)/2){
				if(coords(i,1) <= (max_y+min_y)/2){
					bdrydomain1(subDomBdryLen(0)) = i;
					subDomBdryLen(0)++;
				}else{
					bdrydomain3(subDomBdryLen(2)) = i;
					subDomBdryLen(2)++;
				}
			}else{
				if(coords(i,1) <= (max_y+min_y)/2){
					bdrydomain2(subDomBdryLen(1)) = i;
					subDomBdryLen(1)++;
				}else{
					bdrydomain4(subDomBdryLen(3)) = i;
					subDomBdryLen(3)++;
				}
			}
		}
	}

	subDomains << domain1(Eigen::seqN(0,subDomLen(0))), domain2(Eigen::seqN(0,subDomLen(1))), domain3(Eigen::seqN(0,subDomLen(2))), domain4(Eigen::seqN(0,subDomLen(3)));
	subDomBdry << bdrydomain1(Eigen::seqN(0,subDomBdryLen(0))), bdrydomain2(Eigen::seqN(0,subDomBdryLen(1))),bdrydomain3(Eigen::seqN(0,subDomBdryLen(2))),bdrydomain4(Eigen::seqN(0,subDomBdryLen(3)));

//	std::cout << "subdomains are found" << std::endl;
}

void Mesh::getIntCorNodes(){
	std::cout << "Obtaining the internal nodes that fall within the correction tolerance" << std::endl;
	int nDomains = 4;
	Eigen::ArrayXi subDomains(N_i), subDomBdry(N_m+N_se);
	Eigen::ArrayXi subDomLen, subDomBdryLen;

	subDomLen = Eigen::ArrayXi::Zero(nDomains);
	subDomBdryLen = Eigen::ArrayXi::Zero(nDomains);



	getSubDomains(subDomains, subDomLen, subDomBdry, subDomBdryLen);

	Eigen::ArrayXXd bdryCoords;

	intCorNodes.resize(100);

//	std::cout << subDomLen << std::endl;
	int idxMin,startIdx = 0, startIdxBdry = 0, cnt=0;
	for(int i = 0; i < nDomains; i++){
		bdryCoords = coords(subDomBdry(Eigen::seqN(startIdxBdry,subDomBdryLen(i))), Eigen::all);
		for(int j = startIdx; j < startIdx + subDomLen(i); j++){
//			std::cout << j << std::endl;
			(bdryCoords.rowwise()- coords.row(subDomains(j))).rowwise().squaredNorm().minCoeff(&idxMin);
//			std::cout << "closest bdry node: " << bdryNodes(idxMin) << std::endl;
			if ((bdryCoords.row(idxMin) - coords.row(subDomains(j))).matrix().norm() < tol*gamma){
//				std::cout << "node in question: " << subDomains(j) << " cnt: " << cnt << std::endl;
//				std::cout << "min distance: " << (bdryCoords.row(idxMin) - coords.row(subDomains(j))).matrix().norm() << std::endl;
				intCorNodes(cnt) = subDomains(j);
				cnt++;

				if(cnt == intCorNodes.size()){
//					std::cout << cnt << std::endl;
					intCorNodes.conservativeResize(intCorNodes.size()+5000);
				}
			}
		}
		startIdx += subDomLen(i);
		startIdxBdry += subDomBdryLen(i);

	}


	intCorNodes.conservativeResize(cnt);

	// for debugging puposes
//	std::cout << intCorNodes << std::endl;
//	std::cout << intCorNodes.size() << std::endl;
//	std::cout << "done" << std::endl;
//	std::exit(0);




}

void Mesh::getCharLength(){
//	std::cout << "in the char length function" << std::endl;
	int idxMin;

	int minCoordDir, perDir;
	// only for 2D applications
	if(pDir == "x"){
		minCoordDir = 1;
		perDir = 0;
	}
	else{
		minCoordDir = 0;
		perDir = 1;
	}



	coords.col(minCoordDir).minCoeff(&idxMin);
//	std::cout << coords(idxMin,minCoordDir) << std::endl;

	int size = (coords.col(minCoordDir)==coords(idxMin,minCoordDir)).count();
//	std::cout << "size: " << size << std::endl;

	Eigen::ArrayXi indices(size);

	int cnt = 0;
	for(int i = 0; i < coords.rows(); i++){
		if(coords(i,minCoordDir) == coords(idxMin,minCoordDir)){
			indices(cnt) = i;
			cnt++;
		}
	}

	int min, max;

	Eigen::ArrayXd bdry(size);
	bdry << coords(indices,perDir);

	bdry.minCoeff(&min);
	bdry.maxCoeff(&max);

	lambda = bdry(max) - bdry(min);
}
