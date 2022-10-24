#include "rbf.h"
#include "Mesh.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
#include <math.h>
using Eigen::MatrixXd;

Mesh::Mesh(const std::string& inputFileName,const std::string& outputFileName,const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags,const double& rFac,const int& debugLvl, const std::string& slidingMode)
:ifName(inputFileName), ofName(outputFileName), nNodes(0), nDims(0),nElem(0),lvl(debugLvl), r(0), mode(slidingMode)
{readMeshFile(ibTags,ebTags);
r = rFac*charLength();
}

void Mesh::readMeshFile(const std::vector<std::string>& ibTags,const std::vector<std::string>& ebTags){
	if(lvl>=1){
		std::cout << "Reading mesh file: " << ifName << std::endl;
	}

	int lineNo = 0;								// line number counter
	int nIntBdryElems = 0, nExtBdryElems = 0;	// variables that sum the amount of marker elements
	bool intBdry = false, extBdry = false;		// booleans to acknowledge if on int/ext boundary
	int nIntBdryNodes = 0, nExtBdryNodes = 0;	// counters for int/ext boundary nodes
	int nPnts = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is

	extBdryNodesMat.resize(0,3);
	intBdryNodesMat.resize(0,3);
	int extBdryElemCnt = 0;
	int intBdryElemCnt = 0;// counting the elements of the external boundary
	Eigen::ArrayXi extBdryEndsIdx(2*ebTags.size());

	int MarkerElems;				// locally stores how many elements are in that boundary
	int nExtMarker = 0; 				// Counts the number of external boundary markers

	nrElemsExtBdry.resize(ebTags.size());
	std::string line;	// string containing line obtained by getline() function
	std::ifstream mFile("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" +ifName); 	//opening file name stored in mFile object
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
			else if (line.rfind("NMARK= ",0)==0 && lvl >= 1){
				try{
					if(stoi(line.substr(7)) == int(ibTags.size() + ebTags.size())){
						std::cout << "Number of boundary tags are matched" << std::endl;
					} else{
						throw(stoi(line.substr(7)));
					}
				}
				catch(int nTag){
					std::cout << "Number of tags provided (" << ibTags.size() + ebTags.size() << ") does not match number of tags found in file (" << nTag << ")." << std::endl;
					std::exit(0);
				}
			}

			// Finding tags of the boundaries
			else if (line.rfind("MARKER_TAG= ",0)==0){
				std::string tag =  line.substr(12);

				// Check if found tag is among provided internal boundary tags
				try{
					if(std::find(std::begin(ibTags), std::end(ibTags), tag) != std::end(ibTags)){
						if(lvl >=2){
							std::cout << "Saving nodes of internal boundary: " << tag << std::endl;
						}
						intBdry = true;
						extBdry = false;
					}
					// Check if found tag is among provided external boundary tags
					else if(std::find(std::begin(ebTags), std::end(ebTags), tag) != std::end(ebTags)){
						if(lvl >= 2){
							std::cout << "Saving nodes of external boundary: " << tag << std::endl;
						}
						intBdry = false;
						extBdry = true;
					} else throw(tag);
				}

				catch(std::string& tag){
					std::cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << std::endl;
					std::exit(0);
				}
			}

			// Check whether line corresponds to internal boundary
			else if(intBdry){
				// check if line states the number of elements in the boundary
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					nIntBdryElems += stoi(line.substr(13));	// updating number of int. boundary elements
					intBdryNodesMat.conservativeResize(nIntBdryElems,intBdryNodesMat.cols());

				}
				// check if line contains nodes corresponding to the boundary
				else if(isdigit(line[0])){
					std::istringstream is(line.substr(0)); 	// split line at '\t' and omitting the first number
					int data;
					int lineElem = 0;
					while(is >> data){
						if(lineElem == 0 && data == 5 && intBdryNodesMat.cols() < 4){
							intBdryNodesMat.conservativeResize(nIntBdryElems,4);
						}
						else if(lineElem == 0 && data == 9 && intBdryNodesMat.cols() < 5){
							intBdryNodesMat.conservativeResize(nIntBdryElems,5);
						}
						intBdryNodesMat(intBdryElemCnt,lineElem) = data;
						if(lineElem>0){
							nIntBdryNodes++;
						}
						lineElem++;
					}
					intBdryElemCnt++;
				}
			}

			// Check whether line corresponds to ext. boundary
			else if(extBdry){
				// Check if line states number of elements in boundary
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					MarkerElems = stoi(line.substr(13));
					nExtBdryElems+= MarkerElems; //todo MarkerElems is redundant?
					extBdryEndsIdx(nExtMarker) = nExtBdryElems-1;
					nrElemsExtBdry(nExtMarker) = MarkerElems;
					extBdryNodesMat.conservativeResize(nExtBdryElems,extBdryNodesMat.cols());
					nExtMarker++;
				}
				else if(isdigit(line[0])){
					std::istringstream is(line);
					int data;
					int lineElem = 0;
					while(is >> data){
						if(lineElem == 0 && data == 5 && extBdryNodesMat.cols() < 4){
							extBdryNodesMat.conservativeResize(nExtBdryElems,4);
						}
						else if(lineElem == 0 && data == 9 && extBdryNodesMat.cols() < 5){
							extBdryNodesMat.conservativeResize(nExtBdryElems,5);
						}
						extBdryNodesMat(extBdryElemCnt,lineElem) = data;
						if(lineElem>0){
							nExtBdryNodes++;
						}
						lineElem++;
					}
					extBdryElemCnt++;
				}
			}


			// Check if line corresponds to line containing node coordinates
			if (lineNo > pntsIdx && nPnts < nNodes){
				std::istringstream is(line);
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

		mFile.close();
	}
	else std::cout << "Not able to open input mesh file";

//	Eigen::ArrayXi idxSlidingSurf, idxSlidingEdge, idxExtStatic;


	if(mode=="ds" || mode=="ps"){
		getNodeType(nrElemsExtBdry, extBdryNodesMat);
		N_se = slidingEdgeNodes.size();
		N_ss = slidingSurfNodes.size();
		N_es = extStaticNodes.size();

		getEdgeConnectivity(nrElemsExtBdry);
		if(nDims == 3){
			getSurfConnectivity();
		}
	}

	getBdryNodes(intBdryNodesMat, intBdryNodes, nIntBdryNodes, nIntBdryElems);
	getBdryNodes(extBdryNodesMat, extBdryNodes, nExtBdryNodes, nExtBdryElems);



	midPnts.resize(nExtBdryElems,nDims);
	nVecs.resize(nExtBdryElems,nDims);
	tVecs.resize(nExtBdryElems,nDims);

	getIntNodes();

	N_i = intNodes.size();
	N_ib = intBdryNodes.size();
	N_eb = extBdryNodes.size();

//
//	std::cout << "Internal boundary nodes: \n" << intBdryNodes <<std::endl;
//	std::cout << "External boundary nodes: \n" << extBdryNodes <<std::endl;
//	std::cout << "Static external boundary nodes: \n" << extStaticNodes <<std::endl;
//	std::cout << "Sliding boundary nodes: \n" << slidingNodes << std::endl;
//	std::cout << "Internal nodes: \n " << intNodes << std::endl;
	std::cout << "Mesh file read successfully" << std::endl;

}

void Mesh::getNodeType(Eigen::ArrayXi& nrElemsExtBdry, Eigen::ArrayXXi& extBdryNodesMat){//todo extbdryNodesMat is a public variable
	Eigen::ArrayXi bdryNodesArr;
	Eigen::ArrayXi idxSS, idxSE, idxStatic;
	int cSS,cSE,cStat;

	for(int elem = 0; elem < nrElemsExtBdry.size(); elem++){
		bdryNodesArr.resize(nrElemsExtBdry(elem)*(extBdryNodesMat.cols()-1));
		for(int i =0; i<nrElemsExtBdry(elem);i++){
			for(int j=0; j< extBdryNodesMat.cols()-1; j++ ){
				bdryNodesArr(j+i*(extBdryNodesMat.cols()-1)) = extBdryNodesMat(i+ nrElemsExtBdry(Eigen::seqN(0,elem)).sum(),j+1);
			}
		}
		std::sort(std::begin(bdryNodesArr),std::end(bdryNodesArr));

		if(idxSS.size() != bdryNodesArr.size()){
			idxSS.resize(bdryNodesArr.size());
			idxSE.resize(bdryNodesArr.size());
			idxStatic.resize(bdryNodesArr.size());
		}

		cSS = 0; cSE = 0; cStat = 0;

		if(nDims == 3){
			for (int i= 0; i< bdryNodesArr.size();i++){
				if(i< bdryNodesArr.size()-3 && bdryNodesArr(i) == bdryNodesArr(i+1) && bdryNodesArr(i+1) == bdryNodesArr(i+2) && bdryNodesArr(i+2) == bdryNodesArr(i+3)){
					idxSS(cSS) = bdryNodesArr(i); cSS++;
					i+= 3;
				}else if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){
					idxSE(cSE) = bdryNodesArr(i); cSE++;
					i++;
				}else{
					idxStatic(cStat) = bdryNodesArr(i); cStat++;
				}
			}
		}else if(nDims ==2){
			for (int i= 0; i< bdryNodesArr.size();i++){
				if(i< bdryNodesArr.size()-1 && bdryNodesArr(i) == bdryNodesArr(i+1)){
					idxSE(cSE) = bdryNodesArr(i); cSE++;
					i++;
				}else{
					idxStatic(cStat) = bdryNodesArr(i); cStat++;
				}
			}
		}

		slidingSurfNodes.conservativeResize(slidingSurfNodes.size()+cSS);
		slidingSurfNodes(Eigen::lastN(cSS)) = idxSS(Eigen::seqN(0,cSS));

		slidingEdgeNodes.conservativeResize(slidingEdgeNodes.size()+cSE);
		slidingEdgeNodes(Eigen::lastN(cSE)) = idxSE(Eigen::seqN(0,cSE));

		extStaticNodes.conservativeResize(extStaticNodes.size()+cStat);
		extStaticNodes(Eigen::lastN(cStat)) = idxStatic(Eigen::seqN(0,cStat));

	}

	if(nDims==3){
		removeDuplicates(slidingSurfNodes);
	}

	removeDuplicates(slidingEdgeNodes);
	removeDuplicates(extStaticNodes);
}
void Mesh::getBdryNodes(Eigen::ArrayXXi& bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int& nBdryNodes, int& nBdryElems){
	bdryNodesArr.resize(nBdryNodes);
	int count = 0;
	for(int i = 0; i < nBdryElems; i++){
		int n;
		switch (bdryNodesMat(i,0)){
			case 3:
				n = 2;
				break;
			case 5:
				n = 3;
				break;
			case 9:
				n = 4;
				break;
		}
		int j = 0;
		while(j<n){
			bdryNodesArr(count) = bdryNodesMat(i,j+1);
			j++;count++;
		}
	}

	removeDuplicates(bdryNodesArr);
}

void Mesh::getExtStaticNodes(const std::vector<std::string>& ebTags, Eigen::ArrayXi& extBdryEndsIdx){
	// Adding first element indices
	for(int i = 0; i < int(ebTags.size())-1; i ++){
		extBdryEndsIdx(ebTags.size()+i) = extBdryEndsIdx(i)+1;
	}
	extBdryEndsIdx.tail(1) = 0;

	Eigen::ArrayXi extBdryEnds(extBdryEndsIdx.size()*(extBdryNodesMat.cols()-1));
//	movingNodes.resize(extBdryEndsIdx.size()*(extBdryNodesMat.cols()-1));

	int idx = 0;
	for(int i = 0; i < extBdryEndsIdx.rows(); i++){
		switch(extBdryNodesMat(extBdryEndsIdx(i),0)){
			case 3:
				extBdryEnds(Eigen::seqN(idx,2)) = extBdryNodesMat(extBdryEndsIdx(i),Eigen::seqN(1,2));

				idx = idx+2;
				break;
			case 5:
				extBdryEnds(Eigen::seqN(idx,3)) = extBdryNodesMat(extBdryEndsIdx(i),Eigen::seqN(1,3));
				idx = idx+3;
				break;
		}

	}
	extBdryEnds.conservativeResize(idx+1);
	std::sort(std::begin(extBdryEnds), std::end(extBdryEnds));

	extStaticNodes.conservativeResize(idx+1);
	int extStaticNodeIdx = 0;
	for(int i = 0; i<idx; i++){
		if(extBdryEnds(i+1) == extBdryEnds(i)){
			extStaticNodes(extStaticNodeIdx)=  extBdryEnds(i+1);
			extStaticNodeIdx++;
		}
	}
	extStaticNodes.conservativeResize(extStaticNodeIdx);
}

void Mesh::removeDuplicates(Eigen::ArrayXi& arr){
	/* This function determines the unique elements in the array provided as argument
	 * This is done by initialising a temporary array that will save the unique elements
	 * The provided array is sorted and subsequent values are compared
	 * Since the for loop does not take care of final element, a separate comparison is done
	 */

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

	arr = uniqueElems(Eigen::seq(0,cnt));
}

void Mesh::getSlidingNodes(){
	std::cout << "Determining sliding nodes" << std::endl;
	slidingEdgeNodes.resize(extBdryNodes.size()-extStaticNodes.size());

	int cnt = 0; int i=0, j=0;

	while(j < int(extBdryNodes.size())){							// while loop over all nodes
		if(i < int(extStaticNodes.size())){				// check if j is within the size of the boundary nodes

			if(extBdryNodes(j) == extStaticNodes(i)){				// if node i is equal to the j-th boundary elements then i is not a internal node
				i++; j++;						// going to next node i and element j of the boundary nodes
			}
			else if(extBdryNodes(j) < extStaticNodes(i)){			// if i is smaller than the j-th element of the boundary nodes
				slidingEdgeNodes(cnt) = extBdryNodes(j);				// then i is an internal node
				j++; cnt++;						// update index of internal nodes and i

			}
			else if(extBdryNodes(j) > extStaticNodes(i)){			// if i is larger than the j-th boundary element
				i++;							// then next boundary element should be checked.
			}
		}
		else{									// if all j boundary elements are checked and i < nNodes
			slidingEdgeNodes(cnt) = extBdryNodes(j);					// include remaining points as internal nodes
			cnt++; j++;							// update index and i
		}

	}
}

void Mesh::getIntNodes(){
	/* This function takes the previously found boundary nodes and compares them to an integer array containing
	 * all nodes from 0 to nNodes. If a node is not among the boundary nodes then the node is included in
	 * the intNodes.
	 */
	bdryNodes.resize(intBdryNodes.size()+extBdryNodes.size());
	bdryNodes << intBdryNodes, extBdryNodes;
	std::cout << "Determining internal nodes " << std::endl;

	std::sort(std::begin(bdryNodes), std::end(bdryNodes));		// bdryNodes need to be sorted for this algorithm to work
	intNodes.resize(nNodes-bdryNodes.size());	// Resizing intNodes accordingly

	int cnt = 0; int i=0, j=0;					// cnt keeps track of index of intNodes, i loops over all nodes, j over the boundary nodes

	while(i < nNodes){							// while loop over all nodes
		if(j < (bdryNodes.size())){				// check if j is within the size of the boundary nodes

			if(i == bdryNodes(j)){				// if node i is equal to the j-th boundary elements then i is not a internal node
				i++; j++;						// going to next node i and element j of the boundary nodes
			}
			else if(i < bdryNodes(j)){			// if i is smaller than the j-th element of the boundary nodes
				intNodes(cnt) = i;				// then i is an internal node
				i++; cnt++;						// update index of internal nodes and i

			}
			else if(i > bdryNodes(j)){			// if i is larger than the j-th boundary element
				j++;							// then next boundary element should be checked.
			}
		}
		else{									// if all j boundary elements are checked and i < nNodes
			intNodes(cnt) = i;					// include remaining points as internal nodes
			cnt++; i++;							// update index and i
		}

	}
	bdryNodes << intBdryNodes, extBdryNodes;
}

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

void Mesh::writeMeshFile(){
	std::cout << "Writing output file " << std::endl;
//	std::cout << coords << std::endl;
	std::ofstream outputF;
	outputF.precision(15);		// sets precision of the floats in the file

	outputF.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" + ofName, std::ios::out); // ios::out allows output to file

	std::ifstream inputF("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes\\" + ifName);
	std::string mystring;

	bool printFlag = false;
	int cnt = 0;

	while (getline(inputF, mystring)){
		if(printFlag && cnt < nNodes){
			if(nDims == 2){
				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t'<< cnt << std::endl;
			}
			else if(nDims == 3){
				outputF << coords(cnt,0)<< '\t' << coords(cnt,1) << '\t' << coords(cnt,2) << '\t' << cnt << std::endl;
			}
			cnt++;

		} else outputF << mystring << std::endl;
		if (mystring.rfind("NPOIN= ",0)==0){
			printFlag = true;
		}
	}

	inputF.close();
	outputF.close();
	std::cout << "Done writing mesh file: " << ofName << std::endl;
}

void Mesh::getExtBdryData(){
	double length,dx,dy;
	Eigen::Array2d n,t;
	for(int i=0; i < extBdryNodesMat.rows(); i++){
		dx = coords(extBdryNodesMat(i,2),0) - coords(extBdryNodesMat(i,1),0);
		dy = coords(extBdryNodesMat(i,2),1) - coords(extBdryNodesMat(i,1),1);
		midPnts(i,0) = coords(extBdryNodesMat(i,1),0)+ 0.5*dx;
		midPnts(i,1) = coords(extBdryNodesMat(i,1),1) +0.5*dy;
		t = {dx,dy};
		n = {dy,-dx};
		length = sqrt(pow(dx,2)+pow(dy,2));
		nVecs.row(i) = n/length;
		tVecs.row(i) = t/length;
	}
}



void Mesh::getNodeVecs(Eigen::ArrayXXd& n, Eigen::ArrayXXd& t){

	Eigen::ArrayXd dist;

	// Array with unsorted indices ranging from 0,1...N_midPnts
	Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(midPnts.rows(),0,midPnts.rows()-1);

	// Find for each sliding node the nearest midPnts
	for(int i = 0; i < int(slidingEdgeNodes.size()); i++){
		// Distance from sliding node to all mid-points
		dist = (midPnts.rowwise()-coords.row(slidingEdgeNodes(i))).rowwise().norm();
		// Sorting the indices with a custom comparator based on dist
		std::sort(index.begin(), index.end(),[&](const int& a, const int& b) {
		        return (dist[a] < dist[b]);
		    }
		);


		//todo is weighted average the best option?
		// weighted average of 2 nearest midpoints
//		n.row(i) = ( 1/dist(index(0))*nVecs.row(index(0)) + 1/dist(index(1))*nVecs.row(index(1)))/(1/dist(index(0)) + 1/dist(index(1)));
//		t.row(i) = ( 1/dist(index(0))*tVecs.row(index(0)) + 1/dist(index(1))*tVecs.row(index(1)))/(1/dist(index(0)) + 1/dist(index(1)));

		//simply nearest point
		n.row(i) = nVecs.row(index(0));
		t.row(i) = tVecs.row(index(0));
	}
}

void Mesh::getNormals(Eigen::ArrayXXd& n){

	Eigen::ArrayXd dist;

	// Array with unsorted indices ranging from 0,1...N_midPnts
	Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(midPnts.rows(),0,midPnts.rows()-1);

	// Find for each sliding node the nearest midPnts
	for(int i = 0; i < int(slidingEdgeNodes.size()); i++){
		// Distance from sliding node to all mid-points
		dist = (midPnts.rowwise()-coords.row(slidingEdgeNodes(i))).rowwise().norm();
		// Sorting the indices with a custom comparator based on dist
		std::sort(index.begin(), index.end(),[&](const int& a, const int& b) {
		        return (dist[a] < dist[b]);
		    }
		);
		//todo is weighted average the best option?
		// weighted average of 2 nearest midpoints
//		n.row(i) = ( 1/dist(index(0))*nVecs.row(index(0)) + 1/dist(index(1))*nVecs.row(index(1)))/(1/dist(index(0)) + 1/dist(index(1)));
//		t.row(i) = ( 1/dist(index(0))*tVecs.row(index(0)) + 1/dist(index(1))*tVecs.row(index(1)))/(1/dist(index(0)) + 1/dist(index(1)));

		//simply nearest point
		n.row(i) = nVecs.row(index(0));
	}
}

void Mesh::getEdgeConnectivity(Eigen::ArrayXi& nrElemsExtBdry){
	std::cout << "Obtaining edge connectivity" << std::endl;

	edgeConnectivity.resize(N_se,2);
	int row, idx, col, cnt;
	bool found;

	Eigen::ArrayXi bdryElement;
	// for each sliding edge node
	for(int node = 0; node < N_se; node++){
		col = 1;	// start from second column since first contains the element type
		cnt = 0;	// set count to zero, will go to 2 when both edges corresponding to a node are found

//		 for each column in extBdryNodesMat and while the cnt has not reached 2

		while(col<extBdryNodesMat.cols() && cnt < 2){


			// find the row that contains the slidingEdgeNode
			row = std::distance(std::begin(extBdryNodesMat.col(col)), std::find(std::begin(extBdryNodesMat.col(col)), std::end(extBdryNodesMat.col(col)), slidingEdgeNodes(node)));
//			std::cout << row << std::endl;
//			std::cout << extBdryNodesMat.row(row) << std::endl;
//				bdryElement = extBdryNodesMat(Eigen::seqN(startElem,nrElemsExtBdry(bdry)), col);
//				std::cout << bdryElement << std::endl;
//				row = std::distance(std::begin(bdryElement), std::find(std::begin(bdryElement), std::end(bdryElement), slidingEdgeNodes(node)));
//				std::cout << "row= " << row << std::endl;
			// if row is not equal to the amount of rows in extBdryNodesMat then that column contain the sliding edge node
			if (row!=extBdryNodesMat.rows()){
//					std::cout << "val is found in column" << std::endl;
				found = false; // boolean for checking if one of the edges corresponds to an ext static node
				// first a check is performed if the edge corresponds to an ext static node
//					std::cout << extBdryNodesMat.row(row) << std::endl;
				for (int j=1; j< extBdryNodesMat.cols(); j++){
					// if value j in row is not equal to the sliding edge node

					if(extBdryNodesMat(row,j) != slidingEdgeNodes(node)){
						// check if the value is in the extStaticNodes array
						idx = std::distance(std::begin(extStaticNodes), std::find(std::begin(extStaticNodes), std::end(extStaticNodes),extBdryNodesMat(row,j)));
						// if the value is in the array and either the count is zero or the found value is not equal to a previously found value
						if(idx!=N_es){
							if(cnt == 0 || extBdryNodesMat(row,j) != edgeConnectivity(node,cnt-1)){

//								std::cout << "Coupled Nodes: " << slidingEdgeNodes(node) << '\t' << extBdryNodesMat(row,j) << std::endl;
								edgeConnectivity(node,cnt) = extBdryNodesMat(row,j); // add found node to connectivity array
								cnt++; // update count
								found = true; // if an ext static value is found then there is no need to check the other sliding edge nodes
							} else{

								found = true;
							}
						}
					}
				}

				// if the edge does not correspond to a ext static node, then other sliding edge nodes are checked
				if(found!= true){
					// for each column of the row
					for (int j=1; j< extBdryNodesMat.cols(); j++){
						// if the value is not equal to the sliding edge node
						if(extBdryNodesMat(row,j) != slidingEdgeNodes(node)){
							// obtain index of the row
							idx = std::distance(std::begin(slidingEdgeNodes), std::find(std::begin(slidingEdgeNodes), std::end(slidingEdgeNodes),extBdryNodesMat(row,j)));

							// if the value exists in the array and the cnt is zero or the the value has not been found yet
							if(idx!=N_se && (cnt == 0 || extBdryNodesMat(row,j) != edgeConnectivity(node,cnt-1))){
//								std::cout << "Coupled Nodes: " << slidingEdgeNodes(node) << '\t' << extBdryNodesMat(row,j) << std::endl;
								edgeConnectivity(node,cnt) = extBdryNodesMat(row,j); // add found node to connectivity array
								cnt++;	// update count
							}
						}

					}
				}
			}
		col++;
//			}
//			std::cout << nrElemsExtBdry << std::endl;
//			startElem += nrElemsExtBdry(bdry);
		}

	}
//	for(int i=0; i < edgeConnectivity.rows(); i++){
//		std::cout << slidingEdgeNodes(i) << '\t'  << edgeConnectivity.row(i)<<std::endl;
//	}
}


void Mesh::getSurfConnectivity(){
	surfConnectivity.resize(N_ss, 4);
	int col,idx,cnt;
	for(int i=0; i<N_ss; i++){
		col = 1;
		cnt = 0;
		while(col < extBdryNodesMat.cols()){
			idx = std::distance(std::begin(extBdryNodesMat.col(col)), std::find( std::begin(extBdryNodesMat.col(col)), std::end(extBdryNodesMat.col(col)), slidingSurfNodes(i)));
			if(idx!=extBdryNodesMat.rows()){
//				std::cout << idx << '\t' << extBdryNodesMat.row(idx) << std::endl;
				surfConnectivity(i,cnt) = idx;
				cnt++;
			}
		col++;
		}
	}
//	std::cout << surfConnectivity << std::endl;
}

void Mesh::getVecs(Eigen::ArrayXXd& t_se, Eigen::ArrayXXd& n1_se, Eigen::ArrayXXd& n2_se, Eigen::ArrayXXd& n_ss, Eigen::ArrayXXd& t1_ss, Eigen::ArrayXXd& t2_ss){
//	Eigen::ArrayXXd tVecs_se(N_se,nDims);
//	Eigen::ArrayXXd nVecs1_se(N_se,nDims);
//	Eigen::ArrayXXd nVecs2_se(N_se,nDims);


	Eigen::RowVectorXd t(nDims), n(nDims);
	Eigen::RowVectorXd v1(nDims),v2(nDims), v3(nDims);
	double totalDist, dist;
	Eigen::VectorXd midPntVec;
	std::cout << "Finding sliding edge vectors" << std::endl;
//	for(int i =0;i<N_se;i++){
	for(int i =0;i<N_se;i++){
		//todo 2 terms can possibly be eliminated from the next line
		v1 = coords.row(edgeConnectivity(i,0)) - coords.row(slidingEdgeNodes(i));
		v2 = coords.row(slidingEdgeNodes(i)) - coords.row(edgeConnectivity(i,1));
//		t = (coords.row(edgeConnectivity(i,0)) - coords.row(slidingEdgeNodes(i)) + coords.row(slidingEdgeNodes(i)) - coords.row(edgeConnectivity(i,1)))/2;


//		t = (v1/v1.norm() + v2/v2.norm);
		t = (v1/v1.norm() + v2/v2.norm())/(1/v1.norm()+1/v2.norm());
		t_se.row(i) = t/t.norm();

	}

	getPerpVecs(t_se, n1_se, n2_se);

	for(int j= 0; j<N_ss; j++){
//		std::cout << "sliding surf node in question: " << slidingSurfNodes(j) << std::endl;
		n = Eigen::RowVectorXd::Zero(nDims);
		totalDist = 0;
		for (int k=0;k<4;k++){
//			std::cout << std::endl;
//			std::cout << extBdryNodesMat.row(surfConnectivity(j,k)) << std::endl;
			v1 = coords.row(extBdryNodesMat(surfConnectivity(j,k),2)) - coords.row(extBdryNodesMat(surfConnectivity(j,k),1));
			v2 = coords.row(extBdryNodesMat(surfConnectivity(j,k),4)) - coords.row(extBdryNodesMat(surfConnectivity(j,k),1));

			v3(0) = v1(1)*v2(2) - v1(2)*v2(1);
			v3(1) = v1(2)*v2(0) - v1(0)*v2(2);
			v3(2) = v1(0)*v2(1) - v1(1)*v2(0);
//			std::cout << v1 << '\t' << v2 << '\t' << v3 << std::endl;
			midPntVec = (coords.row(extBdryNodesMat(surfConnectivity(j,k),1)) + coords.row(extBdryNodesMat(surfConnectivity(j,k),2)) + coords.row(extBdryNodesMat(surfConnectivity(j,k),3)) + coords.row(extBdryNodesMat(surfConnectivity(j,k),4)))/4 - coords.row(slidingSurfNodes(j));
			dist = midPntVec.norm();
			totalDist += 1/dist;
			n += v3/dist;
			//
//			/for(int l = 0; l< surfConnectivity.cols();l++){
//				sum = sum + coords.row(extBdryNodesMat(surfConnectivity(j,k),l));
//			}

//			std::cout << totalDist << '\t'<< dist << std::endl;
//			std::cout << n <<  std::endl;

		}

		n = (n/totalDist);
		n_ss.row(j) = n/n.norm();
	}

	getPerpVecs(n_ss, t1_ss, t2_ss);
}

void Mesh::getPerpVecs(Eigen::ArrayXXd& v1,Eigen::ArrayXXd& v2, Eigen::ArrayXXd& v3){

	v2.col(0) = v1.col(1) - v1.col(2);
	v2.col(1) = v1.col(2) - v1.col(0);
	v2.col(2) = v1.col(0) - v1.col(1);

	v3.col(0) = v1.col(1)*v2.col(2) - v1.col(2)*v2.col(1);
	v3.col(1) = v1.col(2)*v2.col(0) - v1.col(0)*v2.col(2);
	v3.col(2) = v1.col(0)*v2.col(1) - v1.col(1)*v2.col(0);

	for(int x = 0; x<v2.rows();x++){
		v2.row(x) = v2.row(x)/ v2.row(x).matrix().norm();
		v3.row(x) = v3.row(x)/ v3.row(x).matrix().norm();
	}


//	for(int i=0; i<v1.rows();i++){
//		Eigen::VectorXd vec1 = v2.row(i);
//		Eigen::VectorXd vec2 = v3.row(i);
//
//		std::cout << vec1.dot(vec2) << std::endl;
//	}


}
