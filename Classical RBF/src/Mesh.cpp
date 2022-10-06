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

Mesh::Mesh(std::string inputFileName,std::string outputFileName,std::vector<std::string> ibTags,std::vector<std::string> ebTags, double rFac,int debugLvl)
:ifName(inputFileName), ofName(outputFileName), nNodes(0), nDims(0),nElem(0),lvl(debugLvl), r(0)
{readMeshFile(ibTags,ebTags,rFac);
r = rFac*charLength();
}

void Mesh::readMeshFile(std::vector<std::string> ibTags,std::vector<std::string> ebTags, double rFac){
	if(lvl>=1){
		std::cout << "Reading mesh file: " << ifName << std::endl;
	}

	int lineNo = 0;								// line number counter
	int nIntBdryElems = 0, nExtBdryElems = 0;	// variables that sum the amount of marker elements
	bool intBdry = false, extBdry = false;		// booleans to acknowledge if on int/ext boundary
	int nIntBdryNodes = 0, nExtBdryNodes = 0;	// counters for int/ext boundary nodes
	int nPnts = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is
	Eigen::ArrayXXi extBdryNodesMat(0,3);
	Eigen::ArrayXXi intBdryNodesMat(0,3);

	int extBdryElemCnt = 0;
	int intBdryElemCnt = 0;// counting the elements of the external boundary
	Eigen::ArrayXi extBdryEndsIdx(2*ebTags.size());

	int MarkerElems;				// locally stores how many elements are in that boundary
	int nExtMarker = 0; 				// Counts the number of external boundary markers

	std::string line;	// string containing line obtained by getline() function
	std::ifstream mFile(ifName); 	//opening file name stored in mFile object
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
					nExtBdryElems+= MarkerElems;
					extBdryEndsIdx(nExtMarker) = nExtBdryElems-1;
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


	getBdryNodes(intBdryNodesMat, intBdryNodes, nIntBdryNodes, nIntBdryElems);
	getBdryNodes(extBdryNodesMat, extBdryNodes, nExtBdryNodes, nExtBdryElems);

	getMovingNodes(ebTags,extBdryNodesMat,extBdryEndsIdx);

	getSlidingNodes();
	std::cout << "Internal boundary nodes: \n" << intBdryNodes <<std::endl;
	std::cout << "External boundary nodes: \n" << extBdryNodes <<std::endl;
	std::cout << "Moving boundary nodes: \n" << movingNodes <<std::endl;
	std::cout << "Sliding boundary nodes: \n" << slidingNodes << std::endl;




//
//	midPntsT.resize(nExtBdryElems,nDims);
//	nVecsT.resize(nExtBdryElems,nDims);
//	tVecsT.resize(nExtBdryElems,nDims);
//	getNormalsTest(nodeData);

	getIntNodes();

	std::cout <<"Internal nodes: \n"<< intNodes << std::endl;
	//TODO if there is only a single element then movingNodes does not exist

//	slidingNodes = extBdryNodes;
	std::cout << "Mesh file read successfully" << std::endl;
}

void Mesh::getBdryNodes(Eigen::ArrayXXi bdryNodesMat, Eigen::ArrayXi& bdryNodesArr, int nBdryNodes, int nBdryElems){
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
		}
		int j = 0;
		while(j<n){
			bdryNodesArr(count) = bdryNodesMat(i,j+1);
			j++;count++;
		}
	}
	bdryNodesArr = removeDuplicates(bdryNodesArr);
}

void Mesh::getMovingNodes(std::vector<std::string> ebTags, Eigen::ArrayXXi extBdryNodesMat, Eigen::ArrayXi extBdryEndsIdx){

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

	movingNodes.conservativeResize(idx+1);
	int movingNodeIdx = 0;
	for(int i = 0; i<idx; i++){
		if(extBdryEnds(i+1) == extBdryEnds(i)){
			movingNodes(movingNodeIdx)=  extBdryEnds(i+1);
			movingNodeIdx++;
		}
	}

	movingNodes.conservativeResize(movingNodeIdx);
}

Eigen::ArrayXi Mesh::removeDuplicates(Eigen::ArrayXi& arr){
	/* This function determines the unique elements in the array provided as argument
	 * This is done by initialising a temporary array that will save the unique elements
	 * The provided array is sorted and subsequent values are compared
	 * Since the for loop does not take care of final element, a separate comparison is done
	 */

	std::sort(std::begin(arr), std::end(arr));

	// Initialise array for the unique elements, could at most contain arr.size() unique elements
	Eigen::ArrayXi uniqueElmnts(arr.size());

	int cnt = 0;											// counter for the amount of unique elements
	for(int x = 0; x < arr.size()-1; x++){					// for loop over all but last element
		if(arr(x+1)!=arr(x)){								// if x+1-th is not equal to x-th element
			uniqueElmnts(cnt) = arr(x);						// include x-th element
			cnt++;											// update count of unique elements
		}
	}

	if(arr(arr.size()-1) != uniqueElmnts(cnt-1)){			// check if last element is already included
		uniqueElmnts(cnt) = arr(arr.size()-1);				// if not then include that element
	}

	return uniqueElmnts(Eigen::seq(0,cnt));
}

void Mesh::getSlidingNodes(){
	std::cout << "Determining sliding nodes" << std::endl;
	slidingNodes.resize(extBdryNodes.size()-movingNodes.size());

	int cnt = 0; int i=0, j=0;

	while(j < int(extBdryNodes.size())){							// while loop over all nodes
		if(i < int(movingNodes.size())){				// check if j is within the size of the boundary nodes

			if(extBdryNodes(j) == movingNodes(i)){				// if node i is equal to the j-th boundary elements then i is not a internal node
				i++; j++;						// going to next node i and element j of the boundary nodes
			}
			else if(extBdryNodes(j) < movingNodes(i)){			// if i is smaller than the j-th element of the boundary nodes
				slidingNodes(cnt) = extBdryNodes(j);				// then i is an internal node
				j++; cnt++;						// update index of internal nodes and i

			}
			else if(extBdryNodes(j) > movingNodes(i)){			// if i is larger than the j-th boundary element
				i++;							// then next boundary element should be checked.
			}
		}
		else{									// if all j boundary elements are checked and i < nNodes
			slidingNodes(cnt) = extBdryNodes(j);					// include remaining points as internal nodes
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

void Mesh::writeMeshFile(Eigen::MatrixXd& newCoords){
	std::ofstream outputF;
	outputF.precision(15);		// sets precision of the floats in the file
	outputF.open(ofName, std::ios::out); // ios::out allows output to file

	std::ifstream inputF(ifName);
	std::string mystring;

	bool printFlag = false;
	int cnt = 0;

	while (getline(inputF, mystring)){
		if(printFlag && cnt < nNodes){
			outputF << newCoords(cnt,0)<< '\t' << newCoords(cnt,1) << '\t'<< cnt << std::endl;
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

void Mesh::getNormals(Eigen::VectorXi nodes, int& cntExtElems){
	Eigen::Vector2d n;
	Eigen::Vector2d t;
	Eigen::Vector2d midPnt;
	double length;
	if(nodes.size()<3){
		double dx = coords(nodes(1),0) - coords(nodes(0),0);
		double dy = coords(nodes(1),1) - coords(nodes(0),1);
		t = {dx,dy};
		n = {dy,-dx};
		length = sqrt(pow(dx,2)+pow(dy,2));
		midPnt = {coords(nodes(0),0)+ 0.5*dx, coords(nodes(0),1) +0.5*dy};
	}
	midPnts.row(cntExtElems) = midPnt;
	nVecs.row(cntExtElems) = n/length;
	tVecs.row(cntExtElems) = t/length;


}

void Mesh::getNormalsTest(Eigen::ArrayXXi nodes){
	double length,dx,dy;
	Eigen::Array2d n,t;
	for(int i=0; i < nodes.rows(); i++){
		dx = coords(nodes(i,1),0) - coords(nodes(i,0),0);
		dy = coords(nodes(i,1),1) - coords(nodes(i,0),1);
		midPntsT(i,0) = coords(nodes(i,0),0)+ 0.5*dx;
		midPntsT(i,1) = coords(nodes(i,0),1) +0.5*dy;
		t = {dx,dy};
		n = {dy,-dx};
		length = sqrt(pow(dx,2)+pow(dy,2));
		nVecsT.row(i) = n/length;
		tVecsT.row(i) = t/length;
	}
}


void Mesh::getNodeVecs(Eigen::ArrayXi& idxs, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t){
	int idx, idx2;
//	struct vecs{
//		Eigen::ArrayXXd n;
//		Eigen::ArrayXXd t;
//	};
//
//	vecs v;
//	v.n.resize(idxs.size(),2);
//	v.t.resize(idxs.size(),2);
	Eigen::ArrayXd dist;
	Eigen::ArrayXd distSort;
	for(int i = 0; i < int(idxs.size()); i++){
		auto d = midPnts.rowwise()-coords.row(idxs(i));
		dist = sqrt(pow(d.col(0),2)+pow(d.col(1),2));
		distSort = dist;
		std::sort(std::begin(distSort),std::end(distSort));

		if(distSort(0)==distSort(1)){
			idx = std::distance(dist.begin(),std::find(dist.begin(), dist.end(),distSort(0)));
			idx2 = std::distance(dist.begin()+idx,std::find(dist.begin()+idx+1, dist.end(),distSort(1)));
			idx2 += idx;
		}
		else{
		idx = std::distance(dist.begin(),std::find(dist.begin(), dist.end(),distSort(0)));

		idx2 = std::distance(dist.begin(),std::find(dist.begin(), dist.end(),distSort(1)));
		}
		n.row(i) = (nVecs.row(idx)+nVecs.row(idx2))/2;
		t.row(i) = (tVecs.row(idx)+tVecs.row(idx2))/2;
	}
}


