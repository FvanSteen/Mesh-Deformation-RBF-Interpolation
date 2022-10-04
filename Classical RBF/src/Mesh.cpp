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

	int lineNo =0;	// keep track of the line number
	int nIntBdry = 0, nExtBdry = 0;		// ints for keeping track of number of internal/ external boundary points
	bool intBound = false, extBound = false;	// booleans to acknowledge if on int/ext boundary
	int cntInt = 0,cntExtNodes = 0;					// counters for int/ext boundary nodes
	int point = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is
	int cntExtElems =0;			// counts the number of elements in the boundaries
	Eigen::ArrayXXi testArr;
	int testCnt = 0;
	//TODO make following statement dynamic
	Eigen::ArrayXi idxMovingNodes(8);	// indices of the moving nodes on the ext boundary
//	idxMovingNodes(0) = 0;			// first element is saved as moving node

	Eigen::ArrayXXi nodeData;	//testing

	int bdryElems;				// locally stores how many elements are in that boundary
	int cntExtBdry = 0;			// counting the number of ext boundaries

	std::string line;	// string containing line obtained by getline() function
	std::ifstream mFile(ifName); 	//opening file name stored in mFile object
	// Check if file is opened
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){

			// Following statements find information on various lines of the input file
			if (line.rfind("NDIME= ",0)==0){
				nDims = stoi(line.substr(7));
			}
			else if (line.rfind("NPOIN= ",0)==0){
				nNodes = stoi(line.substr(7));
				pntsIdx = lineNo;
				coords.resize(nNodes, nDims);
				std::cout << "Saving node coordinates" << std::endl;
			}
			else if (line.rfind("NELEM= ",0)==0){
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
						intBound = true;
						extBound = false;
					}
					// Check if found tag is among provided external boundary tags
					else if(std::find(std::begin(ebTags), std::end(ebTags), tag) != std::end(ebTags)){
						if(lvl >= 2){
							std::cout << "Saving nodes of external boundary: " << tag << std::endl;
						}
						intBound = false;
						extBound = true;
					} else throw(tag);
				}

				catch(std::string& tag){
					std::cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << std::endl;
					std::exit(0);
				}
			}

			// Check whether line corresponds to internal boundary
			else if(intBound){
				// check if line states the number of elements in the boundary
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					nIntBdry += stoi(line.substr(13));	// updating number of int. boundary elements
					intBdryNodes.conservativeResize(4*nIntBdry);	// resizing the array containing the int. boundary nodes

				}
				// check if line contains nodes corresponding to the boundary
				else if(isdigit(line[0])){
					std::istringstream is(line.substr(1)); 	// split line at '\t' and omitting the first number
					int node;
					while(is >> node){					// While the stringstream contains nodes
						intBdryNodes[cntInt]  = node;	// include node in intBdryNodes
						cntInt++;						// update counter
					}
				}
			}

			// Check whether line corresponds to ext. boundary
			else if(extBound){
				// Check if line states number of elements in boundary
				if(line.rfind("MARKER_ELEMS= ",0)==0){

					// If there are more than 2 ext boundaries than the first and last index is saved
					if(cntExtBdry>0 && cntExtBdry<int(ebTags.size())){
						idxMovingNodes(2*cntExtBdry-1) = cntExtNodes-1;
						idxMovingNodes(2*cntExtBdry) = cntExtNodes;
					}

					bdryElems = stoi(line.substr(13));
					nExtBdry += bdryElems;								// updating number of external boundary elements
//					extBdryNodes.conservativeResize(8*nExtBdry); 		// resizing the array containing ext. boundary nodes
					nVecs.conservativeResize(nExtBdry,nDims);
					tVecs.conservativeResize(nExtBdry,nDims);
					midPnts.conservativeResize(nExtBdry,nDims);
					nodeData.conservativeResize(nExtBdry,4);
					nodeData(Eigen::lastN(bdryElems),Eigen::all) = Eigen::ArrayXXi::Zero(bdryElems,4);
//					<< Eigen::ArrayXXi::Zero(nExtBdry,5);

					std::cout << nodeData << std::endl;
					cntExtBdry++;

					testArr.conservativeResize(nExtBdry,8);

				}
				else if(isdigit(line[0])){
//					getExtBdryData(line, cntExtElems, cntExtNodes,ebTags.size(),bdryElems,nodeData);
					std::istringstream is(line);
					int data;
					int nodeCnt = 0;
					while(is >> data){
						testArr(testCnt,nodeCnt) = data;
						nodeCnt++;
					}
					testCnt++;


				}
			}

			// Check if line corresponds to line containing node coordinates
			if (lineNo > pntsIdx && point < nNodes){
				std::istringstream is(line);
				switch(nDims){
					case 2:
						is >> coords(point,0) >> coords(point,1);
						break;
					case 3:
						is >> coords(point,0) >> coords(point,1) >> coords(point,2);
						break;
				}
				point++;
			}
			lineNo++;
		}
		mFile.close();
	}
	else std::cout << "Not able to open input mesh file";
	//
	std::cout << testArr << std::endl;
	std::exit(0);
//	std::cout << nodeData << std::endl;

	int colCnt = 0;
	while(colCnt < nodeData.cols()){
		if(nodeData.col(colCnt).sum() ==0){
			break;
		}
		colCnt++;
	}

	/* Made a better way to collect the boundary nodes. This is done by saving each element on a row of nodeData.
	 * This array is then passed to getNormalsTest in order to obtain the midPnts and vectors.
	 * The array with unique boundary nodes is obtained by conconating the columns of nodeData.
	 * followed by the respective sorting and removing duplicate elements.
	 * Same operation can be performed for the internal nodes.
	 *
	 * When the idx of the movingNodes are saved then the corresponding rows of nodeData can be accessed.
	 * the idx is the cumulative sum of the number of elements per external boundary.
	 * from the resulting rows the duplicate nodes should be saved since these are present on both boundaries and therfore
	 * are moving nodes with zero displacement
	 *
	 */

	nodeData.conservativeResize(nodeData.rows(),colCnt);

	midPntsT.resize(nExtBdry,nDims);
	nVecsT.resize(nExtBdry,nDims);
	tVecsT.resize(nExtBdry,nDims);
	getNormalsTest(nodeData);

	std::cout << nodeData << std::endl;

	extBdryNodes.resize(nodeData.rows()*colCnt);
	extBdryNodes << nodeData.col(0),nodeData.col(1);
	extBdryNodes = removeDuplicates(extBdryNodes);
	std::cout << extBdryNodes << std::endl;
		//
//	idxMovingNodes.tail(1) = cntExtNodes-1;
//	std::cout << idxMovingNodes<< std::endl;

//	movingNodes = extBdryNodes(idxMovingNodes);
//	movingNodes = removeDuplicates(movingNodes);
	intBdryNodes.conservativeResize(cntInt-1);
	intBdryNodes = removeDuplicates(intBdryNodes);

//	extBdryNodes.conservativeResize(cntExtNodes-1);
//	extBdryNodes = removeDuplicates(extBdryNodes);

	bdryNodes.resize(intBdryNodes.size()+extBdryNodes.size());
	bdryNodes << intBdryNodes, extBdryNodes;

	obtainIntNodes();
	//TODO if there is only a single element then movingNodes does not exist
	slidingNodes = uniqueElems(movingNodes, extBdryNodes);
//	slidingNodes = extBdryNodes;
	std::cout << "Mesh file read successfully" << std::endl;
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

Eigen::ArrayXi Mesh::uniqueElems(Eigen::ArrayXi& arr1, Eigen::ArrayXi& arr2){
	/* arr1 is small array
	 * arr2 is large array containing arr1
	 */
	std::cout << "Determining mutually exclusive indices" << std::endl;
	Eigen::ArrayXi nodes(arr2.size());
//	sort(begin(bdryNodes), end(bdryNodes));		// bdryNodes need to be sorted for this algorithm to work
//	intNodes.resize(nNodes-bdryNodes.size());	// Resizing intNodes accordingly

	int cnt = 0; int i=0, j=0;					// cnt keeps track of index of intNodes, i loops over all nodes, j over the boundary nodes

	while(j < int(arr2.size())){							// while loop over all nodes
		if(i < int(arr1.size())){				// check if j is within the size of the boundary nodes

			if(arr2(j) == arr1(i)){				// if node i is equal to the j-th boundary elements then i is not a internal node
				i++; j++;						// going to next node i and element j of the boundary nodes
			}
			else if(arr2(j) < arr1(i)){			// if i is smaller than the j-th element of the boundary nodes
				nodes(cnt) = arr2(j);				// then i is an internal node
				j++; cnt++;						// update index of internal nodes and i

			}
			else if(arr2(j) > arr1(i)){			// if i is larger than the j-th boundary element
				i++;							// then next boundary element should be checked.
			}
		}
		else{									// if all j boundary elements are checked and i < nNodes
			nodes(cnt) = arr2(j);					// include remaining points as internal nodes
			cnt++; j++;							// update index and i
		}

	}
	return nodes(Eigen::seq(0,cnt-1));
}

void Mesh::obtainIntNodes(){
	/* This function takes the previously found boundary nodes and compares them to an integer array containing
	 * all nodes from 0 to nNodes. If a node is not among the boundary nodes then the node is included in
	 * the intNodes.
	 */

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

void Mesh::getExtBdryData(std::string& line, int& cntExtElems, int& cntExtNodes, int nBdry, int bdryElems, Eigen::ArrayXXi& nodeData){


	std::istringstream is(line.substr(0)); 	// split line at '\t' and omitting the first number
	int data;
	Eigen::VectorXi dataVec(4);
	int elemCnt = 0;
	int nNodes;
	while(is >> data){
		if(elemCnt == 0){				// First element indicates the element type
			switch(data){
				case 3:
					nNodes = 2;
					break;
				case 5:
					nNodes = 3;
					break;
				case 9:
					nNodes = 4;
					break;
			}
		}
		else{
			dataVec(elemCnt-1) = data;	// save the remaining nodes
		}
		elemCnt++;
	}
//	std::cout << std::endl;
//	std::cout << cntExtElems << std::endl;
//	std::cout << dataVec << std::endl;


	nodeData(cntExtElems,Eigen::seq(0,nNodes-1)) = dataVec(Eigen::seq(0,nNodes-1));

//	getNormals(dataVec(Eigen::seq(0,nNodes-1)),cntExtElems);
//	extBdryNodes(Eigen::seq(cntExtNodes,cntExtNodes+nNodes-1)) = dataVec(Eigen::seq(0,nNodes-1));
	cntExtNodes+= nNodes;
	cntExtElems++;

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


