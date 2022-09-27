#include "rbf.h"
#include "Mesh.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
using Eigen::MatrixXd;
using namespace std;

Mesh::Mesh(std::string inputFileName,std::string outputFileName,vector<std::string> ibTags,vector<std::string> ebTags, double rFac,int debugLvl)
:ifName(inputFileName), ofName(outputFileName), nNodes(0), nDims(0),nElem(0),lvl(debugLvl), r(0)
{readMeshFile(ibTags,ebTags,rFac);
}

void Mesh::readMeshFile(vector<std::string> ibTags,vector<std::string> ebTags, double rFac){
	if(lvl>=1){
		cout << "Reading mesh file: " << ifName << endl;
	}

	int lineNo =0;	// keep track of the line number
	int nIntBdry = 0, nExtBdry = 0;		// ints for keeping track of number of internal/ external boundary points
	bool intBound = false, extBound = false;	// booleans to acknowledge if on int/ext boundary
	int cntInt = 0,cntExt = 0;					// counters for int/ext boundary nodes
	int point = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is

	string line;	// string containing line obtained by getline() function
	ifstream mFile(ifName); 	//opening file name stored in mFile object
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
				cout << "Saving node coordinates" << endl;
			}
			else if (line.rfind("NELEM= ",0)==0){
				nElem = stoi(line.substr(7));
			}

			// Checking whether provided boundary tags equals the amount in the mesh file
			else if (line.rfind("NMARK= ",0)==0 && lvl >= 1){
				try{
					if(stoi(line.substr(7)) == int(ibTags.size() + ebTags.size())){
						cout << "Number of boundary tags are matched" << endl;
					} else{
						throw(stoi(line.substr(7)));
					}
				}
				catch(int nTag){
					cout << "Number of tags provided (" << ibTags.size() + ebTags.size() << ") does not match number of tags found in file (" << nTag << ")." << endl;
					std::exit(0);
				}
			}

			// Finding tags of the boundaries
			else if (line.rfind("MARKER_TAG= ",0)==0){
				string tag =  line.substr(12);

				// Check if found tag is among provided internal boundary tags
				try{
					if(std::find(std::begin(ibTags), std::end(ibTags), tag) != std::end(ibTags)){
						if(lvl >=2){
							cout << "Saving nodes of internal boundary: " << tag << endl;
						}
						intBound = true;
						extBound = false;
					}
					// Check if found tag is among provided external boundary tags
					else if(std::find(std::begin(ebTags), std::end(ebTags), tag) != std::end(ebTags)){
						if(lvl >= 2){
							cout << "Saving nodes of external boundary: " << tag << endl;
						}
						intBound = false;
						extBound = true;
					} else throw(tag);
				}

				catch(string& tag){
					cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << endl;
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
					istringstream is(line.substr(1)); 	// split line at '\t' and omitting the first number
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
					nExtBdry += stoi(line.substr(13));		// updating number of external boundary elements
					extBdryNodes.conservativeResize(4*nExtBdry); 		// resizing the array containing ext. boundary nodes
				}
				else if(isdigit(line[0])){
					istringstream is(line.substr(1)); 	// split line at '\t' and omitting the first number
					int node;
					while(is >> node){					// While the stringstream contains nodes
						extBdryNodes[cntExt]  = node;	// include node in intBdryNodes
						cntExt++;						// update counter
					}
				}
			}

			// Check if line corresponds to line containing node coordinates
			if (lineNo > pntsIdx && point < nNodes){
				istringstream is(line);
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
	else cout << "Not able to open input mesh file";

	r = rFac*charLength();

	// Finding unique elements in boundary node arrays

	intBdryNodes.conservativeResize(cntInt-1);
	intBdryNodes = UniqueElems(intBdryNodes);

	extBdryNodes.conservativeResize(cntExt-1);
	extBdryNodes = UniqueElems(extBdryNodes);

	bdryNodes.resize(intBdryNodes.size()+extBdryNodes.size());
	bdryNodes << intBdryNodes, extBdryNodes;

	obtainIntNodes();
	cout << "Mesh file read successfully" << endl;
}

Eigen::VectorXi Mesh::UniqueElems(Eigen::ArrayXi& arr){
	/* This function determines the unique elements in the array provided as argument
	 * This is done by initialising a temporary array that will save the unique elements
	 * The provided array is sorted and subsequent values are compared
	 * Since the for loop does not take care of final element, a separate comparison is done
	 */

	sort(begin(arr), end(arr));

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


void Mesh::obtainIntNodes(){
	/* This function takes the previously found boundary nodes and compares them to an integer array containing
	 * all nodes from 0 to nNodes. If a node is not among the boundary nodes then the node is included in
	 * the intNodes.
	 */

	cout << "Determining internal nodes " << endl;

	sort(begin(bdryNodes), end(bdryNodes));		// bdryNodes need to be sorted for this algorithm to work
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

	double charLength = max(maxVal(0) - minVal(0),maxVal(1) - minVal(1));
	return charLength;
}

void Mesh::writeMeshFile(Eigen::MatrixXd& newCoords){
	ofstream outputF;
	outputF.precision(15);		// sets precision of the floats in the file
	outputF.open(ofName, ios::out); // ios::out allows output to file

	ifstream inputF(ifName);
	string mystring;

	bool printFlag = false;
	int cnt = 0;

	while (getline(inputF, mystring)){
		if(printFlag && cnt < nNodes){
			outputF << newCoords(cnt,0)<< '\t' << newCoords(cnt,1) << '\t'<< cnt << endl;
			cnt++;

		} else outputF << mystring << endl;
		if (mystring.rfind("NPOIN= ",0)==0){
			printFlag = true;
		}
	}

	inputF.close();
	outputF.close();
	cout << "Done writing mesh file" << endl;
}
