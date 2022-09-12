//============================================================================
// Name        : ClassicalRBF.cpp
// Author      : F.A. van Steen
// Version     :
// Copyright   : Your copyright notice
// Description : C++ implementation of the standard RBF interpolation as mesh deformation method
//============================================================================

#include <iostream>
#include "Mesh.h"
#include <fstream>

using namespace std;


int main() {
	int nNodes = 5, mNodes =5;
//	int mNodes = 5;5
	double hDomain = 4, wDomain = 4;
//	double wDomain = 4;

	Mesh meshObject(nNodes, mNodes, hDomain, wDomain);
	meshObject.getMeshCoor();
//	meshObject.writeMeshFile();
}
