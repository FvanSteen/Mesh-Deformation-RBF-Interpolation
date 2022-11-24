#include "rbf.h"
#include "probParams.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using namespace std;


int main()
{
	int debugLvl = 2;
	const double rFactor = 2.5;

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "mesh_NACA0012_inv_def.su2";
//	vector<string> intBdryTags = {"airfoil"};
//	vector<string> extBdryTags = {"farfield"};

//	const string ifName = "5x5.su2";
//	const string ofName = "5x5_def.su2";
	string ifName = "25x25.su2";
	string ofName = "25x25_def.su2";
	const vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
	const vector<string> extBdryTags = {"lower","right","upper", "left"};

//	string ifName = "5x5x5.su2";
//	string ofName = "5x5x5_def.su2";
//	const vector<string> intBdryTags = {"BLOC K"};
//	const vector<string> extBdryTags = {"FRONT", "BACK", "LEFT", "RIGHT", "LOWER", "UPPER"};
	const string slidingMode = "ps";
	const bool curved = false; // only relevant in case of ds algorithm
	const string periodicMode = "none";				// none for no periodicity, periodic for making the domain periodic in a to be specified direction ,fixed for allowing periodic boundary displacement with fixed corners, moving for periodic boundary displacement with moving corners
	const vector<string> periodicBdry = {"upper","lower"};
	const string periodicDirection = "y";

	const bool dataReduction = false;

	// initialising class object m, reads mesh input file in constructor.
//	Mesh *meshPtr;
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl, slidingMode, periodicBdry, periodicMode);
//	meshPtr = &meshOb;




	const double xDef = -0.2, yDef = -0.32, zDef= -0.0;
	const int steps = 20;
	struct probParams probParams;
	probParams.dVec.resize(2);
	probParams.dVec << xDef/steps,yDef/steps;
	probParams.rotVec.resize(1);
	probParams.rotVec << 60;
	probParams.steps = 20;
	probParams.rotPnt.resize(2);
	probParams.rotPnt << 0.5,0.5;
	probParams.sMode = slidingMode;
	probParams.pMode = "none";
	probParams.curved = false;
	probParams.pDir = "y";
	probParams.dataRed = true;




//	const double xDef = -0.2, yDef = -0.32, zDef= -0.0;

//	const int steps = 20;
//	Eigen::VectorXd rotVec;
//	Eigen::RowVectorXd rotationPnt(meshOb.nDims);
//	Eigen::VectorXd displacementVector(meshOb.nDims);
//	if(meshOb.nDims == 2){
//		rotVec.resize(1);
//		rotVec << 60;
//		rotationPnt << 0.5,0.5;
//		displacementVector << xDef/steps,yDef/steps;
//	}else{
//		rotationPnt << 0.5,0.5,0.5;
//		rotVec.resize(3);
//		rotVec << 0,0,60;
//		displacementVector << xDef/steps,yDef/steps,zDef/steps;
//	}


	// todo maybe pass problem parameters as a structure later on. Needs definition of the structure in a separate header file.

//	struct problemParams{
//		Eigen::VectorXd dVec;
//		Eigen::RowVectorXd rotPnt;
//		Eigen::VectorXd rotVec;
//		const int steps;
//		const string smode;
//		const bool curved;
//		const string perDir;
//	};

//	struct problemParams pp = {displacementVector, rotationPnt, rotVec, steps, slidingMode, curved, periodicDirection};



	rbf rbf(meshOb, probParams);
	rbf.RBFMain();



}
