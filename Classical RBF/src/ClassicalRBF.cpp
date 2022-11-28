#include "rbf.h"
#include "probParams.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using namespace std;


int main()
{

	// amount of debugLvl messages
	int debugLvl = 2;
	// the radius of influence of the RBFs is the rFactor times the max domain length
	const double rFactor = 2.5;

	// input mesh file, output mesh file, internal and external boundary tags.
	string ifName = "su2mesh.su2";
	string ofName = "su2mesh_def_ps_ref.su2";
	vector<string> intBdryTags = {"wall1","wall2"};
	vector<string> extBdryTags = {"inflow","outflow","periodic1","periodic2"};

//	string ifName = "mesh_original.su2";
//	string ofName = "mesh_original_def.su2";
//	vector<string> intBdryTags = {"tube1","tube2","tube3","tube4","tube5","tube6","tube7","tube8","tube9","tube10",};
//	vector<string> extBdryTags = {"inlet","outlet","periodic_upper","periodic_lower"};



//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "mesh_NACA0012_inv_def_ps_ref.su2";
//	vector<string> intBdryTags = {"airfoil"};
//	vector<string> extBdryTags = {"farfield"};

//	const string ifName = "5x5.su2";
//	const string ofName = "5x5_def.su2";
//	string ifName = "25x25.su2";
//	string ofName = "25x25_def.su2";
//	const vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
//	const vector<string> extBdryTags = {"lower","right","upper", "left"};

//	string ifName = "5x5x5.su2";
//	string ofName = "5x5x5_def.su2";
//	const vector<string> intBdryTags = {"BLOC K"};
//	const vector<string> extBdryTags = {"FRONT", "BACK", "LEFT", "RIGHT", "LOWER", "UPPER"};


	// sliding modes:
	// regular rbf: none
	// pseudo sliding: ps
	// direct sliding: ds
	const string slidingMode = "none";
	// periodic modes:
	// no periodicity: none
	// periodic rbf: periodic
	// moving periodic boundaries with fixed vertices: fixed
	// moving periodic boundaries with moving vertices: moving
	const string periodicMode = "none";				// none for no periodicity, periodic for making the domain periodic in a to be specified direction ,fixed for allowing periodic boundary displacement with fixed corners, moving for periodic boundary displacement with moving corners

	// periodic boundary tags
	const vector<string> periodicBdry = {"upper","lower"};

	// direction of periodicity
	const string periodicDirection = "y";

//	const bool dataReduction = false;

	// initialising class object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl, slidingMode, periodicBdry, periodicMode);



	// constant deformation
	const double xDef = 0.0004, yDef = -0.0004, zDef= -0.0;
//	const double xDef = -0.0, yDef = -0.0, zDef= -0.0;

	// number of deformation steps
	int steps = 10;
	struct probParams probParams;

	probParams.dVec.resize(2);
	probParams.dVec << xDef/steps,yDef/steps;
	probParams.rotVec.resize(1);

	// rotational deformation
	probParams.rotVec << 0;
	probParams.steps = steps;
	probParams.rotPnt.resize(2);
	// point of rotation
	probParams.rotPnt << 0.5,0.5;
	probParams.sMode = slidingMode;
	probParams.pMode = periodicMode;
	// if the boundaries are curved then the normals and tangentials of the boundaries are updated
	probParams.curved = false;
	// periodic direction
	probParams.pDir = periodicDirection;
	// using data reduction techniques (greedy)
	probParams.dataRed = true;
	// tolerance for data reduction techniques
	probParams.tolerance = 5e-6;

	rbf rbf(meshOb, probParams);
	rbf.RBFMain();


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







}
