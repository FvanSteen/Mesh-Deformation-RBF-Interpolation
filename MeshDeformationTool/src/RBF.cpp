#include "getNodeType.h"
#include "rbf.h"

#include "rbfstd.h"
#include "rbfps.h"
#include "rbfds.h"
//#include "greedy.h"
//#include "selectRbf.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <Eigen/Dense>

rbf::rbf(Mesh& meshObject, probParams& probParamsObject)
:m(meshObject), params(probParamsObject)
//:mPtr(meshPtr), dVec(displacementVector), rotPnt(rotationPnt), steps(steps), smode(slidingMode), rotVec(rVec), curved(curved), perDir(periodicDirection), dataRed(dataReduction)
{}

void rbf::RBFMain(){


//	getNodeType n(m, params.dataRed);
//
//	if(params.smode=="none"){
////		 initialising the rbf standard class and performing the interpolation
//		rbf_std rbf(m, params);
//		rbf.perform_rbf(n);
//
//	}
//	else if(params.smode=="ds"){
//		rbf_ds rbf(m, params);
//		rbf.perform_rbf(n);
//	}
//	else if(params.smode=="ps"){
//
//		rbf_ps rbf(m,params);
//		rbf.perform_rbf(n);
//
//	}


}













