#include "WriteResults.h"
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
WriteResults::WriteResults() {

}

void WriteResults::createConvHistFile(std::string& fName){
	std::ofstream convHist;


	convHist.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\convHist\\" + fName, std::ios::out);
	convHist << "step\tlevel\tmax error\tmean error\ttime" << std::endl;
	convHist.close();


}

void WriteResults::setIntResults(int step, int lvl, double maxErr, double meanErr, double time, std::string& fName, int N){
	std::ofstream convHist;


	convHist.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\convHist\\" + fName, std::ios::app);
	convHist.precision(6);
	convHist << step << '\t' << lvl << '\t' << maxErr << '\t' << meanErr <<'\t' << time << '\t' << N << std::endl;

}

void WriteResults::finalResult(probParams& p, long double& time, Eigen::Array3d& quals){

	int greedy = 0;
	int doubleEdged = 0;
	double tol = 0.;
	double tolCrit = 0.;
	if(p.dataRed){
		if(p.multiLvl){
			greedy = 2;
			tolCrit = p.tolCrit;
		}else{
			greedy = 1;
		}
		if(p.doubleEdge){
			doubleEdged = 1;
		}
		tol = p.tol;
	}

	if(!p.generateQuality){
		quals << 0,0,0;
	}


	std::string fName = p.directory + "\\runHistory\\" + p.mesh_ifName.substr(0, p.mesh_ifName.length()-4) + ".txt";
	bool existing = existTest(fName);

	std::ofstream f;
	f.precision(15);
	if(!existing){
		f.open(fName, std::ios::out);
		f << "sliding type\tperiodicity type\tgreedy\tdoubled edge\ttol\ttolCrit\tsteps\tcpu time [ms]\tQmin\tQmean\tQmax"<< std::endl;
	}else{
		f.open(fName, std::ios::app);

	}
	f << p.smode << '\t' << p.pmode << '\t' << greedy << '\t'<< doubleEdged << '\t'<< tol<< '\t'<< tolCrit << '\t'<< p.steps << '\t' << time << '\t' << quals[0] << '\t' << quals[1] << '\t' << quals[2] <<  std::endl;
	f.close();
}

bool WriteResults::existTest(std::string& fName){
	std::ifstream file(fName);
	return file.good();
}

