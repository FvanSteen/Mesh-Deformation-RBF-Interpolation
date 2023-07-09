#include "WriteResults.h"
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
WriteResults::WriteResults() {

}

/* createConvHistFile
 *
 * In case of data reduction methods a convergence history file is created
 *
 */

void WriteResults::createConvHistFile(std::string& dir){
	std::ofstream convHist;

	convHist.open(dir + "\\convergenceHistory.txt", std::ios::out);
	convHist << "step\tlevel\tmax error\ttime\tnr. of control nodes" << std::endl;
	convHist.close();


}


/* setIntResults
 *
 * writing data to the convergence history file
 */

void WriteResults::setIntResults(std::string& dir, int step, int lvl, double& maxErr, long double& time, int N){
	std::ofstream convHist;

	convHist.open(dir + "\\convergenceHistory.txt", std::ios::app);
	convHist.precision(6);
	convHist << step << '\t' << lvl << '\t' << maxErr  <<'\t' << time << '\t' << N << std::endl;

}

/* finalResult
 *
 * 	Writing information on the finalised RBF interpolation to a data file
 */

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


	std::string fName = p.directory + "\\" + p.mesh_ofName.substr(0, p.mesh_ofName.length()-4) + "_data" + ".txt";
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

/* existTest
 *
 * Function that checks if a file exists
 *
 */

bool WriteResults::existTest(std::string& fName){
	std::ifstream file(fName);
	return file.good();
}

