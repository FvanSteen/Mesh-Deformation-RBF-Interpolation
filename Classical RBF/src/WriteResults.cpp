#include "WriteResults.h"
#include <fstream>
#include <string>
WriteResults::WriteResults() {

}

void WriteResults::createConvHistFile(std::string& fName){
	std::ofstream convHist;


	convHist.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\convHist\\" + fName, std::ios::out);
	convHist << "step\tlevel\tmax error\tmean error\ttime" << std::endl;
	convHist.close();


}

void WriteResults::setIntResults(int step, int lvl, double maxErr, double meanErr, double time, std::string& fName, int N){
	std::ofstream convHist;


	convHist.open("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\convHist\\" + fName, std::ios::app);
	convHist.precision(6);
	convHist << step << '\t' << lvl << '\t' << maxErr << '\t' << meanErr <<'\t' << time << '\t' << N << std::endl;

}
