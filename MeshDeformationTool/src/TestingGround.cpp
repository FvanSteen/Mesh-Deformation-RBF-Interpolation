#include "TestingGround.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
#include <chrono>

TestingGround::TestingGround()
{
	testFunction();

}

void TestingGround::testFunction(){
	std::cout << "In the TestingGround source file." << std::endl;
	std::cout << "are changes faster here" << std::endl;

	std::exit(0);
}

