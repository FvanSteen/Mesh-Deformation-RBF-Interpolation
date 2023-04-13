#include "CoordTransform.h"

CoordTransform::CoordTransform()
{}


void CoordTransform::cart_to_polar_spherical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.rowwise().norm();
	for(int row = 0; row < in.rows(); row++){
		out(row,1) = atan2(in(row,1), in(row,0));
	}
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

//	std::cout <<m.coords_polar_spherical << std::endl;
}

void CoordTransform::polar_spherical_to_cart(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.col(0)*cos(in.col(1));
	out.col(1) = in.col(0)*sin(in.col(1));
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

//	std::cout << m.coords << std::endl;
}




