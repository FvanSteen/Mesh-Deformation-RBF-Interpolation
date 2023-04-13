
#ifndef COORDTRANSFORM_H_
#define COORDTRANSFORM_H_
#include <iostream>
#include "Mesh.h"

class CoordTransform {
public:

	CoordTransform();
	void cart_to_polar_spherical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out);
	void polar_spherical_to_cart(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out);

};

#endif /* COORDTRANSFORM_H_ */
