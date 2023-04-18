
#ifndef COORDTRANSFORM_H_
#define COORDTRANSFORM_H_
#include <iostream>
#include "Mesh.h"
#include "getNodeType.h"
class CoordTransform {
public:

	CoordTransform();
	void cart_to_polar_cylindrical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out);
	void polar_cylindrical_to_cart(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out);
	void vector_cart_to_polar_cylindrical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out, Eigen::ArrayXi& idx, Eigen::ArrayXXd& coords);
	void error_to_cart(Eigen::ArrayXXd& error, Mesh* m, getNodeType& n);
};

#endif /* COORDTRANSFORM_H_ */
