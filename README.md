# RBF-Slide a Mesh deformation tool based on the RBF interpolation technique

## Description
 The RBF-SliDe tool is an C++ implementation of radial basis functions interpolation technique for mesh deformation. Special methods are implemented for robust mesh deformation in case of internal flow applications. Sliding boundary nodes methods, periodic boundary conditions and data reduction techniques are implemented to enhance the performance of the mesh deformation method.

# How to use the tool
The tool can be directly run on Windows by running the provided executable. The tool is run with the call $ RBF_SLIDE your_config_file.cfg, where the executable, input mesh file, deformation file and configFile are located in the current working directory. Alternatively, the executable can be placed in a separate folder and its directory can be added to the path environment variable. This allows for calling of the tool from any folder.

Additionally, the source and makefile are in the repository to be able to build from source. The C++ library Eigen is required, since it is used for matrix operations within the RBF interpolations. See https://eigen.tuxfamily.org/index.php?title=Main_Page for more details and how to use the library.

The example configuration file named configFile.cfg shows the various possible inputs for the RBF-SliDe tool and a simple 25x25 square mesh is provided as an example test case.

# Developers
This tool was developed as part of an MSc. thesis project by:
***Floyd van Steen*** MSc. Student, TU Delft.

# Acknowledgements
**Matteo Pini**  Assistant Prof., Power and Propulsion, TU Delft.
**Alexander van Zuijlen** Assistant Prof., Aerodynamics, TU Delft.

# Citation
More details regarding the tool can be found in the following document:
1. F.A. van Steen. Mesh Deformation Using Radial Basis Function Interpolation for Numerical Optimisation of Internal Flow Applications. Master Thesis, TU Delft, 2023.
[![Link paper](https://img.shields.io/badge/MSc%20Thesis-2016-blue.svg)](https://repository.tudelft.nl/islandora/object/uuid%3Acfc2f6aa-5bc8-40b0-b581-2f005239ed0f?collection=education)  
