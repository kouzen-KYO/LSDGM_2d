/*
 * definitions.h
 *
 *  Include all the possible head files might be used in the program and declare namespace for some functions
 *
 */
#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// C libraries
#include <float.h>		// DBL_MAX
#include <math.h> 		// sin, cos
#include <stdio.h>    // printf, fgets
#include <stdlib.h>   // atoi
#include <omp.h>			// open mp stuff
#include <time.h>
#include <sys/stat.h>	// mkdir
#include <unistd.h>		// chdir

// C++ libraries
#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
//using std::isnan;

#include <list>
using std::list;

#include <string>			// strings
using std::string;
using std::to_string;

#include <iostream>		// terminal output
using std::cout;
using std::endl;

#include <fstream>		// file io
using std::stringstream;
using std::ifstream;
using std::getline;
using std::istringstream;

#include <algorithm>		// math stuff
using std::min;
using std::max;

#include <vector>			// standard vector
using std::vector;

#include <ctime>
using std::clock_t;

#include <iomanip>

#include <memory>
using std::shared_ptr;
using std::make_shared;

#include <Eigen/Core>
using Eigen::Vector2d;
using Eigen::Vector2i;
using Eigen::Vector3i;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix2d;

//#include <Eigen/StdVector>




#endif /* DEFINITIONS_H_ */
