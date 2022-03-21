/*
 * readInputFile.h
 *
 *  Read grain information from input file made by MATLAB codes
 *
 */
#ifndef READINPUTFILE_H_
#define READINPUTFILE_H_

#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"

// Creates a vector of grain objects from a file
vector<Grain2d> generateGrainsFromFile(string filename) {

	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// Read first line: type number of particles (not real simulation grain number)
	getline(file, line);
	int numberOfGrains = atoi(line.c_str());

	vector<Grain2d> grainList(numberOfGrains);

	// Temp stuff
	Vector2d point;
	for (int grainidx = 0; grainidx < numberOfGrains; grainidx++) {
		// 01 -- Mass
		getline(file, line);
		double mass = atof(line.c_str());

		// 02 -- Mass center for Current position
		getline(file, line);
		Vector2d position;
		iss.str(line);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		iss.clear();

		// 03 -- Velocity
		getline(file, line);
		Vector2d velocity;
		iss.str(line);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		iss.clear();

		// 04 -- Moment of inertia
		getline(file, line);
		double momentOfInertia = atof(line.c_str());

		// 05 -- Rotation
		getline(file, line);
		double theta = atof(line.c_str());

		// 06 -- Angular speed
		getline(file, line);
		double omega = atof(line.c_str());

		// 07 -- Mass center for Level-set
		getline(file, line);
		Vector2d cmLset;
		iss.str(line);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		iss.clear();

		// 08 -- Number of surface points (INTEGER)
		getline(file, line);
		int npoints = atoi(line.c_str());

		// 09 -- Each surface points positions
		getline(file, line);
		vector<Vector2d> pointList(npoints);
		iss.str(line);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// 10 -- Bounding-box radius
		getline(file, line);
		double bboxRadius = atof(line.c_str());

		// 11 -- level set dimensions (INTEGERS)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		iss.clear();

		// 12 -- level set values in each grid node
		getline(file, line);
		vector<double> lsetvec(xdim*ydim);
		iss.str(line);
		for (int i = 0; i < xdim*ydim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// 13 -- kn
		getline(file, line);
		double kn = atof(line.c_str());

		// 14 -- ks
		getline(file, line);
		double ks = atof(line.c_str());

		// 15 -- mu
		getline(file, line);
		double mu = atof(line.c_str());

		// Create level set object based on above information
		Levelset2d lset(lsetvec, xdim, ydim);
		shared_ptr<Levelset2d> lsetptr = make_shared<Levelset2d>(lset);

		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, velocity, momentOfInertia, theta, omega, cmLset, pointList, bboxRadius, lsetptr, kn, ks, mu, grainidx);
//		grainList[grainidx].updateVariables(mass, position, velocity, momentOfInertia, theta, omega, cmLset, bboxCenterOffset, pointList, bboxRadius, lset, kn, ks, mu);
	}

	return grainList;
}


extern vector<Vector2d> readPositionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector2d> positions;
	positions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		iss.clear();
	}
	return positions;
} // end generateGrainsFromFile

extern vector<double> readRotationFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> rotations;
	rotations.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		rotations[grainidx] = atof(partial.c_str());
		iss.clear();
	}
	return rotations;
}

#endif /* READINPUTFILE_H_ */
