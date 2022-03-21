/*
 * main.cpp
 *
 *  Created on: May 30, 2014
 *      Author: Reid
 */
 
#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
//#include "Wall2dFixed.h"
#include "readInputFile.h"
#include "World2dBiaxial.h"

// drop test for spheres (run after main5)
int main() {
	// file name
	string filename("discShearStart.dat");
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	// grain properties
	double mu = 0.4;
	double kn = 30000;
	double density = 1;
	// assembly properties
//	double width = 350; 
	// world properties
	double gDamping = .1;
	double grav = .1;
	double dt = 0.2*(2*sqrt(.5*113/kn));
	int tdrop = 500001;
	
	vector<Grain2d> grains2 = generateGrainsFromFile(filename);
	vector<Grain2d> grains3 = generateGrainsFromFile(filename);
	grains.insert(grains.end(), grains2.begin(), grains2.end() );
	grains.insert(grains.end(), grains3.begin(), grains3.end() );
	
	size_t ngrains = grains.size();
	
	vector<Vector2d> pos = readPositionFile("positions_disc.dat",ngrains);
	vector<double> rot = readRotationFile("rotations_disc.dat",ngrains);
	
	/* set things up */
	for (size_t i = 0; i < grains.size(); i++) {
		grains[i].changePos(pos[i]);
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(3e4);
		grains[i].changeDensity(density);
	}
	
	// remove grains
	for(size_t i = 0; i < ngrains; i++) {
		if (grains[i].getPosition()(0) > 370 || grains[i].getPosition()(1) > 650) {
			std::swap(grains[i], grains.back());
			grains.pop_back();
			i--;
			ngrains = grains.size();
		}
	}
	
	FILE * IDs;
	IDs = fopen ("IDs_disc.dat","w");
	for (size_t i = 0; i < grains.size(); i++) {
		fprintf(IDs, "%lu\n", grains[i].getId() );
	}
	
	// create walls
	vector<Wall2dFixed> walls; // 0bot 1left 2right
	walls.resize(2);
	// bottom wall
	double y = 0; double dir = -1; string type = "horizontal";
	walls[0] = Wall2dFixed(y, dir, type, kn);
	// left wall
	double x = 0; dir = -1; type = "vertical"; 
	walls[1] = Wall2dFixed(x, dir, type, kn);
	// right wall
//	dir = 1; // 453.904
//	walls[3] = Wall2dFixed(width, dir, type, kn);
	
	// create world
	World2dBiaxial world(grains, walls, gDamping, dt);
	grains.clear();
	
	
	FILE * positions;
	FILE * rotations;
	positions = fopen ("positions_disc_cd.dat","w");
	rotations = fopen ("rotations_disc_cd.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	UtkContactInfo info(world.getGrains().size());
	
	// time integration for biaxial
	
	// time integration for isotropic
	cout << "Starting drop stage" << endl;
	for (int t = 0; t < tdrop; t++) {
		// non output timestep
		world.computeWorldState();
		if (t % 1000 != 0) {
			
		}
		// output timestep
		else {
			// print positions and velocities to file
			grainList = world.getGrains();
			double ke = 0;
			for (size_t i = 0; i < world.getGrains().size(); i++) {
				fprintf(positions, "%4.3f %4.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
				fprintf(rotations, "%4.3f\n", grainList[i].getTheta());
				ke += grainList[i].computeKineticEnergy();
			}
			duration = omp_get_wtime() - start;
			printf( "drop timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tdrop, 100*double(t+1)/double(tdrop), duration/60., duration/60.*tdrop/(t+1) - duration/60.);
			fflush(stdout);
			cout << "total kinetic energy = " << ke << endl;
		}
		world.addGravityForce(grav);
		world.takeTimestep();
	}
	
	
	
	
	
	
	
	fclose(positions);
	fclose(rotations);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %6.2fs (%6.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!";
	return 0;
}
