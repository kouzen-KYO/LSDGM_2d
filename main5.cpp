/*
 * main.cpp
 *
 *  Created on: May 30, 2014
 *      Author: Reid
 */
 
#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
#include "Wall2dFixed.h"
#include "readInputFile.h"
#include "World2dBiaxial.h"

// put the grains in the box
int main() {
	// file name
	string filename("discShearStart.dat");
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	// grain properties
	double mu = 0.0;
	double kn = 30000;
	double density = 1;
	// assembly properties
	double width = 380; 
	// world properties
	double gDamping = 1.0;
	double grav = .1;
	double dt = 0.2*(2*sqrt(.5*113/kn));
	int tgrav = 5e4;
	
	vector<Grain2d> grains2 = generateGrainsFromFile(filename);
	for (size_t i = 0; i < grains2.size(); i++) {
//		grains2[i].changeMu(mu);
//		grains2[i].changeKn(3e4);
//		grains2[i].changeDensity(density);
		grains2[i].changePos( grains2[i].getPosition() + Vector2d(0, 270) );
	}
	vector<Grain2d> grains3 = generateGrainsFromFile(filename);
	for (size_t i = 0; i < grains2.size(); i++) {
//		grains3[i].changeMu(mu);
//		grains3[i].changeKn(3e4);
//		grains3[i].changeDensity(density);
		grains3[i].changePos( grains3[i].getPosition() + Vector2d(0, 2*270) );
	}
	
	grains.insert(grains.end(), grains2.begin(), grains2.end() );
	grains.insert(grains.end(), grains3.begin(), grains3.end() );
	
	/* set things up */
	for (size_t i = 0; i < grains.size(); i++) {
		grains[i].changeMu(mu);
		grains[i].changeKn(3e4);
		grains[i].changeDensity(density);
	}
	
	// create walls
	vector<Wall2dFixed> walls; // 0bot 1left 2right
	walls.resize(4);
	// bottom wall
	double y = 0; double dir = -1; string type = "horizontal";
	walls[0] = Wall2dFixed(y, dir, type, kn);
	// left wall
	double x = 0; dir = -1; type = "vertical"; 
	walls[2] = Wall2dFixed(x, dir, type, kn);
	// right wall
	dir = 1; // 453.904
	walls[3] = Wall2dFixed(width, dir, type, kn);
	
	// create world
	World2dBiaxial world(grains, walls, gDamping, dt);
	grains.clear();
	
	
	FILE * positions;
	FILE * rotations;
	positions = fopen ("positions_disc.dat","w");
	rotations = fopen ("rotations_disc.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	UtkContactInfo info(world.getGrains().size());
	
	// time integration for biaxial
	double bstart = omp_get_wtime();
	
	for (int t = 0; t < tgrav; t++) {
		
		// non output timestep
		world.computeWorldState();
		// output timestep
		if (t % 2000 == 0) {
			duration = omp_get_wtime() - bstart;
			printf( "Drop timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tgrav, 100*double(t+1)/double(tgrav), duration/60., duration/60.*tgrav/(t+1) - duration/60.);
			fflush(stdout);
			
			
			// print positions and velocities to file
			grainList = world.getGrains();
			double ke = 0;
			for (size_t i = 0; i < world.getGrains().size(); i++) {
//				fprintf(positions, "%4.3f %4.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
//				fprintf(rotations, "%4.3f\n", grainList[i].getTheta());
				ke += grainList[i].computeKineticEnergy();
			}
			cout << "kinetic energy: " << ke << endl;
			

		}
		
		world.addGravityForce(grav);
		world.takeTimestep();
		
	}
	
	
	for (size_t i = 0; i < world.getGrains().size(); i++) {
		fprintf(positions, "%4.3f %4.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
		fprintf(rotations, "%4.3f\n", grainList[i].getTheta());
	}
	
	
	
	
	
	
	
	fclose(positions);
	fclose(rotations);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %6.2fs (%6.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!";
	return 0;
}
