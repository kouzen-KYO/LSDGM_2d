/*
 * main.cpp
 *
 *  Created on: May 30, 2014
 *      Author: Reid
 */
 
#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
#include "readInputFile.h"
#include "World2dBiaxial.h"

// drop test for spheres (run after main5)
int main() {
	
	// argv: friction
	double argmu;
	if (argc == 2) {
		argmu = atof(argv[1]);
	}
	else {
		argmu = 0.4;
	}
	
//   omp_set_dynamic(0);     // Explicitly disable dynamic threading
//   omp_set_num_threads(6);
   
	// file name
	string filename("discShearStart.dat");
	vector<Grain2d> grains0 = generateGrainsFromFile(filename); //TODO: fix indices lol
	vector<Grain2d> grains;
	
//	for (size_t i = 0; i < grains0.size(); i++) {
//		cout << grains0[i].getPosition().transpose() << endl;
//	}
	
	for (size_t i = 0; i < 9; i++) {
		grains.insert(grains.end(), grains0.begin(), grains0.end());
	}

	size_t ngrains = grains.size(); 
	// grain properties
	double mu = argmu;
	double wallmu = 0.9;
	double kn = 30000;
	double density = 1;
	// assembly properties
	double width = 322*3; 
	double height = 280*3;
	// world properties
	double gDamping = 0.3;
	double dt = 0.2*(2*sqrt(.5*200/kn));
	int tiso = 20000;
	int tshear = 1;
	
	/* set things up */
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			Vector2d translation(i*322., j*280.);
			size_t idx = i*3+j;
			for (size_t k = 0; k < 384; k++) {
				grains[k+384*idx].changePos( grains[k+384*idx].getPosition() +translation);
				grains[k+384*idx].changeMu(mu);
				grains[k+384*idx].changeKn(kn);
				grains[k+384*idx].changeDensity(density);
			}
		}
	}
	
	
	// create walls
	vector<Wall2d> walls; // 0bot 1left 2right
	walls.resize(2);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, wallmu, INT_MAX - 3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,height), Vector2d(0,-1), kn, 0.9*kn, wallmu, INT_MAX - 3);
	
	Wall2dPeriodicx wallPeriodic(0, width);
	
	// create world
	World2dBiaxial world(grains, walls, wallPeriodic, gDamping, dt);
	grains.clear();
	
	double pressure = 40;
	
	FILE * positions;
	FILE * rotations;
	FILE * wallPos;
	FILE * stress;
	positions = fopen ("positions_disc_3.dat","w");
	rotations = fopen ("rotations_disc_3.dat","w");
	wallPos   = fopen("wallPos_disc_3.dat","w");
	stress    = fopen("stress_disc_3.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	UtkContactInfo info(world.getGrains().size());
	
	// time integration for isotropic compression
	size_t iter = 0;
	cout <<"Starting iso stage" << endl;
	for (int t = 0; t < tiso; t++) {
		// non output timestep
		world.computeWorldState();
		if (t % 1000 != 0) {
			
		}
		// output timestep
		else {
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso, 100*double(t+1)/double(tiso), duration/60., duration/60.*tiso/(t+1) - duration/60.);
			fflush(stdout);
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/world.getWallPeriodic().getWidth()/world.getWalls()[1].getPosition()(1) << endl;
			double ke = 0;
			grainList = world.getGrains();
			for (size_t i = 0; i < world.getGrains().size(); i++) {
				ke += grainList[i].computeKineticEnergy();
			}
			cout << "total kinetic energy = " << ke << endl;
		}
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWallPeriodic().getWidth();
		double outForce = ws.wallForces[1](1);
		if ( .99*fabs(world.getWorldState().stressVoigt(1)) > fabs(world.getWorldState().stressVoigt(0))) {
//			cout << "hi" << endl;
			world.getWallPeriodicNC().moveWidth(-0.0005);
		}
		if ( 1.01*fabs(world.getWorldState().stressVoigt(1)) < fabs(world.getWorldState().stressVoigt(0))) {
//			cout << "hi" << endl;
			world.getWallPeriodicNC().moveWidth(0.0005);
		}
		
		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
		world.takeTimestep();
	}
	
	// set velocities
	world.getWallsNC()[1].changeVelocity(Vector2d(.001/dt,0) );
	world.getWallsNC()[0].changeVelocity(Vector2d(-.001/dt,0) );
	start = omp_get_wtime();
	// time integration for shearing
	cout << "Starting shear stage" << endl;
	for (int t = 0; t < tshear; t++) {
		// non output timestep
		world.computeWorldState();
		if (t % 1000 != 0) {
			
		}
		// output timestep
		else {
			
			CData cState = world.computeCstate();
			stringstream fname;
			fname << "cdata/cstate_disc_3_" << iter << ".dat";
			iter++;
			FILE * cinfo = fopen(fname.str().c_str(), "w");
			
			for (size_t i = 0; i < cState._clocs.size(); i++) {
				fprintf(cinfo, "%d %d ",cState._cpairs[i](0), cState._cpairs[i](1) ); // grains
				fprintf(cinfo, "%.2f %.2f ",cState._forces[i](0), cState._forces[i](1)); // forces
				fprintf(cinfo, "%.3f %.3f ",cState._normals[i](0), cState._normals[i](1)); // normals
				fprintf(cinfo, "%.2f %.2f ",cState._clocs[i](0), cState._clocs[i](1)); // locations
				fprintf(cinfo, "%lu\n", cState._nodes[i]); // node of _cpairs[i](0)
			}
			fclose(cinfo);
			
			
			// print positions and velocities to file
			grainList = world.getGrains();
			double ke = 0;
			for (size_t i = 0; i < world.getGrains().size(); i++) {
				fprintf(positions, "%.3f %.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
				fprintf(rotations, "%.3f\n", grainList[i].getTheta());
				ke += grainList[i].computeKineticEnergy();
			}
			fprintf(stress, "%.3f %.3f %.3f\n", world.getWorldState().stressVoigt(0), world.getWorldState().stressVoigt(1), world.getWorldState().stressVoigt(2));
			fprintf(wallPos,"%.3f %.3f\n", world.getWallPeriodic().getWidth(), world.getWalls()[1].getPosition()(1));
			duration = omp_get_wtime() - start;
			printf( "shear timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tshear, 100*double(t+1)/double(tshear), duration/60., duration/60.*tshear/(t+1) - duration/60.);
			fflush(stdout);
			cout << "total kinetic energy = " << ke << endl;
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/world.getWallPeriodic().getWidth()/world.getWalls()[1].getPosition()(1) << endl;
		}
		
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWallPeriodic().getWidth();
		double outForce = ws.wallForces[1](1);
		
		if ( .99*fabs(world.getWorldState().stressVoigt(1)) > fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(-0.0005);
		}
		if ( 1.01*fabs(world.getWorldState().stressVoigt(1)) < fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(0.0005);
		}

		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
		world.takeTimestep();
	}
	
	
	
	
	
	
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!";
	return 0;
}
