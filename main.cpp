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

int main(int argc, char *argv[]) { 
	
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
   
	cout << "Running m_shear with mu = " << argmu << endl;
	// file name
	string loadfile("martianInitial.dat");
	vector<Grain2d> grains = generateGrainsFromFile(loadfile);
	size_t ngrains = grains.size();
	// grain properties
	double mu = argmu;
	double kn = 30000;
	double density = 1;
	// assembly properties
	double width = 2000; 
	// world properties
	double gDamping = 0.3;
	double dt = 0.2*(2*sqrt(.5*200/kn));
	int tiso = 20001;
	int tshear = 500001;
	
	vector<Vector2d> pos = readPositionFile("positions_martian_2.dat",ngrains);
	vector<double> rot = readRotationFile("rotations_martian_2.dat",ngrains);
	
	/* set things up */
	for (size_t i = 0; i < grains.size(); i++) {
		grains[i].changePos(pos[i]);
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(kn);
		grains[i].changeDensity(density);
	}
	
//	FILE * IDs;
//	IDs = fopen ("IDs_disc.dat","w");
//	for (size_t i = 0; i < grains.size(); i++) {
//		fprintf(IDs, "%lu\n", grains[i].getId() );
//	}
	
	// create walls
	vector<Wall2d> walls; // 0bot 1left 2right
	walls.resize(2);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, 0.9, INT_MAX - 3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,1200), Vector2d(0,-1), kn, 0.9*kn, 0.9, INT_MAX - 3);
	
	Wall2dPeriodicx wallPeriodic(0, width);
	
	// create world
	World2dBiaxial world(grains, walls, wallPeriodic, gDamping, dt);
	grains.clear();
	
	double pressure = 40;
	
	stringstream folder, filename;
	folder << "martian_mu=" << mu << "/";
	mkdir(folder.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	FILE * positions;
	FILE * rotations;
	FILE * wallPos;
	FILE * stress;
	filename << folder.str().c_str() << "positions.dat";
	positions = fopen(filename.str().c_str(),"w");
	filename.str(std::string()); filename << folder.str().c_str() << "rotations.dat";
	rotations = fopen(filename.str().c_str(),"w");
	filename.str(std::string()); filename << folder.str().c_str() << "wallPos.dat";
	wallPos   = fopen(filename.str().c_str(),"w");
	filename.str(std::string()); filename << folder.str().c_str() << "stress.dat";
	stress    = fopen(filename.str().c_str(),"w");
	
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
			double volume = world.getWallPeriodic().getWidth()*world.getWalls()[1].getPosition()(1);
			cout << "volume = " << volume << endl;
			Vector3d stress = world.getWorldState().stressVoigt/volume;
			cout << "stress = " << stress.transpose() << endl;
		}
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWallPeriodic().getWidth();
		double outForce = ws.wallForces[1](1);
		if ( .99*fabs(world.getWorldState().stressVoigt(1)) > fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(-0.001);
		}
		if ( 1.01*fabs(world.getWorldState().stressVoigt(1)) < fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(0.001);
		}
		
		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
		world.takeTimestep(); 
	}
	
	world.getWallsNC()[1].changeVelocity(Vector2d(.002/dt,0) );
	world.getWallsNC()[0].changeVelocity(Vector2d(-.002/dt,0) );
	start = omp_get_wtime(); 
	// time integratio for shearing
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
			fname << folder.str().c_str() << "cstate_" << iter << ".dat";
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
			
			double volume = world.getWallPeriodic().getWidth()*world.getWalls()[1].getPosition()(1);
			Vector3d stress = world.getWorldState().stressVoigt/volume;
			cout << "stress = " << stress.transpose() << endl;
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
	
	
	
	
	
	
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!";
	return 0;
}
