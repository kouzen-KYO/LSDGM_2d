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
 
int main(int argc, char*argv[]) {
//	omp_set_dynamic(0);     // Explicitly disable dynamic threading
//	omp_set_num_threads(6);
	// file name 
	string filename;
	if (argc == 3) {
		filename = string(argv[1]) + ".dat";
	}
	else {
		filename = "caicos.dat";
	}
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	vector<Grain2d> temp = grains;
	for (size_t i = 0; i < 63; i++) {
		grains.insert(grains.end(), temp.begin(), temp.end());
	}
	size_t ngrains = grains.size();
	cout <<"Number of grains: " << ngrains << endl;
	// grain properties
	double mu;
	if (argc == 3) {
		mu = atof(argv[2]);
	}
	else {
		mu = 0.3;
	}
	cout << "Type of grain: " << filename.c_str() << endl;
	cout << "mu: " << mu << endl;
	double kn = 30000;
	double density = 1;
	// assembly properties
//	vector<double> wallp = readRotationFile("walls_19500.dat",2);
//	double widthInit = wallp[0]; 
//	double heightInit = wallp[1];
	
	// world properties
	double gDamping = .7; // .4 or .2
	double dt = 0.2*(2*sqrt(.5*200/kn));
	
	vector<Vector2d> pos = readPositionFile("positions_caicos_2.dat",ngrains);
	vector<double> rot = readRotationFile("rotations_caicos_2.dat",ngrains);
	
	/* set things up */
	for (size_t i = 0; i < ngrains; i++) {
		grains[i].changePos(pos[i]);
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(kn);
		grains[i].changeKs(.9*kn);
		grains[i].changeDensity(density);
		grains[i].changeId(i);
	}
	
	// create walls
	vector<Wall2d> walls; // 0bot 1top
	walls.resize(4);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, 0.0, INT_MAX-3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,1180), Vector2d(0,-1), kn, 0.9*kn, 0.0, INT_MAX-2);
	// left wall
	walls[2] = Wall2d(Vector2d(-10,0), Vector2d(1,0), kn, 0.9*kn, 0.0, INT_MAX-1);
	// right wall
	walls[3] = Wall2d(Vector2d(710,0), Vector2d(-1,0), kn, 0.9*kn, 0.0, INT_MAX-0);
	
	string dir = "biax_" + string(argv[1]) + "_mu=" + string(argv[2]);
	mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (chdir(dir.c_str())) {
		cout << "Invalid folder" << endl;
	}
	string cinfo("cinfo");
	mkdir(cinfo.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//	FILE * props;
	FILE * positions;
	FILE * rotations;
	FILE * wallPos;
	FILE * stress;
//	props 	 = fopen("props.dat","w");
	positions = fopen("positions_biax.dat","w");
	rotations = fopen("rotations_biax.dat","w");
	wallPos   = fopen("wallPos_biax.dat","w");
	stress    = fopen("stress_biax.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	
	int tiso = 50001;
	int tbiax = 30001;
	double pressure = 20;
	double dx = -0.1;
	
	double strain = 0.01; // negative in y direction, positive in x direction
	
	// time integration for isotropic compression
//	size_t iter = 0;
	cout <<"Starting iso stage" << endl;
	World2dBiaxial world(grains, walls, gDamping, dt);
	double width, height, volume;
	for (int t = 0; t < tiso; t++) {
		world.computeWorldState();
		width = (world.getWalls()[3].getPosition()(0)-world.getWalls()[2].getPosition()(0));
		height = (world.getWalls()[1].getPosition()(1)-world.getWalls()[0].getPosition()(1));
		volume = width*height;
							 
		// non output timestep
		if (t % 1000 != 0) {
		}
		// output timestep
		else {
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso, 100*double(t+1)/double(tiso+tbiax), duration/60., duration/60.*(tiso+tbiax)/(t+1) - duration/60.);
			fflush(stdout);
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
		}
		if ( (fabs(world.getWorldState().stressVoigt(1)) + fabs(world.getWorldState().stressVoigt(0)))/volume < pressure*1.5) {
			// move walls
//			world.getWallsNC()[0].moveWall(Vector2d(0, 0.1));
			world.getWallsNC()[1].moveWall(Vector2d(0,dx));
//			world.getWallsNC()[2].moveWall(Vector2d( 0.1,0));
			world.getWallsNC()[3].moveWall(Vector2d(dx,0));
			// move grains accordingly
			world.moveProportional(dx,dx,width,height);
		}
		else {
			WorldState2d ws = world.getWorldState();
			double inForce = pressure*width;
			double outForce = -world.getWorldState().stressVoigt(1)/volume*width;
			world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
			inForce = pressure*height;
			outForce = -world.getWorldState().stressVoigt(0)/volume*height;
			world.getWallsNC()[3].takeForceTimestep(inForce-outForce, ws.wallContacts[3]);
		}
		world.takeTimestep();
	}
	
	width = (world.getWalls()[3].getPosition()(0)-world.getWalls()[2].getPosition()(0));
	height = (world.getWalls()[1].getPosition()(1)-world.getWalls()[0].getPosition()(1));
	
	dx = width*strain/2/(tbiax-1);
	double dy = -height*strain/2/(tbiax-1);
	size_t station = 0;
	cout <<"Starting biax stage" << endl;
	for (int t = 0; t < tbiax; t++) {
		world.computeWorldState();
		
		// output timestep
		if (t % 1000 == 0) {
			duration = omp_get_wtime() - start;
			printf( "biax timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tbiax, 100*double(t+1)/double(tiso+tbiax), duration/60., duration/60.*(tiso+tbiax)/(t+tiso+1) - duration/60.);
			fflush(stdout);
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
		}
		
		if (t%int((tbiax-1)/10) == 0) {
			// output at steps during biax
			// cinfo
			CData cData = world.computeCstate();
			string fname = cinfo + string("/") + to_string(station) + string(".dat");
			FILE * cinfofile;
			cinfofile = fopen(fname.c_str(), "w");
			for (size_t j = 0; j < cData.size(); j++) {
				fprintf(cinfofile, "%d %d %.3f %.3f %.4f %.4f %.3f %.3f %lu\n", cData._cpairs[j](0), cData._cpairs[j](1), cData._forces[j](0), cData._forces[j](1),
																					 cData._normals[j](0), cData._normals[j](1), cData._clocs[j](0), cData._clocs[j](1), cData._nodes[j]);
			}
			station++;
			// print positions and velocities to file
			grainList = world.getGrains();
			for (size_t i = 0; i < world.getGrains().size(); i++) {
				fprintf(positions, "%.8f %.8f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
				fprintf(rotations, "%.8f\n", grainList[i].getTheta());
			}
			for (size_t i = 0; i < world.getWalls().size(); i++){
				fprintf(wallPos,"%.8f %.8f\n", world.getWalls()[i].getPosition()(0), world.getWalls()[i].getPosition()(1) );
			}
			Vector3d stressVoigt = world.getWorldState().stressVoigt;
			fprintf(stress,"%.4f %.4f %.4f \n", stressVoigt(0)/volume, stressVoigt(1)/volume, stressVoigt(2)/volume);
		}
		
		world.getWallsNC()[0].moveWall(Vector2d(0, -dy));
		world.getWallsNC()[1].moveWall(Vector2d(0,dy));
		world.getWallsNC()[2].moveWall(Vector2d(-dx,0));
		world.getWallsNC()[3].moveWall(Vector2d(dx,0));
		
		world.takeTimestep();
	}
	
	
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!\n";
	return 0;
}
