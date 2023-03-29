#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>   

using namespace std;
void check_quads(float point_x,float point_y,float point_z, float shortest_dist, int quad, int which_quad[],float minX, float minY, float minZ, float maxX, float maxY, float maxZ);
int find_quad(int point_x, int point_y, int point_z);
void find_min(double x[],double y[],double z[],  int index1[],int num1,double point_x, double point_y, double point_z, double inflated_sum[], double rad_2_sum[], float min_dist[]);

int main(int argc, char **argv){
	clock_t begin = clock();
	// ./volume_c pdb_name num_residues
	if (argc < 4){
		cout << " Invalid input " << endl ;
		cout << endl;
		return 2;
		
	} //check args
	 
	//Parse file name
	char *file_name_short = 0;
	file_name_short = argv[1];
	char file_name[10];
	strcpy(file_name, file_name_short);
	cout << file_name_short << endl;
	
	//Parse number of residues
	int size_pdb = atoi(argv[2]);
	char *item_num = argv[3];
	//Set random number seed so different every time you run
	srand48(int(time(NULL)/atof(argv[3])/atof(argv[3])));



	ifstream ifs; //stream of input file
	double probe = 1.4; //probe size
	double x[size_pdb+1];  //X coordinates
	double y[size_pdb+1]; // Y coordinates
	double z[size_pdb+1]; // Z coordinates
	double sizes[size_pdb+1]; // residues sizes
	double inflated_sum[size_pdb+1]; // (size + probe)^2
	double rad_2_sum[size_pdb+1]; // size^2

	// reading file variables
	char atom_type[5];
	char Res[4];
	double id;

	
	//protein data
	double minX= 0.0, minY= 0.0, minZ= 0.0, maxX= 0.0, maxY= 0.0, maxZ= 0.0;
	double rangeX= 0.0, rangeY = 0.0, rangeZ = 0.0;
	double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
	
	// mark all atoms as core
	bool is_edge[size_pdb+1];
	for (int a = 0; a <= size_pdb; a++) {
		is_edge[a] = false;
	}
		

	// initialize all counts to 0 and all sizes to 0
	double vol_count[size_pdb+1];
	for (int a = 0; a <= size_pdb; a++) {
		vol_count[a] = 0.0;
		sizes[a] = 0.0;
	}
	
	

	// Open file to read from
	char name_2[150] = "";
	strcat(name_2, file_name_short);
	
	cout << name_2 << endl;
	strcat(name_2,  ".txt");
	cout << name_2 << endl;
	ifs.open(name_2);
	cout << size_pdb << endl;
	cout << ifs.good() << endl;
	//Read each line and store in x, y, z
	//Also find min and max for each dimension
	for(int atoms = 0 ; atoms < size_pdb; atoms ++){
		ifs >> atom_type >> Res >> id >> sizes[atoms] >> x[atoms] >> y[atoms] >> z[atoms];
		if (atoms == 0){
			cout << atom_type << Res << id << sizes[atoms] << x[atoms] << y[atoms] << z[atoms] << endl;
			minZ = z[atoms];
			maxZ = z[atoms];
			minY = y[atoms];
			maxY = y[atoms];
			minX = x[atoms];
			maxX = x[atoms];
			}
		if (x[atoms] > maxX)
			maxX = x[atoms];
		else if (x[atoms] < minX)
			minX = x[atoms];
		if (y[atoms] > maxY)
			maxY = y[atoms];
		else if (y[atoms] < minY)
			minY = y[atoms];
		if (z[atoms] > maxZ)
			maxZ = z[atoms];
		else if (z[atoms] < minZ)
			minZ = z[atoms];	
		
		meanX = meanX+x[atoms];
		meanY = meanY+y[atoms];
		meanZ = meanZ+z[atoms];
		//Calculate distance^2 with and without probe
		rad_2_sum[atoms] = sizes[atoms]*sizes[atoms];
		inflated_sum[atoms] = (sizes[atoms]+probe)*(sizes[atoms]+probe);
	} // atoms < size_pdb
	
	ifs.close();
	
	meanX = meanX/size_pdb;
	meanY = meanY/size_pdb;
	meanZ = meanZ/size_pdb;
	
	//Move everything so center of mass in at origin. Then partition to quadrants
	int num1 = 0, num2 = 0, num3 = 0, num4 = 0, num5 = 0, num6 = 0, num7=0, num8=0;
	for(int atoms = 0 ; atoms < size_pdb; atoms ++){
		x[atoms] = x[atoms]-meanX;
		y[atoms] = y[atoms]-meanY;
		z[atoms] = z[atoms]-meanZ;
		if (x[atoms] <= 0) {
			if (y[atoms] <= 0){
				if (z[atoms] <= 0)
					num7 = num7+1;
				else
					num8 = num8+1;
				}
			else {
				if (z[atoms] <= 0)
					num6 = num6+1;
				else
					num5 = num5+1;
				}
			}
		else {
			if (y[atoms] <= 0){
				if (z[atoms] <= 0)
					num3 = num3+1;
				else
					num4 = num4+1;
				}
			else {
				if (z[atoms] <= 0)
					num2 = num2+1;
				else
					num1 = num1+1;
				}
			}
		}
	
	int index1[num1+1], index2[num2+1], index3[num3+1], index4[num4+1];
	int index5[num5+1], index6[num6+1], index7[num7+1], index8[num8+1];
	int c1=0, c2=0,c3=0, c4=0, c5=0, c6=0, c7=0, c8=0;
	
	for(int atoms = 0 ; atoms < size_pdb; atoms ++){
		if (x[atoms] <= 0) {
			if (y[atoms] <= 0){
				if (z[atoms] <= 0){
					index7[c7]= atoms;
					c7 = c7+1;
					}
				else{
					index8[c8]=atoms;
					c8 = c8+1;
					}
				}
			else {
				if (z[atoms] <= 0){
					index6[c6] = atoms;
					c6 = c6+1;
					}
				else{
					index5[c5] = atoms;
					c5 = c5+1;
					}
				}
			}
		else {
			if (y[atoms] <= 0){
				if (z[atoms] <= 0){
					index3[c3] = atoms;
					c3 = c3+1;
					}
				else{
					index4[c4] = atoms;
					c4 = c4+1;
					}
				}
			else {
				if (z[atoms] <= 0){
					index2[c2] = atoms;
					c2 = c2+1;
					}
				else{
					index1[c1] = atoms;
					c1 = c1+1;
					}
				}
			}
		}
		
		
	
	
	float orig_maxX = maxX;
	float orig_maxY = maxY;
	float orig_maxZ = maxZ;
	float orig_minX = minX;
	float orig_minY = minY;
	float orig_minZ = minZ;
	//Increase and decrease min and max by size of probe and size of largest atom
	//This allows for points outside of protein to get sampled
	minX = minX - meanX -probe*5.0;
	minY = minY - meanY - probe*5.0;
	minZ = minZ - meanZ -probe*5.0;
	maxX = maxX -meanX + probe*5.0;
	maxY = maxY -meanY + probe*5.0;
	maxZ = maxZ -meanZ + probe*5.0;
	rangeX = maxX-minX;
	rangeY = maxY-minY;
	rangeZ = maxZ-minZ;

	double box_vol = rangeX*rangeY*rangeZ;
	cout << rangeX << " " << rangeY << " " << rangeZ << " range" << endl;
	double point_x, point_y, point_z; //random points in box
	double dist_x, dist_y, dist_z; //Distance between random point and an atom in x, y, z
	bool in_sphere = false; //random point is in a sphere
	bool in_range = false; // random point is within probe of sphere surface
	double distemp_old = 100000.1; // distance between random point and atom
	int loc_min = 0; //location of atom that point is closest to
	double clash_min = 100.1;
	int clash_loc = 0;
	double dt, ci, cl;
	double resolution = 0.2;
	int XDIM = rangeX/resolution;
	int YDIM = rangeY/resolution;
	int ZDIM = rangeZ/resolution;
	// int toutalnumofpoints = XDIM * YDIM * ZDIM;
	int toutalnumofpoints = 900000000;
	cout<<resolution<<endl;
	ofstream pointsout;

	//Generate a ton of points and see where they fall
	double num_points = 0.0;
	int quad;
	float min_dist[6];
	while (num_points < toutalnumofpoints){//
		num_points++;
		point_x= drand48()*rangeX+minX;
		point_y= drand48()*rangeY+minY;
		point_z = drand48()*rangeZ+minZ;
	
		//Figure out what quadrant you are in
		quad = 0;
		quad = find_quad(point_x, point_y, point_z);
	// Change to pass combined data for all things it could be near
	
	
		switch (quad){
			case 1:
				find_min(x,y,z,index1,num1,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 2:
				find_min(x,y,z,index2,num2,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 3:
				find_min(x,y,z,index3,num3,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 4:
				find_min(x,y,z,index4,num4,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 5:
				find_min(x,y,z,index5,num5,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 6:
				find_min(x,y,z,index6,num6,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 7:
				find_min(x,y,z,index7,num7,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
			case 8:
				find_min(x,y,z,index8,num8,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist);
				break;
				
			}
		
		
		// Now check to see if you need to send another quadrant
		float shortest_dist = min_dist[4];
		int need_to_run[8];
		
		check_quads(point_x,point_y,point_z, shortest_dist, quad ,need_to_run, minX, minY, minZ, maxX, maxY, maxZ);
		float min_dist1[6];
		min_dist1[4] = 100000;
		for (int run_quad = 0; run_quad < 8; run_quad ++){
			if (need_to_run[run_quad] == 1){
				switch (run_quad+1){
					case 1:
						 find_min(x,y,z,index1,num1,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 2:
						 find_min(x,y,z,index2,num2,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 3:
						 find_min(x,y,z,index3,num3,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 4:
						 find_min(x,y,z,index4,num4,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 5:
						 find_min(x,y,z,index5,num5,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 6:
						 find_min(x,y,z,index6,num6,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 7:
						 find_min(x,y,z,index7,num7,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
					case 8:
						 find_min(x,y,z,index8,num8,point_x, point_y, point_z, inflated_sum, rad_2_sum, min_dist1);
						break;
				
					}
					if (min_dist1[2] > 0){
						min_dist[2] = 1;
						if (min_dist1[5] < min_dist[5]){
							min_dist[5] = min_dist1[5];
							min_dist[3] = min_dist1[3];
							}
						}
					if (min_dist1[0] > 0)
						min_dist[0] = 1;
					if (min_dist1[4] < min_dist[4]){
						min_dist[4] = min_dist1[4];
						min_dist[1] = min_dist1[1];
						}
					
				}
		
		}
		//cout << min_dist[4] << " closest dist" << endl;
	 	// if the point wasn't within 1.5A of any atoms, it is an edge
	 	//cout << static_cast<int>(min_dist[1]) << " dist1" << endl;
		 if (min_dist[0]<0) {
		 	if (min_dist[1] >= 0) {
				is_edge[static_cast<int>(min_dist[1])] = true;
			}
			else {
				//cout << "problem with location" << endl;
			}
		}
		 if (min_dist[2]>0) // If it's in a sphere, increase count of that atoms
			vol_count[static_cast<int>(min_dist[3])] ++;
	}	

cout <<"finished running" << endl;

ofstream outs1; //stream to save to

char name_3[150] = "";
//strcat(name_3, "/Users/jennifergaines/Documents/summer2013/Minimum_energy_algorithm/Core_mutations/");
//strcat(name_3, "/home/fas/ohern/jcg72/volumes/Mutations/");
strcat(name_3, file_name_short);
//strcat(name_3, "1tsr_h.pdb");
strcat(name_3, item_num);
strcat(name_3, "_vol.txt");
cout << name_3 << endl;
outs1.open(name_3);
cout << outs1.is_open() << endl;
//cout << name_3 << endl;
if (outs1.is_open()) {
	cout << name_3 << endl;
	for (int i = 0; i < size_pdb; i ++){
		outs1 << sizes[i] << " " << (is_edge[i] ? 1 : 0) << " " << double(vol_count[i])/toutalnumofpoints*box_vol << endl;
	}
	outs1.close();
	cout << " saved properly" << endl;
}
else {
	cout << "problem with file" << endl;
}

clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "time: " << elapsed_secs << endl;
}

int find_quad(int point_x, int point_y, int point_z){
	int quad = 0;
	if (point_x <= 0) {
		if (point_y <= 0){
			if (point_z <= 0)
				quad = 7;
			else
				quad = 8;
			}
		else {
			if (point_z <= 0)
				quad = 6;
			else
				quad = 5;
			}
		}
	else {
		if (point_y <= 0){
			if (point_z <= 0)
				quad = 3;
			else
				quad = 4;
			}
		else {
			if (point_z <= 0)
				quad = 2;
			else
				quad = 1;
			}
		}
	return quad;
}

void find_min(double x[],double y[],double z[], int index1[],int num1,double point_x, double point_y, double point_z, double inflated_sum[], double rad_2_sum[], float min_dist[]){
	//float min_dist[5];
	bool in_sphere = false; //random point is in a sphere
	bool in_range = false; // random point is within probe of sphere surface
	double distemp_old = 100000.1; // distance between random point and atom
	int loc_min = 0; //location of atom that point is closest to
	double clash_min = 100000000.1;
	int clash_loc = -1;
	double dt, ci, cl;
	double dist_x, dist_y, dist_z; //Distance between random point and an atom in x, y, z

	in_sphere = false;
	in_range = false;
	loc_min = -1;
	int this_ind = 0;
	//loop over all atoms in pdb file and see if the point is in one of them
	for (int i = 0; i < num1; i++){
		this_ind = index1[i];
		dist_x = x[this_ind] - point_x;
		dist_y = y[this_ind] -point_y;
		dist_z = z[this_ind] -point_z;
		//calculate distance between point and atom
		dt = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
		ci = dt -inflated_sum[this_ind]; //if neg, point is inside an inflated sphere
		cl = dt - rad_2_sum[this_ind]; // if neg, point is in an atom

		if (cl < 0) { // If point is in an atom
			in_sphere = true;
			if (cl < clash_min){ // see if closer to center of this atom than to other atoms
				clash_loc = i;
				clash_min = cl;
			} //cl < clash_min
		} // cl<0
		
		if (ci < 0) // if point is inside inflated atom
			in_range = true;
		if (in_sphere && !in_range)
			cout << "Problem" << endl;
			
		// check to if distance to this atom is smaller than to others
		// if so, store so that this can be marked if in_range is false after all atoms
		if (dt < distemp_old) { 
			loc_min = i;
			distemp_old = dt;
			} // dt < distemp_old
	 } // for i < size_pdb
 
	if (in_range)
		min_dist[0] = 1.0;
	else
		min_dist[0] = -1.0;
	min_dist[1] = index1[loc_min];
	if (in_sphere)
		min_dist[2]  = 1.0;
	else
		min_dist[2] = -1.0;
	min_dist[3] = index1[clash_loc];
	min_dist[4] = distemp_old;
	min_dist[5] = clash_min;

	 }
	 
	 
	 
void  check_quads(float point_x,float point_y,float point_z, float shortest_dist, int quad, int which_run[], float minX, float minY, float minZ, float maxX, float maxY, float maxZ){
	for (int i = 0; i< 8; i++){
		which_run[i] = -1;
	}
	
	float to_check_X, to_check_Y, to_check_Z;
	if (maxX > (-minX))
		to_check_X = maxX;
	else
		to_check_X = -minX;
	if (maxY > (-minY))
		to_check_Y = maxY;
	else
		to_check_Y = -minY;
	if (maxZ > (-minZ))
		to_check_Z = maxZ;
	else
		to_check_Z = -minZ;
		
		
	if (point_x < to_check_X)
		to_check_X = point_x;
	if (point_y < to_check_Y)
		to_check_Y = point_y;
	if (point_z < to_check_Z)
		to_check_Z = point_z;
		
	float new_check_val = pow(sqrt(shortest_dist)+1.8,2) ;
	switch (quad){
		// 1.8 is added to shortest_dist (same as subtracted to point*point) to account 
		//for the possibility of an atom extending into the current quadrent, when the center is in the next
		case 1:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val) //
				which_run[5-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[4-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[2-1] = 1;
			if (which_run[5-1]>0 && which_run[4-1]>0)//(point_x*point_x < new_check_val&& point_y*point_y < new_check_val) //which_run[5-1]>0 && which_run[4-1]>0)
				which_run[8-1]=1;
			if (which_run[5-1]>0 && which_run[2-1]>0) //(point_x*point_x < new_check_val&& point_z*point_z < new_check_val)//
				which_run[6-1]=1;
			if (which_run[4-1]>0 && which_run[2-1]>0)//(point_z*point_z < new_check_val&& point_y*point_y < new_check_val)//(which_run[4-1]>0 && which_run[2-1]>0)
				which_run[3-1]=1;
			if (which_run[5-1]>0 && which_run[4-1]>0 && which_run[2-1]>0)//(point_x*point_x < new_check_val && point_y*point_y < new_check_val&& point_z*point_z < new_check_val)//
				which_run[7-1]=1;
			break;
		case 2:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[6-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val)
				which_run[3-1]  = 1;
			if (point_z*point_z  + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val)
				which_run[1-1] = 1;
			if (which_run[6-1]>0 && which_run[3-1]>0)
				which_run[7-1]=1;
			if (which_run[6-1]>0 && which_run[1-1]>0)
				which_run[5-1]=1;
			if (which_run[3-1]>0 && which_run[1-1]>0)
				which_run[4-1]=1;
			if (which_run[6-1]>0 && which_run[3-1]>0 && which_run[1-1]>0)
				which_run[8-1]=1;
			break;
		case 3:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[7-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[2-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[4-1] = 1;
			if (which_run[7-1]>0 && which_run[2-1]>0)
				which_run[6-1]=1;
			if (which_run[7-1]>0 && which_run[4-1]>0)
				which_run[8-1]=1;
			if (which_run[2-1]>0 && which_run[4-1]>0)
				which_run[1-1]=1;
			if (which_run[7-1]>0 && which_run[4-1]>0 && which_run[2-1]>0)
				which_run[5-1]=1;
			break;
		case 4:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[8-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[1-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[3-1] = 1;
			if (which_run[8-1]>0 && which_run[1-1]>0)
				which_run[5-1]=1;
			if (which_run[8-1]>0 && which_run[3-1]>0)
				which_run[7-1]=1;
			if (which_run[1-1]>0 && which_run[3-1]>0)
				which_run[2-1]=1;
			if (which_run[8-1]>0 && which_run[1-1]>0 && which_run[3-1]>0)
				which_run[6-1]=1;
			break;
		case 5:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[1-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[8-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[6-1] = 1;
			if (which_run[1-1]>0 && which_run[8-1]>0)
				which_run[4-1]=1;
			if (which_run[1-1]>0 && which_run[6-1]>0)
				which_run[2-1]=1;
			if (which_run[8-1]>0 && which_run[6-1]>0)
				which_run[7-1]=1;
			if (which_run[1-1]>0 && which_run[8-1]>0 && which_run[6-1]>0)
				which_run[3-1]=1;
			break;
		case 6:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[2-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[7-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[5-1] = 1;
			if (which_run[2-1]>0 && which_run[7-1]>0)
				which_run[3-1]=1;
			if (which_run[2-1]>0 && which_run[5-1]>0)
				which_run[1-1]=1;
			if (which_run[5-1]>0 && which_run[7-1]>0)
				which_run[8-1]=1;
			if (which_run[2-1]>0 && which_run[7-1]>0 && which_run[5-1]>0)
				which_run[4-1]=1;
			break;
		case 7: 
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[3-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[6-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[8-1] = 1;
			if (which_run[3-1]>0 && which_run[6-1]>0)
				which_run[2-1]=1;
			if (which_run[3-1]>0 && which_run[8-1]>0)
				which_run[4-1]=1;
			if (which_run[8-1]>0 && which_run[6-1]>0)
				which_run[5-1]=1;
			if (which_run[8-1]>0 && which_run[6-1]>0 && which_run[3-1]>0)
				which_run[1-1]=1;
			break;
		case 8:
			if (point_x*point_x + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_z-to_check_Z)*(point_z-to_check_Z) < new_check_val)
				which_run[4-1] = 1;
			if (point_y*point_y + (point_z-to_check_Z)*(point_z-to_check_Z) + (point_x-to_check_X)*(point_x-to_check_X)< new_check_val) //
				which_run[5-1]  = 1;
			if (point_z*point_z + (point_y-to_check_Y)*(point_y-to_check_Y) + (point_x-to_check_X)*(point_x-to_check_X) < new_check_val) //
				which_run[7-1] = 1;
			if (which_run[4-1]>0 && which_run[5-1]>0)
				which_run[1-1]=1;
			if (which_run[4-1]>0 && which_run[7-1]>0)
				which_run[3-1]=1;
			if (which_run[7-1]>0 && which_run[5-1]>0)
				which_run[6-1]=1;
			if (which_run[4-1]>0 && which_run[5-1]>0 && which_run[7-1]>0)
				which_run[2-1]=1;
			break;
		}
	//return which_run;
}
				
	
