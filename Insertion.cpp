#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <omp.h>

#include "Functions.h"

using namespace std;

int Functions::insert_particles(int& ins_choice, float& ins_reg_low_x, float& ins_reg_up_x, float& ins_reg_low_y, float& ins_reg_up_y, float& ins_reg_low_z, float& ins_reg_up_z, float& rad_p1, int& no_particles, int& trialcounterlimit, vector<vector<double>>& P)
{

	if (ins_choice == 1) 
	{
		int insx = floor((ins_reg_up_x - ins_reg_low_x) / (2.0 * rad_p1));
		int insy = floor((ins_reg_up_y - ins_reg_low_y) / (2.0 * rad_p1));
		int insz = floor((ins_reg_up_z - ins_reg_low_z) / (2.0 * rad_p1));
		const int insbin = insx * insy * insz;
		int ibin, trialcounter = 0;
		float  ix_pos, iy_pos, iz_pos;
		//cout << insx << " " << insy << " " << insz;
		if (insbin < no_particles)
		{
			cout << "insertion domain is small to accomodate all particles" << endl;
			exit(0);
		}

		vector<vector<int>> ins_B(insbin, vector<int>(1, 0));// array for insertion


		srand(trialcounter);
       // #pragma omp parallel for default(shared)
		for (int i = 0; i < no_particles; i++)
		{
			ibin = int((float(rand()) / (RAND_MAX+1)) * insbin);

			if (ins_B[ibin][0] == 1)
			{

				trialcounter = trialcounter + 1;
				if (trialcounter >= trialcounterlimit)
				{
					cout << "could not insert all particles; only " << (i) << " particles are inserted" << endl;
					break;
				}
				i--;
			}
			else
			{

				trialcounter = 0;
				ins_B[ibin][0] = 1;

				if (ibin == 0 || ibin==1)
				{
					ix_pos = rad_p1;
					iy_pos = rad_p1;
					iz_pos = rad_p1;
					ins_B[0][0] = 1;
					ins_B[1][0] = 1;
				}
				else if ((ibin % (insy * insx)) == 0)
				{
					ix_pos = ins_reg_low_x + (insx - 0.5) * 2.0 * rad_p1;
					iy_pos = ins_reg_low_y + (insy - 0.5) * 2.0 * rad_p1;
					iz_pos = ins_reg_low_z + ((ibin / (insy * insx)) - 0.5) * 2.0 * rad_p1;
				}
				else if (((ibin % (insy * insx)) % insx) == 0)
				{
					ix_pos = ins_reg_low_x + (insx - 0.5) * 2.0 * rad_p1;
					iy_pos = ins_reg_low_y + (((ibin % (insy * insx)) / insx) - 0.5) * 2.0 * rad_p1;
					iz_pos = ins_reg_low_z + (floor(ibin / (insy * insx)) + 0.5) * 2.0 * rad_p1;
				}
				else
				{
					ix_pos = ins_reg_low_x + (((ibin % (insy * insx)) % insx) - 0.5) * 2.0 * rad_p1;
					iy_pos = ins_reg_low_y + (floor((ibin % (insy * insx)) / insx) + 0.5) * 2.0 * rad_p1;
					iz_pos = ins_reg_low_z + (floor(ibin / (insy * insx)) + 0.5) * 2.0 * rad_p1;

				}

	
				P[i][1] = 1;
				P[i][2] = ix_pos;
				P[i][3] = iy_pos;
				P[i][4] = iz_pos;


			}

		}

	}
	else if (ins_choice == 2)
	{
		int insx = floor((ins_reg_up_x - ins_reg_low_x) / (2.0 * rad_p1));
		int insy = floor((ins_reg_up_y - ins_reg_low_y) / (2.0 * rad_p1));
		int insz = floor((ins_reg_up_z - ins_reg_low_z) / (2.0 * rad_p1));
		const int insbin = insx * insy * insz;
		int trialcounter = 0;


		if (insbin < no_particles)
		{
			cout << "insertion domain is small to accomodate all particles" << endl;
			exit(0);
		}
		vector<vector<vector<vector<int>>>> ins_B(insx, vector<vector<vector<int>>>
			(insy, vector<vector<int>>
				(insz, vector<int>
					(1, 0))));
	
		//cout << " vector initialized" << endl;

		srand(no_particles); 
	
        //#pragma omp parallel for default(shared)
		for (int q = 0; q < no_particles; q++)
		{
			//cout << q << endl;

			int i = floor((float(rand()) / (RAND_MAX + 1)) * (insx));
			int j = floor((float(rand()) / (RAND_MAX + 1)) * (insy));
			int k = floor((float(rand()) / (RAND_MAX + 1)) * (insz));
			//cout << i << " " << j << " " << k << endl;

			if (ins_B[i][j][k][0] == 1)
			{

				trialcounter = trialcounter + 1;
				if (trialcounter >= trialcounterlimit)
				{
					cout << "could not insert all particles; only " << (q) << " particles are inserted" << endl;
					break;
				}
				q--;
			}
			else
			{
				trialcounter = 0;
				ins_B[i][j][k][0] = 1;

		
				P[q][1] = 1;
				P[q][2] = ins_reg_low_x + i * 2.0 * rad_p1 + 0.5 * 2.0 * rad_p1;
				P[q][3] = ins_reg_low_y + j * 2.0 * rad_p1 + 0.5 * 2.0 * rad_p1;
				P[q][4] = ins_reg_low_z + k * 2.0 * rad_p1 + 0.5 * 2.0 * rad_p1;


			}

		}


	}
	return(0);
}