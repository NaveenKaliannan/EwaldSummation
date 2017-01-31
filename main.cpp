
// main.cpp computes the Madelung constant for the given crystal structure
// of the material. Ewald summation is used to find the constant. The output 
// can be seen in the output/graph folder

#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include "ewald.h"


using namespace std;

const double Pi = 3.1415927;

// Main implementation
int main(int argc, char* argv[])
{
  // Number of unit cell in x, y and z direction - nx ny nz   
  unsigned int nx = 0, ny = 0, nz = 0;
  // Closest distance between a pair atoms - r0
  double r0 = 0;
  // Number of atoms in the unit cell - N
  unsigned int N = 0;
  // Number of groups - per
  unsigned int per = 0;
  // Length of the basic unit cell - Lx Ly Lz
  double Lx = 0, Ly = 0, Lz = 0, V;
  // Charges of the system - q1 q2
  double q1 = 0 , q2 = 0 ;
  ifstream infile(argv[1]);
  infile >> nx >> ny >> nz ;
  infile >> r0 ;
  infile >> N ;
  infile >> per ;
  per *= nx * ny * nz;
  vector<double> r(N*4*nx*ny*nz);
  infile >> q1 >> q2 ;
  infile >> Lx >> Ly >> Lz;

  for(unsigned int i = 0; i < (N*4); i += 4)
    {
      infile >> r[i+0] >> r[i+1] >> r[i+2] >> r[i+3];
      r[i+0] *= Lx;  r[i+1] *= Ly; r[i+2] *= Lz; 
    }
  infile.close();
  infile.clear();

  // replicates the unit cell
  unsigned int i = 0;
  for(unsigned int x = 0; x < nx; ++x)
    {
      for(unsigned int y = 0; y < ny; ++y)
        {
          for(unsigned int z = 0; z < nz; ++z)
            {
              for(unsigned int j = 0;j < (N*4); j += 4)
                {
                  r[i+0] = Lx * x + r[j+0];
                  r[i+1] = Ly * y + r[j+1];
                  r[i+2] = Lz * z + r[j+2];
                  r[i+3] = r[j+3];
                  //cout << i << " " <<  r[i+0] << " " << r[i+1] << " " << r[i+2] <<" " << r[i+3] << endl;
                  i += 4;
                }
            }
        }
    }
  
  Lx *= nx;
  Ly *= ny;
  Lz *= nz;
  V = Lx * Ly * Lz;

  double U1 = 0, U2 = 0, U3 = 0;
  ofstream outfile (argv[2]);
  for(double k = 0.01;k < 3;k += 0.05)
    {
      U1 = RealandReciprocalSpace(r, Lx, Ly, Lz, k, 2);
      U2 = k * PointEnergy(r) / sqrt(Pi);
      U3 =  2 * Pi * pow(Dipole(r),2) / (3 * V);
      outfile << k << " " << ( U1 + U2 ) * r0 / (-8.0 * per) << "\n" ;
    }
  outfile.close();
  outfile.clear();
  
  return 0;
}
