#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#include "Mie_Functions.h"

#define pi 3.1415926535897932384626433
const complex<double> j(0.0,1.0);


int main(int argc, char **argv)
{
	int n,  i, SEL, npts, plane;
    double r = 0.0, lam0 = 0.0, R, rstart, rstop, wstart, wstop, IN;
    double Lx, Ly, Lz, Cx, Cy, Cz, res;
	char *fname = NULL;
    complex<double>  *F, ers, mrs;
	vector<double> rv, wv;
	vector<complex<double>> CSsca;
    vector<vector<double>> pts;
	vector<vector<complex<double>>> vF;

	
	
	complex<double>  erb(1,0); 		// permittivity of background medium (only real value)
	complex<double>  mrb(1,0);		// permeabilitty of background medium (only real value)
	
	
	
	if (argc == 1)
	{
		printf("Choose what you what to do by specifying the following arguments:\n");
		printf("-----------------------------------------------------------------\n");
		printf("Mie computation on a rectangular grid: 1 lam n r Lx Ly Lz Cx Cy Cz res filename\n");
		printf("Mie computation on the surface of a sphere: 2 lam n r R angular res filename\n");
		printf("Cross sections (fixed wavelength): 3 lam n rstart rstop npts filename\n");
		printf("Cross sections (fixed radius): 4 r n wstart wstop npts filename\n");
		printf("Radiation pattern in the xy-, xz-, or yz-plane: 5 lam n r R angular res plane filename\n");
		printf("-----------------------------------------------------------------\n");
		return 1;
	}
	
	
    // --------------------------------------------------
	// Define parameters
	// --------------------------------------------------
	
	SEL			= int(atof(argv[1]));			// Select computation

	if (SEL == 1)
	{						
		lam0		= atof(argv[2]);   				// Free space wavelength
		n   		= atof(argv[3]);          		// Number of Mie coefficients
		r			= atof(argv[4]);         		// Sphere Radius
		Lx  		= atof(argv[5]);				// Grid length along x
		Ly  		= atof(argv[6]);				// Grid length along y
		Lz  		= atof(argv[7]);				// Grid length along z
		Cx  		= atof(argv[8]);				// Grid center along x
		Cy  		= atof(argv[9]);				// Grid center along y
		Cz  		= atof(argv[10]);				// Grid center along z
		res 		= atof(argv[11]);				// Grid resolution [m]
		fname 	= argv[12];						// File name
	}
	else if (SEL == 2)
	{
		lam0		= atof(argv[2]);   				// Free space wavelength
		n   		= atof(argv[3]);          		// Number of Mie coefficients
		r			= atof(argv[4]);         		// Sphere Radius		
		R  		= atof(argv[5]);				// Points sphere radius
		res  		= atof(argv[6]);				// Angular resolution [°]
		fname 	= argv[7];						// File name
	}
	else if (SEL == 3)
	{
		lam0		= atof(argv[2]);   				// Free space wavelength
		n   		= atof(argv[3]);          		// Number of Mie coefficients
		rstart 	= atof(argv[4]);				// Smallest sphere radius [m]
		rstop  	= atof(argv[5]);				// Largest sphere radius [m]
		npts 		= atof(argv[6]);				// Number of steps between rstart and rstop
		fname 	= argv[7];						// File name
	}
	else if (SEL == 4)
	{
		r			= atof(argv[2]);   				// Radius [m]
		n   		= atof(argv[3]);          		// number of Mie coefficients
		wstart  = atof(argv[4]);				// Smallest free-space wavelength [m]
		wstop  	= atof(argv[5]);				// Largest free-space wavelength [m]
		npts 		= atof(argv[6]);				// Number of steps between wstart and wstop
		fname 	= argv[7];						// File name
	}
	else if (SEL == 5)
	{
		lam0		= atof(argv[2]);   				// Free space wavelength
		n   		= atof(argv[3]);          		// Number of Mie coefficients
		r			= atof(argv[4]);         		// Sphere Radius
		R  		= atof(argv[5]);				// Points sphere radius
		res  		= atof(argv[6]);				// Angular resolution [°]
		plane  	= int(atof(argv[7]));			// Plane (1:xy, 2:xz, 3:yz) 
		fname 	= argv[8];						// File name
	}
	
    // --------------------------------------------------
	// Generate points
	// --------------------------------------------------
    
	// IN: incident + scattered fields (IN = 1), scattered fields only (IN = 0)
	
    if (SEL == 1) {pts = RectPtsGen(Lx, Ly, Lz, Cx, Cy, Cz, res); IN = 1.0;}
	if (SEL == 2) {pts = SpherePtsGen(R, res); IN = 1.0;}
	if (SEL == 3) rv  = linspace(rstart, rstop, npts);
	if (SEL == 4) wv  = linspace(wstart, wstop, npts);
	if (SEL == 5) {pts = CirclePtsGen(plane, R, res); IN = 0.0;}
    
      
	if (SEL == 1 || SEL == 2 || SEL == 5)
	{
		printf("Number of points to be evaluated: %d\n", int(pts.size()));
    
		for(i = 0; i < int(pts.size()); i++) vF.push_back({0,0,0,0,0,0});
		
		// Computes fields (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(pts.size()); i++)
		{
			complex<double> F0[6];
			  
			F = Mie_fieldsPEC(IN, pts[i][0], pts[i][1], pts[i][2], n, r, lam0, erb, mrb, F0);
			
			vF[i][0] = F[0]; vF[i][1] = F[1]; vF[i][2] = F[2];
			vF[i][3] = F[3]; vF[i][4] = F[4]; vF[i][5] = F[5];      
		}
		
		// Save to file
		FieldsToFileXYZ(fname, pts, vF, r, lam0, n, erb, mrb);
	}

	if (SEL == 3)
	{
		for(i = 0; i < int(rv.size()); i++) CSsca.push_back(0);
		
		// Computes scattering cross section (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(rv.size()); i++)	CSsca[i] = CscaPEC(n, rv[i], lam0, erb, mrb);
		
		// Save to file		
		CsToFileR(fname, rv, CSsca, lam0, n, erb, mrb);
	}
	
	if (SEL == 4)
	{	
		for(i = 0; i < int(wv.size()); i++) CSsca.push_back(0); 
		
		// Computes scattering cross section (in parallel)
		#pragma omp parallel for
		for(i = 0; i < int(wv.size()); i++)	CSsca[i] = CscaPEC(n, r, wv[i], erb, mrb);
		
		// Save to file
		CsToFileW(fname, wv, CSsca, r, n, erb, mrb);
	}
	
	return 0;
}


