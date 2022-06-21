//espic1d.c
//A 1-D electrostatic PIC code
//non-relativistic
//electrons are mobile along x-axis
//ions are stationary
//electric field  along x-axis
//constant magnetic field along z-axis
//CGS units everywhere

//Version 1.1, 7 April 2015
//Copyright (c) Kartik Patel
//This program is released under the GNU GPL ver 2 or later
//E-Mail:letapk@gmail.com
//Download: https://letapk.wordpress.com
//See the file COPYING for details

/*
 * Compile:
 * $>gcc espic1d.c -lm -o espic1d
 *
 * The following files are created. All data is human readable ASCII text.
 * The format of the data is specified after the description
 * parameters.dat : run time parameters: see initialize function
 * coordinates.dat : particle coordinates: particle-index x vx
 * energy.dat : energy data at each timestep:
 * 		Efieldenergyone kineticenergy Efieldenergyone+kineticenergy+Bfieldenergy
 *		Bfieldenergy
 *		Efieldenergytwo Efieldenergytwo+kineticenergy+Bfieldenergy
 * 
 * trajectory.dat : trajectory in phase space of particle whose index is PARTICLES/4: timestep x vx
 * thermalspread.dat : thermal velocity spread at different times : timestep vsqmean vmeansq vsqmean-vmeansq
 * 
 * In the files below XXX is an integer denoting the timestep at which the file was created
 * velocityXXX.dat : velocity distribution: vx f(vx)
 * positionXXX.dat : particle position distribution: x f(x)
 * potentialXXX.dat : potential distribution: grid-point potential
 * chargeXXX.dat : charge density distribution: grid-point chargedensity
 * elecricfieldXXX.dat : electric field distribution: grid-point Ex
 * phasespaceXXX.dat : position vs. velocity of all particles: x vx vy
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//if the simulation is noisy uncomment the #define below to 
//activate the low-pass filtering of the potential and charge density
//#define LOWPASSFILTER

//no. of particles
#define PARTICLES 50000
//system length
#define IMAX 100;
//no. of timesteps
#define TMAX 500;

#define PI 3.14159275
//no. of iterations in poisson solver
#define MAXITS 10000
//speed of light, cm per sec
#define cee 2.997825e10
//statcoul
#define ELEMCHG 4.802e-10
//grams
#define ELEMASS 9.109e-28
//convergence factor for poisson solver
#define EPS 1.0e-6

//bins for energy spectrum
#define BINS 100
//these many random nos. averaged for creating a gaussian velocity distribution
#define MAXR 10

void initialize (void);
void allocate (void);
double *vector (int rows);

//particle functions
void init_particles (void);
void advance_position (void);
void advance_velocity (double deltime);
void assign_charge_density (void);

//field functions
void poisson (void);
void compute_field (void);
double get_E_field_at_this_position (double x);

//smoothens the variable passed as the argument
//to be used to smoothen potential and charge density
//when there is too much noise in the system
void low_pass_filter (double *variable);

//diagnostic functions
void save_data (void);
void save_system_data (void);
void find_velocity_distribution (void);
void find_position_distribution (void);

//potential
double *u;
//charge density on mesh point, spare for filtering
double *f, *spare;
//sum of the electron and ion chargedensity deposited on the mesh
double eltotal, iontotal;

//electric field array along x and
//constant magnetic field along z
//at the mesh points
double *Ex, Bz;

//electron number density
double numdens;
//drift velocity
double driftvee;
//fraction of total particles moving at this drift
double fraction;
//electron plasma frequency
double omega;
//Larmor frequency, radius and angle of rotation in one timestep
double larmorfreq, larmorrad, larmorangle;
//sin and cosine of the larmorangle
double coslarmorangle, sinlarmorangle;

//timestep, time, total time,
double dt, t, tmax;
//current time step
int n, nmax;
//width of the mesh, total size of box
double dx, xmax;
//maximum index in the mesh
int imax;

//perturbation amplitude and mode
double pertamp, pertmode;

//mass and charge of each particle
double mass;
double totalcharge;
double elcharge;
//what it says
double ionchargedensity;
//electrons per superparticle
double elperpar;
//electron position and velocity
double *elpos, *elveex, *elveey;
//B-field energy
double Bfieldenergy;
//E-field energy at step n
//Int of (E-sqr/8.pi).dv
double Efieldenergyone;
//Int of 0.5.(rho.phi).dv
double Efieldenergytwo;
//total particle kinetic energy at step n
double kineticenergy;

//bin array used in finding velocity and position distribution
int bin[BINS];
//total particles counted in finding the distributions
int totalpos, totalvee;

//array to hold the name of the data file which is created dynamically at runtime
char filename[20];
//holds a description of the simulation
char description[256];

//============================================================
//code starts here
//============================================================

void initialize (void)
//assign values to the physical parameters
//and compute the values used during the run
{
int i;

	//============================================================
	//these are the parameters usually modified for different runs
	//volume number density
	numdens = 1.0e8;//per cubic centimetres
	//meshwidth
	dx = 1.0e-1;//centimetres

	//perturbation mode as a fraction of total length
	//1.0 = fundamental, 2.0 = second harmonic, etc.
	pertmode = 1.0;
	//perturbation amplitude
	pertamp = 0.01;
	//pertamp = 0.0;

	//Bz in gauss
	Bz = 0.0;

	driftvee = 1.0e8;
	
	//fraction of particles moving at the drift velocity
	fraction = 0.2;
	//============================================================

	//imax is the total length of plasma
	xmax = dx * imax;
	
	//total charge present = numberdensity x e x Volume
	totalcharge = numdens * xmax * ELEMCHG;
	//charge per superparticle
	elcharge = -1.0 * totalcharge / PARTICLES;
	//electrons per superparticle
	elperpar = fabs (elcharge / ELEMCHG);
	//mass per superparticle
	mass = elperpar * ELEMASS;
	
	//neutralizing ion background charge density at each mesh point
	ionchargedensity = numdens * ELEMCHG;
	
	//electron plasma frequency
	omega = sqrt (4.0 * PI * numdens * ELEMCHG * ELEMCHG / ELEMASS);
	//time step is related to omega
	dt = 0.2 / omega;
	
	//drift velocity for maximal growth for the two-stream instability
	driftvee = sqrt (3) * omega * xmax * 0.25 / PI;
	
	//Larmor frequency and radius
	larmorfreq = ELEMCHG * Bz / (ELEMASS * cee);
	larmorangle = larmorfreq * dt;
	coslarmorangle = cos (larmorangle);
	sinlarmorangle = sin (larmorangle);
	
	if (Bz != 0.0)
		//to avoid division by 0
		larmorrad = cee * driftvee * ELEMASS / (Bz * ELEMCHG);
	
	//total run time
	tmax = dt * TMAX;

	Bfieldenergy = Bz * Bz * imax * dx / (8.0 * PI);

	//zero potential boundary conditions
	u[imax] = u[0] = 0.0;

	t = 0.0;
	n = 0;
}

//end initialize

void init_particles (void)
//assign position and velocity to electrons at t=0
//
//various types of initializations are already included
//
//the appropriate init can be activated by removing the comments from that block
//the initialization not needed MUST be commented out
//
//any new initialization added MUST comment out the others
{
int i, i1, i2, j, k;
int par, dpar, paruniform;
double v1, v2;
double x, kx, x0, ddx;
double drho, tch;
FILE *fp;
	
	srand (0);
	
	//====================================================================
	//PARTICLE POSITIONS
	//====================================================================

	//============================================================
	//initialize uniformly spaced particles
	//first particle offset
	x0 = 0.5 * xmax / PARTICLES;
	//particles arranged uniformly
	strcpy (description, "Uniformly spaced particles with a sine perturbation\n");
	for (i = 0; i < PARTICLES; i++) {
		elpos[i] =  x0 + 1.0 * i * xmax / PARTICLES;
		//apply a sine-perturbation to the positions
		//the maximum perturbation is one mesh width
		elpos[i] += dx * sin (2.0 * PI * pertmode * elpos[i] / xmax) ;
	}
	//============================================================

	//============================================================
	/*
	//initialize a perturbed density distribution
	i1 = par = dpar = 0;
	//length of density interval
	ddx = xmax / 100.0;
	//no. of particles present along ddx when distributed uniformly
	paruniform = PARTICLES * ddx / xmax;
	strcpy (description, "Perturbed density distribution\n");
	for (j = 0; j < 100; j++) {
		//beginning of this interval
		x = ddx * j;
		//argument of the sine term
		kx = 2.0 * PI * pertmode * x / xmax;
		//perturbation in particles at this location
		dpar = paruniform * pertamp * sin (kx);
		//total no. of particles at this location
		par = (paruniform - dpar);
		//i1 is the starting and i2 is the ending index of these particles
		i2 = i1 + par;
		//first particle offset
		x0 = 0.5 * ddx / par;
		for (i = i1; i < i2; i++) {
			//particles arranged uniformly on this interval
			elpos[i] = (double) (x0 + x + ddx * (i - i1) / par);
	
			//you can also put the particles randomly on this interval
			//particle coordinate is a random no. between x and x+ddx
			//elpos[i] = (double) (x + ddx * rand () / RAND_MAX);
	
			fprintf (fp, "%i %g\n", i, elpos[i]);
		}
		i1 = i2;
		printf ("x=%g to %g, particles=%i, total=%i\n", x, x+ddx, par, i2);
	}
	printf ("Total particles initialized=%i\n", i2);
	*/
	//============================================================

	//============================================================
	/*
	//initialize particle randomly over the whole mesh
	strcpy (description, "Randomly initialized particles\n");
	for (i = 0; i < PARTICLES; i++) {
		//particle coordinate is a random no. between 0 and xmax
		elpos[i] = xmax * rand () / RAND_MAX;
	}
	*/
	//============================================================
	
	//====================================================================
	//PARTICLE VELOCITIES
	//====================================================================

	//============================================================
	/*
	//initialize two warm counterpropagating streams
	//each with a gaussian velocity distribution
	strcat (description, "Two warm counter-propagating streams\n");
	for (i = 0; i < PARTICLES; i++) {
		v1 = 0;
		for (j = 0; j < MAXR; j++) {
			//generate the sum of MAXR random nos. between 0 and 1
			//v1 is between 0 and MAXR
			v1 += 1.0 * rand () / RAND_MAX;
		}
		//v1 is between 0 and 1
		v1 /= (float)MAXR;
		//assign this to this electron
		//the velocity is between 0 and driftvee
		if (i % 2 == 0) {//even index electrons move in the positive direction
			elveex[i] = v1 * driftvee;
			elveey[i] = 0.0;
		}
		else { //odd index electrons move in the negative direction
			elveex[i] = -1.0 * v1 * driftvee;
			elveey[i] = 0.0;
		}
	}
	*/
	//============================================================

	//============================================================
	/*
	//initialize two cold counterpropagating streams
	strcat (description, "Two cold counterpropagating streams\n");
	for (i = 0; i < PARTICLES; i++) {
		if (i % 2 == 0) {//even index electrons move in the positive direction
			elveex[i] = driftvee;
			elveey[i] = 0.0;
		}
		else {//odd index electrons move in the negative direction
			elveex[i] = -1.0 * driftvee;
			elveey[i] = 0.0;
		}	
	}
	*/
	//============================================================

	//============================================================
	/*
	//initialize a single warm drifting stream with a gaussian
	//velocity distribution consisting of a fraction of the 
	//total number of particles uniformly distributed within the plasma
	strcat (description, "Single warm stream in a background plasma\n");
	for (i = 0; i < PARTICLES; i += (1.0 / fraction)) {
		v1 = 0;
		for (j = 0; j < MAXR; j++) {
			//generate the sum of MAXR random nos.
			v1 += 1.0 * rand () / RAND_MAX;
		}
		//v1 is between 0 and 1
		v1 /= MAXR;
		//the velocity of this electron is between 0 and driftvee
		elveex[i] = v1 * driftvee;
		elveey[i] = 0.0;
	}
	*/
	//============================================================
	
	//============================================================
	/*
	//initialize a single cold drifting stream consisting of a fraction
	//of the total number of particles uniformly distributed within the plasma
	strcat (description, "Single cold stream in a background plasma\n");
	for (i = 0; i < PARTICLES; i += (1.0 / fraction)) {
		elveex[i] = driftvee;
		elveey[i] = 0.0;
	}
	*/
	//============================================================

	//store the particle coordinates and velocities
	fp = fopen ("coordinates.dat", "w");
	for (i = 0; i < PARTICLES; i++) {
		fprintf (fp, "%i %.7g %lg\n", i, elpos[i], elveex[i]);
	}
	fclose (fp);
	fflush (fp);
	
}

//end init_particles	 

void advance_position (void)
//move the particle position from step n to n+1
{
int i;

	for (i = 0; i < PARTICLES; i++) {
		//new position at step n+1
		elpos[i] += elveex[i] * dt;

		//apply periodic boundary conditions
		if (elpos[i] < 0.0){
			//reset position
			elpos[i] = elpos[i] + xmax;
		}

		if (elpos[i] > xmax){
			//reset position
			elpos[i] = elpos[i] - xmax;
		}

		if ((elpos[i] < 0.0) || (elpos[i] > xmax)){
			printf ("Particle out of bounds. Abnormally terminated: i=%i x=%g\n", i, elpos[i]);
			exit(1);
		}
	}
}

//end advance_position

void advance_velocity(double deltime)
//velocity from step n-1/2 to step n+1/2
{
int i;
//elecric field on the particle
double ex;
//particle velocities used in v x B rotation
double elveex1, elveex2, elveexnew;
double elveey1, elveey2, elveeynew;
//half elecric field impulse
double halfEpush;

	kineticenergy = 0.0;
	for (i = 0; i < PARTICLES; i++) {
		//get the E-field acting on this particle
		ex = get_E_field_at_this_position (elpos[i]);

		halfEpush = elcharge * ex * deltime * 0.5 / mass;

		//new velocity at step n+1/2
		//elveexnew = elveex[i] + elcharge * ex * deltime / mass;

		//half-acceleration
		elveex1 = elveex[i] + halfEpush;
		elveey1 = elveey[i];

		//Lorentz (v x B) rotation
		elveex2 = elveex1 * coslarmorangle + elveey1 * sinlarmorangle;
		elveey2 = - elveex1 * sinlarmorangle + elveey1 * coslarmorangle;

		//half-acceleration
		elveexnew = elveex2 + halfEpush;
		elveeynew = elveey2;
		
		//accumulate kineticenergy at step n : sum of 0.5 x m x v-sqr
		//v-sqr is the product of the old and new velocities
		//fabs protects against negative energy!
		kineticenergy += 0.5 * mass * fabs (elveexnew * elveex[i] + elveeynew * elveey[i]);

		//assign the new velocities to this particle
		elveex[i] = elveexnew;
		elveey[i] = elveeynew;
	}
	
}

//end advance_velocity

void assign_charge_density (void)
//deposit chargedensity on the mesh points at step n
//collect electron charge density from particles and assign to mesh points
//and add the constant ion charge density as the background
{
int i, ip1, j;
//x is the particle coordinate
//mx is the mesh coordinate
//lx is the local coordinate - distance from i
//lx1 is the  complementary local coordinates - distance from i+1
double x, mx, lx, lx1;
//chargedensity at mesh points i and i+1
double chgdeni, chgdenip1;

	//clear charge density array
	for (i = 0; i <= imax; i++)
		f[i] = 0.0;
	
	eltotal = 0.0;
	iontotal = 0.0;

	//assign the electron charge density
	for (j = 0; j < PARTICLES; j++) {
		x = elpos[j];
		mx = x / dx;
		//get index of nearest mesh point
		i = floor (mx);
		ip1 = i + 1;
		
		//find local coordinate of particle, in mesh units
		//distance from i
		lx = mx - i;
		//find complementary coordinate of particle, in mesh units
		//distance from ip1
		lx1 = ip1 - mx;
		
		/*
		 *               i       x            i+1
		 *          _____|_______|_____________|_____
		 *               | lx    |     lx1     |
		 */
		
		//interpolate particle charge density to the mesh points i and ip1
		chgdeni = lx1 * elcharge / dx;
		chgdenip1 = lx * elcharge / dx;

		f[i] += chgdeni;
		f[ip1] += chgdenip1;
	}

	//add the neutralizing background ion density to each mesh point
	//accumulate the total electron and ion chargedensity deposited
	for (i = 1; i < imax; i++) {
		eltotal += f[i];
		iontotal += ionchargedensity;
		f[i] += ionchargedensity;
	}

	eltotal += f[0];
	iontotal += ionchargedensity * 0.5;
	f[0] += ionchargedensity * 0.5;

	eltotal += f[imax];
	iontotal += ionchargedensity * 0.5;
	f[imax] += ionchargedensity * 0.5;

	#ifdef LOWPASSFILTER
	//filter the charge density
	low_pass_filter (f);
	#endif
}

//end assign_charge_density

void poisson (void)
//1-D poisson solver
//uses the Gauss-Seidel iteration
//boundary points are fixed at zero potential
{
int ipass, i, ip1, im1, n;
double anorm, anormf=0.0, resid;
double fourpi = 4.0 * PI * dx * dx;

	for (i = 0; i <= imax; i++)
		anormf += fabs(f[i] + u[i]);

     for (n = 1; n <= MAXITS; n++) {
		anorm = 0.0;

		for (i = 1; i < imax; i++) {
			resid = u[i+1] + u[i-1] - 2.0 * u[i] + fourpi * f[i];
			anorm += fabs(resid);
			u[i] += resid * 0.5;
		}

		if (anorm < EPS*anormf) {
               //printf ("n = %i, anorm = %f\n", n, anorm);
               return;
          }
	}
	printf ("Maximum iterations exceeded, anorm = %f, anormf = %f\n", anorm, anormf);
}

//end poisson

void compute_field (void)
//the electric field at each mesh point is
//E= - grad-V
{
int i;

	Efieldenergyone = Efieldenergytwo = 0.0;

	#ifdef LOWPASSFILTER
	//filter the potential
	low_pass_filter (u);
	#endif
	
	Ex[0] = -(u[1] - u[imax-1])/(2.0 * dx);
	for (i = 1; i < imax; i++) {
		Ex[i] = -(u[i+1] - u[i-1])/(2.0 * dx);

		Efieldenergyone += Ex[i] * Ex[i] * dx;
		Efieldenergytwo += f[i] * u[i] * dx;
	}
	Ex[imax] = -(u[1] - u[i-1])/(2.0 * dx);

	Efieldenergyone += Ex[imax] * Ex[imax] * dx * 0.5;
	Efieldenergyone += Ex[0] * Ex[0] * dx * 0.5;
	Efieldenergyone = Efieldenergyone / (8.0 * PI);

	Efieldenergytwo += f[imax] * u[imax] * dx * 0.5;
	Efieldenergytwo += f[0] * u[0]  * 0.5;
	Efieldenergytwo = 0.5 * Efieldenergytwo;
}

//end compute_field

void low_pass_filter (double *variable)
//smmothen out the variable passed to this function
{
int i;

	spare[0] = variable[imax-1] + 4.0 * variable[0] + variable[1];
	spare[0] = spare[0] / 6.0;
	for (i = 1; i < imax; i++) {
		spare[i] = variable[i-1] + 4.0 * variable[i] + variable[i+1];
		spare[i] = spare[i] / 6.0;
	}
	spare[imax] = variable[imax-1] + 4.0 * variable[imax] + variable[1];
	spare[imax] = spare[imax] / 6.0;
	
	for (i = 0; i <= imax; i++) {
		variable[i] = spare[i];
	}	
}

//end low_pass_filter

double get_E_field_at_this_position (double x)
//use linear weighing to interpolate E-field from mesh points to particle position
{
int i, ip1;
//lx and lx1 are the local coordinates in the cell
double lx, lx1;
//mx is the coordinate in mesh units
double mx;
//electric field
double ex;

	mx = x / dx;
	//get index of nearest mesh point
	i = floor (mx);
	ip1 = i + 1;
	
	//find local coordinate of particle, in mesh units
	lx = mx - i;
	//find complementary coordinate of particle, in mesh units
	lx1 = ip1 - mx;
	
	/*
	 *               i       x            i+1
	 *          _____|_______|_____________|_____
	 *               | lx    |    lx1      |
	 */
	
	//interpolate the fields at the mesh points to the particle position
	//to this we can add a constant external electric field, if any
	
	ex = lx * Ex[ip1] + lx1 * Ex[i];

	return ex;
	
}

//end get_E_field_at_this_position

void find_velocity_distribution (void)
//distribution of particle velocities
{
int i, j;
FILE *op;
double v1, v2, vmax, vmin, vmeansq, vsqmean, vsqtot, vtot;
double binrange, binint;

	for (j = 0; j < BINS; j++) {
		bin[j] = 0;
	}
	totalvee = 0;
	
	vmax = vmin = 0.0;
	vtot = vmeansq = vsqtot = vsqmean = 0.0;
	
	for (i = 0; i < PARTICLES; i++) {
		vmin = (elveex[i] <= vmin) ? elveex[i] : vmin;
		vmax = (elveex[i] >= vmax) ? elveex[i] : vmax;
		if (i % 2) {//even index, positive driftvee
			//add up all the velocities
			vtot += elveex[i];
			//add up the square of all the velocities
			vsqtot += elveex[i] * elveex[i];
		}
		
	}
	
	//mean velocity <v>
	vmeansq = 2.0 * vtot / PARTICLES;
	//square of the mean velocity <v>*<v>
	vmeansq = vmeansq * vmeansq;
	
	//mean of square of velocities <v*v>
	vsqmean = 2.0 * vsqtot / PARTICLES;
	
	binrange = vmax - vmin;
	//increase the range by 1% one either side
	vmin = vmin - binrange / 100.0;
	vmax = vmax + binrange / 100.0;

	binrange = vmax - vmin;
	binint = binrange / BINS;
	for (j = 0; j < BINS; j++) {
		v1 = vmin + j * binint;
		for (i = 0; i < PARTICLES; i++) {
			v2 = v1 + binint;
			if ((elveex[i] >= v1) && (elveex[i] < v2))
				bin[j]++;
		}
	}

	sprintf (filename, "velocity%i.dat", n);
	
	op = fopen (filename, "w");
	for (j = 0; j < BINS; j++) {
		v1 = vmin + j * binrange / BINS;
		fprintf (op, "%g %i\n", v1, bin[j]);
		totalvee += bin[j];
	}
	fclose (op);
	fflush (op);
	
	op = fopen ("thermalspread.dat", "a");
	fprintf (op, "%i %g %g %g\n", n, vsqmean, vmeansq, vsqmean-vmeansq);
	fclose (op);
	fflush (op);
	
}

//end find_velocity_distribution

void find_position_distribution (void)
//distribution of particle positions
{
int i, j;
FILE *op;
double x1, x2;

	for (j = 0; j < BINS; j++) {
		bin[j] = 0;
	}
	totalpos = 0;
	
	for (j = 0; j < BINS; j++) {
		//starting position of this bin
		x1 = j * xmax / BINS;
		for (i = 0; i < PARTICLES; i++) {
			//ending positions of this bin
			x2 = x1 + xmax / BINS;
			if ((elpos[i] >= x1) && (elpos[i] < x2))
				bin[j]++;
		}
	}

	sprintf (filename, "position%i.dat", n);
	op = fopen (filename, "w");
	for (j = 0; j < BINS; j++) {
		fprintf (op, "%i %i\n", j, bin[j]);
		totalpos += bin[j];
	}
	fclose (op);
	fflush (op);
}

//end find_position_distribution

void save_system_data (void)
{
FILE *op;

	op = fopen ("parameters.dat", "w");
		fprintf (op, "Description:\n");
		fprintf (op, "%s\n", description);
		fprintf (op, "numdens=%lg per cc\n", numdens);
		fprintf (op, "dx=%lg cm\n", dx);
		fprintf (op, "dt=%lg sec\n", dt);
		fprintf (op, "totalcharge=%lg statcoul\n", totalcharge);
		fprintf (op, "elcharge=%lg statcoul\n", elcharge);
		fprintf (op, "elperpar=%lg\n", elperpar);
		fprintf (op, "mass=%lg g\n", mass);
		fprintf (op, "ionchargedensity=%lg\n", ionchargedensity);
		fprintf (op, "omega=%lg\n", omega);
		fprintf (op, "magnetic field=%lg gauss\n", Bz);
		fprintf (op, "Larmor freq=%lg\nLarmor radius=%lf\n", larmorfreq, larmorrad);
		fprintf (op, "driftvee=%lg cm per sec\n", driftvee);
		fprintf (op, "fraction of particles drifting=%f, increment=%f\n", fraction, (1.0/fraction));
		fprintf (op, "pertmode=%lg\n", pertmode);
		fprintf (op, "pertamp=%lg\n", pertamp);
	fclose (op);
	fflush (op);

}

//end save_system_data

void save_data (void)
//write out the data files at each timestep
{
FILE *op;
int i;

	sprintf (filename, "potential%i.dat", n);
     op = fopen (filename, "w");
     for (i = 0; i <= imax; i++) {
		fprintf (op, "%i %g\n", i, u[i]);
     }
     fflush (op);
     fclose (op);

	sprintf (filename, "charge%i.dat", n);
	op = fopen (filename, "w");
	for (i = 0; i <= imax; i++) {
		fprintf (op, "%i %g\n", i, f[i]);
	}
	fflush (op);
	fclose (op);

	sprintf (filename, "elecricfield%i.dat", n);
	op = fopen (filename, "w");
	for (i = 0; i <= imax; i++) {
		fprintf (op, "%i %g\n", i, Ex[i]);
	}
	fflush (op);
	fclose (op);

	sprintf (filename, "phasespace%i.dat", n);
	op = fopen (filename, "w");
	for (i = 0; i <= PARTICLES; i++) {
		fprintf (op, "%g %g %g\n", elpos[i], elveex[i], elveey[i]);
	}
	fflush (op);
	fclose (op);

	op = fopen ("energy.dat", "a");
	fprintf (op, "%i %g %g %g %g %g %g\n", n,
		Efieldenergyone, kineticenergy, Efieldenergyone+kineticenergy+Bfieldenergy,
		Bfieldenergy,
		Efieldenergytwo, Efieldenergytwo+kineticenergy+Bfieldenergy
	);
	fflush (op);
	fclose (op);

	op = fopen ("trajectory.dat", "a");
	fprintf (op, "%i %g %g\n", n, elpos[PARTICLES/4], elveex[PARTICLES/4]);
	fflush (op);
	fclose (op);

	//at the end of the run
	if (n == nmax) {
		op = fopen ("parameters.dat", "a");
		fprintf (op, "step=%i, eltotal=%g, iontotal=%g, sum=%g, totalpos=%i totalvee=%i\n",
			    n, eltotal, iontotal, iontotal+eltotal, totalpos, totalvee);
		fflush (op);
		fclose (op);
	}
	
}

//end save_data

int main (void)
{
	system ("rm *.dat");

	imax = IMAX;
	nmax  = TMAX;
	
	allocate ();
	
	initialize ();
	
	init_particles ();
	
	save_system_data ();

	//charge density, potential and E-field at t=0
	assign_charge_density ();
	poisson ();
	//field and E-field energy at step 0
	compute_field ();
	
	//move the particle velocity back half a time step, to -1/2
	advance_velocity (-dt*0.5);
	
	//At this stage the following quantities are known at time t = 0:
	//particle positions, electric fields, electric field energy
	//potentials, charge densities
	//and the following are known at time t = -1/2:
	//particle velocities
	//and the following are known at time t = -1/4:
	//kinetic energy
	
	while (n <= nmax) {
		//At this stage the following quantities are known at time t = n:
		//particle positions, electric fields, electric field energy
		//potentials, charge densities
		//and the following are known at time t = n-1/2:
		//particle velocities

		//diagnostics at step n
		find_velocity_distribution ();
		find_position_distribution ();
		
		//particle functions
		//velocity advanced from n-1/2 to n+1/2
		//accumulate KE at step n
		advance_velocity (dt);
		
		//data saved at step n
		save_data ();
		
		//position advanced from n to n+1
		advance_position ();
		//chargedensity assigned at step n+1
		assign_charge_density ();
		
		//field functions
		//potential computed at n+1
		poisson ();
		//field and E-field energy computed at step n+1
		compute_field ();
		
		//increment n
		n++;
		printf ("n=%i steps,\tt=%g sec\n", n, t);
		t += dt;
		
		//At this stage the following quantities are known at time t = n+1:
		//particle positions, electric fields, electric field energy,
		//potentials, charge densities
		//and the following are known at time t = n+1/2:
		//particle velocities
	}
	
	return 0;
}

//end main

void allocate (void)
//assign memory and get the pointers to the variables
{
	//chargedensity
	f = vector (imax);
	//potential
	u = vector (imax);
	//electric field
	Ex = vector (imax);
	//electron position
	elpos = vector (PARTICLES);
	//electron velocity along x and y
	elveex = vector (PARTICLES);
	elveey = vector (PARTICLES);

	#ifdef LOWPASSFILTER
	//spare array used for filtering
	spare = vector (imax);
	#endif
}

//end allocate

double *vector (int rows)
//allocate and return a pointer to the memory block
{
double *v;
	
	v = (double *) calloc (rows+1, sizeof(double));
	
	return v;
}

//end vector

//end file
