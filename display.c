//display.c

//Copyright (c) Kartik Patel
//This program is released under the GNU GPL ver 2 or later
//E-Mail:letapk@gmail.com
//Download: https://letapk.wordpress.com
//See the file COPYING for details
//Last modified 23 May 2014

/*
 * Data display program for the 1-D PIC code
 *
 * Gnuplot is required to be present and accessible
 *
 * Compile:
 * $>gcc display.c -lm -o display
 *
 * Useage : ./display filename
 * where filename is the name of the energy, trajectory, or thermalspread dat file
 * or the prefix of the filenames for the data which is written at each time step
 *
 * eg:
 * or
 * ./display charge : this will display the charge density at each timestep
 * from file charge0.dat onwards
 * 
 * The following data are available for display:
 *
 * velocity : velocity distribution: vx versus f(vx)
 * position : particle position distribution: x versus f(x)
 * potential : potential distribution: x versus grid-point potential
 * charge : charge density distribution: x versus grid-point chargedensity
 * elecricfield : electric field distribution: x versus grid-point Ex
 * phasespace : position vs. velocity of all particles: x versus vx vy
 *
 * In addition:
 * ./display energy.dat : this will display kinetic, electric and total energy versus time
 * ./display thermalspread : this plots the thermalspread versus time
 * ./display trajectory : this plots the trajectory of a sample particle with time
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

char datafile[50];
char plottitle[50];

int main (int argc, char *argv[])
{
int i;
FILE *fp, *fpd;

	printf ("Showing %s data\n", argv[1]);
	
	fp = fopen ("gnuinp.inp", "w");
	
	i = 0;
	for (;;) {
		//three plots for energy
		if (strcmp(argv[1], "energy.dat") == 0) {
			fprintf (fp, "plot '%s' u 1:2 w l title \"Electric field energy\"\n", argv[1]);
			fprintf (fp, "replot '%s' u 1:3 w l title \"Kinetic energy\"\n", argv[1]);
			fprintf (fp, "replot '%s' u 1:4 w l title \"Total energy\"\n", argv[1]);
			//leave
			break;
		}
		//also three plots for thermal velocity spread
		if (strcmp(argv[1], "thermalspread.dat") == 0) {
			fprintf (fp, "plot '%s' u 1:2 w l title \"Mean of velocity squared\"\n", argv[1]);
			fprintf (fp, "replot '%s' u 1:3  title \"Square of mean velocity\"\n", argv[1]);
			fprintf (fp, "replot '%s' u 1:4 title \"Thermal spread\"\n", argv[1]);
			//leave
			break;
		}

		if (strcmp(argv[1], "trajectory.dat") == 0) {
			fprintf (fp, "plot '%s' u 2:3 w l title \"Trajectory of sample particle\"\n", argv[1]);
			//leave
			break;
		}
		
		//append the time step to the file name for the other data files
		sprintf (datafile, "%s%i.dat", argv[1], i);

		//try to open the data file
		fpd = fopen (datafile, "r");
		//file does not exist, so leave
		if (fpd == NULL)
			break;
		else {
			//data file exists, so close it
			fclose (fpd);
			//write gnuplot command
			//phasespace is plotted with dots
			if (strcmp(argv[1], "phasespace") == 0) {
				sprintf (plottitle, "Phasespace, x vs. v(x), step %i", i);
				fprintf (fp, "plot '%s' w d title \"%s\"\n", datafile, plottitle);
			}
			//everything else with lines
			else {
				sprintf (plottitle, "Step %i", i);
				fprintf (fp, "plot '%s' w l title \"%s\"\n", datafile, plottitle);
			}
			//increment file suffix
			i++;
		}
	}
	//close and flush gnuplot command file
	fclose (fp);
	fflush(fp);

	//invoke gnuplot
	system ("gnuplot gnuinp.inp -persist");
	//remove the command file
	i = unlink ("./gnuinp.inp");
	
}

//end file
