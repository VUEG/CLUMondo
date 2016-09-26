/*****************************************************************
 * This file is part of CLUMondo.
 *
 * CLUMondo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * CLUMondois distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CLUMondo.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Development:
 *	Peter Verburg
 *	31.3.2015
 *
 * All rights: Peter Verburg <peter.verburg@vu.nl>
 *
 ********************************************************************/

#include <stdio.h>

#define F_INIT "main.1"				// name of file containing main input parameters
#define F_LOG "log.fil"				// name of log-file
#define CLUMONDO_MAX_FILE_LENGTH 260

#define NCOV 30						// number of land use/cover types
#define NREG 1						// number of regions; default: 1
#define NFACT 40					// maximum number of explanatory factors in one regression equation
#define NATT 100					// total number of grids with explanatory factors
#define NR 4000						// number of rows in grid
#define NC 4000						// number of columns in grid
#define NYEAR 101					// maximum number of years to be simulated
#define NOUTP 10					// maximum number of output types to be included

/* typedefs for memory allocation */
typedef int MATELEMS;
typedef float MATFLOAT;
typedef MATELEMS *vector;
typedef MATFLOAT *fvector;

extern fvector *mat_at[NATT];		// matrix for grids with explanatory factors
extern fvector *mat_co[NCOV];		// matrix for regression results
extern fvector *initoutput[NOUTP];	// matrix with initial values of output types
extern vector *mat_newco;			// matrix with land use results
extern vector *mat_oldco;			// matrix for old land use type grid
extern vector *region;				// matrix for region subdivision
extern vector *mat_age;				// matrix keeping track of time without change
extern vector *luschange;			///matrix keeping track if LUS has changed from initial LUS in order to write output land system characteristics
extern MATFLOAT *param;				// matrix for sorting algorithm

extern int
	start,							// start year of simulation
	end,							// end year of simulation
	year,							// year counter
	allocneigh[NCOV][NCOV],			// indicates those cover types per cover type for which the neighborhood is added to the location
	attr[NREG][NCOV][NFACT],		// explanatory factor numbers for region, land use type and explanatory factor number
	attr2[NREG][NCOV][NFACT],
	no_fact[NREG][NCOV],			// number of explanatory factors in regression equation for region and land use type
	no_fact2[NREG][NCOV],
	loop,							// counter for iteration loop in allocation module
	dem_dir[NREG][NCOV],			// direction of demand: 1 = increase; -1 = decrease; 0 = stable; for land use types
	scale,							// control parameter: default value = 1
	arcview,						// arcview output control parameter: 1 = arcviewheaders 0 = no arcviewheaders (idrisi etc.)
	grchange[NFACT],		        // array containing numbers of changing scgr (explanatory factor) files
	nochange,						// number of scgr files that is changing
	nostable[NREG],					// number of cover types of which demand does not change
	Yearprog,
	diffregregi,					// 1 different regions have different regressions; 0 one regression for all
	lusconvmap,						// switch to turn on reading intitial land use system composition maps, 1 =  output parameter init files read, all other off
	nocells[NREG],
	cropmatrix[NCOV][NCOV],
	bgsmatrix[NCOV][NCOV],
	ppmatrix[NCOV][NCOV],
	builtupmatrix[NCOV][NCOV],
	allowed[NCOV][NCOV],
	allowchange[NCOV][NCOV],
	lusconv[NCOV][NCOV],
	probmaps,
	natcov,
	rncov,
	rnreg,
	rnfact,
	rnatt,
	rnr,
	rnc,
	nooutputmaps,					// number of maps with initial values for output types
	rndemand,						// number of different demand types
	bottomup,						// indicator if bottom-up modification of demand is to be accounted for
	startage,						// maximum age of land use types at start of simulation (randomized)
	agemode,						// binary switch: 1: use random calculation of initial age; 0: read initial age from age.0
	itermode,						// switch for iteration mode: 0: percentage difference between alloc and demand; 1: absolute difference
	influence,						// switch to turn on neighborhood calculations
	locationfactor,					// switch to turn on location specific additions to regression equation
	checkfile,						// switch to turn on module to check driving factor files on missing data; indicator for error in these files
	metamod,						// switch to turn the meta-model mode on or of
	allowmin,						// provides minimum number of conversion matrix that refers to a grid
	allprob,						// switch to turn on the option to create all probability maps in advance
	covdem[NCOV],					// stores the number of the demand array to which this land use type contributes
	maxcovdem,						// highest number of unique land cover type in demand specification
	changeconv,						// switch to turn on dynamic lusmatrix (lusmatrix.year instead of lusmatrix.txt
	writeoutputs;					// switch turning on the writing of output grids value 1 is ON

extern double
	demand[NREG][NCOV][NYEAR],    	// demand for land use type
	constant[NREG][NCOV],
	constant2[NREG][NCOV],	        // constant in regression equation for scale level, region and land use type
	cellsize,                       // surface area of one grid cell (same units as demand)
	xllcorner,                      // x-coordinate of lower left corner (only needed for georeferenced arcview simulations)
	yllcorner;                      // y-coordinate of lower left corner (only needed for georeferenced arcview simulations)

extern double
	areacov[NYEAR][NREG][NCOV];

extern float
	maxdiff,						// maximum difference between demand and allocation allowed at end of iteration
	totdiff,                        // total difference between demand and allocated land use (all cover types)
	demdiff,						// actual difference between demand and allocation
	lusmatrix[NCOV][NCOV],          // matrix containing the quantities of demand per land use system type
	elas[NREG][NCOV],				// iteration parameter, based upon break-off value of logistic regression
	ltstat[NCOV],						// indicator for behaviour allowed for land use type: 0 = increases and decreases always allowed, 1 = only decraese if decreasing demand, -1 = no increase if decreasing demand, 0-1 = relative staticity for competition among increasing demands; larger: more on decreasing covers
	fact[NREG][NCOV][NFACT],        // beta for region, land use type and explanatory factor number
	fact2[NREG][NCOV][NFACT],
	gridsize,						// size of grid as used by arcview
	speed[NREG][NCOV],
	dem_elas[NREG][NCOV],           // elasticity for the respective region and demand group type
	demdiffmax,
	ratio[NCOV],                    // ratio of land use type to total land area
	neighweight[NCOV],
	locationweight[NCOV],           // fraction of location specific addition to be added to probability
	seed,                           // iteration variable
	weightdem[NCOV];                // weight of the different demand types in the modification of the iteration factor

extern FILE*
	flog;							// log-file

extern char    editfile[CLUMONDO_MAX_FILE_LENGTH],
        park[CLUMONDO_MAX_FILE_LENGTH],
        demd[CLUMONDO_MAX_FILE_LENGTH];

extern int g_argc;
extern char** g_argv;

/* sorting rule for qsort command used in void alloc_sort */
inline int sorting_rule( const void* a, const void* b )
{
    printf( "b = %f, a = %f\n", ***(const float***)b, ***(const float***)a );
    if ( ***(const float***)b < ***(const float***)a ) {
        return -1;
    } else {
        return 1;
    }
}
