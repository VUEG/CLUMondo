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
 * CLUMondo globals
 * CLUMondo S.1.1
 *
 * Development:
 *	Peter Verburg
 *	31.3.2015
 *
 * All rights: Peter Verburg <peter.verburg@vu.nl>
 *
 ********************************************************************/
#include "../include/globals.h"

fvector *mat_at[NATT];		// matrix for grids with explanatory factors
fvector *mat_co[NCOV];		// matrix for regression results
fvector *initoutput[NOUTP];	// matrix with initial values of output types
vector *mat_newco;			// matrix with land use results
vector *mat_oldco;			// matrix for old land use type grid
vector *region;				// matrix for region subdivision
vector *mat_age;				// matrix keeping track of time without change
vector *luschange;			///matrix keeping track if LUS has changed from initial LUS in order to write output land system characteristics
//fvector *param;					// matrix for sorting algorithm
MATFLOAT *param;				// matrix for sorting algorithm

int
	start,
	end,
	year,
	allocneigh[NCOV][NCOV],
	attr[NREG][NCOV][NFACT],
	attr2[NREG][NCOV][NFACT],
	no_fact[NREG][NCOV],
	no_fact2[NREG][NCOV],
	loop,
	dem_dir[NREG][NCOV],
	scale,
	arcview,
	grchange[NFACT],
	nochange,
	nostable[NREG],
	Yearprog,
	diffregregi,
	lusconvmap,
	nocells[NREG],
	cropmatrix[NCOV][NCOV],
	bgsmatrix[NCOV][NCOV],
	ppmatrix[NCOV][NCOV],
	builtupmatrix[NCOV][NCOV],
	allowed[NCOV][NCOV],
	allowchange[NCOV][NCOV],
	lusconv[NCOV][NCOV],
	probmaps = 0,
	natcov,
	rncov,
	rnreg,
	rnfact,
	rnatt,
	rnr,
	rnc,
	nooutputmaps,
	rndemand,
	bottomup,
	startage,
	agemode,
	itermode,
	influence,
	locationfactor,
	checkfile,
	metamod,
	allowmin,
	allprob,
	covdem[NCOV],
	maxcovdem,
	changeconv,
	writeoutputs;

double
	demand[NREG][NCOV][NYEAR],
	constant[NREG][NCOV],
	constant2[NREG][NCOV],
	cellsize,
	xllcorner,
	yllcorner;

double
	areacov[NYEAR][NREG][NCOV];

float
	maxdiff,
	totdiff,
	demdiff,
	lusmatrix[NCOV][NCOV],
	elas[NREG][NCOV],
	stat[NCOV],
	fact[NREG][NCOV][NFACT],
	fact2[NREG][NCOV][NFACT],
	gridsize,
	speed[NREG][NCOV],
	dem_elas[NREG][NCOV],
	demdiffmax,
	ratio[NCOV],
	neighweight[NCOV],
	locationweight[NCOV],
	seed,
	weightdem[NCOV];

FILE*
	flog;

char
	editfile[CLUMONDO_MAX_FILE_LENGTH],
	park[CLUMONDO_MAX_FILE_LENGTH],
	demd[CLUMONDO_MAX_FILE_LENGTH];

int g_argc;
char** g_argv;