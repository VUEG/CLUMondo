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
 * The CLUMondo model calculation modules
 * CLU-mondo S.1.1
 *
 * Development:
 *      Peter Verburg
 *      Sanneke van Asselen
 *
 * 31.3.2015
 *
 * peter.verburg@vu.nl
 ********************************************************************/

#include "../include/globals.h"
#include "../include/util.h"

#include <cmath>
#include <cstring>
#include <time.h>
#include <fstream>

/*
 * Read main simulation parameters from input file.
 */
void all_init()
{
	int i, dimensionerror;
	dimensionerror = 0;

	time_t now = time(0);
	struct tm time_struct = *localtime(&now);

	fprintf(flog, "-- Start of simulation at: %2d:%02d:%02d -- \n\n\n",
		time_struct.tm_hour, time_struct.tm_min, time_struct.tm_sec);

	fprintf(flog, " demand file used: %s\n region file used: %s\n", demd, park);

	std::ifstream f1;
	f1.open(F_INIT);
	if (!f1) {
		show_error();
		exit(0);
	}

	f1 >> rncov >> rnreg >> rnfact >> rnatt >> rndemand;
	fprintf(flog, " number of cover types: %d\n number of regions: %d\n number of factors in 1 equation: %d\n total number of factors: %d\n number of demand types: %d\n", rncov, rnreg, rnfact, rnatt, rndemand);
	f1 >> rnr >> rnc;
	fprintf(flog, " grid specification (rows, columns): %d, %d\n", rnr, rnc);
	f1 >> cellsize;
	fprintf(flog, " cellsize: %15.10lf\n", cellsize);
	f1 >> xllcorner >> yllcorner;
	fprintf(flog, " xll and yll coordinates: %f, %f", xllcorner, yllcorner);
	fprintf(flog, "\n land use system types: ");
	maxcovdem = 0;
	for (i = 0; i<rncov; i++) {
		f1 >> covdem[i];
		fprintf(flog, "%d ", covdem[i]);
		if (covdem[i] > maxcovdem) {
			maxcovdem = covdem[i];
		}
	}

	fprintf(flog, "\n elasticity of land use types: ");
	for (i = 0; i < rncov; i++) {
		f1 >> stat[i];                               /* elasticity of land use types*/
		fprintf(flog, "%3.2f ", stat[i]);
	}
	fprintf(flog, "\n number of demand types: %d \n weight of demand types: ", rndemand);
	for (i = 0; i < rndemand; i++) {
		f1 >> weightdem[i];
		fprintf(flog, "%f ", weightdem[i]);
	}
	f1 >> itermode >> demdiff >> demdiffmax;
	fprintf(flog, "\n iteration mode: ");
	if (itermode == 0) {
		fprintf(flog, "percentage");
	} else {
		fprintf(flog, "absolute");
	}
	fprintf(flog, "\n maximum average difference in demand allowed: %f", demdiff);
	fprintf(flog, "\n maximum individual difference in demand allowed: %f", demdiffmax);
	f1 >> start >> end;
	fprintf(flog, "\n start: %d end: %d \n", start, end);
	f1 >> nochange;
	fprintf(flog, " changing driver numbers: ");
	for (i = 0; i < nochange; i++) {
		f1 >> grchange[i];
		fprintf(flog, "%d ", grchange[i]);
	}
	f1 >> arcview;
	fprintf(flog, "\n output switch: %d\n", arcview);
	f1 >> diffregregi;
	fprintf(flog, " region specific regression switch: %d\n", diffregregi);
	f1 >> agemode;
	fprintf(flog, " random calculation of land use history: ");
	if (agemode > 0) {
		f1 >> startage;
		if (agemode == 2) {
			//randomize();
			srand((unsigned int)time(NULL));
			fprintf(flog, "randomizer ");
		}
		fprintf(flog, "on   - maximum age: %d\n", startage);
	} else {
		fprintf(flog, "off\n");
	}
	f1 >> influence;
	fprintf(flog, " neighborhood functions: ");
	if (influence > 0) {
		fprintf(flog, "on\n");
		if ((rnatt + rncov) > NATT) {
			dimensionerror = 1;
		}
	} else {
		fprintf(flog, "off\n");
	}
	fprintf(flog, " location specific addition: ");
	f1 >> locationfactor;
	if (locationfactor > 0) {
		fprintf(flog, "on, weights: ");
		for (i = 0; i < rncov; i++) {
			f1 >> locationweight[i];
			fprintf(flog, "%3.2f ", locationweight[i]);
			if ((rnatt + rncov) > NATT)
				dimensionerror = 1;
		}
	} else {
		fprintf(flog, "off");
	}
	f1 >> changeconv;
	if (changeconv == 1) {
		fprintf(flog, "\n dynamic lusmatrix");
	}
	f1 >> writeoutputs;
	if (writeoutputs == 1) {
		fprintf(flog, "\n write output grids if outputmatrix.txt is present");
	}
	f1 >> seed;
	if ((seed < 0.001) || (seed > 0.9)) {
		seed = 0.05f;
	}
	fprintf(flog, "\n seed for iteration: %f ", seed);
	f1 >> lusconvmap;
	if (lusconvmap == 1) {
		f1 >> nooutputmaps;
		fprintf(flog, "\n initial composition of land use systems will be read from %d maps called lusoutputmap.# ");
	}
	if ((rncov > NCOV) || (rnreg > NREG) || (rnfact > NFACT) || (rnatt > NATT) || (rnr > NR) || (rnc > NC) || (rndemand > NCOV) || (nooutputmaps > NOUTP)) {
		dimensionerror = 1;
	}
	if ((influence > 0) && (locationfactor > 0) && ((rnatt + rncov + rncov) > NATT)) {
		dimensionerror = 1;
	}
	if (dimensionerror > 0) {
		fprintf(flog, "\n MAXIMUM DIMENSIONS of attribute file EXCEEDED -- Exit");
		fclose(flog);
		show_error();
		exit(0);
	}
	f1.close();
	//gridsize = ((sqrt(cellsize))* 100);
	fprintf(flog, "\n cellsize: %f \n\n", cellsize);
}

void set_oldco()
{
	FILE *f7;
	char f_oldco[20];
	int i, j, k, nl, l, dum;
	float totalno, cum[NCOV];

	for (i = 0; i < rncov; i++) {
		cum[i] = 0;
	}
	totalno = 0;
	if (arcview == 3) {
		sprintf(f_oldco, "cov_all.%d.asc", (year - 1));
	} else {
		sprintf(f_oldco, "cov_all.%d", (year - 1));
	}
	if ((f7 = fopen(f_oldco, "r")) == NULL) {
		fprintf(flog, "\n %s not found \n", f_oldco);
		fclose(flog);
		show_error();
		exit(0);
	}
	rewind(f7);
	if (arcview > 0) {
		nl = 0;
		while (nl < 6) {
			dum = fgetc(f7);
			if (dum == '\n') nl++;
		}
	}
	for (j = 0; j < rnr; j++) {
		for (k = 0; k < rnc; k++) {
			fscanf(f7, "%d", &mat_oldco[j][k]);
			for (i = 0; i < rncov; i++) {
				if (i == mat_oldco[j][k]) {
					cum[i]++;
					totalno++;
				}
			}
		}
	}
	fclose(f7);
	fprintf(flog, "\n");

	for (i = 0; i < rncov; i++) {

		if (year == 1) {
			for (l = 0; l < rnreg; l++) {
				areacov[0][l][i] = 0;
			}
			areacov[0][0][i] = (cum[i]);
		}
		ratio[i] = cum[i] / totalno;
		fprintf(flog, "fraction occupied by land use type %d: %f\n", i, ratio[i]);
		if (ratio[i] == 0) {
			ratio[i] = 9;
		}
	}
}

void load_reg()
{
	FILE *f2;
	char F_ALREG[11];
	int i,j,k,l,tempcov;

	fprintf(flog,"\n regression parameters: \n");
	sprintf(F_ALREG,"alloc1.reg");
	if ((f2 = fopen(F_ALREG, "r")) == NULL) {
		fprintf(flog,"no file: alloc1.reg");
		fclose(flog); show_error(); exit(0);
	}
	if (diffregregi >= 1) {
		l=rnreg;
	} else {
		l=1;
	}
	for (i=0; i<l; i++) {
		for (j=0; j<rncov; j++) {
			fscanf(f2,"%d",&tempcov);
			fscanf(f2,"%lf",&constant[i][tempcov]);
			fscanf(f2,"%d",&no_fact[i][tempcov]);
			for (k=0; k<no_fact[i][tempcov]; k++) {
				fscanf(f2,"%f%d",&fact[i][tempcov][k],&attr[i][tempcov][k]);
			}
			fprintf(flog,"loc:  %d  %lf  %d \n",tempcov,constant[i][tempcov],no_fact[i][tempcov]);
		}
	}
	fclose(f2);
}

void load_reg2()
{
	FILE *f2;
	int i,j,k,l,tempcov;

    fprintf(flog,"\n regression parameters neighborhood functions: \n");
	if ((f2 = fopen("alloc2.reg", "r")) == NULL) {
		fprintf(flog,"no file: alloc2.reg");
		fclose(flog);
		show_error();
		exit(0);
	}
	if (diffregregi >= 1) {
		l=rnreg;
	} else {
		l=1;
	}

	for (i=0; i<l; i++) {
		for (j=0; j<rncov; j++) {
			fscanf(f2,"%d",&tempcov);
			fscanf(f2,"%lf",&constant2[i][tempcov]);
			fscanf(f2,"%d",&no_fact2[i][tempcov]);
			for (k=0; k<no_fact2[i][tempcov]; k++) {
				fscanf(f2,"%f%d",&fact2[i][tempcov][k],&attr2[i][tempcov][k]);
                                attr2[i][tempcov][k] += rnatt;
            }
			fprintf(flog,"nbh:  %d  %lf  %d  \n",tempcov,constant2[i][tempcov],no_fact2[i][tempcov]);
		}
	}
	fclose(f2);
}

void load_locationfactor()
{
	FILE *f3;
	int i, j, k, nl, dum;
	float l;
	char filename[20];

	for (i = 0; i < rncov; i++) {
		if (arcview == 3) {
			sprintf(filename, "locspec%d.fil.asc", i);
		} else {
			sprintf(filename, "locspec%d.fil", i);
		}
		if ((f3 = fopen(filename, "r")) == NULL) {
			fprintf(flog, " locspec* file missing \n");
			fclose(flog);
			show_error();
			exit(0);
		}
		if (arcview > 0) {
			nl = 0;
			while (nl < 6) {
				dum = fgetc(f3);
				if (dum == '\n') nl++;
			}
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				fscanf(f3, "%f ", &l);
				if (influence == 0) {
					mat_at[rnatt + i][j][k] = l;
				} else {
					mat_at[rnatt + rncov + i][j][k] = l;
				}
			}
		}
		fclose(f3);
	}
}

void read_allowed()
{
	int i, j, countbottom;
	FILE *f9;

	bottomup = 0;
	if ((f9 = fopen("allow.txt", "r")) == NULL) {
		fprintf(flog, "no file: allow.txt");
		fclose(flog);
		show_error();
		exit(0);
	}
	fprintf(flog, "\n allowed changes matrix\n");
	for (i = 0; i < rncov; i++) {
		countbottom = 0;
		for (j = 0; j < rncov; j++) {
			fscanf(f9, "%d", &allowed[i][j]);
			if (allowed[i][j] > 1000) {
				bottomup = 1;
				countbottom++;
				if (countbottom > 1) {
					fprintf(flog, "\n\n error reading allow.txt: double specification of autonomous change");
					fclose(flog);
					show_error();
					exit(0);
				}
			}
			fprintf(flog, "%d ", allowed[i][j]);
		}
		fprintf(flog, "\n");
	}
	if (bottomup == 1) {
		fprintf(flog, "\n - demand may be modified by bottom up interactions as specified in the change matrix - \n");
	}

	fclose(f9);
}

void make_mat()
{
	int i, j, result, totnatt;
	result = 1;

	totnatt = rnatt - allowmin;
	if ((locationfactor > 0) && (metamod == 0)) totnatt += rncov;
	if (influence > 0) totnatt += rncov;
	for (j = 0; j < totnatt; j++) {
		mat_at[j] = (fvector *)calloc(rnr, sizeof(fvector));
		if (mat_at == NULL) result = 0;
		for (i = 0; i < rnr; i++) {
			mat_at[j][i] = (fvector)calloc(rnc, sizeof(MATFLOAT));
			if (mat_at[j][i] == NULL) {
				result = 0;
			}
		}
	}
	if (lusconvmap == 1) {
		for (j = 0; j < nooutputmaps; j++) {
			initoutput[j] = (fvector *)calloc(rnr, sizeof(fvector));
			if (initoutput == NULL) result = 0;
			for (i = 0; i < rnr; i++) {
				initoutput[j][i] = (fvector)calloc(rnc, sizeof(MATFLOAT));
				if (initoutput[j][i] == NULL) {
					result = 0;
				}
			}
		}
		luschange = (vector *)calloc(rnr, sizeof(vector));
		if (luschange == NULL) result = 0;
		for (i = 0; i < rnr; i++) {
			luschange[i] = (vector)calloc(rnc, sizeof(MATELEMS));
			if (luschange[i] == NULL) {
				result = 0;
			}
		}
	}

	for (j = 0; j < rncov; j++) {
		mat_co[j] = (fvector *)calloc(rnr, sizeof(fvector));
		if (mat_co == NULL) result = 0;
		for (i = 0; i < rnr; i++) {
			mat_co[j][i] = (fvector)calloc(rnc, sizeof(MATFLOAT));
			if (mat_co[j][i] == NULL) {
				result = 0;
			}
		}
	}

	mat_newco = (vector *)calloc(rnr, sizeof(vector));
	if (mat_newco == NULL) {
		result = 0;
	}
	for (i = 0; i < rnr; i++) {
		mat_newco[i] = (vector)calloc(rnc, sizeof(MATELEMS));
		if (mat_newco[i] == NULL) {
			result = 0;
		}
	}

	region = (vector *)calloc(rnr, sizeof(vector));
	if (region == NULL) result = 0;
	for (i = 0; i < rnr; i++) {
		region[i] = (vector)calloc(rnc, sizeof(MATELEMS));
		if (region[i] == NULL) {
			result = 0;
		}
	}
	mat_oldco = (vector *)calloc(rnr, sizeof(vector));
	if (mat_oldco == NULL) {
		result = 0;
	}
	for (i = 0; i < rnr; i++) {
		mat_oldco[i] = (vector)calloc(rnc, sizeof(MATELEMS));
		if (mat_oldco[i] == NULL) {
			result = 0;
		}
	}

	mat_age = (vector *)calloc(rnr, sizeof(vector));
	if (mat_age == NULL) {
		result = 0;
	}
	for (i = 0; i < rnr; i++) {
		mat_age[i] = (vector)calloc(rnc, sizeof(MATELEMS));
		if (mat_age[i] == NULL) {
			result = 0;
		}
	}

	if (result == 0) {
		fprintf(flog, "Make_mat out of memory\n");
		fclose(flog);
		show_error(); exit(0);
	}
}

void scgr_change()
{
	FILE *pf1, *pf1a;
	int i, j, k, nl, dum;
	float /*c,*/l;
	char inname[20], outname[20], c;
	fprintf(flog, "\n update of driving factor files for drivers: ");

	for (i = 0; i < nochange; i++) {
		if ((grchange[i] >= allowmin) || (metamod == 0)) {
			if (arcview == 3) {
				sprintf(inname, "sc1gr%d.%d.asc", grchange[i], year);
				sprintf(outname, "sc1gr%d.fil.asc", grchange[i]);
			} else {
				sprintf(inname, "sc1gr%d.%d", grchange[i], year);
				sprintf(outname, "sc1gr%d.fil", grchange[i]);
			}
			fprintf(flog, "%d ", grchange[i]);
			if ((pf1 = fopen(inname, "r")) == NULL) {
				fprintf(flog, "\n\n\n\n NO NEW FILE %s \n\n", inname);
				fclose(flog);
				show_error();
				exit(0);
			}
			if ((pf1a = fopen(outname, "w")) == NULL) {
				fprintf(flog, "\n\n\n\n\n ERROR: %s \n\n\n\n\n", outname);
				fclose(flog);
				show_error();
				exit(0);
			}
			while ((c = getc(pf1)) != EOF) {
				putc(c, pf1a);
			}
			fclose(pf1);
			fclose(pf1a);
			if ((pf1a = fopen(outname, "r")) == NULL) {
				fprintf(flog, "\n\n\n\n\n ERROR: %s \n\n\n\n\n", outname);
				fclose(flog);
				show_error();
				exit(0);
			}
			if (arcview > 0) {
				nl = 0;
				while (nl < 6) {
					dum = fgetc(pf1a);
					if (dum == '\n') {
						nl++;
					}
				}
			}
			for (j = 0; j < rnr; j++) {
				for (k = 0; k < rnc; k++) {
					fscanf(pf1a, "%f ", &l);
					mat_at[(grchange[i] - allowmin)][j][k] = l;
				}
			}
			fclose(pf1a);
		}
	}
	fprintf(flog, "\n");
}

void load_lusconv()
{
	FILE* f9;
	int i, j;

	if ((f9 = fopen("lusconv.txt", "r")) == NULL) {
		fprintf(flog, "no file: lusconv.txt");
		fclose(flog);
		show_error();
		exit(0);
	}
	fprintf(flog, "\n\n land use system conversion matrix\n");
	for (i = 0; i < rncov; i++) {
		for (j = 0; j < rndemand; j++) {
			fscanf(f9, "%d", &lusconv[i][j]);
			fprintf(flog, "%d ", lusconv[i][j]);
		}
		fprintf(flog, "\n");
	}
	fprintf(flog, "\n\n");
}

void load_lusmatrix()
{
	FILE* f9;
	int i, j;
	char inname[20];

	if (changeconv == 1) {
		sprintf(inname, "lusmatrix.%d", year);
		if ((f9 = fopen(inname, "r")) == NULL) {
			fprintf(flog, "no file: lusmatrix.year");
			fclose(flog);
			show_error();
			exit(0);
		}
	} else {
		if ((f9 = fopen("lusmatrix.txt", "r")) == NULL) {
			fprintf(flog, "no file: lusmatrix.txt");
			fclose(flog);
			show_error();
			exit(0);
		}
	}

	fprintf(flog, "\n land use system matrix\n");
	for (i = 0; i < rncov; i++) {
		for (j = 0; j < rndemand; j++) {
			fscanf(f9, "%f", &lusmatrix[i][j]);
			fprintf(flog, "%f ", lusmatrix[i][j]);
		}
		fprintf(flog, "\n");
	}

	fprintf(flog, "\n\n");
}

void load_grid()
{
	FILE *f3;
	int i, j, k, nl, dum;
	float l;
	char filename[20];

	for (i = (allowmin); i < rnatt; i++) {
		if (arcview == 3) {
			sprintf(filename, "sc1gr%d.fil.asc", i);
		} else {
			sprintf(filename, "sc1gr%d.fil", i);
		}
		if ((f3 = fopen(filename, "r")) == NULL) {
			fprintf(flog, " sc1gr* file missing \n");
			fclose(flog);
			show_error();
			exit(0);
		}
		if (arcview > 0) {
			nl = 0;
			while (nl < 6) {
				dum = fgetc(f3);
				if (dum == '\n') nl++;
			}
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				fscanf(f3, "%f ", &l);
				mat_at[(i - allowmin)][j][k] = l;
			}
		}
		fclose(f3);
	}
}

void load_region()
{
	//char F_REGION[MAX_PATH];
	FILE *f4;
	int i,j,m,tempregion, nl,dum;

    for (m=0;m<rnreg;m++) {
		nocells[m] = 0;
	}

    //sprintf(F_REGION,"%s", park);

	if ((f4 = fopen(park, "r")) == NULL) {
		fprintf(flog,"no file: %s", park);
		fclose(flog);
		show_error();
		exit(0);
	}
	if (arcview > 0) {
		nl = 0;
		while (nl < 6) {
			dum = fgetc(f4);
			if (dum == '\n') nl++;
		}
	}
	for (i = 0; i < rnr; i++) {
		for (j = 0; j < rnc; j++) {
			fscanf(f4, "%d", &region[i][j]);
			if (region[i][j] < -9900) {
				tempregion = region[i][j] + 9998;
			} else {
				tempregion = region[i][j];
			}
			if (tempregion > -1) {
				nocells[tempregion]++;
			}
		}
	}
	for (m = 0; m < rnreg; m++) {
		fprintf(flog, "\n number of cells in region %d: %d\n", m, nocells[m]);
	}
	fclose(f4);
}

void demand_read()
{
	int i, j, l, nyear, tempreg;
	FILE *in;

	if (diffregregi < 2) {
		tempreg = rnreg;
	} else {
		tempreg = 1;
	}
	if ((in = fopen(demd, "r")) == NULL) {
		show_error();
		fprintf(flog, "/n no demand file");
		fclose(flog);
		exit(0);
	}
	fscanf(in, "%d", &nyear);
	for (l = 0; l < tempreg; l++) {
		for (i = 0; i < nyear; i++) {
			for (j = 0; j < rndemand; j++) {
				fscanf(in, "%lf", &demand[l][j][i]);
			}
		}
	}
	fclose(in);
}

void demand_dir()
{
	int i, j, tempreg;

	/* determine if demand increases or decreases */
	if (diffregregi < 2) {
		tempreg = rnreg;
	} else {
		tempreg = 1;
	}

	for (j = 0; j < tempreg; j++) {
		nostable[j] = 0;
	}

	for (i = 0; i < rndemand; i++) {
		fprintf(flog, "demand type: %d\\ ", i);
		for (j = 0; j < tempreg; j++) {
			if (demand[j][i][year - 1] < demand[j][i][year]) {
				dem_dir[j][i] = 1;
			}
			if (demand[j][i][year - 1] == demand[j][i][year]) {
				dem_dir[j][i] = 0;
				nostable[j]++;
			}
			if (demand[j][i][year - 1] > demand[j][i][year]) {
				dem_dir[j][i] = -1;
			}
			fprintf(flog, " demand direction for region %d is %d; demand:\t%8.1lf\n", j, dem_dir[j][i], demand[j][i][year]);
		}
		fprintf(flog, "\n");
	}
	fprintf(flog, "\n");
	for (j = 0; j < tempreg; j++) {
		for (i = 0; i < rndemand; i++) {
			speed[j][i] = seed;
		}
	}
}

void init_iter()
{
   int i,j,tempreg;

   if (diffregregi < 2) {
        tempreg = rnreg;
   } else {
        tempreg = 1;
   }
   if (year == 1) {
	   for (i = 0; i < rndemand; i++) {
		   for (j = 0; j < tempreg; j++) {
			   dem_elas[j][i] = 0;
			   speed[j][i] = seed;
		   }
	   }
   }
   loop = 0;
}

/*
 * Calculates logistic regressions.
 */
void calc_reg()
{
	int i, j, k, m, reg, nl, dum, temploc;
	float temp, temp2;
	FILE *fa7;
	char f_prob[20];

	if (influence > 0) {
		temploc = rncov;
	} else {
		temploc = 0;
	}

	if (metamod == 1) {
		for (i = 0; i < rncov; i++) {
			if (arcview == 3) {
				sprintf(f_prob, "prob1_%d.%d.asc", i, year);
			} else {
				sprintf(f_prob, "prob1_%d.%d", i, year);
			}
			if ((fa7 = fopen(f_prob, "r")) == NULL) {
				fprintf(flog, "\n %s not found \n", f_prob);
				fclose(flog);
				show_error();
				exit(0);
			}
			rewind(fa7);
			if (arcview > 0) {
				nl = 0;
				while (nl < 6) {
					dum = fgetc(fa7);
					if (dum == '\n') {
						nl++;
					}
				}
			}
			for (j = 0; j < rnr; j++) {
				for (k = 0; k < rnc; k++) {
					fscanf(fa7, "%f", &mat_co[i][j][k]);
				}
			}
			fclose(fa7);
		}
	}

	/* Just calculate the probability as determined by the regression equations, put in mat_co */

	if ((metamod == 0) || (influence > 0)) {
		for (i = 0; i < rncov; i++) {
			for (j = 0; j < rnr; j++) {
				for (k = 0; k < rnc; k++) {
					if (region[j][k] < 0) {
						mat_co[i][j][k] = -9999;
					} else {
						if (diffregregi >= 1) {
							reg = region[j][k];
						} else {
							reg = 0;
						}
						if (metamod == 0) {
							temp = constant[reg][i];
							for (m = 0; m < no_fact[reg][i]; m++) {
								temp += (fact[reg][i][m] * mat_at[attr[reg][i][m]][j][k]);          /*Preference = x0+x1*fac1+...Xn*facn*/
							}
							if (temp > 87) {
								fprintf(flog, "\n\n error: regression cannot be calculated due to large value in cell %d,%d for land cover %d \n", j, k, i);
								fclose(flog);
								show_error();
								exit(0);
							}
							mat_co[i][j][k] = (exp(temp) / (1 + exp(temp)));                            /*Probability = (exp(pref)/(1+exp(pref)*/
							if (loop == 5000) {
								mat_co[i][j][k] += ((rand() % 10000) / 1000000);
							}
						}

						if ((influence == 1) && (allprob == 0) && (neighweight[i] > 0) && (allocneigh[i][mat_oldco[j][k]] == 1)) {
							// neighborhood calculation on; only for land systems with weight larger than 0 //
							temp2 = constant2[reg][i];
							for (m = 0; m < no_fact2[reg][i]; m++) {
								if ((attr2[reg][i][m] - rnatt) < 99) {
									temp2 += (fact2[reg][i][m] * mat_at[(attr2[reg][i][m] - allowmin)][j][k]);
								} else {
									temp2 += (fact2[reg][i][m] * mat_at[(attr2[reg][i][m] - allowmin - rnatt - 100)][j][k]);
								}
							}
							//    if (temp2 > 700) {fprintf(flog,"\n\n error: nbh regression cannot be calculated due to large value in cell %d,%d for land cover %d \n",j,k,i);fclose(flog); show_error(); exit(0);}

							mat_co[i][j][k] = ((mat_co[i][j][k]) + (neighweight[i] * temp2));
						}

						if (metamod == 0) {
							// location specific addition
							if (locationfactor > 0) {
								mat_co[i][j][k] += (mat_at[(rnatt + temploc + i)][j][k] * locationweight[i]);
							}
						}
					}
				}
			}
		}
	}
}

void free_mat()
{
	int i, j, totnatt;

	totnatt = rnatt - allowmin;
	if ((locationfactor > 0) && (metamod == 0)) totnatt += rncov;
	if (influence > 0) totnatt += rncov;
	for (i = 0; i < totnatt; i++) {
		for (j = 0; j < rnr; j++) {
			free(mat_at[i][j]);
		}
		free(mat_at[i]);
	}
	for (i = 0; i < rncov; i++) {
		for (j = 0; j < rnr; j++) {
			free(mat_co[i][j]);
		}
		free(mat_co[i]);
	}
	if (lusconvmap == 1) {
		for (i = 0; i < nooutputmaps; i++) {
			for (j = 0; j < rnr; j++) {
				free(initoutput[i][j]);
			}
			free(initoutput[i]);
		}
	}
	for (j = 0; j < rnr; j++) {
		free(mat_newco[j]);
		free(region[j]);
		free(mat_age[j]);
		free(mat_oldco[j]);
		if (lusconvmap == 1) {
			free(luschange[j]);
		}
	}
	free(mat_newco);
	free(region);
	free(mat_age);
	free(mat_oldco);
	if (lusconvmap == 1) {
		free(luschange);
	}
}

/*
 * Calculates the actual changes in land use systems.
 */
void calc_change_ch()
{
	int i, j, k, maxcov, oldcov, tempreg, cdeman, tempallow;
	float max, tempmax, sumelas;

	for (j = 0; j < rnr; j++) {
		for (k = 0; k < rnc; k++) {
			maxcov = -1;
			// loop over all rows and columns
			if ((region[j][k] > -9999) && (region[j][k] < -9900)) {
				// If pixel belongs to restricted area allocate the old land use system
				mat_newco[j][k] = mat_oldco[j][k];
			}
			if (region[j][k] > -9900) {
				// If pixel is outside protected area or outside nodata region
				oldcov = mat_oldco[j][k];             // assign variable oldcov the old land use system at t-1
				max = -30;                             // set initial value for maximum total transition probability

				if (diffregregi > 1) {
					// In case of different demands by region but allocation for the whole area set the region number of the pixel to region 0
					tempreg = 0;
				} else {
					tempreg = region[j][k];
				}
				for (i = 0; i < rncov; i++) {
					// loop all land use systems to determine which land use system has the highest transition probability
					sumelas = 0;
					// set the iteration variable to be used for this land use system to 0
					for (cdeman = 0; cdeman < rndemand; cdeman++) {
						// loop over all demand_types to calculate the iteration variable value for the land use system; only include those land use systems that have a value > -1 in the lusconv.txt file
						if (lusconv[i][cdeman] != -1) {
							if (lusconv[oldcov][cdeman] < lusconv[i][cdeman]) {
								sumelas += (weightdem[cdeman] * dem_elas[tempreg][cdeman]);
								//in case the new potential land use system provides more for this demand type increase the iteration variable (corrected for the weight of the demand type)
							}
							if (lusconv[oldcov][cdeman] > lusconv[i][cdeman]) {
								sumelas -= (weightdem[cdeman] * dem_elas[tempreg][cdeman]);
								//in case the new potential land use system provides less for this demand type decrease the iteration variable (corrected for the weight of the demand type)
							}
						}
					}
					if (oldcov == i) {
						// In case the land use system equals the current land use system and changes are allowed
						// the total probability is the location suitability, iteration variable (=elas) and resistance (=stat)
						tempmax = (mat_co[i][j][k] + sumelas + stat[i]); 
					} else {
						// In all other cases do not add the resistance because it does not apply
						tempmax = (mat_co[i][j][k] + sumelas);
					}
					//if ((dem_dir[tempreg][i] == -1) && (oldcov != i) && (stat[i] == 1))     // in case of a decreasing demand and 'stability' allocate a very low value to tempmax for all other land use systems
					//        tempmax = -30;


					if (allowed[oldcov][i] != 1) {
						// in case the conversion is not fully allowed
						tempallow = 1;
						if (allowed[oldcov][i] == 0) {
				           // in case the conversion matrix value is 0 give a very low tempmax is given to this land use system
							tempmax = -30;
						}
						if ((allowed[oldcov][i] > 1) && (allowed[oldcov][i] < 100)) {
							// In case the conversion matrix refers to a map, read the value from that map
							tempallow = mat_at[(allowed[oldcov][i] - allowmin)][j][k];
						}
						if (tempallow == 0) {
							// In case the conversion matrix value of the map is 0 give  a very low tempmax to this land use system
							tempmax = -30;
						}
						if (abs(allowed[oldcov][i]) > 100) {
							tempallow = allowed[oldcov][i];
						}
						if ((tempallow > 100) && (tempallow < 1000) && (mat_age[j][k] < (tempallow - 100))) {
							// if age of the land use system does not fulfill the conversion setting requirements:
							// low value of tempmax
							tempmax = -30;
						}
						if ((tempallow < -100) && (mat_age[j][k] >= ((abs(tempallow)) - 100))) {
							// if age of the land use system does not fulfill the conversion setting requirements:
							// low value of tempmax
							tempmax = -30;
						}
						if ((tempallow > 1000) && (mat_age[j][k] < (tempallow - 1000))) {
							// In case of succession/autonomous conversions give low value tempmax
							// if conditions for autonomous conversion are not yet fulfilled
							tempmax = -30;
						}
						if (tempallow > 1000) {
							// flag bottom-up conversions for succession/autonomous change calculations
							bottomup = 1;
						}
					}

					if (tempmax > max) {
						// In case the resulting value of tempmax is higher than of previous land use systems
						// assign the value to max and assign the new land use system.
						max = tempmax;
						// The current land use system is (for now) the one with the now highest tempmax.
						maxcov = i;
					}

				}
				if (max == -30) {
					// Give error in case no change is single change is allowed (tempmax for all land systems at minimal value) while stability is not requested
					fprintf(flog, "\n too many constraints defined in model input files"); fclose(flog); show_error(); exit(0);
				}

				mat_newco[j][k] = maxcov;
			}
			if (region[j][k] <= -9999) {
				// if nodata, allocate nodata
				mat_newco[j][k] = -9999;
			}
		}
	}
}

/*
 * not tested in clumondo: superimposes change for those land uses that are assigned for autonomous change in the allow matrix
 */
void autonomous_change()
{
	int i, j, k, oldcov, tempallow, countbottom;
	for (j = 0; j < rnr; j++) {
		for (k = 0; k < rnc; k++) {
			if (region[j][k] > -9999) {
				oldcov = mat_oldco[j][k];
				countbottom = 0;

				for (i = 0; i < rncov; i++) {
					if (allowed[oldcov][i] != 1) {
						tempallow = 1;
						if ((allowed[oldcov][i] > 1) && (allowed[oldcov][i] < 100)) {
							tempallow = mat_at[(allowed[oldcov][i] - allowmin)][j][k];
						}
						if (allowed[oldcov][i] > 1000) {
							tempallow = allowed[oldcov][i];
						}
						if (tempallow > 1000) {
							countbottom++;
							if (countbottom > 1) { fprintf(flog, "\n error: same land cover cannot change autonomous in two different new land covers"); fclose(flog); show_error(); exit(0); }
						}
						if (mat_newco[j][k] == oldcov) {
							if ((tempallow > 1000) && (mat_age[j][k] >= (tempallow - 1000))) {
								mat_newco[j][k] = i;
							}
						}
					}
				}
			}
		}
	}
}

/*
 * Compares the demand with the allocated area.
 */
void comp_demand()
{
	int i, j, k, l, m, tempreg;
	float ran[NREG][NCOV];
	double area[NREG][NCOV];

	/* calculate allocated cover */
	if (abs(arcview) != 2) {
		// output to log file control
		fprintf(flog, "\n%d\t", loop);
	}
	tempreg = 1;
	if (diffregregi < 2) {
        // flag in case different regions have different demands
		tempreg = rnreg;
	}

	// loop this for all demand_types
	for (i = 0; i < rndemand; i++) {
        // loop this for all regions
		for (l = 0; l < tempreg; l++) {
			area[l][i] = 0;             // set allocated quantity to 0
			speed[l][i] += 0.0002f;     // increase step-size iteration variable

			// initialize random number, add to avoid seed*ran being < 1
			srand(101 - (1 / speed[l][i]));
			ran[l][i] = rand() + 1 / speed[l][i];
		}
	}

    // loop all land use systems
	for (i = 0; i < rncov; i++) {
		// loop all regions
		for (l = 0; l < tempreg; l++) {
			// put area of each land use system to 0
			areacov[year][l][i] = 0;
		}
	}

	// loop all rows and columns
	for (j = 0; j < rnr; j++) {
		for (k = 0; k < rnc; k++) {
			if (diffregregi > 1) {
				m = 0;
			} else {
                // determine region number (=m) for restricted areas
				if (region[j][k] < -9900) {
					m = region[j][k] + 9998;
				} else {
					m = region[j][k];
				}
			}
			if (m > -1) {
				// loop all demand_types
				for (l = 0; l < rndemand; l++) {
					// calculate the allocated quantity of the demand types by
					// cellsize * demand contribution of LUS
					area[m][l] += lusmatrix[mat_newco[j][k]][l];
				}
				// add to area of land use system
				areacov[year][m][mat_newco[j][k]] += 1;
			}
		}
	}

    // log-file control, only print if requested
	if (abs(arcview) != 2) {
		for (i = 0; i < rndemand; i++) {
			for (l = 0; l < tempreg; l++) {
				// print in log-file the demand allocated in iteration step of each demand_type
				fprintf(flog, "%8.1f\t", area[l][i]);
			}
			fprintf(flog, "-\t");
		}
		fprintf(flog, "*\t");
	}

	/* adapt iteration factor */
	maxdiff = 0;
	totdiff = 0;

	for (i = 0; i < rndemand; i++)                                          // loop demand_types
	{
		for (l = 0; l < tempreg; l++)                                      // loop regions
		{
			dem_elas[l][i] -= (((area[l][i] - demand[l][i][year]) / demand[l][i][year]) / (speed[l][i] * ran[l][i]));   // adapt iteration varaible based on difference and adapt step size and add small random component
			if (fabs(dem_elas[l][i]) > 2)                     // if iteration variable gets too big set to 0
				dem_elas[l][i] = 0;
			if (abs(arcview) != 2)
				fprintf(flog, "%4.3f\t", dem_elas[l][i]);       // write iteration variables to log file
		}
	}


	for (i = 0; i<rndemand; i++)                  //loop demand types
	{
		for (l = 0; l < tempreg; l++)              // loop regions
		{
			if (itermode == 0)
			{
				if (demand[l][i][year] > 0)        // calculate the average deviation between demand and allocated and the maximum difference
				{
					if (fabs(((area[l][i] - demand[l][i][year]) / demand[l][i][year] * 100)) > maxdiff)
						maxdiff = fabs(((area[l][i] - demand[l][i][year]) / demand[l][i][year] * 100));
					totdiff += (fabs(((area[l][i] - demand[l][i][year]) / demand[l][i][year] * 100)));
				} else {
					if (fabs(((area[l][i] - demand[l][i][year]) / (rnr*rnc) * 100)) > maxdiff) {
						maxdiff = fabs(((area[l][i] - demand[l][i][year]) / (rnr*rnc) * 100));
					}
					totdiff += (fabs(((area[l][i] - demand[l][i][year]) / (rnr*rnc) * 100)));
				}
			}
			if (itermode == 1) {
				if (fabs(area[l][i] - demand[l][i][year]) > maxdiff) {
					maxdiff = fabs(area[l][i] - demand[l][i][year]);
				}
				totdiff += fabs(area[l][i] - demand[l][i][year]);
			}
		}

	}
	if (abs(arcview) != 2) {
		// print the calculated differences to log-file
		fprintf(flog, "*\t%4.2f\t%4.2f\t", maxdiff, (totdiff / (rndemand)));
	}

	/*
	if (itermode == 0) {
        // calculate bar position in user interface
		Barpos = ((10 * (log(((totdiff + 0.00001) * 100) / (rndemand)))) - 10);
	}
	if (itermode == 1) {
		i = (totdiff / (rndemand)) - demdiff;
		if (i <= 0) i = 1;
		Barpos = (5 * log(i));
	}
	*/
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void write_grid()
{
	FILE *f7, *f8, *f9, *f10, *f11;
	int i, j, k, tm, l, tempreg, nooutputs;
	char filenam4[20], filenam5[20], filenam6[20], filenam7[20];
	double tempdemand, tempoutput;
	static double outputmatrix[NCOV][NCOV][100];

	/* read outputmatrix input file */
	nooutputs = 0;
	sprintf(filenam7, "outputmatrix.%d", year);
	if (((f11 = fopen("outputmatrix.txt", "r")) != NULL) || ((f11 = fopen(filenam7, "r")) != NULL)) {
		fscanf(f11, "%d", &nooutputs);
		if ((lusconvmap == 1) && (nooutputs > nooutputmaps))  { fprintf(flog, "\n more output types included than output maps \n"); fclose(flog); show_error(); exit(0); }
		for (i = 0; i < rncov; i++) {
			for (j = 0; j < nooutputs; j++) {
				fscanf(f11, "%lf", &outputmatrix[i][j][year]);
			}
		}
		fclose(f11);
	}

	if (year == 1) {
		if (((f11 = fopen("outputmatrix.txt", "r")) != NULL) || ((f11 = fopen("outputmatrix.0", "r")) != NULL)) {
			fscanf(f11, "%d", &nooutputs);
			for (i = 0; i < rncov; i++) {
				for (j = 0; j < nooutputs; j++) {
					fscanf(f11, "%lf", &outputmatrix[i][j][0]);
				}
			}
			fclose(f11);
		}
	}

	/* write summary tables if in end year */
	if (year == (end - start)) {

		/* write summary table of allocated land use system areas */

		if ((f10 = fopen("landarea.txt", "w")) == NULL) { fprintf(flog, "\n not able to make file landarea.txt"); fclose(flog); show_error(); exit(0); }
		fprintf(f10, "area allocated to each of the land use systems\n");
		if (diffregregi < 2) {
			tempreg = rnreg;
		} else {
			tempreg = 1;
		}
		for (i = 0; i <= (end - start); i++) {
			fprintf(f10, "%d\t", i);
			for (j = 0; j < rncov; j++) {
				tempdemand = 0;
				for (l = 0; l < tempreg; l++) {
					tempdemand += areacov[i][l][j];
				}
				fprintf(f10, "%12.4lf\t", (tempdemand));
			}
			fprintf(f10, "\n");
		}
		/* also write outputs for other outputvariables if a file outputmatrix.txt exists */
		if (nooutputs > 0) {
			fprintf(f10, "\n\nquantity allocated to the output types as specified in outputmatrix.txt\n");
			for (tm = 0; tm <= (end - start); tm++) {
				fprintf(f10, "%d\t", tm);
				for (j = 0; j < nooutputs; j++) {
					tempoutput = 0;
					for (i = 0; i < rncov; i++) {
						tempdemand = 0;
						for (l = 0; l < tempreg; l++) {
							tempdemand += areacov[tm][l][i];
						}
						tempoutput += (tempdemand * outputmatrix[i][j][tm]);
					}
					fprintf(f10, "%12.4lf\t", tempoutput);
				}
				fprintf(f10, "\n");
			}
		}
		fclose(f10);
	}

	/* export output land system maps, age file and influence*/

	if ((probmaps != 1) && (allprob == 0)) {
		/* write land system maps */
		if (arcview == 3) {
			sprintf(filenam5, "cov_all.%d.asc", year);
		} else {
			sprintf(filenam5, "cov_all.%d", year);
		}
		if ((f8 = fopen(filenam5, "w")) == NULL)  {
			fprintf(flog, "\n not able to write file \n");
			fclose(flog);
			show_error();
			exit(0);
		}
		if (arcview > 0) {
			fprintf(f8, "ncols    %d\n", rnc);
			fprintf(f8, "nrows    %d\n", rnr);
			fprintf(f8, "xllcorner  %f\n", xllcorner);
			fprintf(f8, "yllcorner  %f\n", yllcorner);
			fprintf(f8, "cellsize   %f\n", cellsize);
			fprintf(f8, "NODATA_value   -9999\n");
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				fprintf(f8, "%d\n", mat_newco[j][k]);
			}
		}
		fclose(f8);

		/* write outputmatrix maps for year 0 */
		if ((writeoutputs > 0) && (year == 1) && (nooutputs > 0)) {
			for (l = 0; l < nooutputs; l++) {
				if (arcview == 3) {
					sprintf(filenam6, "output%d.0.asc", l);
				} else {
					sprintf(filenam6, "output%d.0", l);
				}
				if ((f9 = fopen(filenam6, "w")) == NULL) {
					fprintf(flog, "\n not able to write file %s\n", filenam6);
					fclose(flog);
					show_error();
					exit(0);
				}
				if (arcview > 0) {
					fprintf(f9, "ncols    %d\n", rnc);
					fprintf(f9, "nrows    %d\n", rnr);
					fprintf(f9, "xllcorner  %f\n", xllcorner);
					fprintf(f9, "yllcorner  %f\n", yllcorner);
					fprintf(f9, "cellsize   %f\n", cellsize);
					fprintf(f9, "NODATA_value   -9999\n");
				}
				for (j = 0; j < rnr; j++) {
					for (k = 0; k < rnc; k++) {
						if (mat_newco[j][k] >= 0) {
							if (lusconvmap == 1) {
								// If option to read initial maps is on write the value of the initial map
								fprintf(f9, "%10.4f\n", initoutput[l][j][k]);
							} else {
								fprintf(f9, "%10.4f\n", outputmatrix[mat_newco[j][k]][l][year]);
							}
						} else {
							fprintf(f9, "-9999\n");
						}
					}
				}
				fclose(f9);
			}
		}
		/* write outputmatrix maps if requested*/
		if ((nooutputs > 0) && (writeoutputs == 1) && (year == (end - start)) || (writeoutputs == 2)) {
			for (l = 0; l < nooutputs; l++) {
				if (arcview == 3) {
					sprintf(filenam6, "output%d.%d.asc", l, year);
				} else {
					sprintf(filenam6, "output%d.%d", l, year);
				}
				if ((f9 = fopen(filenam6, "w")) == NULL) {
					fprintf(flog, "\n not able to write file %s\n", filenam6);
					fclose(flog);
					show_error();
					exit(0);
				}
				if (arcview > 0) {
					fprintf(f9, "ncols    %d\n", rnc);
					fprintf(f9, "nrows    %d\n", rnr);
					fprintf(f9, "xllcorner  %f\n", xllcorner);
					fprintf(f9, "yllcorner  %f\n", yllcorner);
					fprintf(f9, "cellsize   %f\n", cellsize);
					fprintf(f9, "NODATA_value   -9999\n");
				}
				for (j = 0; j < rnr; j++) {
					for (k = 0; k < rnc; k++) {
						if (mat_newco[j][k] >= 0) {
							if ((lusconvmap == 1) && (luschange[j][k] == 0)) {
								// If the land use system has not changed from the intitial system,
								// write the value of the initial map
								fprintf(f9, "%10.4f\n", initoutput[l][j][k]);
							} else {
								fprintf(f9, "%10.4f\n", outputmatrix[mat_newco[j][k]][l][year]);
							}
						} else {
							fprintf(f9, "-9999\n");
						}
					}
				}
				fclose(f9);
			}
		}

		/* write age file */
		if (arcview == 3) {
			sprintf(filenam6, "age.%d.asc", year);
		} else {
			sprintf(filenam6, "age.%d", year);
		}
		if ((f9 = fopen(filenam6, "w")) == NULL) {
			fprintf(flog, "\n not able to write file %s\n", filenam6);
			fclose(flog);
			show_error();
			exit(0);
		}
		if (arcview > 0) {
			fprintf(f9, "ncols    %d\n", rnc);
			fprintf(f9, "nrows    %d\n", rnr);
			fprintf(f9, "xllcorner  %f\n", xllcorner);
			fprintf(f9, "yllcorner  %f\n", yllcorner);
			fprintf(f9, "cellsize   %f\n", cellsize);
			fprintf(f9, "NODATA_value   -9999\n");
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				fprintf(f9, "%d\n ", mat_age[j][k]);
			}
		}
		fclose(f9);

		/* influence */
		if (influence > 1) {
			for (i = 0; i < rncov; i++) {
				if (arcview == 3) {
					sprintf(filenam4, "infl_%d.%d.asc", i, (year - 1));
				} else {
					sprintf(filenam4, "infl_%d.%d", i, (year - 1));
				}
				if ((f7 = fopen(filenam4, "w")) == NULL) {
					fprintf(flog, "\n not able to write file influence file \n");
					fclose(flog);
					show_error();
					exit(0);
				}
				if (arcview > 0) {
					fprintf(f7, "ncols    %d\n", rnc);
					fprintf(f7, "nrows    %d\n", rnr);
					fprintf(f7, "xllcorner  %f\n", xllcorner);
					fprintf(f7, "yllcorner  %f\n", yllcorner);
					fprintf(f7, "cellsize   %f\n", cellsize);
					fprintf(f7, "NODATA_value   -9999\n");
				}
				for (j = 0; j < rnr; j++) {
					for (k = 0; k < rnc; k++) {
						fprintf(f7, "%6.2f\n ", (mat_at[(rnatt + i - allowmin)][j][k]));
					}
				}
				fclose(f7);
			}
		}
	}
	/* probability maps */
	if ((probmaps == 1) || (allprob == 1)) {
		for (i = 0; i < rncov; i++) {
			if (arcview == 3) {
				sprintf(filenam4, "prob1_%d.%d.asc", i, year);
			} else {
				sprintf(filenam4, "prob1_%d.%d", i, year);
			}
			if ((f7 = fopen(filenam4, "w")) == NULL) {
				fprintf(flog, "\n not able to write file \n");
				fclose(flog);
				show_error();
				exit(0);
			}
			if (arcview > 0) {
				fprintf(f7, "ncols    %d\n", rnc);
				fprintf(f7, "nrows    %d\n", rnr);
				fprintf(f7, "xllcorner  %f\n", xllcorner);
				fprintf(f7, "yllcorner  %f\n", yllcorner);
				fprintf(f7, "cellsize   %f\n", cellsize);
				fprintf(f7, "NODATA_value   -9999\n");
			}
			for (j = 0; j < rnr; j++) {
				for (k = 0; k < rnc; k++)
					fprintf(f7, "%9.7f\n ", (mat_co[i][j][k]));
			}
			fclose(f7);
		}
	}
}

void calc_age()
{
	int i, j, k, nl, dum;
	FILE *f11, *f19;
	bool same;
	char filenam6[20];

	if ((year == 0) && (agemode > 0)) {
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				mat_age[j][k] = (rand() % (startage)) + 1;
			}
		}

		if (arcview == 3) {
			sprintf(filenam6, "age.%d.asc", year);
		} else {
			sprintf(filenam6, "age.%d", year);
		}
		if ((f19 = fopen(filenam6, "w")) == NULL) {
			fprintf(flog, "\n not able to write file %s\n", filenam6);
			fclose(flog);
			show_error();
			exit(0);
		}
		if (arcview > 0) {
			fprintf(f19, "ncols    %d\n", rnc);
			fprintf(f19, "nrows    %d\n", rnr);
			fprintf(f19, "xllcorner  %f\n", xllcorner);
			fprintf(f19, "yllcorner  %f\n", yllcorner);
			fprintf(f19, "cellsize   %f\n", cellsize);
			fprintf(f19, "NODATA_value   -9999\n");
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				fprintf(f19, "%d\n ", mat_age[j][k]);
			}
		}
		fclose(f19);

	}
	if ((year == 0) && (agemode == 0)) {
		if (arcview == 3) {
			if ((f11 = fopen("age.0.asc", "r")) == NULL) {
				fprintf(flog, "no file: age.0.asc");
				fclose(flog);
				show_error(); exit(0);
			}
		} else {
			if ((f11 = fopen("age.0", "r")) == NULL) {
				fprintf(flog, "no file: age.0");
				fclose(flog);
				show_error();
				exit(0);
			}
		}
		if (arcview > 0) {
			nl = 0;
			while (nl < 6) {
				dum = fgetc(f11);
				if (dum == '\n') nl++;
			}
		}
		for (i = 0; i < rnr; i++) {
			for (j = 0; j < rnc; j++) {
				fscanf(f11, "%d", &mat_age[i][j]);
			}
		}
		fclose(f11);
	}

	if (year > 0) {
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				if (region[j][k] > -9900) {
					same = false;
					if (mat_newco[j][k] == mat_oldco[j][k]) {
						mat_age[j][k] += 1;
						same = true;
					}
					if (same == false) {
						mat_age[j][k] = 1;
						if (lusconvmap == 1)
							luschange[j][k] = 1;  /* indicate in luschange matrix that the lus is no longer the initial one in order to steer output variable writing */
					}
				} else {
					mat_age[j][k] = -9999;
				}
			}
		}
	}
}

void calc_neigh()
{
	int i, j, k, m, n, tel, tempval, teltypes, s[NCOV], t[NCOV], u, radius[NCOV];
	float cumvalue[NCOV], novalue, neighmat[NCOV][50][50];
	FILE *f13;

	//read kernel definition matrix//

	if ((f13 = fopen("neighmat.txt", "r")) == NULL) {
		fprintf(flog, "no file: neighmat.txt");
		fclose(flog);
		show_error();
		exit(0);
	}
	fprintf(flog, "\n NEIGHBORHOOD CONFIGURATION\n");
	fprintf(flog, "\n weights for neighborhood function: ");

	for (i = 0; i < rncov; i++) {
		fscanf(f13, "%f", &neighweight[i]);
		fprintf(flog, " %3.2f ", neighweight[i]);
	}
	for (u = 0; u < rncov; u++) {
		for (tel = 0; tel < rncov; tel++) {
			allocneigh[u][tel] = 0;
		}
	}

	for (u = 0; u < rncov; u++) {
		if (neighweight[u] > 0) {
			fscanf(f13, "%d", &teltypes);
			fprintf(flog, "\n land system %d: %d different land use systems receive neighborhood addition", u, teltypes);
			for (tel = 0; tel < teltypes; tel++) {
				fscanf(f13, "%d", &tempval);
				allocneigh[u][tempval] = 1;
			}
			fscanf(f13, "%d", &radius[u]);
			fprintf(flog, "\n kernel size: %d\n", radius[u]);
			s[u] = -1 * radius[u];
			t[u] = radius[u] + 1;
			for (i = s[u]; i < t[u]; i++) {
				for (j = s[u]; j < t[u]; j++) {
					fscanf(f13, "%f", &neighmat[u][radius[u] + i][radius[u] + j]);
					fprintf(flog, "%3.2f ", neighmat[u][radius[u] + i][radius[u] + j]);
				}
				fprintf(flog, "\n");
			}
		}
	}
	fclose(f13);

	//process all grid cells for all covers that have neighweight larger than 0
	for (u = 0; u < rncov; u++) {
		if (neighweight[u] > 0) {
			for (j = 0; j < rnr; j++) {
				for (k = 0; k < rnc; k++) {
					if (region[j][k] == -9999) {
						// exclude no-data region
						mat_at[rnatt + u - allowmin][j][k] = -9999;
					} else {
						// process all cells with a value
						novalue = 0;    // initialize neighborhood counters
						cumvalue[u] = 0;
						for (m = s[u]; m < t[u]; m++) {
							for (n = s[u]; n < t[u]; n++) {
								if (((j + m) >= 0) && ((j + m) < rnr) && ((k + n) >= 0) && ((k + n) < rnc) && (region[j + m][k + n] > -9999)) {
									// kernel inside grid with value
									novalue += neighmat[u][radius[u] + m][radius[u] + n];
									if (mat_oldco[j + m][k + n] == u) {
										cumvalue[u] += neighmat[u][radius[u] + m][radius[u] + n];
									}
								}
							}
						}
						if (novalue > 0) {
							mat_at[(rnatt + u - allowmin)][j][k] = (cumvalue[u] / novalue);   //calculates the fraction of the neighborhood occupied by the land system types
						} else {
							mat_at[(rnatt + u - allowmin)][j][k] = 0;
						}
					}
				}
			}
		}
	}
}

void check_file()
{
	int i, j, k, ercount;

	fprintf(flog, "\n\n-- Results of file check --");
	ercount = 0;
	checkfile = 0;
	for (j = 0; j < rnr; j++) {
		for (k = 0; k < rnc; k++) {
			if (region[j][k] > -9999) {
				ercount = 0;
				for (i = 0; i < rncov; i++) {
					if (mat_co[i][j][k] == -9999) {
						ercount += 1;
					}
				}
				if (ercount > 0) {
					fprintf(flog, "=>  %d errors in file cov_all.0\n", ercount);
					checkfile = 1;
				}
				for (i = 0; i < rnatt; i++) {
					if (mat_at[i][j][k] == -9999) {
						fprintf(flog, "  error in file sc1gr%d.fil\n", i);
						checkfile = 1;
					}
				}

			}
		}
	}

	if (checkfile == 0) {
		fprintf(flog, " => no errors found in input files -\n\n");
	}
	if (checkfile == 1) {
		fprintf(flog, "\n\n application terminated because of file errors");
		fclose(flog);
		show_error();
		exit(0);
	}
}

/*
 * This module determines the lowest number of the sc1gr grids used in the conversion matrix.
 * Result is stored in global variable 'allowmin'.
 */
void init_allow()
{
	int i, j, k;
	FILE* fa9;

	allowmin = rnatt;
	if ((fa9 = fopen("allow.txt", "r")) == NULL) {
		fprintf(flog, "no file: allow.txt");
		fclose(flog);
		show_error();
		exit(0);
	}

	for (i = 0; i < rncov; i++) {
		for (j = 0; j < rncov; j++) {
			fscanf(fa9, "%d", &k);
			if ((k > 1) && (k < 100)) {
				if (k < allowmin) {
					allowmin = k;
				}
			}
		}
	}
	fprintf(flog, "\n - lowest sc1gr number in conversion matrix is: %d\n", allowmin);
	fclose(fa9);
}

void unfinished()
{
	FILE *ferror;
	char ername[29];

	fprintf(flog, "\nERROR: no solution; program terminated");
	year = (end - start);
	fclose(flog);
	if (g_argc < 1) {
		show_error();
	} else {
		time_t now = time(0);
		struct tm time_struct = *localtime(&now);
		sprintf(ername, "error%2d%02d.txt", time_struct.tm_hour, time_struct.tm_min);
		if ((ferror = fopen(ername, "w")) == NULL) {
			show_error();
			exit(0);
		}
		fprintf(ferror, "ERROR REPORT\n\n");
		fprintf(ferror, "error at: %2d:%02d\n", time_struct.tm_hour, time_struct.tm_min);
		std::strcpy(ername, g_argv[0]);
		fprintf(ferror, "directory of error: %s", ername);
		fclose(ferror);
		//Beep(1800, 2800);
	}
	exit(0);
}

/*
 * Reads the initial composition of the land use systems in terms of
 * land use system output characteristics.
 */
void load_ini_output()
{
	FILE *f3;
	int i, j, k, nl, dum;
	float l;
	char filename[20];

	for (i = 0; i < nooutputmaps; i++) {
		if (arcview == 3) {
			sprintf(filename, "initoutputmap.%d.asc", i);
		} else {
			sprintf(filename, "initoutputmap.%d", i);
		}
		if ((f3 = fopen(filename, "r")) == NULL) { fprintf(flog, " initoutputmap* file missing \n"); fclose(flog); show_error(); exit(0); }
		if (arcview > 0) {
			nl = 0;
			while (nl < 6) {
				dum = fgetc(f3);
				if (dum == '\n') nl++;
			}
		}
		for (j = 0; j < rnr; j++) {
			for (k = 0; k < rnc; k++) {
				luschange[j][k] = 0;    /* set before the first year the luschange matrix to 0, meaning that there has not yet been a land system change */
				fscanf(f3, "%f ", &l);
				initoutput[i][j][k] = l;
			}
		}
		fclose(f3);
	}
}