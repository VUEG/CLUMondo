/*****************************************************************
 *
 * CLUMondo main calculation loop
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
#include "../include/clumondomodules.h"
#include "../include/util.h"

#include <stdlib.h>
#include <iostream>
#include <time.h>

int main(int argc, char* argv[])
{
	if (argc <= 2) {
		show_usage();
		exit(1);
	}

	g_argc = argc;
	g_argv = argv;

	// Convert command line switches to global state variables
	strcpy(demd, argv[1]);
	strcpy(park, argv[2]);

	if (argc > 3) {
		if (argv[3] == std::string("1")) {
			metamod = 1;
		}
		if (argv[3] == std::string("2")) {
			allprob = 1;
		}
		if (argv[3] == std::string("3")) {
			probmaps = 1;
		}
		if (argv[3] == std::string("4")) {
			probmaps = 2;
		}
	} else {
		probmaps = 0; metamod = 0; allprob = 0;
	}

	//bool selectcheck;
	int /*i,*/ stepper, stepsize;

	/*
	pausebutton->Visible = true;
	Image3->Visible = false;
	selectcheck = false;
	if (FileCheckon1->Checked == false)
		checkfile = 0;
	for (int i = 0; i < FileListBox2->Items->Count; i++)
	{
		if (FileListBox2->Selected[i])
			selectcheck = true;
	}
	if ((selectcheck == false) && (ParamCount() < 2))
		bbar->SimpleText = "select demand";
	if ((selectcheck == false) && (ParamCount() > 1))
	{
		strcpy(demd, ParamStr(1).c_str());
		selectcheck = true;
	}
	*/
	/*
	if (g_argc > 1) {
		strcpy(demd, argv[1]);
		//selectcheck = true;
	}
	*/

	/*
	if (selectcheck == true)
	{
		selectcheck = false;
		for (i = 0; i < FileListBox1->Items->Count; i++)
		{
			if (FileListBox1->Selected[i])
				selectcheck = true;
		}
		if ((selectcheck == false) && (ParamCount() < 2))
			bbar->SimpleText = "select region";
		if ((selectcheck == false) && (ParamCount() > 1))
		{
			strcpy(park, ParamStr(2).c_str());
			if (ParamCount() > 2)
			{
				if (ParamStr(3) == "1")
					metamod = 1;
				if (ParamStr(3) == "2")
					allprob = 1;
				if (ParamStr(3) == "4")
					probmaps = 2;
				if (ParamStr(3) == "3")
					probmaps = 1;
				if (ParamStr(4) == "beep")
					Beep(1000, 500);
			}
			selectcheck = true;
			checkfile = 0;

		}
	}

	if (selectcheck == true)
	{
		BitBtn1->Visible = true;
		for (int i = 0; i < FileListBox1->Items->Count; i++)	 /* show list of variables in box 1: Area restrictions* /
		{
			if (FileListBox1->Selected[i])
				strcpy(park, FileListBox1->Items->Strings[i].c_str());
		}
		for (i = 0; i < FileListBox2->Items->Count; i++)		 /* show list of demand files in box 2* /
		{
			if (FileListBox2->Selected[i])
				strcpy(demd, FileListBox2->Items->Strings[i].c_str());
		}
	*/

	year = 0;
	allowmin = 0;
	//YearBar->Position = 0;
	if ((flog = fopen(F_LOG, "w")) == NULL) {
		/* open log-file */
		show_error();
		fclose(flog);
		exit(0);
	}
	//bbar->SimpleText = "Loading...";
	std::cout << "Loading..." << std::endl;
	//Application->ProcessMessages();

	//READING OF ALL INPUT FILES AT START SIMULATIONS
	all_init();				/* read main simulation parameters from input file */
	stepsize = (4000000 / (rncov*rnr*rnc));
	if (stepsize < 1) {
		stepsize = 1;
	}
	if (stepsize > 30) {
		stepsize = 30;
	}
	if (metamod == 1) {
		// In case of a meta-model run determine the minimum value of the sc1gr grids used
		init_allow();
	}
	make_mat();				/* allocates memory for matrices and checks success */
	calc_age();				/* initializes age of land use */
	load_grid();			/* loads explanatory factor grids into memory */
	demand_read();			/* read the demand file for all years */
	load_lusconv();			/* read the lusconv matrix with hierarchy of land use systems in terms of demand types */
	if (changeconv != 1) {
		// read the lusmatrix with values of all land use systems in terms of demand types
		load_lusmatrix();
	}
	load_reg();				/* read for the different land use types the regression equations */
	if (influence > 0) {	/* read the regression equations for neighborhoods if option is on */
		load_reg2();
	}
	read_allowed();			/* read the allowed changes between land use types */
	load_region();			/* loads region distribution from regionfile */
	if (lusconvmap == 1) {
		// loads the initial composition of the lusconv output maps
		load_ini_output();
	}
	//PerformanceGraph1->Visible = true;

	// Year loop; needs to be connected to start and end
	for (year = 1; year <= (end - start); year++) {
		std::cout << " " << std::endl;
		std::cout << "Now simulating year: " << year << std::endl;
		std::cout << "Updating input information..." << std::endl;

		fprintf(flog, "\n\n\n______   year %d  ______\n", year);

		if (allprob == 0) {
			// Read coverage of land use types for year - 1 from file
			set_oldco();
		}
		if ((locationfactor > 0) && (metamod == 0)) {
			// Read the locationspecific addition if requested
			load_locationfactor();
		}
		if (nochange > 0) {
			// Replace scgr (explanatory factor) files that have changed by current files
			scgr_change();
		}
		if (checkfile == 1) {
			check_file();
		}
		if (changeconv == 1) {
			load_lusmatrix();	   /* read the lusmatrix with values of all land use systems in terms of demand types */
		}
		demand_dir();			  /* loads demands for land use types and determines if demand increases or decreases */
		if (influence > 0) {
			// Calculate neighborhood attributes
			calc_neigh();
			if (influence == 2) {
				write_grid();
				fprintf(flog, "\n influence calculations made");
				fclose(flog);
				//finishbox->ShowModal();
				show_finished();
				exit(0);
			}
		}

		// Initializes iteration parameters
		init_iter();

		std::cout << "Calculating preference maps...." << std::endl;

		calc_reg();			  /* calculate probabilities for each land use type for each grid from regression equations */

		if (probmaps == 1) {
			write_grid();
			loop = 20001;
			year = (end - start);
		}
		std::cout << "Iterating....";
		// Start of iteration for numerical solution
		if ((allprob == 0) && (probmaps < 2)) {
			stepper = 0;
			do {
				loop++;
				std::cout << ".";
				calc_change_ch();	// calculates change in land use types
				comp_demand();		// compares demand with allocated land use, adjust elas (iteration parameter)
				stepper++;
				if (stepper == stepsize) {
					/*
					if ((maxdiff>demdiffmax) && ((totdiff / (rncov)) < demdiff)) {
						maxdevlabel->Visible = true;
					} else {
						maxdevlabel->Visible = false;
					}
					*/
					stepper = 0;
					//Iterbar->SelEnd = Barpos;
					//PerformanceGraph1->DataPoint(clYellow, (Barpos + 1));
					//PerformanceGraph1->Update();
					//bbar->SimpleText = loop;
					//Application->ProcessMessages();
					if (loop == 5000) {
						calc_reg();
					}
				}
			} while ((((totdiff / (rndemand)) > demdiff) && (loop < 20000)) || ((loop < 20000) && (maxdiff > demdiffmax)));   /* end of loop for iteration of allocation on demand */

			if ((loop < 20000) && (bottomup == 1)) {
				// Calculate autonomous change */
				autonomous_change();
			}

			if (loop == 20000) {
				unfinished();
			}

			if (loop < 20000) {
				// Update number of time-steps with stable land use type
				calc_age();
			}

			if ((probmaps == 0) && (loop < 20000)) {
				// Write results in output files
				write_grid();
			}
		}
		// END OF ITERATION
		if (allprob == 1) {
			write_grid();
		}
		/*
		Yearprog = (((year)* 100) / (end - start));
		if (loop < 20000)
			YearBar->Position = Yearprog;
			*/
	}
	// END OF YEAR LOOP

	free_mat();				 /* empty memory */

	time_t now = time(NULL);
	struct tm time_struct = *localtime(&now);
	fprintf(flog, "\n\n-- End of simulation at: %2d:%02d:%02d -- \n\n\n", time_struct.tm_hour, time_struct.tm_min, time_struct.tm_sec);

	fclose(flog);
}
