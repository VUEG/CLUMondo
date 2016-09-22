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
#pragma once

// I/O functions
void all_init();
void load_reg();
void load_grid();
void load_region();
void write_grid();
void load_reg2();
void load_locationfactor();
void read_allowed();
void demand_read();
void init_allow();
void load_lusconv();
void load_lusmatrix();
void load_ini_output();
void check_file();

// Storage management
void make_mat();
void free_mat();

// Model functions
void calc_reg();
void calc_change_ch();
void comp_demand();
void scgr_change();
void demand_dir();
void set_oldco();
void init_iter();
void calc_age();
void calc_neigh();
void autonomous_change();
void unfinished();