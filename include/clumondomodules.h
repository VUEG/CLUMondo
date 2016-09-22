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