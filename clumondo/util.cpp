#include "../include/util.h"

#include <iostream>

void show_error()
{
	std::cout << "Error - see LOG.FIL for details" << std::endl;
}

void show_finished()
{
	std::cout << "Finished" << std::endl;
}

void show_usage()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "CLUMondo.exe demand_file region_file [model_mode]" << std::endl;
	std::cout << "  Valid values for the optional parameter [model mode] are:" << std::endl;
	std::cout << "    - 1: use meta model mode" << std::endl;
	std::cout << "    - 2: create all probability maps in advance" << std::endl;
	std::cout << "    - 3: write yearly probability map" << std::endl;
	std::cout << "    - 4: only write initial probability map" << std::endl;
}